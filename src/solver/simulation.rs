use rayon::prelude::*;

use simplelog::warn;

use crate::error::SolverError;
use crate::solver::hydraulicsolver::HydraulicSolver;
use crate::solver::state::SolverState;
use crate::solver::result::SolverResult;
use crate::model::network::Network;
use crate::model::node::NodeType;
use crate::model::units::{Cfs, Ft3};
use crate::model::options::DemandModel;
use crate::model::control::ControlCondition;

pub struct Simulation<'a> {
  pub network: &'a Network,
  pub solver: HydraulicSolver,
  pub state: SolverState,
  pub time: usize,
  pub skip_timesteps: bool,
}

impl<'a> Simulation<'a> {

  /// Creates a new simulation, allocating the hydraulic solver and initializing state.
  /// Functionally equivalent to EN_openH followed by EN_initH.
  pub fn new(network: &'a Network) -> Result<Self, SolverError> {
    let solver = HydraulicSolver::new(network)?;
    let state = SolverState::new_with_initial_values(network);
    Ok(Self { network, solver, state, time: 0, skip_timesteps: true })
  }

  /// Resets the simulation to initial conditions (time = 0, fresh state).
  /// Equivalent to EN_initH.
  pub fn initialize_hydraulics(&mut self) {
    self.state = SolverState::new_with_initial_values(self.network);
    self.time = 0;
  }

  /// Applies patterns and controls, then solves hydraulics at the given time.
  /// Returns the solved state at the given time.
  /// Equivalent to EN_runH.
  pub fn run_hydraulics(&mut self, time: usize) -> Result<usize, SolverError> {
    self.time = time;
    Self::apply_patterns(self.network, &mut self.state, self.time);
    Self::apply_controls(self.network, &mut self.state, self.time);
    self.state = self.solver.solve(self.network, &self.state)?;
    Ok(self.time)
  }

  /// Advances to the next hydraulic time step (updating tank levels).
  /// Returns the time step length in seconds. Returns 0 when the simulation is complete.
  /// Equivalent to EN_nextH.
  pub fn next_hydraulic_timestep(&mut self) -> usize {
    let duration = self.network.options.time_options.duration;
    if self.time >= duration {
      return 0;
    }
    let remaining = duration - self.time;
    let timestep = Self::calculate_time_step(self.network, &self.state, self.time, self.skip_timesteps).min(remaining);
    Self::update_tanks(self.network, &mut self.state, timestep);
    self.time += timestep;
    timestep
  }

  /// Runs a complete hydraulic simulation, collecting results at each report step.
  /// Equivalent to EN_solveH.
  /// Simulations can be run in parallel or sequentially. Parallel mode is only supported for networks without tanks or pressure controls.
  /// Returns a SolverResult containing the flows and heads at each report step.
  pub fn solve_hydraulics(&mut self, parallel: bool) -> Result<SolverResult, SolverError> {
    self.initialize_hydraulics();

    if parallel && !self.network.has_tanks() && !self.network.has_pressure_controls() {
      return self.solve_parallel();
    } else if parallel {
      warn!("Networks with tanks or pressure controls cannot be solved in parallel, running sequentially");
    }

    let report_timestep = self.network.options.time_options.report_timestep;
    let report_steps = (self.network.options.time_options.duration / report_timestep) + 1;

    // Initialize the results vector
    let mut results = SolverResult::new(self.network.links.len(), self.network.nodes.len(), report_steps);

    // Run the simulation step by step
    let mut time = 0;
    loop {
      let t = self.run_hydraulics(time)?;
      // Append the results to the SolverResult if the time step is a report step
      if t % report_timestep == 0 {
        results.append(&self.state, t / report_timestep);
      }
      // Advance the time step, returns 0 if the simulation is complete
      let dt = self.next_hydraulic_timestep();
      if dt == 0 { break; }
      time += dt;
    }

    // Convert the units of the results to the original units
    results.convert_units(&self.network.options.flow_units, &self.network.options.unit_system);
    Ok(results)
  }

  /// Solves the hydraulics in parallel, only supported for networks without tanks or pressure controls.
  fn solve_parallel(&mut self) -> Result<SolverResult, SolverError> {

    // parallel solve is only supported for networks without tanks or pressure controls, so we only use the report timestep to determine the number of steps
    let report_timestep = self.network.options.time_options.report_timestep;
    let report_steps = (self.network.options.time_options.duration / report_timestep) + 1;

    // Initialize the results vector
    let mut results = SolverResult::new(self.network.links.len(), self.network.nodes.len(), report_steps);

    // First solve the first time step, and use that as the initial state for the parallel solve
    // allowing for faster convergence.
    Self::apply_patterns(self.network, &mut self.state, 0);
    self.state = self.solver.solve(self.network, &self.state)?;
    results.append(&self.state, 0);

    let initial_state = self.state.clone();

    // Solve the remaining time steps in parallel
    let par_results: Vec<Result<SolverState, SolverError>> = (1..report_steps).into_par_iter().map(|step| {
      let mut state = initial_state.clone();
      Self::apply_patterns(self.network, &mut state, step * report_timestep);
      self.solver.solve(self.network, &state)
    }).collect();

    // update the results vector with the parallel solve results
    for (step, step_result) in par_results.into_iter().enumerate() {
      results.append(&step_result?, step + 1);
    }

    // Convert the units of the results to the original units
    results.convert_units(&self.network.options.flow_units, &self.network.options.unit_system);
    Ok(results)
  }

  /// Applies demand and head patterns to the state at the given time.
  /// TODO: Add support for multiple patterns
  fn apply_patterns(network: &Network, state: &mut SolverState, time: usize) {
    let time_options = &network.options.time_options;
    let pattern_time = time_options.pattern_start + time;
    let pattern_index = pattern_time / time_options.pattern_timestep;

    for (i, node) in network.nodes.iter().enumerate() {
      let Some(head_pattern) = node.head_pattern() else { continue };
      let pattern = &network.patterns[head_pattern];
      state.heads[i] = node.elevation * pattern.multipliers[pattern_index % pattern.multipliers.len()];
    }

    // Resolve the default demand pattern: explicit PATTERN option, or EPANET's default "1"
    let default_pattern_id: Box<str> = network.options.pattern.clone()
      .unwrap_or_else(|| "1".into());
    let default_pattern = network.patterns.get(&default_pattern_id);

    state.demands = network.nodes.iter().map(|n| {
      let NodeType::Junction(junction) = &n.node_type else { return 0.0 };
      let pattern = junction.pattern.as_ref()
        .and_then(|id| network.patterns.get(id))
        .or(default_pattern);
      let multiplier = match pattern {
        Some(p) => p.multipliers[pattern_index % p.multipliers.len()],
        None => 1.0,
      };
      junction.basedemand * multiplier * network.options.demand_multiplier
    }).collect::<Vec<Cfs>>();

    if network.options.demand_model == DemandModel::PDA {
      state.demand_flows = state.demands.clone();
    }
  }

  /// Applies controls to the state at the given time.
  fn apply_controls(network: &Network, state: &mut SolverState, time: usize) {
    let clocktime = (time + network.options.time_options.start_clocktime) % (24 * 3600);

    for control in &network.controls {
      // skip controls that are not pressure or level controls
      if matches!(control.condition, ControlCondition::LowPressure { .. } | ControlCondition::HighPressure { .. }) {
        continue;
      }
      if matches!(control.condition, ControlCondition::LowLevel { .. } | ControlCondition::HighLevel { .. }) && time == 0 {
        continue;
      }
      if control.is_active(state, network, time, clocktime) {
        control.activate(state, network);
      }
    }
  }

  /// Calculates the time step for the next hydraulic time step.
  /// Returns the time step length in seconds. Returns 0 when the simulation is complete.
  fn calculate_time_step(network: &Network, state: &SolverState, current_time: usize, skip_timesteps: bool) -> usize {
    let times = &network.options.time_options;
    let clocktime = (current_time + times.start_clocktime) % (24 * 3600);

    // if there are no controls, tanks, or quality, and the skip_timesteps flag is set, return the report timestep
    if network.controls.is_empty() && !network.has_tanks() && !network.has_quality() && skip_timesteps {
      return times.report_timestep;
    }

    // calculate the time to reach the next report step
    let t_report = times.report_timestep - (current_time % times.report_timestep);
    // calculate the time to reach the next pattern step
    let t_pattern = times.pattern_timestep - (current_time % times.pattern_timestep);

    // calculate the time to reach the next tank level
    let t_tanks = network.nodes.iter()
      .zip(state.heads.iter())
      .zip(state.demands.iter())
      .map(|((node, head), demand)| {
        let NodeType::Tank(tank) = &node.node_type else { return usize::MAX };
        let level = head - node.elevation;
        tank.time_to_fill_or_drain(level, *demand)
      })
      .min().unwrap_or(usize::MAX);

    // calculate the time to reach the next control event
    let t_controls = network.controls.iter()
      .map(|control| {
        match &control.condition {
          ControlCondition::Time { seconds } => {
            if *seconds < current_time { return usize::MAX; }
            seconds - current_time
          }
          ControlCondition::ClockTime { seconds } => {
            if *seconds < clocktime {
              return ((*seconds as i32 - clocktime as i32) + (24 * 3600)) as usize
            }
            seconds - clocktime
          }
          ControlCondition::LowLevel { tank_index, target } | ControlCondition::HighLevel { tank_index, target } => {
            let node = &network.nodes[*tank_index];
            if let NodeType::Tank(tank) = &node.node_type {
              let level = state.heads[*tank_index] - node.elevation;
              tank.time_to_reach_level(level, *target, state.demands[*tank_index])
            } else {
              usize::MAX
            }
          }
          _ => usize::MAX,
        }
      })
      .filter(|&x| x > 0)
      .min().unwrap_or(usize::MAX);

    // return the minimum of the calculated time steps
    *[t_report, t_pattern, t_tanks, t_controls, times.hydraulic_timestep].iter().filter(|&x| *x > 0).min().unwrap()
  }

  /// Updates the tank levels for the given time step.
  fn update_tanks(network: &Network, state: &mut SolverState, timestep: usize) {
    for (tank_index, node) in network.nodes.iter().enumerate() {
      if let NodeType::Tank(tank) = &node.node_type {
        // calculate the delta volume and new head for the tank
        let delta_volume = state.demands[tank_index] * timestep as Ft3;
        let new_head = tank.new_head(delta_volume, state.heads[tank_index]);
        state.heads[tank_index] = new_head;
      }
    }
  }
}
