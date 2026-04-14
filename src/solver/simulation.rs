use rayon::prelude::*;

use simplelog::warn;

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
  pub solver: HydraulicSolver<'a>,
  pub state: SolverState,
  pub time: usize,
  pub skip_timesteps: bool,
}

impl<'a> Simulation<'a> {

  /// Creates a new simulation, allocating the hydraulic solver and initializing state.
  /// Functionally equivalent to EN_openH followed by EN_initH.
  pub fn new(network: &'a Network) -> Self {
    let solver = HydraulicSolver::new(network);
    let state = SolverState::new_with_initial_values(network);
    Self { network, solver, state, time: 0, skip_timesteps: true }
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
  pub fn run_hydraulics(&mut self, time: usize) -> Result<usize, String> {
    self.time = time;
    Self::apply_patterns(self.network, &mut self.state, self.time);
    Self::apply_controls(self.network, &mut self.state, self.time);
    self.state = self.solver.solve(&self.state)?;
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
  pub fn solve_hydraulics(&mut self, parallel: bool) -> SolverResult {
    self.initialize_hydraulics();

    if parallel && (self.network.has_tanks() || self.network.has_pressure_controls()) {
      warn!("Networks with tanks or pressure controls cannot be solved in parallel, running sequentially");
    }

    if parallel && !self.network.has_tanks() && !self.network.has_pressure_controls() {
      return self.solve_parallel();
    }

    let report_timestep = self.network.options.time_options.report_timestep;
    let report_steps = (self.network.options.time_options.duration / report_timestep) + 1;
    let mut results = SolverResult::new(self.network.links.len(), self.network.nodes.len(), report_steps);

    let mut time = 0;
    loop {
      let t = self.run_hydraulics(time).unwrap();
      if t % report_timestep == 0 {
        results.append(&self.state, t / report_timestep);
      }
      let dt = self.next_hydraulic_timestep();
      if dt == 0 { break; }
      time += dt;
    }

    results.convert_units(&self.network.options.flow_units, &self.network.options.unit_system);
    results
  }

  fn solve_parallel(&mut self) -> SolverResult {
    let report_timestep = self.network.options.time_options.report_timestep;
    let report_steps = (self.network.options.time_options.duration / report_timestep) + 1;
    let mut results = SolverResult::new(self.network.links.len(), self.network.nodes.len(), report_steps);

    Self::apply_patterns(self.network, &mut self.state, 0);
    self.state = self.solver.solve(&self.state).unwrap();
    results.append(&self.state, 0);

    let initial_state = self.state.clone();
    let par_results: Vec<SolverState> = (1..report_steps).into_par_iter().map(|step| {
      let mut state = initial_state.clone();
      Self::apply_patterns(self.network, &mut state, step * report_timestep);
      self.solver.solve(&state).unwrap()
    }).collect();

    for (step, step_result) in par_results.iter().enumerate() {
      results.append(step_result, step + 1);
    }

    results.convert_units(&self.network.options.flow_units, &self.network.options.unit_system);
    results
  }

  fn apply_patterns(network: &Network, state: &mut SolverState, time: usize) {
    let time_options = &network.options.time_options;
    let pattern_time = time_options.pattern_start + time;
    let pattern_index = pattern_time / time_options.pattern_timestep;

    for (i, node) in network.nodes.iter().enumerate() {
      let Some(head_pattern) = node.head_pattern() else { continue };
      let pattern = &network.patterns[head_pattern];
      state.heads[i] = node.elevation * pattern.multipliers[pattern_index % pattern.multipliers.len()];
    }

    state.demands = network.nodes.iter().map(|n| {
      let NodeType::Junction(junction) = &n.node_type else { return 0.0 };
      if let Some(pattern_id) = &junction.pattern {
        let pattern = &network.patterns[pattern_id];
        let multiplier = pattern.multipliers[pattern_index % pattern.multipliers.len()];
        return junction.basedemand * multiplier * network.options.demand_multiplier;
      }
      junction.basedemand * network.options.demand_multiplier
    }).collect::<Vec<Cfs>>();

    if network.options.demand_model == DemandModel::PDA {
      state.demand_flows = state.demands.clone();
    }
  }

  fn apply_controls(network: &Network, state: &mut SolverState, time: usize) {
    let clocktime = (time + network.options.time_options.start_clocktime) % (24 * 3600);

    for control in &network.controls {
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

  fn calculate_time_step(network: &Network, state: &SolverState, current_time: usize, skip_timesteps: bool) -> usize {
    let times = &network.options.time_options;
    let clocktime = (current_time + times.start_clocktime) % (24 * 3600);

    if network.controls.is_empty() && !network.has_tanks() && !network.has_quality() && skip_timesteps {
      return times.report_timestep;
    }

    let t_report = times.report_timestep - (current_time % times.report_timestep);
    let t_pattern = times.pattern_timestep - (current_time % times.pattern_timestep);

    let t_tanks = network.nodes.iter()
      .zip(state.heads.iter())
      .zip(state.demands.iter())
      .map(|((node, head), demand)| {
        let NodeType::Tank(tank) = &node.node_type else { return usize::MAX };
        let level = head - node.elevation;
        tank.time_to_fill_or_drain(level, *demand)
      })
      .min().unwrap_or(usize::MAX);

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

    *[t_report, t_pattern, t_tanks, t_controls, times.hydraulic_timestep].iter().filter(|&x| *x > 0).min().unwrap()
  }

  fn update_tanks(network: &Network, state: &mut SolverState, timestep: usize) {
    for (tank_index, node) in network.nodes.iter().enumerate() {
      if let NodeType::Tank(tank) = &node.node_type {
        let delta_volume = state.demands[tank_index] * timestep as Ft3;
        let new_head = tank.new_head(delta_volume, state.heads[tank_index]);
        state.heads[tank_index] = new_head;
      }
    }
  }
}
