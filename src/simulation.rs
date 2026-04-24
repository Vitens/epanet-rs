//! High-level extended-period simulation driver for sequential and parallel hydraulic solves.

/*!
The simulation driver is responsible for managing the simulation workflow for sequential and parallel solving of hydraulic networks and collecting results at each report step.

The simulation driver updates the network state with the changes to the network, applies patterns and controls, updates tank levels and advances the time step until the simulation is complete.

Example:
```rust no_run
use epanet_rs::simulation::Simulation;
use epanet_rs::model::network::modify::PipeUpdate;

# fn main() -> Result<(), Box<dyn std::error::Error>> {
// load the network from a file
let mut simulation = Simulation::from_file("tests/pump.inp")?;
// run a complete simulation
simulation.solve_hydraulics(false)?;
// modify the network
simulation.network.update_pipe("P1", &PipeUpdate {
  roughness: Some(0.2),
  ..Default::default()
})?;
// run the simulation again
simulation.solve_hydraulics(false)?;
# Ok(())
# }
```
 */

use rayon::prelude::*;

use simplelog::{debug, warn};

use crate::error::{InputError, SolverError};
use crate::solver::hydraulicsolver::HydraulicSolver;
use crate::solver::state::SolverState;
use crate::solver::result::SolverResult;
use crate::model::network::Network;
use crate::model::node::NodeType;
use crate::model::units::FlowUnits;
use crate::model::options::HeadlossFormula;
use crate::model::options::SimulationOptions;

use crate::model::control::ControlCondition;


pub struct Simulation {
  pub network: Network,
  pub solver: Option<HydraulicSolver>,
  pub state: Option<SolverState>,
  pub time: usize,
  pub skip_timesteps: bool,
  /// flag to indicate if the simulation is solved
  pub solved: bool,
}

impl Simulation {

  /// Creates a new simulation from an input file.
  pub fn from_file(path: &str) -> Result<Self, InputError> {
    let network = Network::from_file(path)?;
    Ok(Self::new(network))
  }
  /// Creates a new simulation with the given network and options.
  pub fn init(flow_units: FlowUnits, headloss_formula: HeadlossFormula) -> Self {
    let mut network = Network::default();
    network.options = SimulationOptions::new(flow_units, headloss_formula);
    Self::new(network)
  }

  /// Creates a new simulation, allocating the hydraulic solver and initializing state.
  /// Functionally equivalent to EN_openH followed by EN_initH.
  pub fn new(network: Network) -> Self {
    Self { network, solver: None, state: None, time: 0, skip_timesteps: true, solved: false }
  }

  /// Resets the simulation to initial conditions (time = 0, fresh state).
  /// Equivalent to EN_initH.
  pub fn initialize_hydraulics(&mut self) -> Result<(), SolverError> {
    // update the solver if the network version has changed
    self.solver = Some(HydraulicSolver::new(&self.network)?);
    self.state = Some(SolverState::new_with_initial_values(&self.network));
    self.network.reset_changes();
    self.time = 0;
    self.solved = false;
    Ok(())
  }

  /// Applies patterns and controls, then solves hydraulics at the given time.
  /// Returns the solved state at the given time.
  /// Equivalent to EN_runH.
  pub fn run_hydraulics(&mut self) -> Result<usize, SolverError> {
    // check if the solver is initialized and that its topology matches the network
    let needs_reinit = match self.solver.as_ref() {
      Some(solver) => self.network.topology_version != solver.topology_version,
      None => return Err(SolverError::NotInitialized),
    };
    if needs_reinit {
      self.initialize_hydraulics()?;
    }

    let solver = self.solver.as_ref().unwrap();
    let state = self.state.as_mut().unwrap();

    // update the state with the changes to the network
    state.update_with_network_changes(&mut self.network);

    state.apply_patterns(&self.network, self.time);
    state.apply_controls(&self.network, self.time);
    self.state = Some(solver.solve(&self.network, state)?);
    // set the solved flag to true
    self.solved = true;
    Ok(self.time)
  }

  /// Advances to the next hydraulic time step (updating tank levels).
  /// Returns the time step length in seconds. Returns 0 when the simulation is complete or the state is not initialized.
  /// Equivalent to EN_nextH.
  pub fn next_hydraulic_timestep(&mut self) -> usize {

    // unwrap the state, if it is none, return 0
    let Some(state) = self.state.as_mut() else { return 0 };

    let duration = self.network.options.time_options.duration;
    if self.time >= duration {
      return 0;
    }
    let remaining = duration - self.time;
    let timestep = Self::calculate_time_step(&self.network, &state, self.time, self.skip_timesteps).min(remaining);
    state.update_tanks(&self.network, timestep);
    self.time += timestep;
    timestep
  }

  /// Runs a complete hydraulic simulation, collecting results at each report step.
  /// Equivalent to EN_solveH.
  /// Simulations can be run in parallel or sequentially. Parallel mode is only supported for networks without tanks or pressure controls.
  /// Returns a SolverResult containing the flows and heads at each report step.
  pub fn solve_hydraulics(&mut self, parallel: bool) -> Result<SolverResult, SolverError> {
    // initialize the solver if it is not initialized
    if self.solver.is_none() {
      self.initialize_hydraulics()?;
    }

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
    loop {
      let t = self.time;
      debug!("Running hydraulics at hour: {}:{:02}:{:02}", t / 3600, (t % 3600) / 60, t % 60);
      let t = self.run_hydraulics()?;
      // Append the results to the SolverResult if the time step is a report step
      if t % report_timestep == 0 {
        results.append(self.state.as_ref().unwrap(), t / report_timestep);
      }
      // Advance the time step, returns 0 if the simulation is complete
      let dt = self.next_hydraulic_timestep();
      if dt == 0 { break; }
    }
    self.solved = true;

    // Convert the units of the results to the original units
    results.convert_units(&self.network.options.flow_units, &self.network.options.unit_system);
    Ok(results)
  }

  /// Solves the hydraulics in parallel, only supported for networks without tanks or pressure controls.
  fn solve_parallel(&mut self) -> Result<SolverResult, SolverError> {
    let solver = self.solver.as_ref().ok_or(SolverError::NotInitialized)?;
    let state = self.state.as_mut().unwrap();

    // parallel solve is only supported for networks without tanks or pressure controls, so we only use the report timestep to determine the number of steps
    let report_timestep = self.network.options.time_options.report_timestep;
    let report_steps = (self.network.options.time_options.duration / report_timestep) + 1;

    // Initialize the results vector
    let mut results = SolverResult::new(self.network.links.len(), self.network.nodes.len(), report_steps);

    // First solve the first time step, and use that as the initial state for the parallel solve
    // allowing for faster convergence.
    state.apply_patterns(&self.network, 0);
    let initial_state = solver.solve(&self.network, state)?;
    results.append(&initial_state, 0);

    // Solve the remaining time steps in parallel
    let par_results: Vec<Result<SolverState, SolverError>> = (1..report_steps).into_par_iter().map(|step| {
      let mut state = initial_state.clone();
      state.apply_patterns(&self.network, step * report_timestep);
      solver.solve(&self.network, &state)
    }).collect();

    // update the results vector with the parallel solve results
    for (step, step_result) in par_results.into_iter().enumerate() {
      results.append(&step_result?, step + 1);
    }
    self.solved = true;

    // Convert the units of the results to the original units
    results.convert_units(&self.network.options.flow_units, &self.network.options.unit_system);
    Ok(results)
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

  /// get the solved state if it is available, otherwise return None
  pub(crate) fn solved_state(&self) -> Option<&SolverState> {
    if self.solved {
      self.state.as_ref()
    } else {
      None
    }
  }
}
