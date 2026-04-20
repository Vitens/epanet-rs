use rayon::prelude::*;

use simplelog::warn;

use crate::error::{InputError, SolverError};
use crate::solver::hydraulicsolver::HydraulicSolver;
use crate::solver::state::SolverState;
use crate::solver::result::SolverResult;
use crate::model::network::Network;
use crate::model::node::NodeType;
use crate::model::units::{UnitConversion, FlowUnits};
use crate::model::options::HeadlossFormula;
use crate::model::node::Node;
use crate::model::junction::Junction;
use crate::model::options::SimulationOptions;

use crate::model::control::ControlCondition;

use crate::constants::*;

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
    self.time = 0;
    self.solved = false;
    Ok(())
  }

  /// Applies patterns and controls, then solves hydraulics at the given time.
  /// Returns the solved state at the given time.
  /// Equivalent to EN_runH.
  pub fn run_hydraulics(&mut self) -> Result<usize, SolverError> {
    // check if the solver is initialized
    let solver = self.solver.as_ref().ok_or(SolverError::NotInitialized)?;
    let state = self.state.as_mut().unwrap();

    state.apply_patterns(&self.network, self.time);
    state.apply_controls(&self.network, self.time);
    self.state = Some(solver.solve(&self.network, &state)?);
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

/// Methods to modify the network and solverstate
impl Simulation {

  /// Add a new junction node to the network.
  pub fn add_junction(&mut self, id: &str, elevation: f64, basedemand: f64, pattern: Option<&str>, emitter_coefficient: f64, coordinates: Option<(f64, f64)>) -> Result<(), InputError> {

    if self.solver.is_some() {
      return Err(InputError::TopologyChangeWhileSolverOpen);
    }

    // get the pattern index if it is provided
    let pattern_index = if let Some(pattern) = pattern {
      Some(*self.network.pattern_map.get(pattern).ok_or(InputError::PatternNotFound { pattern_id: pattern.into() })?)
    } else {
      None
    };

    let mut junction = Node {
      id: id.into(),
      elevation: elevation,
      node_type: NodeType::Junction(Junction { basedemand: basedemand, pattern: pattern.map(|s| s.into()), pattern_index: pattern_index, emitter_coefficient: emitter_coefficient }),
      coordinates: coordinates,
    };

    junction.convert_to_standard(&self.network.options);

    self.network.add_node(junction)?;

    Ok(())
  }



  /// Sets the elevation of a node and updates the state heads if the hydraulic simulation is open.
  pub fn set_node_elevation(&mut self, node_id: &str, elevation: f64) {
    let index = self.network.node_map.get(node_id).unwrap();
    let node = &mut self.network.nodes[*index];
    // update elevation
    let value = elevation / self.network.options.unit_system.per_feet();
    
    // update the state heads
    if let Some(state) = self.state.as_mut() {
      match &mut node.node_type {
        NodeType::Tank(_) => state.heads[*index] = (node.elevation - state.heads[*index]) + value,
        NodeType::Reservoir(_) => state.heads[*index] = value,
        _ => (),
      }
    }

    node.elevation = value;
  }

}


#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_network_addition() {
    let network = Network::default();
    let mut simulation = Simulation::init(FlowUnits::CFS, HeadlossFormula::HazenWilliams);
    simulation.add_junction("J1", 100.0, 10.0, None, 0.0, None).unwrap();

    assert_eq!(simulation.network.nodes.len(), 1);
    assert_eq!(simulation.network.nodes[0].id, "J1".into());
    assert_eq!(simulation.network.nodes[0].elevation, 100.0);

  }

  #[test]
  fn test_network_addition_with_units() {
    let mut simulation = Simulation::init(FlowUnits::LPS, HeadlossFormula::DarcyWeisbach);
    simulation.add_junction("J1", 100.0, 10.0, None, 0.0, None).unwrap();

    dbg!(&simulation.network.options);

    assert_eq!(simulation.network.nodes.len(), 1);
    assert_eq!(simulation.network.nodes[0].id, "J1".into());
    assert_eq!(simulation.network.nodes[0].elevation, 100.0 / MperFT); // test conversion to standard units

    let NodeType::Junction(junction) = &simulation.network.nodes[0].node_type else {
      panic!("Node is not a junction");
    };
    assert_eq!(junction.basedemand, 10.0 / LPSperCFS); // test conversion to standard units
  }



  // #[test]
  // fn test_set_node_elevation() {
  //   let network = test_network();
  //   let mut simulation = Simulation::new(network);
  //   simulation.initialize_hydraulics().unwrap();
  //   simulation.set_node_elevation("R1", 123.0);
  //   assert_eq!(simulation.network.nodes[0].elevation, 123.0);
  //   assert_eq!(simulation.state.as_ref().unwrap().heads[0], 123.0);
  // }
}