use faer::sparse::{SparseColMat, Triplet, SymbolicSparseColMat};
use faer::sparse::linalg::LltError;
use faer::linalg::cholesky::llt::factor::LltError::NonPositivePivot;
use faer::sparse::linalg::solvers::{SymbolicLlt, Llt};
use faer::{Mat, Side};
use faer::prelude::*;

use rayon::prelude::*;

use crate::solver::state::SolverState;
use crate::solver::result::SolverResult;
use crate::solver::matrix::{CSCIndex, ResistanceCoefficients, find_csc_index};

use crate::model::units::{Cfs, Ft3};


use simplelog::{warn, debug, error};

use crate::constants::{BIG_VALUE, Q_ZERO, PDA_MIN_DIFF};
use crate::model::node::NodeType;
use crate::model::link::{LinkType, LinkTrait, LinkStatus};
use crate::model::network::Network;
use crate::model::valve::ValveType;
use crate::model::control::ControlCondition;
use crate::model::options::{SimulationOptions, DemandModel};

pub struct FlowBalance {
  pub total_demand: Cfs,
  pub total_supply: Cfs,
  pub error: Cfs,
}


#[derive(Default)]
pub struct IterationStatistics {
  pub sum_dq: Cfs,
  pub sum_q: Cfs,
  pub max_dq: Cfs,
  pub max_dq_index: usize,
  pub status_changed: bool,
}
impl IterationStatistics {
  pub fn relative_change(&self) -> f64 {
    self.sum_dq / (self.sum_q + Q_ZERO)
  }
  pub fn max_dq_converted(&self, options: &SimulationOptions) -> Cfs {
    self.max_dq * options.flow_units.per_cfs()
  }
}

pub struct HydraulicSolver<'a> {
  network: &'a Network, 
  node_to_unknown: Vec<Option<usize>>,
  sparsity_pattern: SymbolicSparseColMat<usize>,
  symbolic_llt: SymbolicLlt<usize>,
  /// precomputed Jacobian matrix
  jac: SparseColMat<usize, f64>,
  /// precomputed CSC indices for the links
  csc_indices: Vec<CSCIndex>,
  /// precomputed indices for the rows of the Jacobian matrix for each node
  node_rows: Vec<Option<usize>>,

  pub skip_timesteps: bool, // set to false to match epanet timestep behaviour
}

impl<'a> HydraulicSolver<'a> {

  pub fn new(network: &'a Network) -> Self {
    // build global unknown-numbering map
    let node_to_unknown = Self::build_unknown_numbering_map(network);
    // generate sparsity pattern
    let sparsity_pattern = Self::build_sparsity_pattern(network, &node_to_unknown);
    // map each link to its CSC (Compressed Sparse Column) indices
    let csc_indices = Self::map_links_to_csc_indices(network, &sparsity_pattern, &node_to_unknown);
    // precompute the indices for the rows of the Jacobian matrix for each node
    let node_rows = Self::map_nodes_to_rows(network, &sparsity_pattern, &node_to_unknown);

    // compute the Jacobian matrix
    let values = vec![0.0; sparsity_pattern.as_ref().row_idx().len()]; // Jacobian matrix values
    let jac = SparseColMat::new(sparsity_pattern.clone(), values.clone());

    // compute the symbolic Cholesky factorization
    let symbolic_llt = SymbolicLlt::try_new(jac.symbolic(), Side::Lower).expect("Failed to compute symbolic Cholesky factorization");

    Self { network, node_to_unknown, sparsity_pattern, symbolic_llt, jac, csc_indices, node_rows, skip_timesteps: true }
  }

  /// Run the hydraulic solver
  pub fn run(self, mut parallel: bool) -> SolverResult {
    
    // calculate number of report steps to run the solver for
    let report_steps = (self.network.options.time_options.duration / self.network.options.time_options.report_timestep) + 1;

    // initialize the results struct
    let mut results = SolverResult::new(self.network.links.len(), self.network.nodes.len(), report_steps);

    // calculate the initial state
    let mut state = SolverState::new_with_initial_values(self.network);

    // check network is elegible for parallel solving if requested
    if (self.network.has_tanks() || self.network.has_pressure_controls()) && parallel {
      warn!("Networks with tanks or pressure controls cannot be solved in parallel, running sequentially");
      parallel = false;
    }

    // run the solver in parallel using Rayon if enabled
    if parallel {
      // solve the first step to use as initial values for the next, parallel computed steps
      self.apply_patterns(&mut state, 0);
      let step_result = self.solve(&mut state, 0).unwrap();

      // store the results of step 0
      results.append(&step_result, 0);

      // do parallel solves using Rayon
      let par_results: Vec<SolverState> = (1..report_steps).into_par_iter().map(|step| {
        let mut state = state.clone();
        // apply the head/demand patterns to the state
        self.apply_patterns(&mut state, step * self.network.options.time_options.report_timestep);
        self.solve(&mut state, step).unwrap()
      }).collect();
      for (step,step_result) in par_results.iter().enumerate() {
        results.append(step_result, step+1);
      }
    } else {
      // do sequential solves
      let mut time = 0;
      while time <= self.network.options.time_options.duration {
        // apply the head pattern to reservoirs with a head pattern
        self.apply_patterns(&mut state, time);
        // apply controls to the state
        self.apply_controls(&mut state, time);
        // solve the step, update the state
        state = self.solve(&mut state, time).unwrap();
        // update the heads of the tanks (= elevation + level) before running the next step
        if time % self.network.options.time_options.report_timestep == 0 {
          results.append(&state, time / self.network.options.time_options.report_timestep);
        }
        let timestep = self.next_time_step(time, &mut state);
        time += timestep;
      }
    }
    // convert the units back to the original units
    results.convert_units(&self.network.options.flow_units, &self.network.options.unit_system);

    results
  }

  // apply pressure controls to the state
  fn apply_pressure_controls(&self, state: &mut SolverState) -> bool {

    let mut changed = false;

    for control in &self.network.controls {
      if matches!(control.condition, ControlCondition::LowPressure { .. } | ControlCondition::HighPressure { .. }) {
        let active = control.is_active(state, self.network, 0, 0);
        if active {
          changed = changed || control.activate(state, self.network);
        }
      }
    }
    changed
  }

  // apply controls to the state
  fn apply_controls(&self, state: &mut SolverState, time: usize) {
    let clocktime = (time + self.network.options.time_options.start_clocktime) % (24 * 3600);

    // evaluate the controls
    for control in &self.network.controls {
      // skip pressure controls
      if matches!(control.condition, ControlCondition::LowPressure { .. } | ControlCondition::HighPressure { .. }) {
        continue;
      }
      // skip level controls at time 0 (no levels computed yet)
      if matches!(control.condition, ControlCondition::LowLevel { .. } | ControlCondition::HighLevel { .. }) && time == 0 {
        continue;
      }
      // evaluate the control
      if control.is_active(state, self.network, time, clocktime) {
        // activate the control
        control.activate(state, self.network);
      }
    }
  }


  // Calculate the next time step. The next time step is:
  fn next_time_step(&self, current_time: usize, state: &mut SolverState) -> usize {

    let times = &self.network.options.time_options;
    let clocktime = (current_time + times.start_clocktime) % (24 * 3600); // get clocktime

    // if no quality, tanks, rules and controls, simply return report time, but dont do this when skipping timesteps is set to false, 
    // as skipping timesteps will result in slight mismatches in results due to different initial values

    if self.network.controls.len() == 0 && !self.network.has_tanks() && !self.network.has_quality() && self.skip_timesteps {
      return times.report_timestep;
    }

    // time to next report
    let t_report = times.report_timestep - (current_time % times.report_timestep);
    // time to next pattern
    let t_pattern = times.pattern_timestep - (current_time % times.pattern_timestep);

    // time for the next tank to fill or drain
    let t_tanks = self.network.nodes.iter()
        .zip(state.heads.iter())
        .zip(state.demands.iter())
        .map(|((node, head), demand)| {
          let NodeType::Tank(tank) = &node.node_type else { return usize::MAX };
          let level = head - node.elevation;
          tank.time_to_fill_or_drain(level, *demand)
        })
        .min().unwrap_or(usize::MAX);
    
    // time for the next control to activate
    let t_controls = self.network.controls.iter()
        .map(|control| {
          match &control.condition {
            ControlCondition::Time { seconds } => {
              if *seconds < current_time {
                return usize::MAX;
              }
              seconds - current_time
            }
            ControlCondition::ClockTime { seconds } => {
              if *seconds < clocktime {
                // wrap around to the next day
                return ((*seconds as i32 - clocktime as i32) + (24 * 3600)) as usize
              }
              seconds - clocktime
            }
            ControlCondition::LowLevel { tank_index, target } | ControlCondition::HighLevel { tank_index, target } => {
              let node = &self.network.nodes[*tank_index];
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

    let timestep = *[t_report, t_pattern, t_tanks, t_controls, times.hydraulic_timestep].iter().filter(|&x| *x > 0).min().unwrap();
    // update the heads of the tanks
    self.update_tanks(state, timestep);

    return timestep;
  }

  /// Solve the network using the Global Gradient Algorithm (Todini & Pilati, 1987) for a single step
  /// Returns the flows and heads
  fn solve(&self, state: &mut SolverState, step: usize) -> Result<SolverState, String> {

    debug!("Solving step {}...", step);

    // perform GGA iterations
    // solve the system of equations: A * h = rhs
    // where A is the Jacobian matrix, h is the vector of heads, and rhs is the vector of right-hand side values
    // A is a sparse matrix, so we use the Compressed Sparse Column (CSC) format to store it
    // h is the vector of heads, and rhs is the vector of right-hand side values
    let unknown_nodes = self.node_to_unknown.iter().filter(|&x| x.is_some()).count();

    // pre-allocate the Jacobian matrix values and the right-hand side values
    let mut values = vec![0.0; self.sparsity_pattern.as_ref().row_idx().len()]; // Jacobian matrix values
    let mut rhs = vec![0.0; unknown_nodes]; // (unknown nodes only)

    // pre-allocate the link coefficients and the Jacobian matrix
    let mut link_coefficients = ResistanceCoefficients::new(self.network.links.len());
    let mut jac = self.jac.clone();

    // pre-allocate the excess flows vector
    let mut excess_flows = vec![0.0; self.network.nodes.len()];

    'gga: for iteration in 1..=self.network.options.max_trials {

      // reset values and rhs
      values.fill(0.0);
      rhs.fill(0.0);

      // set RHS to -demand (unknown nodes only)
      // for (global, &head_id) in self.node_to_unknown.iter().enumerate() {
      //   if let Some(i) = head_id {
      //     rhs[i] = -state.demands[global];
      //   }
      // }

      // calculate excess flows at each node (needed for PSV/PRV valves)
      if self.network.contains_pressure_control_valve {
        self.calculate_excess_flows(state, &mut excess_flows);
      }

      // assemble Jacobian and RHS contributions from links
      self.assemble_jacobian(state, &mut values, &mut rhs, &mut link_coefficients, &excess_flows);

      // solve the system of equations: J * dh = rhs
      jac.val_mut().copy_from_slice(&values);

      // Perform numerical factorization using pre-computed symbolic factorization
      let llt = match Llt::try_new_with_symbolic(self.symbolic_llt.clone(), jac.as_ref(), Side::Lower) {
        Ok(llt) => llt,
        // if the pivot is non-positive, try to fix the bad valve
        Err(LltError::Numeric(NonPositivePivot { index })) => {
          // fix the bad valve
          if self.fix_bad_valve(index-1, &mut state.statuses) {
            continue 'gga;
          }
          // if no valves found, panic
          let node_index = self.node_to_unknown.iter().position(|&x| x.is_some() && x.unwrap() == index-1).unwrap();
          let error_message = format!("Singular matrix: check connectivity at node '{}'", self.network.nodes[node_index+1].id);
          error!("{}", error_message);
          return Err(error_message);

        }
        Err(e) => {
          // if other error, panic
          panic!("{}", e);
        }
      };

      // solve the system of equations: J * dh = rhs
      let dh = llt.solve(&Mat::from_fn(unknown_nodes, 1, |r, _| rhs[r]));

      // update the heads of the nodes (unknown nodes only)
      for (global, &head_id) in self.node_to_unknown.iter().enumerate() {
        if let Some(i) = head_id {
          state.heads[global] = dh[(i, 0)];
        }
      }

      // update flows and statuses and get the iteration statistics
      let mut stats = self.update_links(state, &link_coefficients);
      // update the emitter flows and iteration statistics
      self.update_emitter_flows(state, &mut stats);

      // if PDA, update the demand flows and get the iteration statistics
      if self.network.options.demand_model == DemandModel::PDA {
        self.update_demand_flows(state, &mut stats);
      }


      // update links connected to tanks (close/open them based on tank level) and gather flow balance into/out of tanks
      self.update_tank_links(state);


      debug!(">> Iteration {}: Relative change: {:.6}, Status changed: {}", iteration, stats.relative_change(), stats.status_changed);
      debug!(">>>> Max flow change: {:.6} for link {}", stats.max_dq_converted(&self.network.options), self.network.links[stats.max_dq_index].id);

      let max_flow_change = self.network.options.max_flow_change.unwrap_or(BIG_VALUE);

      if stats.relative_change() < self.network.options.accuracy && !stats.status_changed && stats.max_dq_converted(&self.network.options) < max_flow_change {

          // apply pressure controls to the state, if any pressure controls are active, return the state
          if self.apply_pressure_controls(state) {
            continue 'gga;
          }

          let flow_balance = self.flow_balance(&state.demands, &state.flows);
          // update the demands of the junctions
          for i in 0..state.emitter_flows.len() {
            state.demands[i] += state.emitter_flows[i];
          }
          debug!("Converged in {} iterations: Error = {:.4}, Supply = {:.4}, Demand = {:.4}", iteration, flow_balance.error, flow_balance.total_supply, flow_balance.total_demand);

        return Ok(state.clone());
      }
    }
    Err(format!("Maximum number of iterations reached: {}", self.network.options.max_trials))
  }

  // Assemble the Jacobian matrix and the right-hand side vector for the system of equations
  // First assemble the contributions from the links
  // Then assemble the contributions from the emitters (virtual links)
  fn assemble_jacobian(&self, state: &mut SolverState, values: &mut Vec<f64>, rhs: &mut Vec<f64>, link_coefficients: &mut ResistanceCoefficients, excess_flows: &Vec<f64>) {

    // assemble the contributions from the links
    self.link_contributions(state, values, rhs, link_coefficients, excess_flows);
    dbg!(&rhs[0]);
    // assemble the contributions from the emitters
    self.emitter_contributions(state, values, rhs);
    dbg!(&rhs[0]);
    // if Pressure Dependent Analysis, assemble the contributions from the pressure dependent demand
    if self.network.options.demand_model == DemandModel::PDA {
      self.demand_contributions(state, values, rhs);
    }
    dbg!(&rhs[0]);
  }

  fn link_contributions(&self, state: &mut SolverState, values: &mut Vec<f64>, rhs: &mut Vec<f64>, link_coefficients: &mut ResistanceCoefficients, excess_flows: &Vec<f64>) {
    // iterate over the links
    for (i, link) in self.network.links.iter().enumerate() {
      let q = state.flows[i];
      let csc_index = &self.csc_indices[i];
      let coefficients = link.coefficients(q, state.resistances[i], state.settings[i], state.statuses[i], excess_flows[link.start_node], excess_flows[link.end_node]);

      link_coefficients.g_inv[i] = coefficients.g_inv;
      link_coefficients.y[i] = coefficients.y;

      // update the status of the link if it has a new status
      if let Some(status) = coefficients.new_status {
        state.statuses[i] = status;
      }

      // Get the CSC indices for the start and end nodes
      let u = self.node_to_unknown[link.start_node];
      let v = self.node_to_unknown[link.end_node];

      if let Some(i) = u {
          values[csc_index.diag_u.unwrap()] += coefficients.g_inv;
          rhs[i] -= q - coefficients.y;
          if self.network.nodes[link.end_node].is_fixed() {
            rhs[i] += coefficients.g_inv * state.heads[link.end_node];
          }
      }
      if let Some(j) = v {
          values[csc_index.diag_v.unwrap()] += coefficients.g_inv;
          rhs[j] += q - coefficients.y;
          if self.network.nodes[link.start_node].is_fixed() {
            rhs[j] += coefficients.g_inv * state.heads[link.start_node];
          }
      }
      if let (Some(_i), Some(_j)) = (u, v) {
          values[csc_index.off_diag.unwrap()] -= coefficients.g_inv;
      }

      // apply the upstream/downstream modifications to the nodes (for PSV/PRV valves)
      if let Some(upstream_modification) = coefficients.upstream_modification {
        rhs[u.unwrap()] += upstream_modification.rhs_add;
        values[csc_index.diag_u.unwrap()] += upstream_modification.diagonal_add;
      }

      if let Some(downstream_modification) = coefficients.downstream_modification {
        rhs[v.unwrap()] += downstream_modification.rhs_add;
        values[csc_index.diag_v.unwrap()] += downstream_modification.diagonal_add;
      }
    }
  }

  fn emitter_contributions(&self, state: &mut SolverState, values: &mut Vec<f64>, rhs: &mut Vec<f64>) {
    // iterate over emitters
    for (i, node) in self.network.nodes.iter().enumerate() {
      if let NodeType::Junction(junction) = &node.node_type {
        if junction.emitter_coefficient > 0.0 {
          // get the index for the diagonal entry in the Jacobian matrix
          let row = self.node_rows[i].unwrap();
          // get the index for the unknown node in the RHS vector
          let idx = self.node_to_unknown[i].unwrap();

          let (g_inv, y) = junction.emitter_coefficients(state.emitter_flows[i], self.network.options.emitter_exponent);
          // update RHS
          rhs[idx] += (y + node.elevation) * g_inv - state.emitter_flows[i];
          // update matrix diagonal
          values[row] += g_inv;

        }

      }
    }
  }

  fn demand_contributions(&self, state: &mut SolverState, values: &mut Vec<f64>, rhs: &mut Vec<f64>) {

    let options = &self.network.options;

    // Get demand function parameters
    let dp = (options.required_pressure - options.minimum_pressure).max(PDA_MIN_DIFF);
    let n = 1.0 / options.pressure_exponent;

    // Iterate over all junctions
    for (i, node) in self.network.nodes.iter().enumerate() {
      if let NodeType::Junction(junction) = &node.node_type {
        // only consider junctions with a positive demand
        if state.demands[i] > 0.0 {
          // get the index for the diagonal entry in the Jacobian matrix
          let row = self.node_rows[i].unwrap();
          // get the index for the unknown node in the RHS vector
          let idx = self.node_to_unknown[i].unwrap();
          // get coefficients for the demand function
          let (g_inv, y) = junction.demand_coefficients(state.demand_flows[i], state.demands[i], dp, n);

          if (1.0 / g_inv) > 0.0 {
            // update RHS
            rhs[idx] += (y + node.elevation) * g_inv - state.demand_flows[i];
            // update matrix diagonal
            values[row] += g_inv;
          }

        }
      }
    }

  }

  fn calculate_excess_flows(&self, state: &SolverState, excess_flows: &mut Vec<Cfs>) {
      // calculate excess flows at each node (needed for PSV/PRV valves)
      for (i, emitter_flow) in state.emitter_flows.iter().enumerate() {
        excess_flows[i] = -emitter_flow;
      }
      for (i, demand) in state.demands.iter().enumerate() {
        excess_flows[i] -= demand;
      }
      for (i, link) in self.network.links.iter().enumerate() {
        let q = state.flows[i];
        excess_flows[link.start_node] -= q;
        excess_flows[link.end_node] += q;
      }
  }

  fn update_links(&self, state: &mut SolverState, coefficients: &ResistanceCoefficients) -> IterationStatistics {

    let mut stats = IterationStatistics::default();

    for (i, link) in self.network.links.iter().enumerate() {
      // calculate the head difference between the start and end nodes
      let dh = state.heads[link.start_node] - state.heads[link.end_node];
      // calculate the 1/G_ij and Y_ij coefficients
      let g_inv = coefficients.g_inv[i];
      let y = coefficients.y[i];
    
      // Flow update: dq = y - g_inv * dh
      let dq = y - g_inv * dh;

      if dq.abs() > stats.max_dq {
        stats.max_dq = dq.abs();
        stats.max_dq_index = i;
      }

      // update the flow of the link
      state.flows[i] -= dq;

      dbg!(state.flows[i]);

      // update the link status and check for status changes
      let new_status = link.update_status(state.settings[i], state.statuses[i], state.flows[i], state.heads[link.start_node], state.heads[link.end_node]);
      if let Some(status) = new_status {
        // ignore temporary closed status changes (Check valve) and pump status changes
        if state.statuses[i] != LinkStatus::TempClosed && state.statuses[i] != LinkStatus::Xhead {
          stats.status_changed = true;
        }
        debug!("<yellow>Status changed for link {} from {:?} to {:?}</>", link.id, state.statuses[i], status);
        state.statuses[i] = status;
      }

      // update the sum of the absolute changes in flow and the sum of the absolute flows
      stats.sum_dq += dq.abs();
      stats.sum_q  += state.flows[i].abs();
    }
    // return the iteration statistics
    stats
  }
  fn update_emitter_flows(&self, state: &mut SolverState, stats: &mut IterationStatistics) {
    // iterate over all junctions
    for (i, node) in self.network.nodes.iter().enumerate() {
      if let NodeType::Junction(junction) = &node.node_type {
        if junction.emitter_coefficient > 0.0 {
          // get the driving head difference
          let dh = state.heads[i] - node.elevation;
          let (g_inv, y) = junction.emitter_coefficients(state.emitter_flows[i], self.network.options.emitter_exponent);
          // update the emitter flow
          let dq = (y - dh) * g_inv;
          // update the emitter flow
          state.emitter_flows[i] -= dq;
          // update the iteration statistics
          stats.sum_dq += dq.abs();
          stats.sum_q += state.emitter_flows[i].abs();
          if dq.abs() > stats.max_dq {
            stats.max_dq = dq.abs();
          }
        }
      }
    }

  }

  fn update_demand_flows(&self, state: &mut SolverState, stats: &mut IterationStatistics) {

    let options = &self.network.options;

    let dp = (options.required_pressure - options.minimum_pressure).max(PDA_MIN_DIFF);
    let n = 1.0 / options.pressure_exponent;

    // iterate over all junctions
    for (i, node) in self.network.nodes.iter().enumerate() {
      if let NodeType::Junction(junction) = &node.node_type {
        if state.demands[i] > 0.0 {

          let (g_inv, y) = junction.demand_coefficients(state.demand_flows[i], state.demands[i], dp, n);

          let dh = state.heads[i] - node.elevation - options.minimum_pressure;
          let dq = (y - dh) * g_inv;

          // update the demand flow
          state.demand_flows[i] -= dq;
          dbg!(state.demand_flows[i]);
          // update the iteration statistics
          stats.sum_dq += dq.abs();
          stats.sum_q += state.demand_flows[i].abs();
          if dq.abs() > stats.max_dq {
            stats.max_dq = dq.abs();
          }
        }
      }
    }


  }


  fn apply_patterns(&self, state: &mut SolverState, time: usize) {

    let time_options = &self.network.options.time_options;
    // get the pattern time
    let pattern_time = time_options.pattern_start + time;
    // get pattern index
    let pattern_index = pattern_time / time_options.pattern_timestep;

    // apply the head pattern to reservoirs with a head pattern
    for (i, node) in self.network.nodes.iter().enumerate() {
      let Some(head_pattern) = node.head_pattern() else { continue };
      let pattern = &self.network.patterns[head_pattern];
      state.heads[i] = node.elevation * pattern.multipliers[pattern_index % pattern.multipliers.len()];
    }

    // apply the demand pattern to junctions
    state.demands = self.network.nodes.iter().map(|n| {
        let NodeType::Junction(junction) = &n.node_type else { return 0.0 };

        if let Some(pattern_id) = &junction.pattern {
          let pattern = &self.network.patterns[pattern_id];

          // get the multiplier for the pattern index (wrap around if needed)
          let multiplier = pattern.multipliers[pattern_index % pattern.multipliers.len()];
          return junction.basedemand * multiplier * self.network.options.demand_multiplier;
        }
        // if no pattern, return the basedemand times the global demand multiplier
        return junction.basedemand * self.network.options.demand_multiplier;

      }).collect::<Vec<Cfs>>();

    // if PDA, set the demand flows to the demands
    if self.network.options.demand_model == DemandModel::PDA {
      state.demand_flows = state.demands.clone();
    }

  }

  /// Update the links connected to tanks and gather flow balance into/out of tanks
  fn update_tank_links(&self, state: &mut SolverState) {
    // iterate over the tanks
    for (tank_index, node) in self.network.nodes.iter().enumerate() {

      if let NodeType::Tank(tank) = &node.node_type {
        // reset the demand of the tank
        state.demands[tank_index] = 0.0;
        // check if the tank is closed for filling
        let fill_closed = state.heads[tank_index] >= tank.elevation + tank.max_level && !tank.overflow;
        let empty_closed = state.heads[tank_index] <= tank.elevation + tank.min_level;

        // iterate over the links connected to the tank
        for link_index in &tank.links_to {
          // add the flow to the demand of the tank (positive flow is into the tank)
          state.demands[tank_index] += state.flows[*link_index];
          if fill_closed && state.flows[*link_index] > 0.0 { state.statuses[*link_index] = LinkStatus::TempClosed; }
          if empty_closed && state.flows[*link_index] < 0.0 { state.statuses[*link_index] = LinkStatus::TempClosed; }
        }
        for link_index in &tank.links_from {
          // subtract the flow from the demand of the tank (negative flow is out of the tank)
          state.demands[tank_index] -= state.flows[*link_index];
          if fill_closed && state.flows[*link_index] < 0.0 { state.statuses[*link_index] = LinkStatus::TempClosed; }
          if empty_closed && state.flows[*link_index] > 0.0 { state.statuses[*link_index] = LinkStatus::TempClosed; }
        }
      }
    }
    

  }
  
  /// Update the heads of tanks, and the statuses of the links connected to tanks
  fn update_tanks(&self, state: &mut SolverState, timestep: usize) {

    for (tank_index, node) in self.network.nodes.iter().enumerate() {
      if let NodeType::Tank(tank) = &node.node_type {
        // calculate flow balance of the tank
        let delta_volume = state.demands[tank_index] * timestep as Ft3; // in ft^3
        // update the head of the tank
        let new_head = tank.new_head(delta_volume, state.heads[tank_index]);
        // update the head of the tank
        state.heads[tank_index] = new_head;
      }
    }

  }

  /// Fix a bad valve by setting its status to Closed
  fn fix_bad_valve(&self, unknown_node_index: usize, statuses: &mut Vec<LinkStatus>) -> bool {
    // get the node index
    let node_index = self.node_to_unknown.iter().position(|&x| x.is_some() && x.unwrap() == unknown_node_index).unwrap();
    // find links connected to the node
    for (i, link) in self.network.links.iter().enumerate() {
      // skip if the link is not connected to the node
      if link.start_node != node_index && link.end_node != node_index {
        continue;
      }
      // check if the link is a PSV or PRV valve
      if let LinkType::Valve(valve) = &link.link_type {
        if valve.valve_type == ValveType::PSV || valve.valve_type == ValveType::PRV {
          if statuses[i] == LinkStatus::Active {
            // set the status to XPressure to indicate a bad valve
            statuses[i] = LinkStatus::XPressure;
            return true;
          }
        }
      }
    }
    // if no valves found, return false
    false
  }

  /// Calculate the flow balance error
  fn flow_balance(&self, demands: &Vec<Cfs>, flows: &Vec<Cfs>) -> FlowBalance {

    let sum_demand: Cfs = demands.iter().sum();

    let mut sum_supply: Cfs = 0.0;
    for (i, link) in self.network.links.iter().enumerate() {
      if !matches!(self.network.nodes[link.end_node].node_type, NodeType::Junction { .. }) {
        sum_supply -= flows[i];
      }
      if !matches!(self.network.nodes[link.start_node].node_type, NodeType::Junction { .. }) {
        sum_supply += flows[i];
      }
    }
    let error = sum_demand - sum_supply;
    return FlowBalance { total_demand: sum_demand, total_supply: sum_supply, error: error };
  }
  fn map_nodes_to_rows(network: &Network, sparsity_pattern: &SymbolicSparseColMat<usize>, node_to_unknown: &Vec<Option<usize>>) -> Vec<Option<usize>> {

    let mut node_rows = Vec::with_capacity(network.nodes.len());

    for (i, _) in network.nodes.iter().enumerate() {
      if let Some(idx) = node_to_unknown[i] {
        let row = find_csc_index(sparsity_pattern.as_ref(), idx, idx).unwrap();
        node_rows.push(Some(row));
      } else {
        node_rows.push(None);
      }
    }
    node_rows
  }

  /// Map each link to its CSC (Compressed Sparse Column) indices
  fn map_links_to_csc_indices(network: &Network, sparsity_pattern: &SymbolicSparseColMat<usize>, node_to_unknown: &Vec<Option<usize>>) -> Vec<CSCIndex> {

    let mut csc_indices = Vec::with_capacity(network.links.len());
    for link in network.links.iter() {
      let mut csc_index = CSCIndex::default();
      let u = node_to_unknown[link.start_node];
      let v = node_to_unknown[link.end_node];

      let sym = sparsity_pattern.as_ref();
      if let Some(i) = u { csc_index.diag_u = find_csc_index(sym, i, i); }
      if let Some(j) = v { csc_index.diag_v = find_csc_index(sym, j, j); }
      if let (Some(i), Some(j)) = (u, v) {
          let (row, col) = if i >= j { (i, j) } else { (j, i) };
          csc_index.off_diag = find_csc_index(sym, row, col);
      }
      csc_indices.push(csc_index);
    }
    csc_indices
  }

  /// Build sparsity pattern
  fn build_sparsity_pattern(network: &Network, node_to_unknown: &Vec<Option<usize>>) -> SymbolicSparseColMat<usize> {
    let n_unknowns = node_to_unknown.iter().filter(|&x| x.is_some()).count();

    // Pre-allocate: at most 3 triplets per link (2 diagonal + 1 lower off-diagonal)
    let mut triplets = Vec::with_capacity(3 * network.links.len());
    for link in network.links.iter() {
      let u = node_to_unknown[link.start_node]; // u is the index of the start node (None if unknown)
      let v = node_to_unknown[link.end_node]; // v is the index of the end node (None if unknown)
      // diagonal elements (self-connectivity)
      if let Some(i) = u { triplets.push(Triplet::new(i, i, 0.0)); } // add diagonal element for u
      if let Some(j) = v { triplets.push(Triplet::new(j, j, 0.0)); } // add diagonal element for v
      // lower-triangular off diagonal element (connectivity)
      if let (Some(i), Some(j)) = (u, v) {
        let (row, col) = if i >= j { (i, j) } else { (j, i) };
        triplets.push(Triplet::new(row, col, 0.0));
      }
    }
    // convert triplets to sparse matrix
    let sparsity_matrix = SparseColMat::try_new_from_triplets(n_unknowns, n_unknowns, &triplets).unwrap();
    let symbolic = sparsity_matrix.symbolic().to_owned().unwrap();

    symbolic
  }

  /// Build global unknown-numbering map
  fn build_unknown_numbering_map(network: &Network) -> Vec<Option<usize>> {
    let mut unknown_id = 0;
    let node_to_unknown: Vec<Option<usize>> = network.nodes
        .iter()
        .map(|n| if matches!(n.node_type, NodeType::Junction { .. }) { let id = unknown_id; unknown_id += 1; Some(id) } else { None })
        .collect();
    node_to_unknown
  }

}
