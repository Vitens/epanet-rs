use faer::sparse::{SparseColMat, Triplet, SymbolicSparseColMat};
use faer::sparse::linalg::LltError;
use faer::linalg::cholesky::llt::factor::LltError::NonPositivePivot;
use faer::sparse::linalg::solvers::{SymbolicLlt, Llt};
use faer::{Mat, Side};
use faer::prelude::*;
use serde::Serialize;
use rayon::prelude::*;


use simplelog::{warn, debug, error};

use crate::constants::{BIG_VALUE, Q_ZERO};
use crate::model::node::NodeType;
use crate::model::link::{LinkType, LinkTrait, LinkStatus};
use crate::model::network::Network;
use crate::model::valve::ValveType;
use crate::model::units::{FlowUnits, UnitSystem};

#[derive(Serialize)]
pub struct SolverResult {
  pub flows: Vec<Vec<f64>>,
  pub heads: Vec<Vec<f64>>,
  pub demands: Vec<Vec<f64>>,
}

impl SolverResult {
  pub fn new(n_links: usize, n_nodes: usize, n_steps: usize) -> Self {
    Self { flows: vec![vec![0.0; n_links]; n_steps], heads: vec![vec![0.0; n_nodes]; n_steps], demands: vec![vec![0.0; n_nodes]; n_steps] }
  }

  fn append(&mut self, state: &SolverState, step: usize) {
    self.flows[step] = state.flows.clone();
    self.heads[step] = state.heads.clone();
    self.demands[step] = state.demands.clone();
  }

  // convert the solver units back to the original units
  pub fn convert_units(&mut self, flow_units: &FlowUnits, unit_system: &UnitSystem) {

    let flow_scale = flow_units.per_cfs();
    let head_scale = unit_system.per_feet();

    self.flows.iter_mut().flatten().for_each(|flow| *flow *= flow_scale);
    self.heads.iter_mut().flatten().for_each(|head| *head *= head_scale);
    self.demands.iter_mut().flatten().for_each(|demand| *demand *= flow_scale);
  }
}

pub struct FlowBalance {
  pub total_demand: f64,
  pub total_supply: f64,
  pub error: f64,
}

pub struct ResistanceCoefficients {
  pub g_inv: Vec<f64>,
  pub y: Vec<f64>,
}

impl ResistanceCoefficients {
  pub fn new(size: usize) -> Self {
    Self { g_inv: vec![0.0; size], y: vec![0.0; size] }
  }
}

#[derive(Default)]
pub struct IterationStatistics {
  pub sum_dq: f64,
  pub sum_q: f64,
  pub max_dq: f64,
  pub max_dq_index: usize,
  pub status_changed: bool,
  pub relative_change: f64,
}

/// The solver state is the initial/final state of the solver for a single step
#[derive(Debug, Clone)]
pub struct SolverState {
  pub flows: Vec<f64>,
  pub heads: Vec<f64>,
  pub demands: Vec<f64>,
  pub statuses: Vec<LinkStatus>,
  pub resistances: Vec<f64>,
}

impl SolverState {
  /// Create a new solver state with the initial values for the flows, heads, demands and statuses and calculate resistances
  pub fn new_with_initial_values(network: &Network) -> Self {
    Self { flows: network.links.iter().map(|l| l.initial_flow()).collect::<Vec<f64>>(), 
           heads: network.nodes.iter().map(|n| n.initial_head()).collect::<Vec<f64>>(), 
           demands: vec![0.0; network.nodes.len()], 
           statuses: network.links.iter().map(|l| l.initial_status).collect::<Vec<LinkStatus>>(),
           resistances: network.links.iter().map(|l| l.resistance()).collect::<Vec<f64>>(),
         }
  }
}

/// CSC (Compressed Sparse Column) indices for the Jacobian matrix used in the Global Gradient Algorithm
#[derive(Default)]
pub struct CSCIndex {
  pub diag_u: Option<usize>,      // CSC index for J[u,u]
  pub diag_v: Option<usize>,      // CSC index for J[v,v]
  pub off_diag: Option<usize>, // CSC index for lower triangular off-diagonal entry
}

pub struct HydraulicSolver<'a> {
  network: &'a Network, 
  node_to_unknown: Vec<Option<usize>>,
  sparsity_pattern: SymbolicSparseColMat<usize>,
  symbolic_llt: SymbolicLlt<usize>,
  jac: SparseColMat<usize, f64>,
  csc_indices: Vec<CSCIndex>,
}

impl<'a> HydraulicSolver<'a> {

  pub fn new(network: &'a Network) -> Self {
    // build global unknown-numbering map
    let node_to_unknown = Self::build_unknown_numbering_map(network);
    // generate sparsity pattern
    let sparsity_pattern = Self::build_sparsity_pattern(network, &node_to_unknown);
    // map each link to its CSC (Compressed Sparse Column) indices
    let csc_indices = Self::map_links_to_csc_indices(network, &sparsity_pattern, &node_to_unknown);

    // compute the Jacobian matrix
    let values = vec![0.0; sparsity_pattern.as_ref().row_idx().len()]; // Jacobian matrix values
    let jac = SparseColMat::new(sparsity_pattern.clone(), values.clone());

    // compute the symbolic Cholesky factorization
    let symbolic_llt = SymbolicLlt::try_new(jac.symbolic(), Side::Lower).expect("Failed to compute symbolic Cholesky factorization");

    Self { network, node_to_unknown, sparsity_pattern, symbolic_llt, jac, csc_indices }
  }

  /// Run the hydraulic solver
  pub fn run(self, mut parallel: bool) -> SolverResult {
    
    // calculate number of steps to run the solver for
    let steps = (self.network.options.time_options.duration / self.network.options.time_options.report_timestep) + 1;

    // initialize the results struct
    let mut results = SolverResult::new(self.network.links.len(), self.network.nodes.len(), steps);

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
      let par_results: Vec<SolverState> = (1..steps).into_par_iter().map(|step| {
        let mut state = state.clone();
        // apply the head/demand patterns to the state
        self.apply_patterns(&mut state, step);
        self.solve(&mut state, step).unwrap()
      }).collect();
      for (step,step_result) in par_results.iter().enumerate() {
        results.append(step_result, step);
      }
    } else {
      // do sequential solves
      for step in 0..steps {
        // apply the head pattern to reservoirs with a head pattern
        self.apply_patterns(&mut state, step);
        // solve the step, update the state
        state = self.solve(&mut state, step).unwrap();
        // update the heads of the tanks (= elevation + level) before running the next step
        results.append(&state, step);
        if steps > 1 {
          self.update_tanks(&state.flows, &mut state.heads);
        }
      }
    }
    // convert the units back to the original units
    results.convert_units(&self.network.options.flow_units, &self.network.options.unit_system);

    results
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
      for (global, &head_id) in self.node_to_unknown.iter().enumerate() {
        if let Some(i) = head_id {
          rhs[i] = -state.demands[global];
        }
      }

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
      let stats = self.update_links(state, &link_coefficients);

      // update links connected to tanks (close/open them based on tank level)
      self.update_tank_links(&state.flows, &state.heads, &mut state.statuses);


      debug!(">> Iteration {}: Relative change: {:.6}, Status changed: {}", iteration, stats.relative_change, stats.status_changed);
      debug!(">>>> Max flow change: {:.6} for link {}", stats.max_dq, self.network.links[stats.max_dq_index].id);

      let max_flow_change = self.network.options.max_flow_change.unwrap_or(BIG_VALUE);

      if stats.relative_change < self.network.options.accuracy && !stats.status_changed && stats.max_dq < max_flow_change {

          let flow_balance = self.flow_balance(&state.demands, &state.flows);
          debug!("Converged in {} iterations: Error = {:.4}, Supply = {:.4}, Demand = {:.4}", iteration, flow_balance.error, flow_balance.total_supply, flow_balance.total_demand);

        return Ok(state.clone());
      }
    }
    Err(format!("Maximum number of iterations reached: {}", self.network.options.max_trials))
  }
  fn assemble_jacobian(&self, state: &mut SolverState, values: &mut Vec<f64>, rhs: &mut Vec<f64>, link_coefficients: &mut ResistanceCoefficients, excess_flows: &Vec<f64>) {
    for (i, link) in self.network.links.iter().enumerate() {
      let q = state.flows[i];
      let csc_index = &self.csc_indices[i];
      let coefficients = link.coefficients(q, state.resistances[i], state.statuses[i], excess_flows[link.start_node], excess_flows[link.end_node]);

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

  fn calculate_excess_flows(&self, state: &SolverState, excess_flows: &mut Vec<f64>) {
      // calculate excess flows at each node (needed for PSV/PRV valves)
      for (i, demand) in state.demands.iter().enumerate() {
        excess_flows[i] = -demand;
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

      // update the link status and check for status changes
      let new_status = link.update_status(state.statuses[i], state.flows[i], state.heads[link.start_node], state.heads[link.end_node]);
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
    // update the relative change
    stats.relative_change = stats.sum_dq / (stats.sum_q + Q_ZERO);
    // convert max_dq and max_dh to correct units
    stats.max_dq *= self.network.options.flow_units.per_cfs();
    // return the iteration statistics
    stats
  }

  fn apply_patterns(&self, state: &mut SolverState, step: usize) {

    let time_options = &self.network.options.time_options;
    // get the time
    let time = step * time_options.report_timestep;
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

      }).collect::<Vec<f64>>()

  }

  fn update_tank_links(&self, flows: &Vec<f64>, heads: &Vec<f64>, statuses: &mut Vec<LinkStatus>) {
    // iterate over the tanks
    for (tank_index, node) in self.network.nodes.iter().enumerate() {
      if let NodeType::Tank(tank) = &node.node_type {
        // check if the tank is closed for filling
        let fill_closed = heads[tank_index] >= tank.elevation + tank.max_level && !tank.overflow;
        let empty_closed = heads[tank_index] <= tank.elevation + tank.min_level;

        // iterate over the links connected to the tank
        for link_index in &tank.links_to {
          if fill_closed && flows[*link_index] > 0.0 { statuses[*link_index] = LinkStatus::TempClosed; }
          if empty_closed && flows[*link_index] < 0.0 { statuses[*link_index] = LinkStatus::TempClosed; }
        }
        for link_index in &tank.links_from {
          if fill_closed && flows[*link_index] < 0.0 { statuses[*link_index] = LinkStatus::TempClosed; }
          if empty_closed && flows[*link_index] > 0.0 { statuses[*link_index] = LinkStatus::TempClosed; }
        }
      }
    }
    

  }
  
  /// Update the heads of tanks, and the statuses of the links connected to tanks
  fn update_tanks(&self, flows: &Vec<f64>, heads: &mut Vec<f64>) {

    for (tank_index, node) in self.network.nodes.iter().enumerate() {
      if let NodeType::Tank(tank) = &node.node_type {
        // calculate flow balance of the tank
        let mut flow_balance = 0.0;

        for link_index in &tank.links_to {
          flow_balance += flows[*link_index];
        }
        for link_index in &tank.links_from {
          flow_balance -= flows[*link_index];
        }

        let delta_volume = flow_balance * self.network.options.time_options.hydraulic_timestep as f64; // in ft^3
        // update the head of the tank
        let new_head = tank.new_head(delta_volume, heads[tank_index]);
        // update the head of the tank
        heads[tank_index] = new_head;
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
  fn flow_balance(&self, demands: &Vec<f64>, flows: &Vec<f64>) -> FlowBalance {

    let sum_demand: f64 = demands.iter().sum();

    let mut sum_supply: f64 = 0.0;
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

/// Helper function to find the CSC index for a given row and column
#[inline]
fn find_csc_index(
    sym: faer::sparse::SymbolicSparseColMatRef<usize>,
    row: usize,
    col: usize,
    ) -> Option<usize> {
    let col_start = sym.col_ptr()[col];
    let col_end = sym.col_ptr()[col + 1];
    sym.row_idx()[col_start..col_end]
        .iter()
        .position(|&r| r == row)
        .map(|pos| col_start + pos)}