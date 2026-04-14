use faer::sparse::{SparseColMat, Triplet, SymbolicSparseColMat};
use faer::sparse::linalg::LltError;
use faer::linalg::cholesky::llt::factor::LltError::NonPositivePivot;
use faer::sparse::linalg::solvers::{SymbolicLlt, Llt};
use faer::sparse::linalg::amd;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::{Mat, Side};
use faer::prelude::*;

use crate::solver::state::SolverState;
use crate::solver::matrix::{CSCIndex, ResistanceCoefficients, find_csc_index};

use crate::model::units::Cfs;


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
  pub fn relative_change(&self, options: &SimulationOptions) -> f64 {
    if self.sum_q > options.accuracy {
      self.sum_dq / (self.sum_q + Q_ZERO)
    } else {
      self.sum_dq
    }
  }
  pub fn max_dq_converted(&self, options: &SimulationOptions) -> Cfs {
    self.max_dq * options.flow_units.per_cfs()
  }
}

pub struct HydraulicSolver<'a> {
  /// network
  pub network: &'a Network, 
  /// global unknown-numbering map: node_to_unknown[node_index] = unknown_index
  pub node_to_unknown: Vec<Option<usize>>,
  /// symbolic sparsity pattern
  pub sparsity_pattern: SymbolicSparseColMat<usize>,
  /// symbolic Cholesky factorization
  pub symbolic_llt: SymbolicLlt<usize>,
  /// AMD fill-reducing permutation: perm_fwd[permuted] = original
  pub perm_fwd: Vec<usize>,
  /// precomputed Jacobian matrix
  pub jac: SparseColMat<usize, f64>,
  /// precomputed CSC indices for the links
  pub csc_indices: Vec<CSCIndex>,
  /// precomputed indices for the rows of the Jacobian matrix for each node
  pub node_rows: Vec<Option<usize>>,
}

impl<'a> HydraulicSolver<'a> {

  pub fn new(network: &'a Network) -> Self {
    let node_to_unknown = Self::build_unknown_numbering_map(network);
    let sparsity_pattern = Self::build_sparsity_pattern(network, &node_to_unknown);
    let csc_indices = Self::map_links_to_csc_indices(network, &sparsity_pattern, &node_to_unknown);
    let node_rows = Self::map_nodes_to_rows(network, &sparsity_pattern, &node_to_unknown);

    let values = vec![0.0; sparsity_pattern.as_ref().row_idx().len()];
    let jac = SparseColMat::new(sparsity_pattern.clone(), values.clone());

    let perm_fwd = Self::compute_amd_permutation(&sparsity_pattern, &node_to_unknown);
    let symbolic_llt = SymbolicLlt::try_new(jac.symbolic(), Side::Lower).expect("Failed to compute symbolic Cholesky factorization");

    Self { network, node_to_unknown, sparsity_pattern, symbolic_llt, perm_fwd, jac, csc_indices, node_rows }
  }

  fn compute_amd_permutation(sparsity_pattern: &SymbolicSparseColMat<usize>, node_to_unknown: &Vec<Option<usize>>) -> Vec<usize> {
    let n_unknowns = node_to_unknown.iter().filter(|x| x.is_some()).count();
    let a_nnz = sparsity_pattern.as_ref().compute_nnz();
    let mut perm_fwd = vec![0usize; n_unknowns];
    let mut perm_inv = vec![0usize; n_unknowns];
    {
      let scratch_size = amd::order_maybe_unsorted_scratch::<usize>(n_unknowns, a_nnz);
      let mut mem = MemBuffer::new(scratch_size);
      amd::order_maybe_unsorted(
        &mut perm_fwd,
        &mut perm_inv,
        sparsity_pattern.as_ref(),
        amd::Control::default(),
        MemStack::new(&mut mem),
      ).expect("Failed to compute AMD ordering");
    }

    perm_fwd
  }

  /// Solve the network for a single timestep using the Global Gradient Algorithm (Todini & Pilati, 1987).
  /// Takes an immutable solver state and returns a new state after convergence.
  pub fn solve(&self, state: &SolverState) -> Result<SolverState, String> {

    let mut state = state.clone();

    let unknown_nodes = self.node_to_unknown.iter().filter(|&x| x.is_some()).count();

    let mut values = vec![0.0; self.sparsity_pattern.as_ref().row_idx().len()];
    let mut rhs = vec![0.0; unknown_nodes];

    let mut link_coefficients = ResistanceCoefficients::new(self.network.links.len());
    let mut jac = self.jac.clone();

    let mut excess_flows = vec![0.0; self.network.nodes.len()];

    let mut grounded_nodes = vec![false; self.network.nodes.len()];

    'gga: for iteration in 1..=self.network.options.max_trials {

      values.fill(0.0);
      rhs.fill(0.0);

      if self.network.contains_pressure_control_valve {
        self.calculate_excess_flows(&state, &mut excess_flows);
      }

      self.assemble_jacobian(&mut state, &mut values, &mut rhs, &mut link_coefficients, &excess_flows, &grounded_nodes);

      jac.val_mut().copy_from_slice(&values);

      let llt = match Llt::try_new_with_symbolic(self.symbolic_llt.clone(), jac.as_ref(), Side::Lower) {
        Ok(llt) => llt,
        Err(LltError::Numeric(NonPositivePivot { index })) => {
          let original_unknown = self.perm_fwd[index-1];
          if self.fix_bad_valve(original_unknown, &mut state.statuses) {
            continue 'gga;
          }
          let node_index = self.node_to_unknown.iter().position(|&x| x.is_some() && x.unwrap() == original_unknown).unwrap();
          if !grounded_nodes[node_index] {
            warn!("Grounding node '{}' with virtual reservoir (elevation 0) to fix singular matrix", self.network.nodes[node_index].id);
            grounded_nodes[node_index] = true;
            continue 'gga;
          }
          let error_message = format!("Singular matrix: check connectivity at node '{}'", self.network.nodes[node_index].id);
          error!("{}", error_message);
          return Err(error_message);
        }
        Err(e) => {
          return Err(e.to_string());
        }
      };

      let dh = llt.solve(&Mat::from_fn(unknown_nodes, 1, |r, _| rhs[r]));

      for (global, &head_id) in self.node_to_unknown.iter().enumerate() {
        if let Some(i) = head_id {
          state.heads[global] = dh[(i, 0)];
        }
      }

      let mut stats = self.update_links(&mut state, &link_coefficients);
      self.update_emitter_flows(&mut state, &mut stats);

      if self.network.options.demand_model == DemandModel::PDA {
        self.update_demand_flows(&mut state, &mut stats);
      }

      self.update_tank_links(&mut state);

      debug!(">> Iteration {}: Relative change: {:.6}, Status changed: {}", iteration, stats.relative_change(&self.network.options), stats.status_changed);
      debug!(">>>> Max flow change: {:.6} for link {}", stats.max_dq_converted(&self.network.options), self.network.links[stats.max_dq_index].id);

      let max_flow_change = self.network.options.max_flow_change.unwrap_or(BIG_VALUE);

      if stats.relative_change(&self.network.options) < self.network.options.accuracy && !stats.status_changed && stats.max_dq_converted(&self.network.options) < max_flow_change {

          if self.apply_pressure_controls(&mut state) {
            continue 'gga;
          }

          let flow_balance = self.flow_balance(&state.demands, &state.flows);
          if self.network.options.demand_model == DemandModel::PDA {
            state.demands = state.demand_flows.clone();
          }
          for i in 0..state.emitter_flows.len() {
            state.demands[i] += state.emitter_flows[i];
          }
          debug!("Converged in {} iterations: Error = {:.4}, Supply = {:.4}, Demand = {:.4}", iteration, flow_balance.error, flow_balance.total_supply, flow_balance.total_demand);

        return Ok(state);
      }
    }
    Err(format!("Maximum number of iterations reached: {}", self.network.options.max_trials))
  }


  fn calculate_excess_flows(&self, state: &SolverState, excess_flows: &mut Vec<Cfs>) {
      for (i, emitter_flow) in state.emitter_flows.iter().enumerate() {
        excess_flows[i] = -emitter_flow;
      }
      for (i, demand) in state.demands.iter().enumerate() {
        if self.network.options.demand_model == DemandModel::PDA {
          excess_flows[i] -= state.demand_flows[i];
        } else {
          excess_flows[i] -= demand;
        }
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
      let dh = state.heads[link.start_node] - state.heads[link.end_node];
      let g_inv = coefficients.g_inv[i];
      let y = coefficients.y[i];
    
      let dq = y - g_inv * dh;

      if dq.abs() > stats.max_dq {
        stats.max_dq = dq.abs();
        stats.max_dq_index = i;
      }

      state.flows[i] -= dq;

      let new_status = link.update_status(state.settings[i], state.statuses[i], state.flows[i], state.heads[link.start_node], state.heads[link.end_node]);
      if let Some(status) = new_status {
        if state.statuses[i] != LinkStatus::TempClosed && state.statuses[i] != LinkStatus::Xhead {
          stats.status_changed = true;
        }
        debug!("<yellow>Status changed for link {} from {:?} to {:?}</>", link.id, state.statuses[i], status);
        state.statuses[i] = status;
      }

      stats.sum_dq += dq.abs();
      stats.sum_q  += state.flows[i].abs();
    }
    stats
  }
  fn update_emitter_flows(&self, state: &mut SolverState, stats: &mut IterationStatistics) {
    for (i, node) in self.network.nodes.iter().enumerate() {
      if let NodeType::Junction(junction) = &node.node_type {
        if junction.emitter_coefficient > 0.0 {
          let dh = state.heads[i] - node.elevation;
          let (g_inv, y) = junction.emitter_coefficients(state.emitter_flows[i], self.network.options.emitter_exponent);
          let dq = (y - dh) * g_inv;
          state.emitter_flows[i] -= dq;
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

    for (i, node) in self.network.nodes.iter().enumerate() {
      if let NodeType::Junction(junction) = &node.node_type {
        if state.demands[i] > 0.0 {

          let (g_inv, y) = junction.demand_coefficients(state.demand_flows[i], state.demands[i], dp, n);

          let dh = state.heads[i] - node.elevation - options.minimum_pressure;
          let dq = (y - dh) * g_inv;

          state.demand_flows[i] -= dq;
          stats.sum_dq += dq.abs();
          stats.sum_q += state.demand_flows[i].abs();
          if dq.abs() > stats.max_dq {
            stats.max_dq = dq.abs();
          }
        }
      }
    }


  }

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

  /// Update the links connected to tanks and gather flow balance into/out of tanks
  fn update_tank_links(&self, state: &mut SolverState) {
    for (tank_index, node) in self.network.nodes.iter().enumerate() {

      if let NodeType::Tank(tank) = &node.node_type {
        state.demands[tank_index] = 0.0;
        let fill_closed = state.heads[tank_index] >= tank.elevation + tank.max_level && !tank.overflow;
        let empty_closed = state.heads[tank_index] <= tank.elevation + tank.min_level;

        for link_index in &tank.links_to {
          state.demands[tank_index] += state.flows[*link_index];
          if fill_closed && state.flows[*link_index] > 0.0 { state.statuses[*link_index] = LinkStatus::TempClosed; }
          if empty_closed && state.flows[*link_index] < 0.0 { state.statuses[*link_index] = LinkStatus::TempClosed; }
        }
        for link_index in &tank.links_from {
          state.demands[tank_index] -= state.flows[*link_index];
          if fill_closed && state.flows[*link_index] < 0.0 { state.statuses[*link_index] = LinkStatus::TempClosed; }
          if empty_closed && state.flows[*link_index] > 0.0 { state.statuses[*link_index] = LinkStatus::TempClosed; }
        }
      }
    }
    

  }

  /// Fix a bad valve by setting its status to Closed
  fn fix_bad_valve(&self, unknown_node_index: usize, statuses: &mut Vec<LinkStatus>) -> bool {
    let node_index = self.node_to_unknown.iter().position(|&x| x.is_some() && x.unwrap() == unknown_node_index).unwrap();
    for (i, link) in self.network.links.iter().enumerate() {
      if link.start_node != node_index && link.end_node != node_index {
        continue;
      }
      if let LinkType::Valve(valve) = &link.link_type {
        if valve.valve_type == ValveType::PSV || valve.valve_type == ValveType::PRV {
          if statuses[i] == LinkStatus::Active {
            debug!("Fixing bad valve for node index: {}", node_index);
            statuses[i] = LinkStatus::XPressure;
            return true;
          }
        }
      }
    }
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

    let mut triplets = Vec::with_capacity(3 * network.links.len());

    for link in network.links.iter() {
      let u = node_to_unknown[link.start_node];
      let v = node_to_unknown[link.end_node];
      if let Some(i) = u { triplets.push(Triplet::new(i, i, 0.0)); }
      if let Some(j) = v { triplets.push(Triplet::new(j, j, 0.0)); }
      if let (Some(i), Some(j)) = (u, v) {
        let (row, col) = if i >= j { (i, j) } else { (j, i) };
        triplets.push(Triplet::new(row, col, 0.0));
      }
    }
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
