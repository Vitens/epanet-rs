use faer::sparse::{SparseColMat, Triplet};
use faer::sparse::SymbolicSparseColMat;
use faer::sparse::linalg::solvers::{SymbolicLlt, Llt};
use faer::{Mat, Side};
use faer::prelude::*;
use serde::Serialize;
use rayon::prelude::*;

use crate::model::node::NodeType;
use crate::model::link::{LinkType, LinkTrait, LinkStatus};
use crate::model::network::Network;
use crate::model::units::{FlowUnits, PressureUnits};

#[derive(Serialize)]
pub struct SolverResult {
  pub flows: Vec<Vec<f64>>,
  pub heads: Vec<Vec<f64>>
}
impl SolverResult {
  // convert the solver units back to the original units
  pub fn convert_units(&mut self, flow_units: &FlowUnits, pressure_units: &PressureUnits) {
    for flow in self.flows.iter_mut() {
      for flow in flow.iter_mut() {
        *flow = *flow * flow_units.per_cfs();
      }
    }
    for head in self.heads.iter_mut() {
      for head in head.iter_mut() {
        *head = *head * pressure_units.per_feet();
      }
    }
  }
}

pub struct FlowBalance {
  pub total_demand: f64,
  pub total_supply: f64,
  pub error: f64,
}

/// CSC (Compressed Sparse Column) indices for the Jacobian matrix used in the Global Gradient Algorithm
#[derive(Default)]
pub struct CSCIndex {
  pub diag_u: Option<usize>,      // CSC index for J[u,u]
  pub diag_v: Option<usize>,      // CSC index for J[v,v]
  pub off_diag_uv: Option<usize>, // CSC index for J[u,v]
  pub off_diag_vu: Option<usize>, // CSC index for J[v,u]
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
  pub fn run(self, parallel: bool, verbose: bool) -> SolverResult {
    
    // calculate setps
    let steps = (self.network.options.time_options.duration / self.network.options.time_options.hydraulic_timestep) + 1;

    // initialize the flows and heads result vectors
    let mut flows: Vec<Vec<f64>> = vec![vec![0.0; self.network.links.len()]; steps];
    let mut heads: Vec<Vec<f64>> = vec![vec![0.0; self.network.nodes.len()]; steps];

    // set the initial flows
    let mut initial_flows: Vec<f64> = self.get_initial_flows();

    // set the initial heads
    let mut initial_heads: Vec<f64> = self.get_initial_heads();

    // run the solver in parallel using Rayon if enabled
    if parallel {

      // solve the first step to use as initial values for the next, parallel computed steps
      let (flow, head) = self.solve(&initial_flows, &initial_heads, 0, verbose).unwrap();

      // store the results
      flows[0] = flow.clone();
      heads[0] = head.clone();

      // update the initial flows and heads for the next step
      initial_flows = flow;
      initial_heads = head;

      // do parallel solves using Rayon
      let results: Vec<(Vec<f64>, Vec<f64>)> = (1..steps).into_par_iter().map(|step| {
        self.solve(&initial_flows, &initial_heads, step, verbose).unwrap()
      }).collect();
      for (step, (flow, head)) in results.iter().enumerate() {
        flows[step+1] = flow.clone();
        heads[step+1] = head.clone();
      }

    } else {
      // do sequential solves
      for step in 0..steps {
        let (flow, head) = self.solve(&initial_flows, &initial_heads, step, verbose).unwrap();
        flows[step] = flow.clone();
        heads[step] = head.clone();
        initial_flows = flow;
        initial_heads = head;
      }
    }
    let mut result = SolverResult { flows, heads };
    // convert the units back to the original units
    result.convert_units(&self.network.options.flow_units, &self.network.options.pressure_units);
    result
  }

  /// Solve the network using the Global Gradient Algorithm (Todini & Pilati, 1987) for a single step
  /// Returns the flows and heads
  fn solve(&self, initial_flows: &Vec<f64>, initial_heads: &Vec<f64>, step: usize, verbose: bool) -> Result<(Vec<f64>, Vec<f64>), String> {

    let time_options = &self.network.options.time_options;
    let mut flows = initial_flows.clone();
    let mut heads = initial_heads.clone();


    // get the time
    let time = step * time_options.hydraulic_timestep;
    // get the pattern time
    let pattern_time = time_options.pattern_start + time;
    // get pattern index
    let pattern_index = pattern_time / time_options.pattern_timestep;

    // apply the head pattern to reservoirs with a head pattern
    for (i, node) in self.network.nodes.iter().enumerate() {
      let Some(head_pattern) = node.head_pattern() else { continue };
      let pattern = &self.network.patterns[head_pattern];
      // TODO: FIX THIS TO USE THE CORRECT UNIT CONVERSION
      heads[i] = pattern.multipliers[pattern_index % pattern.multipliers.len()] / 0.3048;
      // println!("{}:{}", node.id, heads[i]);

    }

    // gather demands
    let demands: Vec<f64> = self.get_demands(pattern_index);
    // gather link initial statuses
    let mut statuses: Vec<LinkStatus> = self.get_initial_status();

    // calculate the resistances of the links
    let resistances: Vec<f64> = self.network.links.iter().map(|l| l.resistance()).collect::<Vec<f64>>();

    // perform GGA iterations
    // solve the system of equations: A * h = rhs
    // where A is the Jacobian matrix, h is the vector of heads, and rhs is the vector of right-hand side values
    // A is a sparse matrix, so we use the Compressed Sparse Column (CSC) format to store it
    // h is the vector of heads, and rhs is the vector of right-hand side values

    let unknown_nodes = self.node_to_unknown.iter().filter(|&x| x.is_some()).count();

    let mut values = vec![0.0; self.sparsity_pattern.as_ref().row_idx().len()]; // Jacobian matrix values
    let mut rhs = vec![0.0; unknown_nodes]; // (unknown nodes only)

    let mut g_invs: Vec<f64> = vec![0.0; self.network.links.len()];
    let mut ys: Vec<f64> = vec![0.0; self.network.links.len()];
    let mut jac = self.jac.clone();

    for iteration in 1..=self.network.options.max_trials {
      // reset values and rhs
      values.fill(0.0);
      rhs.fill(0.0);

      // set RHS to -demand (unknown nodes only)
      for (global, &head_id) in self.node_to_unknown.iter().enumerate() {
        if let Some(i) = head_id {
          rhs[i] = -demands[global];
        }
      }

      // clone the symbolic LLT to avoid borrowing issues
      let symbolic_llt = self.symbolic_llt.clone();

      // assemble Jacobian and RHS contributions from links
      for (i, link) in self.network.links.iter().enumerate() {
        let q = flows[i];
        let csc_index = &self.csc_indices[i];
        let (g_inv, y, status) = link.coefficients(q, resistances[i], statuses[i]);

        g_invs[i] = g_inv;
        ys[i] = y;
        // update the status of the link
        statuses[i] = status;

        // Get the CSC indices for the start and end nodes
        let u = self.node_to_unknown[link.start_node];
        let v = self.node_to_unknown[link.end_node];

        if let Some(i) = u {
            values[csc_index.diag_u.unwrap()] += g_inv;
            rhs[i] -= q - y;
            if self.network.nodes[link.end_node].is_fixed() {
              rhs[i] += g_inv * heads[link.end_node];
            }
        }
        if let Some(j) = v {
            values[csc_index.diag_v.unwrap()] += g_inv;
            rhs[j] += q - y;
            if self.network.nodes[link.start_node].is_fixed() {
              rhs[j] += g_inv * heads[link.start_node];
            }
        }
        if let (Some(_i), Some(_j)) = (u, v) {
            values[csc_index.off_diag_uv.unwrap()] -= g_inv;
            values[csc_index.off_diag_vu.unwrap()] -= g_inv;
        }
      }

      // solve the system of equations: J * dh = rhs
      jac.val_mut().copy_from_slice(&values);
      
      // Perform numerical factorization using pre-computed symbolic factorization
      let llt = Llt::try_new_with_symbolic(symbolic_llt, jac.as_ref(), Side::Lower)
        .expect("Singular matrix – check connectivity");

      
      let dh = llt.solve(&Mat::from_fn(unknown_nodes, 1, |r, _| rhs[r]));

      // update the heads of the nodes
      for (global, &head_id) in self.node_to_unknown.iter().enumerate() {
        if let Some(i) = head_id {
          heads[global] = dh[(i, 0)];
        }
      }

      // update the flows of the links (Equation 12.12)
      let mut sum_dq = 0.0;
      let mut sum_q  = 0.0;

      for (i, link) in self.network.links.iter().enumerate() {
          // calculate the head difference between the start and end nodes
          let dh = heads[link.start_node] - heads[link.end_node];

          // calculate the 1/G_ij and Y_ij coefficients
          let g_inv = g_invs[i];
          let y = ys[i];
        
          // Flow update: dq = y - g_inv * dh
          let dq = y - g_inv * dh;

          // update the flow of the link
          flows[i] -= dq;
          
          // update the sum of the absolute changes in flow and the sum of the absolute flows
          sum_dq += dq.abs();
          sum_q  += flows[i].abs();
      }

      let rel_change = sum_dq / (sum_q + 1e-6);

      if rel_change < self.network.options.accuracy {


        if verbose {
          let flow_balance = self.flow_balance(&demands, &flows);
          println!("Converged in {} iterations: Error = {:.4}, Supply = {:.4}, Demand = {:.4}", iteration, flow_balance.error, flow_balance.total_supply, flow_balance.total_demand);
        }

        return Ok((flows, heads));
      }
    }
    Err(format!("Maximum number of iterations reached: {}", self.network.options.max_trials))
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
          csc_index.off_diag_uv = find_csc_index(sym, i, j);
          csc_index.off_diag_vu = find_csc_index(sym, j, i);
      }
      csc_indices.push(csc_index);
    }
    csc_indices
  }

  /// Build sparsity pattern
  /// 
  fn build_sparsity_pattern(network: &Network, node_to_unknown: &Vec<Option<usize>>) -> SymbolicSparseColMat<usize> {
    let n_unknowns = node_to_unknown.iter().filter(|&x| x.is_some()).count();

    // Pre-allocate: at most 4 triplets per link (2 diagonal + 2 off-diagonal)
    let mut triplets = Vec::with_capacity(4 * network.links.len());
    for link in network.links.iter() {
      let u = node_to_unknown[link.start_node]; // u is the index of the start node (None if unknown)
      let v = node_to_unknown[link.end_node]; // v is the index of the end node (None if unknown)
      // diagonal elements (self-connectivity)
      if let Some(i) = u { triplets.push(Triplet::new(i, i, 0.0)); } // add diagonal element for u
      if let Some(j) = v { triplets.push(Triplet::new(j, j, 0.0)); } // add diagonal element for v
      // off diagonal elements (connectivity)
      if let (Some(i), Some(j)) = (u, v) {
        triplets.push(Triplet::new(i, j, 0.0)); // add off-diagonal element for u,v
        triplets.push(Triplet::new(j, i, 0.0)); // add off-diagonal element for v,u
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

  fn get_initial_status(&self) -> Vec<LinkStatus> {
    self.network.links.iter().map(|l| l.initial_status).collect::<Vec<LinkStatus>>()
  }
  /// Compute the initial flows ofthe links (1 ft/s velocity)
  fn get_initial_flows(&self) -> Vec<f64> {
    self.network.links.iter().map(|l| {
      if let LinkType::Pipe(pipe) = &l.link_type {
        return 0.25 * std::f64::consts::PI * pipe.diameter.powi(2);
      } else if let LinkType::Pump(pump) = &l.link_type {
        if let Some(head_curve) = &pump.head_curve {
          return head_curve.statistics.q_initial;
        } else {
          return 0.0;
        }
      } else if let LinkType::Valve(valve) = &l.link_type {
        return 0.25 * std::f64::consts::PI * valve.diameter.powi(2);
      } else {
        return 0.0;
      }
    }).collect::<Vec<f64>>()
  }

  /// Get the demands of the junctions for a given pattern index
  fn get_demands(&self, pattern_index: usize) -> Vec<f64> {
    let demands = self.network.nodes.iter().map(|n| {
          if let NodeType::Junction(junction) = &n.node_type {
            if let Some(pattern_id) = &junction.pattern {
              // TODO: Can be improved to reduce hashmap lookups by storing the pattern index in the junction
              let pattern = &self.network.patterns[pattern_id];
              // get the multiplier for the pattern index (wrap around if needed)
              let multiplier = pattern.multipliers[pattern_index % pattern.multipliers.len()];
              return junction.basedemand * multiplier;
            }
            // if no pattern, return the basedemand
            return junction.basedemand
          } else {
            return 0.0;
          }
        }).collect::<Vec<f64>>();
    demands
  }

  /// Set the initial heads of the nodes
  fn get_initial_heads(&self) -> Vec<f64> {
    self.network.nodes.iter().map(|n| {
      if n.is_fixed() {
        return n.elevation;
      } else {
        return 0.0;
      }
    }).collect::<Vec<f64>>()
  }

}

/// Helper function to find the CSC index for a given row and column
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
        .map(|pos| col_start + pos)
}