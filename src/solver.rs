use faer::sparse::{SparseColMat, Triplet};
use faer::traits::num_traits::ToPrimitive;
use faer::sparse::SymbolicSparseColMat;
use faer::{Mat, Side};
use faer::prelude::*;

use crate::network::*;

const MAX_ITER: usize = 50;
const CONVERGENCE_TOL: f64 = 1e-3;

const A1: f64 = 2.7662548476475583e-06; // Hazen-Williams conversion for m3/h to m3/s (1/3600)**H_EXPONENT * 10.67
const H_EXPONENT: f64 = 1.852; // Hazen-Williams exponent

impl Link {
  pub fn update_resistance(&mut self) {
    match self.link_type {
      LinkType::Pipe { diameter, length, roughness } => {
        self.resistance = A1 * length / ((diameter/1000.0).powf(4.87) * roughness.powf(H_EXPONENT));
      }
      _ => (),
    }
  }
}

impl Network {
  /// Solve the network using the Global Gradient Algorithm (Todini & Pilati, 1987)
  /// 
  /// The following steps are performed:
  /// 1. Build global unknown-numbering map
  /// 2. Build sparsity pattern (symbolic phase)
  /// 3. Map each pipe to its CSC indices
  /// 4. Update the resistance of all links
  /// 5. Set demand for unknown nodes
  /// 
  pub fn solve(&mut self) -> Result<(), String> {

    // build global unknown-numbering map
    let node_to_unknown = self.build_unknown_numbering_map();

    // generate sparsity pattern
    let sparsity_pattern = self.build_sparsity_pattern(&node_to_unknown);

    // map each link to its CSC indices
    self.map_links_to_csc_indices(&sparsity_pattern, &node_to_unknown);

    // update the resistance of all links
    self.update_link_resistances();

    // set demand for unknown nodes
    for node in self.nodes.iter_mut() {
      if let NodeType::Junction { basedemand } = node.node_type {
        node.demand = basedemand;
      }
    }

    // perform GGA iterations
    // solve the system of equations: J * dh = rhs

    let mut values = vec![0.0; sparsity_pattern.as_ref().row_idx().len()]; // Jacobian matrix values
    let mut rhs = vec![0.0; node_to_unknown.len()]; // RHS = –demand  (unknown nodes only)

    for _ in 1..=MAX_ITER {
      // reset values and rhs
      values.fill(0.0);
      rhs.fill(0.0);

      // set RHS to -demand (unknown nodes only)
      for (global, &head_id) in node_to_unknown.iter().enumerate() {
        if let Some(i) = head_id {
          rhs[i] = -self.nodes[global].demand;
        }
      }

      // assemble Jacobian and RHS contributions from links
      for link in self.links.iter() {
        let q_abs = link.result.flow.abs().max(1e-8); // get absolute flow, avoid division by zero by setting a small value
        let r_term = link.resistance * q_abs.powf(H_EXPONENT - 1.0); // resistance term
        let gradient = H_EXPONENT * r_term; // gradient term
        let h_loss = r_term * link.result.flow; // head loss term
        let g = 1.0 / gradient; // conductance term
        let y = link.result.flow - (h_loss / gradient); // head loss term

        // get the CSC indices for the start and end nodes
        let u = node_to_unknown[link.start_node];
        let v = node_to_unknown[link.end_node];

        if let Some(i) = u {
          values[link.csc_index.diag_u.unwrap()] += g;
          rhs[i] -= y;
          // if the end node is not a junction (i.e. a reservoir or tank), add the head of the end node to the RHS
          if !matches!(self.nodes[link.end_node].node_type, NodeType::Junction { .. }) {
            rhs[i] += g * self.nodes[link.end_node].result.head;
          }
        }
        if let Some(j) = v {
          values[link.csc_index.diag_v.unwrap()] += g;
          rhs[j] += y;
          // if the start node is not a junction (i.e. a reservoir or tank), add the head of the start node to the RHS
          if !matches!(self.nodes[link.start_node].node_type, NodeType::Junction { .. }) {
            rhs[j] += g * self.nodes[link.start_node].result.head;
          }
        }
        if let (Some(i), Some(j)) = (u, v) {
          values[link.csc_index.off_diag_uv.unwrap()] -= g;
          values[link.csc_index.off_diag_vu.unwrap()] -= g;
        }
      }
    
      // solve the system of equations: J * dh = rhs
      let jac = SparseColMat::new(sparsity_pattern.clone(), values.clone());
      let solver = jac.sp_cholesky(Side::Lower).expect("Singular matrix – check connectivity");
      let dh = solver.solve(&Mat::from_fn(node_to_unknown.len(), 1, |r, _| rhs[r]));

      // update the heads of the nodes
      for (global, &head_id) in node_to_unknown.iter().enumerate() {
        if let Some(i) = head_id {
          self.nodes[global].result.head = dh[(i, 0)];
        }
      }
      // update the flows of the links
      let mut sum_dq = 0.0;
      let mut sum_q  = 0.0;

      for link in self.links.iter_mut() {
        // calculate the head difference between the start and end nodes
        let dh = self.nodes[link.start_node].result.head - self.nodes[link.end_node].result.head;
        let q_abs = link.result.flow.abs().max(1e-8);
        let r_coef = link.resistance * q_abs.powf(H_EXPONENT - 1.0);
        let h_loss = r_coef * link.result.flow;
        let grad = H_EXPONENT * r_coef;
        let dq = (dh - h_loss) / grad;
        link.result.flow += dq;
        // update the sum of the absolute changes in flow and the sum of the absolute flows
        sum_dq += dq.abs();
        sum_q  += link.result.flow.abs();
      }
      let rel_change = sum_dq / (sum_q + 1e-6);
      if rel_change < CONVERGENCE_TOL {
        return Ok(());
      }
    }
    Err(format!("Maximum number of iterations reached: {}", MAX_ITER))
  }

  /// Update the resistance of all links
  fn update_link_resistances(&mut self) {
    for link in self.links.iter_mut() {
      link.update_resistance();
    }
  }

  /// Map each link to its CSC (Compressed Sparse Column) indices
  fn map_links_to_csc_indices(&mut self, sparsity_pattern: &SymbolicSparseColMat<usize>, node_to_unknown: &Vec<Option<usize>>) {
    for link in self.links.iter_mut() {
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
      link.csc_index = csc_index;
    }
  }

  /// Build sparsity pattern
  /// 
  fn build_sparsity_pattern(&self, node_to_unknown: &Vec<Option<usize>>) -> SymbolicSparseColMat<usize> {
    let n_unknowns = node_to_unknown.iter().filter(|&x| x.is_some()).count();

    let mut triplets = Vec::new();
    for link in self.links.iter() {
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
  fn build_unknown_numbering_map(&self) -> Vec<Option<usize>> {
    let mut unknown_id = 0;
    let node_to_unknown: Vec<Option<usize>> = self.nodes
        .iter()
        .map(|n| if matches!(n.node_type, NodeType::Junction { .. }) { let id = unknown_id; unknown_id += 1; Some(id) } else { None })
        .collect();
    node_to_unknown
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