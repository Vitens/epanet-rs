use faer::sparse::{SparseColMat, Triplet};
use faer::sparse::SymbolicSparseColMat;
use faer::{Mat, Side};
use faer::prelude::*;

use std::time::Instant;

use crate::network::*;

const MAX_ITER: usize = 10;
const CONVERGENCE_TOL: f64 = 1e-3;

const H_EXPONENT: f64 = 1.852; // Hazen-Williams exponent

// impl Link {
//   pub fn update_resistance(&mut self) {
//     match self.link_type {
//       LinkType::Pipe { diameter, length, roughness } => {
//         self.resistance = A1 * length / ((diameter/1000.0).powf(4.87) * roughness.powf(H_EXPONENT));
//       }
//       _ => (),
//     }
//   }
// }

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

    // map each link to its CSC (Compressed Sparse Column) indices
    self.map_links_to_csc_indices(&sparsity_pattern, &node_to_unknown);

    // update the resistance of all links
    self.update_link_resistances();

    // set demand for unknown nodes
    for node in self.nodes.iter_mut() {
      if let NodeType::Junction { basedemand } = node.node_type {
        node.demand = basedemand;
      }
      if node.is_fixed() {
        node.result.head = node.elevation;
      }
    }
    for link in self.links.iter_mut() {
      if let LinkType::Pipe { diameter, .. } = link.link_type {
        // initialize the flow to correspond with 1.0 f/s
        let area = 0.25 * std::f64::consts::PI * (diameter).powf(2.0);
        let velocity = 1.0;
        let flow = area * velocity;
        link.result.flow = flow;
      }
    }

    // perform GGA iterations
    // solve the system of equations: A * h = rhs
    // where A is the Jacobian matrix, h is the vector of heads, and rhs is the vector of right-hand side values
    // A is a sparse matrix, so we use the Compressed Sparse Column (CSC) format to store it
    // h is the vector of heads, and rhs is the vector of right-hand side values

    let unknown_nodes = node_to_unknown.iter().filter(|&x| x.is_some()).count();

    let mut values = vec![0.0; sparsity_pattern.as_ref().row_idx().len()]; // Jacobian matrix values
    let mut rhs = vec![0.0; unknown_nodes]; // (unknown nodes only)

    let mut jac = SparseColMat::new(sparsity_pattern.clone(), values.clone());

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
        let q = link.result.flow;
        let q_abs = q.abs().max(1e-8);
        let r = link.resistance;
        let m = 0.0; // minor loss coefficient (0 for now)
        let n = H_EXPONENT;

        // Calculate head loss gradient (g_ij) - EPANET Eq. 12.8
        let g = n * r * q_abs.powf(n - 1.0) + 2.0 * m * q_abs;
        let g_inv = 1.0 / g;
        
        // Calculate head loss (h_Lij) - EPANET Eq. 12.1
        let y = (r * q_abs.powf(n) + m * q_abs.powf(2.0)) * q.signum();

        // Get the CSC indices for the start and end nodes
        let u = node_to_unknown[link.start_node];
        let v = node_to_unknown[link.end_node];

        if let Some(i) = u {
            values[link.csc_index.diag_u.unwrap()] += g_inv;
            rhs[i] -= q - (y * g_inv);
            if self.nodes[link.end_node].is_fixed() {
              rhs[i] += g_inv * self.nodes[link.end_node].result.head;
            }
        }
        if let Some(j) = v {
            values[link.csc_index.diag_v.unwrap()] += g_inv;
            rhs[j] += q - (y * g_inv);
            if self.nodes[link.start_node].is_fixed() {
              rhs[j] += g_inv * self.nodes[link.start_node].result.head;
            }
        }
        if let (Some(_i), Some(_j)) = (u, v) {
            values[link.csc_index.off_diag_uv.unwrap()] -= g_inv;
            values[link.csc_index.off_diag_vu.unwrap()] -= g_inv;
        }
      }

      // solve the system of equations: J * dh = rhs
      jac.val_mut().copy_from_slice(&values);
      let solver = jac.sp_cholesky(Side::Lower).expect("Singular matrix – check connectivity");
      let dh = solver.solve(&Mat::from_fn(unknown_nodes, 1, |r, _| rhs[r]));

      // update the heads of the nodes
      for (global, &head_id) in node_to_unknown.iter().enumerate() {
        if let Some(i) = head_id {
          self.nodes[global].result.head = dh[(i, 0)];
        }
      }
      // update the flows of the links - EPANET Eq. 12.12
      let mut sum_dq = 0.0;
      let mut sum_q  = 0.0;

      for link in self.links.iter_mut() {
          // calculate the head difference between the start and end nodes
          let dh = self.nodes[link.start_node].result.head - self.nodes[link.end_node].result.head;
          let q = link.result.flow;
          let q_abs = q.abs().max(1e-8);
          let r = link.resistance;
          let m = 0.0; // minor loss coefficient
          let n = H_EXPONENT;
          
          // Calculate head loss gradient (g_ij) - EPANET Eq. 12.8
          let g = n * r * q_abs.powf(n - 1.0) + 2.0 * m * q_abs;
  
          // Calculate head loss (h_Lij) - EPANET Eq. 12.1
          let y = (r * q_abs.powf(n) + m * q_abs.powf(2.0)) * q.signum();
        
          // Flow update: dq = (h_L - dh) / g
          let dq = 1.0/g*(y - dh);

          link.result.flow -= dq;
          
          // update the sum of the absolute changes in flow and the sum of the absolute flows
          sum_dq += dq.abs();
          sum_q  += link.result.flow.abs();
      }

      let rel_change = sum_dq / (sum_q + 1e-6);

      if rel_change < CONVERGENCE_TOL {
        return Ok(());
      }
    }
    // Err(format!("Maximum number of iterations reached: {}", MAX_ITER))
    Ok(())
  }

  /// Update the resistance of all links
  fn update_link_resistances(&mut self) {
    for link in self.links.iter_mut() {
      if let LinkType::Pipe { diameter, length, roughness } = link.link_type {
        link.resistance = 4.727 * roughness.powf(-1.852) * diameter.powf(-4.87) * length;
      }
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