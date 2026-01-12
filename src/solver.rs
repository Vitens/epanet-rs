use faer::sparse::{SparseColMat, Triplet};
use faer::sparse::SymbolicSparseColMat;
use faer::sparse::linalg::solvers::{SymbolicLlt, Llt};
use faer::{Mat, Side};
use faer::prelude::*;

use rayon::prelude::*;

use crate::network::*;

const MAX_ITER: usize = 200;
const CONVERGENCE_TOL: f64 = 0.01;

const H_EXPONENT: f64 = 1.852; // Hazen-Williams exponent

impl Network {
  /// Solve the network using the Global Gradient Algorithm (Todini & Pilati, 1987) for a single time step
  /// 
  /// The following steps are performed:
  /// 1. Build global unknown-numbering map
  /// 2. Build sparsity pattern (symbolic phase)
  /// 3. Map each pipe to its CSC indices
  /// 4. Update the resistance of all links
  /// 5. Set demand for unknown nodes
  /// 
  pub fn run(&mut self, parallel: bool) {
    // build global unknown-numbering map
    let node_to_unknown = self.build_unknown_numbering_map();

    // generate sparsity pattern
    let sparsity_pattern = self.build_sparsity_pattern(&node_to_unknown);

    // map each link to its CSC (Compressed Sparse Column) indices
    self.map_links_to_csc_indices(&sparsity_pattern, &node_to_unknown);

    let values = vec![0.0; sparsity_pattern.as_ref().row_idx().len()]; // Jacobian matrix values
    let jac = SparseColMat::new(sparsity_pattern.clone(), values.clone());

    let symbolic_llt = SymbolicLlt::try_new(jac.symbolic(), Side::Lower)
      .expect("Failed to compute symbolic Cholesky factorization");

    let steps = 24*4;

    let mut flows: Vec<Vec<f64>> = vec![vec![0.0; self.links.len()]; steps];
    let mut heads: Vec<Vec<f64>> = vec![vec![0.0; self.nodes.len()]; steps];


    if parallel {
      // do parallel solve with Rayon
      let results: Vec<(Vec<f64>, Vec<f64>)> = (0..steps).into_par_iter().map(|step| {
        self.solve(&node_to_unknown, &sparsity_pattern, &symbolic_llt, &jac).unwrap()
      }).collect();
      for (step, (flow, head)) in results.iter().enumerate() {
        flows[step] = flow.clone();
        heads[step] = head.clone();
      }

    } else {
      // do sequential solves
      for step in 0..steps {
        let (flow, head) = self.solve(&node_to_unknown, &sparsity_pattern, &symbolic_llt, &jac).unwrap();
        flows[step] = flow;
        heads[step] = head;
      }
    }

    // assign the flows and heads to the links and nodes
    for (i, link) in self.links.iter_mut().enumerate() {
      link.result.flow = flows[0][i];
    }
    for (i, node) in self.nodes.iter_mut().enumerate() {
      node.result.head = heads[0][i];
    }


  }

  pub fn solve(&self, node_to_unknown: &Vec<Option<usize>>, sparsity_pattern: &SymbolicSparseColMat<usize>, symbolic_llt: &SymbolicLlt<usize>, jac: &SparseColMat<usize, f64>) -> Result<(Vec<f64>, Vec<f64>), String> {

    // initialize the flows and heads
    let mut flows: Vec<f64> = self.links.iter().map(|l| {
      if let LinkType::Pipe { diameter, .. } = l.link_type {
        let area = 0.25 * std::f64::consts::PI * diameter.powi(2);
        let velocity = 1.0;
        let flow = area * velocity;
        return flow;
      } else {
        return 0.0;
      }
    }).collect::<Vec<f64>>();

    // initialize the heads
    let mut heads: Vec<f64> = self.nodes.iter().map(|n| {
      if n.is_fixed() {
        return n.elevation;
      } else {
        return 0.0;
      }
    }).collect::<Vec<f64>>();

    // gather demands
    let demands: Vec<f64> = self.nodes.iter().map(|n| {
      if let NodeType::Junction { basedemand } = n.node_type {
        return basedemand;
      } else {
        return 0.0;
      }
    }).collect::<Vec<f64>>();

    let resistances: Vec<f64> = self.links.iter().map(|l| {
      if let LinkType::Pipe { diameter, length, roughness } = l.link_type {
        return 4.727 * roughness.powf(-1.852) * diameter.powf(-4.87) * length;
      } else {
        return 0.0;
      }
    }).collect::<Vec<f64>>();

    // perform GGA iterations
    // solve the system of equations: A * h = rhs
    // where A is the Jacobian matrix, h is the vector of heads, and rhs is the vector of right-hand side values
    // A is a sparse matrix, so we use the Compressed Sparse Column (CSC) format to store it
    // h is the vector of heads, and rhs is the vector of right-hand side values

    let unknown_nodes = node_to_unknown.iter().filter(|&x| x.is_some()).count();

    let mut values = vec![0.0; sparsity_pattern.as_ref().row_idx().len()]; // Jacobian matrix values
    let mut rhs = vec![0.0; unknown_nodes]; // (unknown nodes only)

    let mut g_invs: Vec<f64> = vec![0.0; self.links.len()];
    let mut ys: Vec<f64> = vec![0.0; self.links.len()];
    let mut jac = jac.clone();

    for iteration in 1..=MAX_ITER {
      // reset values and rhs
      values.fill(0.0);
      rhs.fill(0.0);

      // set RHS to -demand (unknown nodes only)
      for (global, &head_id) in node_to_unknown.iter().enumerate() {
        if let Some(i) = head_id {
          rhs[i] = -demands[global];
        }
      }

      // clone the symbolic LLT to avoid borrowing issues
      let symbolic_llt = symbolic_llt.clone();

      // assemble Jacobian and RHS contributions from links
      for (i, link) in self.links.iter().enumerate() {
        let q = flows[i];
        let (g_inv, y) = link.coefficients(q, resistances[i]);

        g_invs[i] = g_inv;
        ys[i] = y;

        // Get the CSC indices for the start and end nodes
        let u = node_to_unknown[link.start_node];
        let v = node_to_unknown[link.end_node];

        if let Some(i) = u {
            values[link.csc_index.diag_u.unwrap()] += g_inv;
            rhs[i] -= q - (y * g_inv);
            if self.nodes[link.end_node].is_fixed() {
              rhs[i] += g_inv * heads[link.end_node];
            }
        }
        if let Some(j) = v {
            values[link.csc_index.diag_v.unwrap()] += g_inv;
            rhs[j] += q - (y * g_inv);
            if self.nodes[link.start_node].is_fixed() {
              rhs[j] += g_inv * heads[link.start_node];
            }
        }
        if let (Some(_i), Some(_j)) = (u, v) {
            values[link.csc_index.off_diag_uv.unwrap()] -= g_inv;
            values[link.csc_index.off_diag_vu.unwrap()] -= g_inv;
        }
      }

      // solve the system of equations: J * dh = rhs
      jac.val_mut().copy_from_slice(&values);
      
      // Perform numerical factorization using pre-computed symbolic factorization
      let llt = Llt::try_new_with_symbolic(symbolic_llt, jac.as_ref(), Side::Lower)
        .expect("Singular matrix – check connectivity");

      
      let dh = llt.solve(&Mat::from_fn(unknown_nodes, 1, |r, _| rhs[r]));

      // update the heads of the nodes
      for (global, &head_id) in node_to_unknown.iter().enumerate() {
        if let Some(i) = head_id {
          heads[global] = dh[(i, 0)];
        }
      }

      // update the flows of the links (Equation 12.12)
      let mut sum_dq = 0.0;
      let mut sum_q  = 0.0;

      for (i, link) in self.links.iter().enumerate() {
          // calculate the head difference between the start and end nodes
          let dh = heads[link.start_node] - heads[link.end_node];

          // calculate the 1/G_ij and Y_ij coefficients
          let g_inv = g_invs[i];
          let y = ys[i];
        
          // Flow update: dq = (h_L - dh) / g
          let dq = g_inv*(y - dh);

          flows[i] -= dq;
          
          // update the sum of the absolute changes in flow and the sum of the absolute flows
          sum_dq += dq.abs();
          sum_q  += flows[i].abs();
      }

      let rel_change = sum_dq / (sum_q + 1e-6);
      // println!("Iteration {} relative change = {:.4}", iteration, rel_change);

      if rel_change < CONVERGENCE_TOL {
        return Ok((flows, heads));
      }
    }
    Err(format!("Maximum number of iterations reached: {}", MAX_ITER))
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

    // Pre-allocate: at most 4 triplets per link (2 diagonal + 2 off-diagonal)
    let mut triplets = Vec::with_capacity(4 * self.links.len());
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