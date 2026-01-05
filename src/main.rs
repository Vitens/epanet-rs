use faer::sparse::{SparseColMat, Triplet};
use faer::{Mat, Side};
use faer::prelude::*;
use rand::Rng;
use std::time::Instant;

const MAX_ITER: usize = 20;
const CONVERGENCE_TOL: f64 = 1e-5;
const H_EXPONENT: f64 = 2.0; 

#[derive(Clone, Debug)]
struct Pipe {
    start: usize, end: usize,
    resistance: f64, flow: f64,
    diag_u: Option<usize>, diag_v: Option<usize>,
    off_diag_uv: Option<usize>, off_diag_vu: Option<usize>,
}

#[derive(Clone, Debug)]
struct Node { head: f64, demand: f64, is_fixed: bool }

fn find_csc_index(symbolic: faer::sparse::SymbolicSparseColMatRef<usize>, row: usize, col: usize) -> Option<usize> {
    let col_ptr = symbolic.col_ptr();
    let row_idx = symbolic.row_idx();
    let start = col_ptr[col];
    let end = col_ptr[col + 1];
    row_idx[start..end].iter().position(|&r| r == row).map(|pos| start + pos)
}

fn main() {
    let dim = 100; 
    let (mut nodes, mut pipes) = generate_grid(dim);
    let n_unknowns = nodes.len() - 1;

    // --- 1. SYMBOLIC PHASE ---
    let mut tmp_triplets = Vec::new();
    for p in &pipes {
        let u = if nodes[p.start].is_fixed { None } else { Some(p.start - 1) };
        let v = if nodes[p.end].is_fixed { None } else { Some(p.end - 1) };
        if let Some(i) = u { tmp_triplets.push(Triplet::new(i, i, 0.0)); }
        if let Some(j) = v { tmp_triplets.push(Triplet::new(j, j, 0.0)); }
        if let (Some(i), Some(j)) = (u, v) {
            tmp_triplets.push(Triplet::new(i, j, 0.0));
            tmp_triplets.push(Triplet::new(j, i, 0.0));
        }
    }
    let matrix_tmp = SparseColMat::<usize, f64>::try_new_from_triplets(n_unknowns, n_unknowns, &tmp_triplets).unwrap();
    let symbolic = matrix_tmp.symbolic().to_owned().unwrap();

    for pipe in pipes.iter_mut() {
        let u = if nodes[pipe.start].is_fixed { None } else { Some(pipe.start - 1) };
        let v = if nodes[pipe.end].is_fixed { None } else { Some(pipe.end - 1) };
        let sym = symbolic.as_ref();
        if let Some(i) = u { pipe.diag_u = find_csc_index(sym, i, i); }
        if let Some(j) = v { pipe.diag_v = find_csc_index(sym, j, j); }
        if let (Some(i), Some(j)) = (u, v) {
            pipe.off_diag_uv = find_csc_index(sym, i, j);
            pipe.off_diag_vu = find_csc_index(sym, j, i);
        }
    }

    // --- 2. SOLVE PHASE ---
    let nnz = symbolic.as_ref().row_idx().len();
    let mut matrix_values = vec![0.0f64; nnz];
    let mut rhs_vec = vec![0.0f64; n_unknowns];

    println!("Starting EPANET-style GGA Solve...");

    for iter in 0..MAX_ITER {
        matrix_values.fill(0.0);
        rhs_vec.fill(0.0);

        // Load demands into RHS
        for i in 0..n_unknowns { rhs_vec[i] = -nodes[i+1].demand; }

        for pipe in &pipes {
            let q_abs = pipe.flow.abs().max(1e-8);
            // Gradient of headloss: n * R * Q^(n-1)
            let f_prime = H_EXPONENT * pipe.resistance * q_abs.powf(H_EXPONENT - 1.0);
            let head_loss = pipe.resistance * pipe.flow * q_abs.powf(H_EXPONENT - 1.0);
            
            let g = 1.0 / f_prime;
            let y = head_loss / f_prime - pipe.flow; // The GGA Flow Correction

            if let Some(idx) = pipe.diag_u { 
                matrix_values[idx] += g; 
                rhs_vec[pipe.start - 1] -= y;
                if nodes[pipe.end].is_fixed { rhs_vec[pipe.start - 1] += g * nodes[pipe.end].head; }
            }
            if let Some(idx) = pipe.diag_v { 
                matrix_values[idx] += g; 
                rhs_vec[pipe.end - 1] += y;
                if nodes[pipe.start].is_fixed { rhs_vec[pipe.end - 1] += g * nodes[pipe.start].head; }
            }
            if let Some(idx) = pipe.off_diag_uv { matrix_values[idx] -= g; }
            if let Some(idx) = pipe.off_diag_vu { matrix_values[idx] -= g; }
        }

        let matrix_a = SparseColMat::new(symbolic.clone(), matrix_values.clone());
        let solver = matrix_a.sp_cholesky(Side::Lower).expect("Singular Matrix");
        let new_heads = solver.solve(&Mat::from_fn(n_unknowns, 1, |i, _| rhs_vec[i]));

        let mut max_dh = 0.0f64;
        for i in 0..n_unknowns {
            let diff = (new_heads[(i, 0)] - nodes[i+1].head).abs();
            max_dh = max_dh.max(diff);
            nodes[i+1].head = new_heads[(i, 0)];
        }

        // Update flows using the new heads
        for pipe in pipes.iter_mut() {
            let h_diff = nodes[pipe.start].head - nodes[pipe.end].head;
            let q_abs = pipe.flow.abs().max(1e-8);
            let f_prime = H_EXPONENT * pipe.resistance * q_abs.powf(H_EXPONENT - 1.0);
            let head_loss = pipe.resistance * pipe.flow * q_abs.powf(H_EXPONENT - 1.0);
            pipe.flow -= (head_loss - h_diff) / f_prime;
        }

        println!("Iter {}: max_dH = {:.2e}", iter + 1, max_dh);
        if max_dh < CONVERGENCE_TOL { break; }
    }
}

fn generate_grid(dim: usize) -> (Vec<Node>, Vec<Pipe>) {
    let mut rng = rand::rng(); 
    let mut nodes = Vec::new();
    let mut pipes = Vec::new();
    for r in 0..dim {
        for c in 0..dim {
            let is_res = r == 0 && c == 0;
            nodes.push(Node {
                head: 100.0,
                demand: if is_res { 0.0 } else { rng.random_range(0.001..0.005) },
                is_fixed: is_res,
            });
        }
    }
    for r in 0..dim {
        for c in 0..dim {
            let u = r * dim + c;
            if c < dim - 1 { pipes.push(create_pipe(u, u + 1, &mut rng)); }
            if r < dim - 1 { pipes.push(create_pipe(u, u + dim, &mut rng)); }
        }
    }
    (nodes, pipes)
}

fn create_pipe(s: usize, e: usize, rng: &mut rand::rngs::ThreadRng) -> Pipe {
    Pipe { start: s, end: e, resistance: rng.random_range(1.0..10.0), flow: 0.01,
           diag_u: None, diag_v: None, off_diag_uv: None, off_diag_vu: None }
}