// main.rs
use faer::sparse::{SparseColMat, Triplet};
use faer::{Mat, Side};
use faer::prelude::*;
use std::time::Instant;
use std::env;
use rand::Rng;

const MAX_ITER: usize = 50;
const CONVERGENCE_TOL: f64 = 1e-3;
const H_EXPONENT: f64 = 2.0; // Hazen-Williams exponent

#[derive(Clone, Debug)]
struct Pipe {
    start: usize,
    end: usize,
    resistance: f64,
    flow: f64,
    diag_u: Option<usize>,      // CSC index for J[u,u]
    diag_v: Option<usize>,      // CSC index for J[v,v]
    off_diag_uv: Option<usize>, // CSC index for J[u,v]
    off_diag_vu: Option<usize>, // CSC index for J[v,u]
}

#[derive(Clone, Debug)]
struct Node {
    head: f64,
    demand: f64,
    is_fixed: bool,
}

fn main() {
    //------------------------------------------------------------------
    // 1.  Build a small test grid
    //------------------------------------------------------------------
    let dim = env::args().nth(1).unwrap_or("100".to_string()).parse::<usize>().unwrap();
    let (mut nodes, mut pipes) = generate_grid(dim);

    //------------------------------------------------------------------
    // 2.  Build global unknown-numbering map
    //------------------------------------------------------------------
    let start_time = Instant::now();

    let mut unknown_id = 0;
    let node_to_unknown: Vec<Option<usize>> = nodes
        .iter()
        .map(|n| if n.is_fixed { None } else { let id = unknown_id; unknown_id += 1; Some(id) })
        .collect();
    let n_unknowns = unknown_id;

    //------------------------------------------------------------------
    // 3.  Symbolic phase – build sparsity pattern
    //------------------------------------------------------------------
    let mut triplets = Vec::new();
    for p in &pipes {
        let u = node_to_unknown[p.start];
        let v = node_to_unknown[p.end];
        if let Some(i) = u { triplets.push(Triplet::new(i, i, 0.0)); }
        if let Some(j) = v { triplets.push(Triplet::new(j, j, 0.0)); }
        if let (Some(i), Some(j)) = (u, v) {
            triplets.push(Triplet::new(i, j, 0.0));
            triplets.push(Triplet::new(j, i, 0.0));
        }
    }
    let tmp = SparseColMat::<usize, f64>::try_new_from_triplets(n_unknowns, n_unknowns, &triplets)
        .unwrap();
    let symbolic = tmp.symbolic().to_owned().unwrap();

    //------------------------------------------------------------------
    // 4.  Map each pipe to its CSC indices
    //------------------------------------------------------------------
    for p in pipes.iter_mut() {
        let u = node_to_unknown[p.start];
        let v = node_to_unknown[p.end];
        let sym = symbolic.as_ref();
        if let Some(i) = u { p.diag_u = find_csc_index(sym, i, i); }
        if let Some(j) = v { p.diag_v = find_csc_index(sym, j, j); }
        if let (Some(i), Some(j)) = (u, v) {
            p.off_diag_uv = find_csc_index(sym, i, j);
            p.off_diag_vu = find_csc_index(sym, j, i);
        }
    }

    //------------------------------------------------------------------
    // 5.  GGA iterations
    //------------------------------------------------------------------
    println!("Starting GGA solve ({} unknowns)", n_unknowns);
    for iter in 1..=MAX_ITER {
        let mut values = vec![0.0; symbolic.as_ref().row_idx().len()];
        let mut rhs = vec![0.0; n_unknowns];

        // RHS = –demand  (unknown nodes only)
        for (global, &head_id) in node_to_unknown.iter().enumerate() {
            if let Some(i) = head_id {
                rhs[i] = -nodes[global].demand;
            }
          }

        // Assemble Jacobian and RHS contributions from pipes
        for p in &pipes {
            let q_abs = p.flow.abs().max(1e-8);
            let r_term = p.resistance * q_abs.powf(H_EXPONENT - 1.0);
            let gradient = H_EXPONENT * r_term;
            let h_loss = r_term * p.flow;
            let g = 1.0 / gradient;
            // let y = (h_loss / gradient) - p.flow;
            let y = p.flow - (h_loss / gradient);

            let u = node_to_unknown[p.start];
            let v = node_to_unknown[p.end];

            if let Some(i) = u {
                values[p.diag_u.unwrap()] += g;
                rhs[i] -= y;
                if nodes[p.end].is_fixed {
                    rhs[i] += g * nodes[p.end].head;
                }
            }
            if let Some(j) = v {
                values[p.diag_v.unwrap()] += g;
                rhs[j] += y;
                if nodes[p.start].is_fixed {
                    rhs[j] += g * nodes[p.start].head;
                }
            }
            if let (Some(i), Some(j)) = (u, v) {
                values[p.off_diag_uv.unwrap()] -= g;
                values[p.off_diag_vu.unwrap()] -= g;
            }
        }

        //------------------------------------------------------------------
        // Solve linear system
        //------------------------------------------------------------------
        let jac = SparseColMat::new(symbolic.clone(), values.clone());

        let solver = jac.sp_cholesky(Side::Lower).expect("Singular matrix – check connectivity");
        let dh = solver.solve(&Mat::from_fn(n_unknowns, 1, |r, _| rhs[r]));

        //------------------------------------------------------------------
        // Update heads
        //------------------------------------------------------------------
        for (global, &head_id) in node_to_unknown.iter().enumerate() {
            if let Some(i) = head_id {
                // nodes[global].head += dh[(i, 0)];
                nodes[global].head = dh[(i, 0)];
            }
        }

        let mut sum_dq = 0.0;
        let mut sum_q  = 0.0;

        for (idx, p) in pipes.iter_mut().enumerate() {
            let dh = nodes[p.start].head - nodes[p.end].head;
            let q_abs  = p.flow.abs().max(1e-8);
            let r_coef = p.resistance * q_abs.powf(H_EXPONENT - 1.0);
            let h_loss = r_coef * p.flow;
            let grad   = H_EXPONENT * r_coef;
            let dq = (dh - h_loss) / grad;
            p.flow += dq;

            sum_dq += dq.abs();
            sum_q  += p.flow.abs();
        }
        let rel_change = sum_dq / (sum_q + 1e-6);
        println!("Iter {:2}: relative change = {:.3e}", iter, rel_change);
        if rel_change < CONVERGENCE_TOL {
            println!("Converged!");
            break;
        }
    }

    //------------------------------------------------------------------
    // 6.  Simple mass-balance report
    //------------------------------------------------------------------
    let total_demand: f64 = nodes.iter().map(|n| n.demand).sum();
    let reservoir_supply: f64 = pipes
        .iter()
        .map(|p| {
            if nodes[p.start].is_fixed { p.flow }
            else if nodes[p.end].is_fixed { -p.flow }
            else { 0.0 }
        })
        .sum();

    println!(
        "\nMass balance: demand = {:.4}, supply = {:.4}, error = {:.2e}",
        total_demand,
        reservoir_supply,
        (total_demand - reservoir_supply).abs()
    );
    // print how long the solve took
    let end_time = Instant::now();
    println!("Solve took {:?}", end_time.duration_since(start_time));

    //------------------------------------------------------------------
    // Optional: print heads/flows
    //------------------------------------------------------------------
    // for (i,n) in nodes.iter().enumerate() {
    //     println!("Node {} head = {:.3}", i, n.head);
    // }
    // for (i,p) in pipes.iter().enumerate() {
    //     println!("Pipe {} flow = {:.3}", i, p.flow);
    // }
}

fn generate_grid(dim: usize) -> (Vec<Node>, Vec<Pipe>) {
    let mut nodes = Vec::new();
    for r in 0..dim {
        for c in 0..dim {
            let is_res = r == 0 && c == 0; // top-left corner is reservoir
            // random demand between 0 and 10
            let demand = rand::rng().random_range(0.0..1.0);
            nodes.push(Node {
                head: if is_res { 100.0 } else { 30.0 },
                demand: if is_res { 0.0 } else { demand },
                is_fixed: is_res,
            });
        }
    }
    let mut pipes = Vec::new();
    for r in 0..dim {
        for c in 0..dim {
            let u = r * dim + c;
            if c + 1 < dim { pipes.push(create_pipe(u, u + 1)); }
            if r + 1 < dim { pipes.push(create_pipe(u, u + dim)); }
        }
    }
    (nodes, pipes)
}

fn create_pipe(s: usize, e: usize) -> Pipe {
  // random resistance between 0.01 and 0.02
    let resistance = rand::rng().random_range(0.01..0.02);
    Pipe {
        start: s,
        end: e,
        resistance: resistance,
        flow: 0.1,
        diag_u: None,
        diag_v: None,
        off_diag_uv: None,
        off_diag_vu: None,
    }
}

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