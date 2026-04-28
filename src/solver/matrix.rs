//! Sparsity pattern, CSC indexing and AMD ordering helpers for the Jacobian matrix.

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::linalg::amd;
use faer::sparse::{SparseColMat, SymbolicSparseColMat, Triplet};

use crate::model::network::Network;
use crate::model::node::NodeType;

pub struct ResistanceCoefficients {
    pub g_inv: Vec<f64>,
    pub y: Vec<f64>,
}

impl ResistanceCoefficients {
    pub fn new(size: usize) -> Self {
        Self {
            g_inv: vec![0.0; size],
            y: vec![0.0; size],
        }
    }
}

/// CSC (Compressed Sparse Column) indices for the Jacobian matrix used in the Global Gradient Algorithm
#[derive(Default)]
pub struct CSCIndex {
    pub diag_u: Option<usize>,   // CSC index for J[u,u]
    pub diag_v: Option<usize>,   // CSC index for J[v,v]
    pub off_diag: Option<usize>, // CSC index for lower triangular off-diagonal entry
}

/// Helper function to find the CSC index for a given row and column
#[inline]
pub fn find_csc_index(
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

/// Build sparsity pattern
pub fn build_sparsity_pattern(
    network: &Network,
    node_to_unknown: &[Option<usize>],
) -> SymbolicSparseColMat<usize> {
    let n_unknowns = node_to_unknown.iter().filter(|&x| x.is_some()).count();

    let mut triplets = Vec::with_capacity(3 * network.links.len());

    for link in network.links.iter() {
        let u = node_to_unknown[link.start_node];
        let v = node_to_unknown[link.end_node];
        if let Some(i) = u {
            triplets.push(Triplet::new(i, i, 0.0));
        }
        if let Some(j) = v {
            triplets.push(Triplet::new(j, j, 0.0));
        }
        if let (Some(i), Some(j)) = (u, v) {
            let (row, col) = if i >= j { (i, j) } else { (j, i) };
            triplets.push(Triplet::new(row, col, 0.0));
        }
    }
    let sparsity_matrix =
        SparseColMat::try_new_from_triplets(n_unknowns, n_unknowns, &triplets).unwrap();

    sparsity_matrix.symbolic().to_owned().unwrap()
}

/// Build global unknown-numbering map
pub fn build_unknown_numbering_map(network: &Network) -> Vec<Option<usize>> {
    let mut unknown_id = 0;
    let node_to_unknown: Vec<Option<usize>> = network
        .nodes
        .iter()
        .map(|n| {
            if matches!(n.node_type, NodeType::Junction { .. }) {
                let id = unknown_id;
                unknown_id += 1;
                Some(id)
            } else {
                None
            }
        })
        .collect();
    node_to_unknown
}

/// Map each link to its CSC (Compressed Sparse Column) indices
pub fn map_links_to_csc_indices(
    network: &Network,
    sparsity_pattern: &SymbolicSparseColMat<usize>,
    node_to_unknown: &[Option<usize>],
) -> Vec<CSCIndex> {
    let mut csc_indices = Vec::with_capacity(network.links.len());
    for link in network.links.iter() {
        let mut csc_index = CSCIndex::default();
        let u = node_to_unknown[link.start_node];
        let v = node_to_unknown[link.end_node];

        let sym = sparsity_pattern.as_ref();
        if let Some(i) = u {
            csc_index.diag_u = find_csc_index(sym, i, i);
        }
        if let Some(j) = v {
            csc_index.diag_v = find_csc_index(sym, j, j);
        }
        if let (Some(i), Some(j)) = (u, v) {
            let (row, col) = if i >= j { (i, j) } else { (j, i) };
            csc_index.off_diag = find_csc_index(sym, row, col);
        }
        csc_indices.push(csc_index);
    }
    csc_indices
}

/// Map each node to its row in the Jacobian matrix
pub fn map_nodes_to_rows(
    network: &Network,
    sparsity_pattern: &SymbolicSparseColMat<usize>,
    node_to_unknown: &[Option<usize>],
) -> Vec<Option<usize>> {
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

/// Compute AMD fill-reducing permutation
/// Returns the forward permutation, used to map LLT errors back to the original unknowns
pub fn compute_amd_permutation(
    sparsity_pattern: &SymbolicSparseColMat<usize>,
    node_to_unknown: &[Option<usize>],
) -> Vec<usize> {
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
        )
        .expect("Failed to compute AMD ordering");
    }

    perm_fwd
}
