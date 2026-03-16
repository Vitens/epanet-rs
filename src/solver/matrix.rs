pub struct ResistanceCoefficients {
  pub g_inv: Vec<f64>,
  pub y: Vec<f64>,
}

impl ResistanceCoefficients {
  pub fn new(size: usize) -> Self {
    Self { g_inv: vec![0.0; size], y: vec![0.0; size] }
  }
}

/// CSC (Compressed Sparse Column) indices for the Jacobian matrix used in the Global Gradient Algorithm
#[derive(Default)]
pub struct CSCIndex {
  pub diag_u: Option<usize>,      // CSC index for J[u,u]
  pub diag_v: Option<usize>,      // CSC index for J[v,v]
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
        .map(|pos| col_start + pos)}