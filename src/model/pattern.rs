use serde::{Deserialize, Serialize};

/// Pattern struct
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Pattern {
  pub multipliers: Vec<f64>,
}