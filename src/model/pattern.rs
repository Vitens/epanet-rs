use serde::{Deserialize, Serialize};

/// Pattern struct
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Pattern {
  pub id: Box<str>,
  pub multipliers: Vec<f64>,
}