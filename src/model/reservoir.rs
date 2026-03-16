use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Reservoir {
  pub head_pattern: Option<Box<str>>,
}