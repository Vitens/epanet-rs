//! `Reservoir` node: a fixed-head source with an optional head pattern.

use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Reservoir {
    pub head_pattern: Option<Box<str>>,
    #[serde(skip)]
    pub head_pattern_index: Option<usize>,
}
