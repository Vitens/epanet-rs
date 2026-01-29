use serde::{Deserialize, Serialize};
use crate::model::link::LinkStatus;

#[derive(Debug, Deserialize, Serialize)]
pub enum ControlCondition {
  Pressure { node_id: Box<str>, above: bool, below: f64 },
  Time { seconds: usize },
  ClockTime { seconds: usize },
}


#[derive(Debug, Deserialize, Serialize)]
pub struct Control {
  pub condition: ControlCondition,
  pub link_id: Box<str>,
  pub setting: Option<f64>,
  pub status: Option<LinkStatus>
}