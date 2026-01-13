use crate::model::pipe::Pipe;
use crate::model::pump::Pump;
use crate::model::valve::Valve;

use crate::model::options::HeadlossFormula;

/// Link struct
pub struct Link {
  pub id: Box<str>,
  pub link_type: LinkType,
  pub minor_loss: f64,
  pub start_node: usize,
  pub end_node: usize,
}

pub enum LinkType {
    Pipe(Pipe),
    Pump(Pump),
    Valve(Valve)
}

pub trait LinkTrait {
  /// Calculate the 1/G_ij and Y_ij coefficients for the link
  fn coefficients(&self, q: f64, resistance: f64) -> (f64, f64);
  /// Calculate the resistance of the link
  fn resistance(&self) -> f64;
}

impl LinkTrait for Link {
  fn coefficients(&self, q: f64, resistance: f64) -> (f64, f64) {
    match &self.link_type {
      LinkType::Pipe(pipe) => pipe.coefficients(q, resistance),
      LinkType::Pump(pump) => pump.coefficients(q, resistance),
      LinkType::Valve(valve) => valve.coefficients(q, resistance),
    }
  }
  fn resistance(&self) -> f64 {
    match &self.link_type {
      LinkType::Pipe(pipe) => pipe.resistance(),
      LinkType::Pump(pump) => pump.resistance(),
      LinkType::Valve(valve) => valve.resistance(),
    }
  }
}