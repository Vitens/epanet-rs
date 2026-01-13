use crate::model::link::LinkTrait;
use crate::model::options::HeadlossFormula;

pub enum ValveType {
  PRV, // Pressure Reducing Valve
  PSV, // Pressure Sensing Valve
  PBV, // Pressure Breaking Valve
  FCV, // Flow Control Valve
  TCV, // Throttle Control Valve
  PCV, // Positional Control Valve
  GPV, // General Purpose Valve
}

pub struct Valve {
  pub diameter: f64,
  pub setting: f64,
  pub curve: Option<Box<str>>,
  pub valve_type: ValveType,
}

impl LinkTrait for Valve {
  fn coefficients(&self, _q: f64, _resistance: f64) -> (f64, f64) {
    (0.0, 0.0)
  }
  fn resistance(&self) -> f64 {
    0.0
  }
}