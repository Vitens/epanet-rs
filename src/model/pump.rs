use crate::model::link::LinkTrait;
use crate::model::options::HeadlossFormula;

pub struct Pump {
  pub speed: f64,
  pub head_curve: Box<str>,
  pub power: f64,
}

impl LinkTrait for Pump {
  fn coefficients(&self, _q: f64, _resistance: f64) -> (f64, f64) {
    (0.0, 0.0)
  }
  fn resistance(&self) -> f64 {
    0.0
  }
}