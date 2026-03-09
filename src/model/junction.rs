use crate::model::units::{Cfs, UnitConversion, UnitSystem};
use crate::model::options::SimulationOptions;
use crate::constants::{SMALL_VALUE, RQ_TOL, PSIperFT, MperFT};


use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize)]
pub struct Junction {
  pub basedemand: Cfs,
  pub pattern: Option<Box<str>>,
  pub emitter_coefficient: f64,
}

impl Junction {
  pub fn emitter_coefficients(&self, q: Cfs) -> (f64, f64) {

    let q_exp = 2.0;
    // get the emitter coefficient
    let ke = self.emitter_coefficient.max(SMALL_VALUE);

    // compute the gradient of headloss through emitter
    let hgrad = q_exp * ke * q.abs().powf(q_exp - 1.0);

    if hgrad < RQ_TOL {
      let hgrad = RQ_TOL;
      return (1.0 / hgrad, q * hgrad);
    }

    let y = hgrad * q / q_exp;

    return (1.0 / hgrad, y);

  }
}


impl UnitConversion for Junction {
  fn convert_to_standard(&mut self, options: &SimulationOptions) {
    self.basedemand = self.basedemand / options.flow_units.per_cfs();
    // convert the emitter coefficient to the standard unit system
    if self.emitter_coefficient > 0.0 {
      let pressure_factor = if options.unit_system == UnitSystem::US { PSIperFT } else { MperFT };
      let conversion_factor = options.flow_units.per_cfs().powf(options.emitter_exponent) / pressure_factor;

      self.emitter_coefficient = conversion_factor / self.emitter_coefficient.powf(options.emitter_exponent);
    }

  }
  fn convert_from_standard(&mut self, options: &SimulationOptions) {
    self.basedemand = self.basedemand * options.flow_units.per_cfs();
  }
}