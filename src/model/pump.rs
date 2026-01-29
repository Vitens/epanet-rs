use crate::model::link::{LinkTrait, LinkStatus, LinkCoefficients};
use crate::model::curve::{HeadCurve};
use crate::model::units::{FlowUnits, UnitSystem, UnitConversion};
use crate::constants::*;
use serde::{Deserialize, Serialize};



#[derive(Debug, Deserialize, Serialize)]
pub struct Pump {
  pub speed: f64,
  pub head_curve_id: Option<Box<str>>,
  pub power: f64,

  // Pump curve for custom curves 
  #[serde(skip)]
  pub head_curve: Option<HeadCurve>,
}

impl LinkTrait for Pump {
  #[inline]
  fn coefficients(&self, q: f64, _resistance: f64, status: LinkStatus, _:f64, _:f64) -> LinkCoefficients {

    // for closed pumps, stalled pumps, or pumps with speed, act as closed pipe
    if status == LinkStatus::Closed || status == LinkStatus::Xhead || status == LinkStatus::TempClosed || self.speed == 0.0 {
      return LinkCoefficients::simple(1.0 / BIG_VALUE, q);
    }

    // get the maximum head from the pump curve or use BIG_VALUE if no curve (constant power pump)
    let h_max = if let Some(curve) = self.head_curve.as_ref() {
      curve.statistics.h_max
    }
    else {
      BIG_VALUE
    };


    // Prevent negative flow
    if q < 0.0 {
      let hloss = -(self.speed.powi(2) * h_max) + BIG_VALUE * q;
      let hgrad = BIG_VALUE;
      return LinkCoefficients::new_status(1.0/hgrad, hloss/hgrad, LinkStatus::Xhead);
    }
    // if no pump curve, treat pump as open valve
    if self.head_curve.is_none() && self.power == 0.0 {
      return LinkCoefficients::simple(1.0 / SMALL_VALUE, q);
    }

    // if no pump curve, and power is non-zero, treat as constant power (HP) pump
    if self.head_curve.is_none() && self.power > 0.0 {

      let r = -8.814 * self.power;
      let hgrad = -r / q.powi(2);
      if hgrad > BIG_VALUE {
        let hloss = -hgrad * q;
        return LinkCoefficients::simple(1.0 / BIG_VALUE, hloss / hgrad);
      }
      else if hgrad < RQ_TOL {
        let hloss = -hgrad * q;
        return LinkCoefficients::simple(1.0 / RQ_TOL, hloss / hgrad);
      }
      else {
        let hloss = r / q;
        return LinkCoefficients::simple(1.0 / hgrad, hloss / hgrad);
      }
    }

   // get the curve coefficients (save to unwrap now)
    let curve = self.head_curve.as_ref().unwrap();

    let q_abs = q.abs();
    // get the curve coefficients
    let (hgrad, hloss) = curve.curve_coefficients(q_abs, self.speed);

    LinkCoefficients::simple(1.0 / hgrad, hloss/hgrad)
  }
  fn resistance(&self) -> f64 {
    BIG_VALUE
  }

  fn update_status(&self, _: LinkStatus, _: f64, _: f64, _: f64) -> Option<LinkStatus> {
    // Reopen the pump if it was temporarily closed
    None
  }

  fn initial_flow(&self) -> f64 {
    if self.speed == 0.0 {
      return Q_ZERO;
    }
    if let Some(head_curve) = &self.head_curve {
      return head_curve.statistics.q_initial;
    }
    else {
      return 1.0; // constant power pump
    }
  }
}

impl UnitConversion for Pump {
  fn convert_units(&mut self, _flow: &FlowUnits, system: &UnitSystem, reverse: bool) {
    // convert the pump units
    if system == &UnitSystem::SI {
      if reverse {
        self.power = self.power * KWperHP;
      }
      else {
        self.power = self.power / KWperHP;
      }
    }
  }
}