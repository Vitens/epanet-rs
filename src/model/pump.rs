use crate::model::link::{LinkTrait, LinkStatus, LinkCoefficients};
use crate::model::curve::{HeadCurve};
use crate::constants::*;
use serde::{Deserialize, Serialize};



#[derive(Debug, Deserialize, Serialize)]
pub struct Pump {
  pub speed: f64,
  pub head_curve_id: Box<str>,
  pub power: f64,

  // Pump curve for custom curves 
  #[serde(skip)]
  pub head_curve: Option<HeadCurve>,
}

impl LinkTrait for Pump {
  fn coefficients(&self, q: f64, _resistance: f64, status: LinkStatus, _:f64, _:f64) -> LinkCoefficients {

    // for closed pumps, stalled pumps, or pumps with speed, act as closed pipe
    if status == LinkStatus::Closed || status == LinkStatus::Xhead || self.speed == 0.0 {
      return LinkCoefficients::simple(1.0 / BIG_VALUE, q);
    }

    let curve = self.head_curve.as_ref().unwrap();

    // Prevent negative flow
    if q < 0.0 {
      let hloss = -(self.speed.powi(2) * curve.statistics.h_max) + BIG_VALUE * q;
      let hgrad = BIG_VALUE;
      return LinkCoefficients::new_status(1.0/hgrad, hloss/hgrad, LinkStatus::Xhead);
    }
    // if no pump curve, treat pump as open valve
    if self.head_curve.is_none() {
      return LinkCoefficients::simple(1.0 / SMALL_VALUE, q);
    }

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
}