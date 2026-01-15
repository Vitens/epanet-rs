use crate::model::link::{LinkTrait, LinkStatus};
use crate::model::curve::HeadCurveStatistics;
use crate::constants::*;

#[derive(Debug)]
pub struct Pump {
  pub speed: f64,
  pub head_curve: Box<str>,
  pub power: f64,
  pub head_curve_statistics: Option<HeadCurveStatistics>,
}

impl LinkTrait for Pump {
  fn coefficients(&self, q: f64, _resistance: f64, status: LinkStatus) -> (f64, f64, LinkStatus) {

    // for closed pumps, stalled pumps, or pumps with speed, act as closed pipe
    if status == LinkStatus::Closed || status == LinkStatus::Xhead || self.speed == 0.0 {
      return (1.0 / BIG_VALUE, q, status);
    }
    let curve = self.head_curve_statistics.as_ref().unwrap();

    // Prevent negative flow
    if q < 0.0 {
      let hloss = -(self.speed.powi(2) * curve.h_max) + BIG_VALUE * q;
      let hgrad = BIG_VALUE;
      return (1.0/hgrad, hloss/hgrad, LinkStatus::Xhead);
    }
    // if no pump curve, treat pump as open valve
    if self.head_curve.is_empty() {
      return (1.0 / SMALL_VALUE, q, status);
    }

    let q_abs = q.abs();
    // if custom curve type
    // TODO: Implement custom curve type
    // TODO: Implement constant HP pump

    // shutoff head is negative to represent head gain
    let h0 = self.speed.powi(2) * -curve.h_shutoff;
    let mut n = curve.n;
    if (curve.n-1.0) < TINY { n = 1.0; }
    let r = curve.r * self.speed.powf(n-1.0);

    // curve is nonlinear
    let (hgrad, hloss) = if n != 1.0 {
      // compute curve gradient
      let hgrad = n * r * q_abs.powf(n - 1.0);
      // ... otherwise compute head loss from pump curve
      let hloss = h0 + hgrad * q/ n;
      (hgrad, hloss)
    }
    // curve is linear
    else {
      let hgrad = r;
      let hloss = h0 + hgrad * q;
      (hgrad, hloss)
    };

    (1.0 / hgrad, hloss/hgrad, status)
  }
  fn resistance(&self) -> f64 {
    BIG_VALUE
  }
}