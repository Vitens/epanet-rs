use crate::model::link::{LinkTrait, LinkStatus, LinkCoefficients};
use crate::model::units::{FlowUnits, UnitSystem, UnitConversion};
use crate::constants::*;
use serde::{Deserialize, Serialize};

#[derive(Debug, PartialEq, Eq, Deserialize, Serialize)]
pub enum ValveType {
  PRV, // Pressure Reducing Valve
  PSV, // Pressure Sensing Valve
  PBV, // Pressure Breaking Valve
  FCV, // Flow Control Valve
  TCV, // Throttle Control Valve
  PCV, // Positional Control Valve
  GPV, // General Purpose Valve
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Valve {
  pub diameter: f64,
  pub setting: f64,
  pub curve: Option<Box<str>>,
  pub valve_type: ValveType,
  pub minor_loss: f64,
}


impl LinkTrait for Valve {
  fn coefficients(&self, q: f64, _resistance: f64, status: LinkStatus) -> LinkCoefficients {
    if status == LinkStatus::Closed {
      return LinkCoefficients::simple(1.0/BIG_VALUE, q, status);
    }
    match self.valve_type {
      ValveType::TCV => {
        // Minor loss coefficient is the setting of the valve
        let km = 0.02517 * self.setting / self.diameter.powi(4);
        let (g_inv, y) = self.valve_coefficients(q, km);
        return LinkCoefficients::simple(g_inv, y, status);
      }
      // Positional Control Valve (PCV)
      ValveType::PCV => {
        let km = self.pcv_minor_loss();
        let (g_inv, y) = self.valve_coefficients(q, km);
        return LinkCoefficients::simple(g_inv, y, status);
      }
      // Flow Control Valve (FCV)
      ValveType::FCV => {
        let (g_inv, y) = self.valve_coefficients(q, SMALL_VALUE);
        // if flow is less than the setting, treat as a regular valve (no flow control/restrictions)
        if q < self.setting {
          return LinkCoefficients::simple(g_inv, y, status);
        }
        else {
          let hloss = y / g_inv + BIG_VALUE * (q - self.setting);
          let hgrad = BIG_VALUE;
          return LinkCoefficients::simple(1.0 / hgrad, hloss / hgrad, status);
        }
      }
      // Pressure Reducing Valve (PRV)
      _ => {
        return LinkCoefficients::simple(1.0/SMALL_VALUE, q, status);
      }
    }
  }
  /// Return the resistance of the valve
  fn resistance(&self) -> f64 {
    SMALL_VALUE
  }


}
impl Valve {


  fn pcv_minor_loss(&self) -> f64 {
    // Minor loss coefficient for a completely open valve
    let k_open = 0.02517 * self.minor_loss / self.diameter.powi(4);

    // Valve is completely closed
    if self.setting <= 0.0 {
      return BIG_VALUE;
    }
    // Valve is completely open
    if self.setting >= 100.0 {
      return k_open;
    }
    // Valve is partially open
    let ratio = self.setting / 100.0;

    // clamp the ratio to SMALL_VALUE and 1.0 to avoid division by zero
    let ratio = ratio.clamp(SMALL_VALUE, 1.0);

    // convert the ratio to a minor loss coefficient
    let km = k_open / ratio.powi(2);

    return km.min(BIG_VALUE);
  }
  /// Compute the coefficients for throttle control valve with a minor loss coefficient km and flow q
  fn valve_coefficients(&self, q: f64, km: f64) -> (f64, f64) {

    if km > 0.0 {
      let q_abs = q.abs();
      let hgrad = 2.0 * km * q_abs;

      // guard against too small a head loss gradient
      if hgrad < RQ_TOL {
        let hgrad = RQ_TOL;
        let hloss = q * hgrad;
        return (1.0/hgrad, hloss/hgrad)
      }
      else {
        let hloss = q * hgrad / 2.0;
        return (1.0/hgrad, hloss/hgrad)
      }
    }
    // if no minor loss coefficient, use a low resistance linear head loss relation
    else {
      (1.0/SMALL_VALUE, q)
    }
  }
}

impl UnitConversion for Valve {
  fn convert_units(&mut self, _flow: &FlowUnits, system: &UnitSystem, reverse: bool) {
    if system == &UnitSystem::SI {
      if reverse {
        self.diameter = self.diameter * MperFT * 1e3; // convert in to mm
      }
      else {
        self.diameter = self.diameter / 1e3 / MperFT; // convert mm to in
      }
    } else {
      self.diameter = self.diameter / 12.0; // convert in to ft
    }
  }
}