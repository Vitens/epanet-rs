use crate::model::link::{LinkTrait, LinkStatus, LinkCoefficients, NodeModification};
use crate::model::curve::Curve;
use crate::model::units::{FlowUnits, UnitSystem, UnitConversion};
use crate::constants::*;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

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
  pub curve_id: Option<Box<str>>,
  pub valve_type: ValveType,
  pub minor_loss: f64,
  #[serde(skip)]
  pub valve_curve: Option<Arc<Curve>>
}


impl LinkTrait for Valve {
  fn coefficients(&self, q: f64, _resistance: f64, status: LinkStatus, excess_flow_upstream: f64, excess_flow_downstream: f64) -> LinkCoefficients {

    // if the valve is closed or XPressure, return a high resistance valve
    if status == LinkStatus::Closed || status == LinkStatus::XPressure {
      return LinkCoefficients::simple(1.0/BIG_VALUE, q);
    }

    match self.valve_type {
      // Throttle Control Valve (TCV)
      ValveType::TCV => {
        // Minor loss coefficient is the setting of the valve
        let km = 0.02517 * self.setting / self.diameter.powi(4);
        let (g_inv, y) = self.valve_coefficients(q, km);
        return LinkCoefficients::simple(g_inv, y);
      }
      // Pressure Breaking Valve (PBV)
      ValveType::PBV => {
        return LinkCoefficients::simple(BIG_VALUE, self.setting * BIG_VALUE);
      }
      // Positional Control Valve (PCV)
      ValveType::PCV => {
        let km = self.pcv_minor_loss();
        let (g_inv, y) = self.valve_coefficients(q, km);
        return LinkCoefficients::simple(g_inv, y);
      }
      // Flow Control Valve (FCV)
      ValveType::FCV => {
        let (g_inv, y) = self.valve_coefficients(q, SMALL_VALUE);
        // if flow is less than the setting, treat as a regular valve (no flow control/restrictions)
        if q < self.setting {
          return LinkCoefficients::simple(g_inv, y);
        }
        else {
          let hloss = y / g_inv + BIG_VALUE * (q - self.setting);
          let hgrad = BIG_VALUE;
          return LinkCoefficients::simple(1.0 / hgrad, hloss / hgrad);
        }
      }
      // Pressure Reducing Valve (PRV)
      ValveType::PRV => {
        if status == LinkStatus::Active {
          return self.prv_coefficients(excess_flow_downstream);
        }
        else {
          return LinkCoefficients::simple(1.0/SMALL_VALUE, q);
        }
      }
      ValveType::PSV => {
        return self.psv_coefficients(excess_flow_upstream);
      }
      ValveType::GPV => {
        return self.gpv_coefficients(q);
      }
    }
  }
  /// Return the resistance of the valve
  fn resistance(&self) -> f64 {
    SMALL_VALUE
  }

  fn update_status(&self, status: LinkStatus, q: f64, head_upstream: f64, head_downstream: f64) -> Option<LinkStatus> {
    match self.valve_type {
      ValveType::PRV => {
        self.prv_status(status, q, head_upstream, head_downstream)
      }
      ValveType::PSV => {
        self.psv_status(status, q, head_upstream, head_downstream)
      }
      _ => None,
    }
  }


}
impl Valve {

  /// Compute the updated status of the pressure reducing valve
  fn prv_status(&self, status: LinkStatus, q: f64, head_upstream: f64, head_downstream: f64) -> Option<LinkStatus> {
    match status {
      LinkStatus::Active => {
        if q < -Q_TOL { Some(LinkStatus::Closed) } // closed if flow is negative
        else if head_upstream < self.setting - H_TOL { Some(LinkStatus::Open) } // open if head upstream is less than the setting
        else { None } // no status change
      }
      LinkStatus::Open => {
        if q < -Q_TOL { Some(LinkStatus::Closed) } // closed if flow is negative
        else if head_downstream > self.setting + H_TOL { Some (LinkStatus::Active)} // active if head downstream is greater than the setting
        else { None } // no status change
      }
      LinkStatus::Closed => {
        if head_upstream >= self.setting + H_TOL && head_downstream < self.setting - H_TOL {
          Some(LinkStatus::Active) // active if head upstream is greater than the setting and head downstream is less than the setting
        } else if head_upstream < self.setting - H_TOL && head_upstream > head_downstream + H_TOL {
          Some(LinkStatus::Open) // open if head upstream is less than the setting and head upstream is greater than head downstream plus the tolerance
        } else {
          None // no status change
        }
      }
      LinkStatus::XPressure => {
        if q < -Q_TOL { Some(LinkStatus::Closed)}
        else { None }
      }
      _ => None
    }
  }

  /// Compute the updated status of the pressure sustaining valve
  fn psv_status(&self, status: LinkStatus, q: f64, head_upstream: f64, head_downstream: f64) -> Option<LinkStatus> {
    match status {
      LinkStatus::Active => {
        if q < -Q_TOL { Some(LinkStatus::Closed)} // closed if flow is negative
        else if head_downstream > self.setting + H_TOL { Some(LinkStatus::Open)} //  open if head downstream is less then the setting
        else { None } // no status change
      }
      LinkStatus::Open => {
        if q < -Q_TOL { Some(LinkStatus::Closed)} // closed if flow is negative
        else if head_upstream < self.setting - H_TOL { Some(LinkStatus::Active) } // Active when upstream head is less than the setting
        else { None } // no status change
      }
      LinkStatus::Closed => {
        if head_downstream > self.setting + H_TOL && head_upstream > head_downstream + H_TOL { Some(LinkStatus::Open) } 
        else if head_upstream > self.setting + H_TOL && head_upstream > head_downstream + H_TOL { Some(LinkStatus::Active) }
        else { None } // no status change
      }
      LinkStatus::XPressure => {
        if q < -Q_TOL { Some(LinkStatus::Closed)}
        else { None }
      }
      _ => None
    }
  }


  /// Compute the coefficients for general purpose valve with a flow q
  fn gpv_coefficients(&self, q: f64) -> LinkCoefficients {
    let curve = self.valve_curve.as_ref().unwrap();

    let q_abs = q.abs().max(TINY);

    let (h0, r) = curve.coefficients(q_abs);
    let r = r.max(TINY);

    return LinkCoefficients::simple(1.0/r, (h0 / r + q_abs) * q.signum());
  }

  /// Compute the coefficients for pressure sustaining valve with a flow q and excess flow upstream
  fn psv_coefficients(&self, excess_flow_upstream: f64) -> LinkCoefficients {

    let mut rhs_add = self.setting * BIG_VALUE;

    if excess_flow_upstream > 0.0 {
      rhs_add += excess_flow_upstream;
    }

    return LinkCoefficients {
      g_inv: 0.0,
      y: -excess_flow_upstream,
      new_status: None,
      upstream_modification: Some(NodeModification {
        diagonal_add: BIG_VALUE,
        rhs_add: rhs_add,
      }),
      downstream_modification: None,
    }
  }

  /// Compute the coefficients for pressure reducing valve with a flow q and excess flow downstream
  fn prv_coefficients(&self, excess_flow_downstream: f64) -> LinkCoefficients {

    let mut rhs_add = self.setting * BIG_VALUE;
    if excess_flow_downstream < 0.0 {
      rhs_add += excess_flow_downstream;
    }

    return LinkCoefficients {
      g_inv: 0.0,
      y: excess_flow_downstream,
      new_status: None,
      upstream_modification: None,
      downstream_modification: Some(NodeModification {
        diagonal_add: BIG_VALUE,
        rhs_add: rhs_add,
      })
    }

  }

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
    let ratio = if self.valve_curve.is_some() {
      // Use the valve curve to compute the ratio
      let curve = self.valve_curve.as_ref().unwrap();
      let (h0, r) = curve.coefficients(self.setting);
      let ratio = h0 + r * self.setting;
      ratio / 100.0
    } else {
      // Assume a linear relationship between the setting and the minor loss coefficient
      self.setting / 100.0
    };

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
  fn convert_units(&mut self, flow: &FlowUnits, system: &UnitSystem, reverse: bool) {
    if system == &UnitSystem::SI {
      if reverse {
        self.diameter = self.diameter * MperFT * 1e3; // convert in to mm
      }
      else {
        self.diameter = self.diameter / 1e3 / MperFT; // convert mm to in
      }
    } else {

      self.diameter = self.diameter / 12.0; // convert in to ft

      if self.valve_type == ValveType::PRV || self.valve_type == ValveType::PSV || self.valve_type == ValveType::PBV {
        self.setting = self.setting / PSIperFT; // convert PSI to feet
      }
    }

    // for FCV, convert the setting to CFS
    if self.valve_type == ValveType::FCV {
      self.setting = self.setting * flow.per_cfs();
    }



  }
}