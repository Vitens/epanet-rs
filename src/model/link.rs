use crate::model::pipe::Pipe;
use crate::model::pump::Pump;
use crate::model::valve::Valve;
use crate::model::units::{FlowUnits, UnitSystem, UnitConversion};

use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize)]
/// Link struct
pub struct Link {
  /// Link ID
  pub id: Box<str>,
  /// Link type (pipe, pump, valve)
  pub link_type: LinkType,
  /// Start node ID
  pub start_node_id: Box<str>,
  /// End node ID
  pub end_node_id: Box<str>,
  /// Initial status (open, closed, active)
  pub initial_status: LinkStatus,

  /// Cached start and end node indices to avoid looking up the node map every time
  #[serde(skip)]
  pub start_node: usize,
  #[serde(skip)]
  pub end_node: usize,
}

#[derive(Deserialize, Serialize)]
pub enum LinkType {
    Pipe(Pipe),
    Pump(Pump),
    Valve(Valve)
}

// Source: EPANET 2.3 types.h
#[derive(PartialEq, Eq, Debug, Clone, Copy, Deserialize, Serialize)]
pub enum LinkStatus {
  Xhead,         // pump cannot deliver head (closed)
  TempClosed,    // temporarily closed
  Closed,        // closed
  Open,          // open
  Active,        // valve active (partially open)
  Xflow,         // pump exceeds maximum flow
  XFCV,          // FCV cannot supply flow
  XPressure,     // valve cannot supply pressure
  Filling,       // tank filling
  Emptying,      // tank emptying
  Overflowing    // tank overflowing
}

impl LinkStatus {
  pub fn from_str(status: &str) -> LinkStatus {
    match status.to_uppercase().as_str() {
      "CLOSED" => LinkStatus::Closed,
      "OPEN" => LinkStatus::Open,
      "ACTIVE" => LinkStatus::Active,
      _ => panic!("Invalid link status {}", status)
    }
  }
}
pub struct NodeModification {
  pub diagonal_add: f64,
  pub rhs_add: f64
}

pub struct LinkCoefficients {
  pub g_inv: f64,
  pub y: f64,
  /// New status of the link
  pub new_status: Option<LinkStatus>,
  /// Optional modification to upstream node (for PSV valves)
  pub upstream_modification: Option<NodeModification>,
  /// Optional modification to downstream node (for PSV valves)
  pub downstream_modification: Option<NodeModification>,
}

impl LinkCoefficients {
  /// Create a simple link coefficients struct with no matrix modifications (all links except PRV/PSV valves)
  pub fn simple(g_inv: f64, y: f64) -> Self {
    Self { g_inv, y, new_status: None, upstream_modification: None, downstream_modification: None }
  }
  pub fn new_status(g_inv: f64, y: f64, new_status: LinkStatus) -> Self {
    Self { g_inv, y, new_status: Some(new_status), upstream_modification: None, downstream_modification: None }
  }
}

pub trait LinkTrait {
  /// Calculate the 1/G_ij and Y_ij coefficients for the link
  fn coefficients(&self, q: f64, resistance: f64, status: LinkStatus, excess_flow_upstream: f64, excess_flow_downstream: f64) -> LinkCoefficients;
  /// Calculate the resistance of the link
  fn resistance(&self) -> f64;
  /// Update the status of the link
  fn update_status(&self, status: LinkStatus, flow: f64, head_upstream: f64, head_downstream: f64) -> Option<LinkStatus>;
}

impl LinkTrait for Link {
  fn coefficients(&self, q: f64, resistance: f64, status: LinkStatus, excess_flow_upstream: f64, excess_flow_downstream: f64) -> LinkCoefficients {
    match &self.link_type {
      LinkType::Pipe(pipe) => pipe.coefficients(q, resistance, status, excess_flow_upstream, excess_flow_downstream),
      LinkType::Pump(pump) => pump.coefficients(q, resistance, status, excess_flow_upstream, excess_flow_downstream),
      LinkType::Valve(valve) => valve.coefficients(q, resistance, status, excess_flow_upstream, excess_flow_downstream),
    }
  }
  fn resistance(&self) -> f64 {
    match &self.link_type {
      LinkType::Pipe(pipe) => pipe.resistance(),
      LinkType::Pump(pump) => pump.resistance(),
      LinkType::Valve(valve) => valve.resistance(),
    }
  }
  fn update_status(&self, status: LinkStatus, flow: f64, head_upstream: f64, head_downstream: f64) -> Option<LinkStatus> {
    match &self.link_type {
      LinkType::Pipe(pipe) => pipe.update_status(status, flow, head_upstream, head_downstream),
      LinkType::Pump(pump) => pump.update_status(status, flow, head_upstream, head_downstream),
      LinkType::Valve(valve) => valve.update_status(status, flow, head_upstream, head_downstream),
    }
  }
}


impl UnitConversion for Link {
  fn convert_units(&mut self, flow: &FlowUnits, system: &UnitSystem, reverse: bool) {
    match &mut self.link_type {
      LinkType::Pipe(pipe) => pipe.convert_units(flow, system, reverse),
      LinkType::Pump(pump) => pump.convert_units(flow, system, reverse),
      LinkType::Valve(valve) => valve.convert_units(flow, system, reverse),
    }
  }
}