//! Generic `Link` type, `LinkTrait` and the `LinkType` variant for pipes, pumps and valves.

use crate::constants::Q_ZERO;
use crate::error::InputError;
use crate::model::options::SimulationOptions;
use crate::model::pipe::Pipe;
use crate::model::pump::Pump;
use crate::model::units::UnitConversion;
use crate::model::valve::Valve;
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize, Clone)]
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

    /// Optional vertices for the link (only relevant for graphical display)
    pub vertices: Option<Vec<(f64, f64)>>,

    /// Cached start and end node indices to avoid looking up the node map every time
    #[serde(skip)]
    pub start_node: usize,
    #[serde(skip)]
    pub end_node: usize,
}

#[derive(Deserialize, Serialize, Clone)]
pub enum LinkType {
    Pipe(Pipe),
    Pump(Pump),
    Valve(Valve),
}

// Source: EPANET 2.3 types.h
#[derive(PartialEq, Eq, Debug, Clone, Copy, Deserialize, Serialize)]
#[derive(Default)]
pub enum LinkStatus {
    Xhead,       // pump cannot deliver head (closed)
    TempClosed,  // temporarily closed
    Closed,      // closed
    #[default]
    Open,        // open
    Active,      // valve active (partially open)
    Xflow,       // pump exceeds maximum flow
    XFCV,        // FCV cannot supply flow
    XPressure,   // valve cannot supply pressure
    FixedOpen,   // fixed open
    FixedClosed, // fixed closed
}


impl std::fmt::Display for LinkStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl LinkStatus {
    pub fn from_str(status: &str, is_valve: bool) -> Result<LinkStatus, InputError> {
        match status.to_uppercase().as_str() {
            "CLOSED" => {
                if is_valve {
                    Ok(LinkStatus::FixedClosed)
                } else {
                    Ok(LinkStatus::Closed)
                }
            }
            "OPEN" => {
                if is_valve {
                    Ok(LinkStatus::FixedOpen)
                } else {
                    Ok(LinkStatus::Open)
                }
            }
            "ACTIVE" => Ok(LinkStatus::Active),
            _ => Err(InputError::new(format!("Invalid link status '{}'", status))),
        }
    }
}
pub struct NodeModification {
    pub diagonal_add: f64,
    pub rhs_add: f64,
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
    #[inline]
    pub fn simple(g_inv: f64, y: f64) -> Self {
        Self {
            g_inv,
            y,
            new_status: None,
            upstream_modification: None,
            downstream_modification: None,
        }
    }
    #[inline]
    pub fn new_status(g_inv: f64, y: f64, new_status: LinkStatus) -> Self {
        Self {
            g_inv,
            y,
            new_status: Some(new_status),
            upstream_modification: None,
            downstream_modification: None,
        }
    }
}

pub trait LinkTrait {
    /// Calculate the 1/G_ij and Y_ij coefficients for the link
    fn coefficients(
        &self,
        q: f64,
        resistance: f64,
        setting: f64,
        status: LinkStatus,
        excess_flow_upstream: f64,
        excess_flow_downstream: f64,
    ) -> LinkCoefficients;
    /// Calculate the resistance of the link
    fn resistance(&self) -> f64;
    /// Update the status of the link
    fn update_status(
        &self,
        setting: f64,
        status: LinkStatus,
        flow: f64,
        head_upstream: f64,
        head_downstream: f64,
    ) -> Option<LinkStatus>;
    /// Get the initial flow of the link
    fn initial_flow(&self) -> f64;
}

impl LinkTrait for Link {
    #[inline]
    fn coefficients(
        &self,
        q: f64,
        resistance: f64,
        setting: f64,
        status: LinkStatus,
        excess_flow_upstream: f64,
        excess_flow_downstream: f64,
    ) -> LinkCoefficients {
        match &self.link_type {
            LinkType::Pipe(pipe) => pipe.coefficients(
                q,
                resistance,
                setting,
                status,
                excess_flow_upstream,
                excess_flow_downstream,
            ),
            LinkType::Pump(pump) => pump.coefficients(
                q,
                resistance,
                setting,
                status,
                excess_flow_upstream,
                excess_flow_downstream,
            ),
            LinkType::Valve(valve) => valve.coefficients(
                q,
                resistance,
                setting,
                status,
                excess_flow_upstream,
                excess_flow_downstream,
            ),
        }
    }
    #[inline]
    fn resistance(&self) -> f64 {
        match &self.link_type {
            LinkType::Pipe(pipe) => pipe.resistance(),
            LinkType::Pump(pump) => pump.resistance(),
            LinkType::Valve(valve) => valve.resistance(),
        }
    }
    #[inline]
    fn update_status(
        &self,
        setting: f64,
        status: LinkStatus,
        flow: f64,
        head_upstream: f64,
        head_downstream: f64,
    ) -> Option<LinkStatus> {
        match &self.link_type {
            LinkType::Pipe(pipe) => {
                pipe.update_status(setting, status, flow, head_upstream, head_downstream)
            }
            LinkType::Pump(pump) => {
                pump.update_status(setting, status, flow, head_upstream, head_downstream)
            }
            LinkType::Valve(valve) => {
                valve.update_status(setting, status, flow, head_upstream, head_downstream)
            }
        }
    }
    fn initial_flow(&self) -> f64 {
        // if the link is fixed closed or closed, return 0 flow
        if self.initial_status == LinkStatus::FixedClosed
            || self.initial_status == LinkStatus::Closed
        {
            return Q_ZERO;
        }

        match &self.link_type {
            LinkType::Pipe(pipe) => pipe.initial_flow(),
            LinkType::Pump(pump) => pump.initial_flow(),
            LinkType::Valve(valve) => valve.initial_flow(),
        }
    }
}

impl UnitConversion for Link {
    fn convert_to_standard(&mut self, options: &SimulationOptions) {
        match &mut self.link_type {
            LinkType::Pipe(pipe) => pipe.convert_to_standard(options),
            LinkType::Pump(pump) => pump.convert_to_standard(options),
            LinkType::Valve(valve) => valve.convert_to_standard(options),
        }
    }
    fn convert_from_standard(&mut self, options: &SimulationOptions) {
        match &mut self.link_type {
            LinkType::Pipe(pipe) => pipe.convert_from_standard(options),
            LinkType::Pump(pump) => pump.convert_from_standard(options),
            LinkType::Valve(valve) => valve.convert_from_standard(options),
        }
    }
}

impl Link {
    /// Get the initial setting of the link
    pub fn initial_setting(&self) -> f64 {
        match &self.link_type {
            LinkType::Pipe(_) => 0.0,
            LinkType::Pump(pump) => pump.speed,
            LinkType::Valve(valve) => valve.setting,
        }
    }
}
