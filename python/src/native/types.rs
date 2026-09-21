//! Read-only snapshot types mirroring the native Rust model structs
//! (`Node`, `Link`, `Pattern`, `Curve`, ...). Values are converted to the
//! network's configured unit system, matching what `add_*`/`update_*` accept.
//!
//! These are plain data snapshots, not live references: per
//! `epanet-rs`'s own contract, mutation always goes through the
//! `Network::update_*`/`add_*`/`remove_*` methods on [`crate::native::network::Network`].

use pyo3::prelude::*;

use epanet::model::control::{Control as RsControl, ControlCondition as RsControlCondition};
use epanet::model::curve::Curve as RsCurve;
use epanet::model::demand::Demand as RsDemand;
use epanet::model::junction::Junction as RsJunction;
use epanet::model::link::{Link as RsLink, LinkType as RsLinkType};
use epanet::model::network::Network as RsNetwork;
use epanet::model::node::{Node as RsNode, NodeType as RsNodeType};
use epanet::model::options::SimulationOptions;
use epanet::model::pattern::Pattern as RsPattern;
use epanet::model::pipe::Pipe as RsPipe;
use epanet::model::pump::Pump as RsPump;
use epanet::model::reservoir::Reservoir as RsReservoir;
use epanet::model::tank::Tank as RsTank;
use epanet::model::units::UnitConversion;
use epanet::model::valve::Valve as RsValve;

use super::enums::{HeadlossFormula, LinkStatus, LinkType, NodeType, ValveType};

/// Mirrors `epanet_rs::model::demand::Demand`.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Demand {
    pub basedemand: f64,
    pub pattern: Option<String>,
    pub name: Option<String>,
}

/// Mirrors `epanet_rs::model::junction::Junction`.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Junction {
    pub emitter_coefficient: f64,
    pub demands: Vec<Demand>,
}

/// Mirrors `epanet_rs::model::tank::Tank`.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Tank {
    pub elevation: f64,
    pub initial_level: f64,
    pub min_level: f64,
    pub max_level: f64,
    pub diameter: f64,
    pub min_volume: f64,
    pub volume_curve_id: Option<String>,
    pub overflow: bool,
}

/// Mirrors `epanet_rs::model::reservoir::Reservoir`.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Reservoir {
    pub head_pattern: Option<String>,
}

/// Mirrors `epanet_rs::model::node::Node`. Exactly one of `junction`, `tank`,
/// `reservoir` is set, matching `node_type`.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Node {
    pub id: String,
    pub node_type: NodeType,
    pub elevation: f64,
    pub coordinates: Option<(f64, f64)>,
    pub junction: Option<Junction>,
    pub tank: Option<Tank>,
    pub reservoir: Option<Reservoir>,
}

/// Mirrors `epanet_rs::model::pipe::Pipe`. `minor_loss` is the dimensionless
/// user-facing coefficient (K), not the internally re-normalized value.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Pipe {
    pub diameter: f64,
    pub length: f64,
    pub roughness: f64,
    pub minor_loss: f64,
    pub check_valve: bool,
    pub headloss_formula: HeadlossFormula,
}

/// Mirrors `epanet_rs::model::pump::Pump`.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Pump {
    pub speed: f64,
    pub power: f64,
    pub head_curve_id: Option<String>,
}

/// Mirrors `epanet_rs::model::valve::Valve`. `minor_loss` is the
/// dimensionless user-facing coefficient (K).
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Valve {
    pub diameter: f64,
    pub setting: f64,
    pub valve_type: ValveType,
    pub minor_loss: f64,
    pub curve_id: Option<String>,
}

/// Mirrors `epanet_rs::model::link::Link`. Exactly one of `pipe`, `pump`,
/// `valve` is set, matching `link_type`.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Link {
    pub id: String,
    pub link_type: LinkType,
    pub start_node_id: String,
    pub end_node_id: String,
    pub initial_status: LinkStatus,
    pub vertices: Option<Vec<(f64, f64)>>,
    pub pipe: Option<Pipe>,
    pub pump: Option<Pump>,
    pub valve: Option<Valve>,
}

/// Mirrors `epanet_rs::model::pattern::Pattern`.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Pattern {
    pub id: String,
    pub multipliers: Vec<f64>,
}

/// Mirrors `epanet_rs::model::curve::Curve`. Curves are stored in the
/// network in whatever units they were entered with (they are not
/// unit-converted like nodes/links), so no conversion is applied here either.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Curve {
    pub id: String,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

/// Mirrors `epanet_rs::model::control::ControlCondition`. Node/tank indices
/// are resolved to ids for stable, human-readable access.
#[pyclass(module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub enum ControlCondition {
    HighPressure { node_id: String, target: f64 },
    LowPressure { node_id: String, target: f64 },
    HighLevel { tank_id: String, target: f64 },
    LowLevel { tank_id: String, target: f64 },
    Time { seconds: usize },
    ClockTime { seconds: usize },
}

/// Mirrors `epanet_rs::model::control::Control`.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct Control {
    pub link_id: String,
    pub condition: ControlCondition,
    pub setting: Option<f64>,
    pub status: Option<LinkStatus>,
}

// --- snapshot conversion helpers (standard units -> user units) ---

pub fn demand_snapshot(demand: &RsDemand) -> Demand {
    Demand {
        basedemand: demand.basedemand,
        pattern: demand.pattern.as_deref().map(String::from),
        name: demand.name.as_deref().map(String::from),
    }
}

pub fn junction_snapshot(junction: &RsJunction) -> Junction {
    Junction {
        emitter_coefficient: junction.emitter_coefficient,
        demands: junction.demands.iter().map(demand_snapshot).collect(),
    }
}

pub fn tank_snapshot(tank: &RsTank) -> Tank {
    Tank {
        elevation: tank.elevation,
        initial_level: tank.initial_level,
        min_level: tank.min_level,
        max_level: tank.max_level,
        diameter: tank.diameter,
        min_volume: tank.min_volume,
        volume_curve_id: tank.volume_curve_id.as_deref().map(String::from),
        overflow: tank.overflow,
    }
}

pub fn reservoir_snapshot(reservoir: &RsReservoir) -> Reservoir {
    Reservoir {
        head_pattern: reservoir.head_pattern.as_deref().map(String::from),
    }
}

/// Convert a node to user units (clones + calls `convert_from_standard`,
/// mirroring the same conversion the C-toolkit-style getters apply
/// per-property).
pub fn node_snapshot(node: &RsNode, options: &SimulationOptions) -> Node {
    let mut converted = node.clone();
    converted.convert_from_standard(options);

    let (node_type, junction, tank, reservoir) = match &converted.node_type {
        RsNodeType::Junction(j) => (NodeType::Junction, Some(junction_snapshot(j)), None, None),
        RsNodeType::Tank(t) => (NodeType::Tank, None, Some(tank_snapshot(t)), None),
        RsNodeType::Reservoir(r) => (NodeType::Reservoir, None, None, Some(reservoir_snapshot(r))),
    };

    Node {
        id: converted.id.to_string(),
        node_type,
        elevation: converted.elevation,
        coordinates: converted.coordinates,
        junction,
        tank,
        reservoir,
    }
}

/// The minor-loss coefficient (K) is dimensionless; it is stored internally
/// re-normalized against the (standard-units) diameter, so it must be
/// recovered from the *pre-conversion* diameter rather than converted
/// directly. See `Network::update_pipe`/`update_valve`.
fn pipe_minor_loss_k(pipe: &RsPipe) -> f64 {
    pipe.minor_loss * pipe.diameter.powi(4) / 0.02517
}

fn valve_minor_loss_k(valve: &RsValve) -> f64 {
    valve.minor_loss * valve.diameter.powi(4) / 0.02517
}

pub fn pipe_snapshot(pipe: &RsPipe, options: &SimulationOptions) -> Pipe {
    let minor_loss = pipe_minor_loss_k(pipe);
    let mut converted = pipe.clone();
    converted.convert_from_standard(options);
    Pipe {
        diameter: converted.diameter,
        length: converted.length,
        roughness: converted.roughness,
        minor_loss,
        check_valve: converted.check_valve,
        headloss_formula: converted.headloss_formula.into(),
    }
}

pub fn pump_snapshot(pump: &RsPump, options: &SimulationOptions) -> Pump {
    let mut converted = pump.clone();
    converted.convert_from_standard(options);
    Pump {
        speed: converted.speed,
        power: converted.power,
        head_curve_id: converted.head_curve_id.as_deref().map(String::from),
    }
}

pub fn valve_snapshot(valve: &RsValve, options: &SimulationOptions) -> Valve {
    let minor_loss = valve_minor_loss_k(valve);
    let mut converted = valve.clone();
    converted.convert_from_standard(options);
    Valve {
        diameter: converted.diameter,
        setting: converted.setting,
        valve_type: converted.valve_type.into(),
        minor_loss,
        curve_id: converted.curve_id.as_deref().map(String::from),
    }
}

pub fn link_snapshot(link: &RsLink, options: &SimulationOptions) -> Link {
    let (link_type, pipe, pump, valve) = match &link.link_type {
        RsLinkType::Pipe(p) => (LinkType::Pipe, Some(pipe_snapshot(p, options)), None, None),
        RsLinkType::Pump(p) => (LinkType::Pump, None, Some(pump_snapshot(p, options)), None),
        RsLinkType::Valve(v) => (
            LinkType::Valve,
            None,
            None,
            Some(valve_snapshot(v, options)),
        ),
    };

    Link {
        id: link.id.to_string(),
        link_type,
        start_node_id: link.start_node_id.to_string(),
        end_node_id: link.end_node_id.to_string(),
        initial_status: link.initial_status.into(),
        vertices: link.vertices.clone(),
        pipe,
        pump,
        valve,
    }
}

pub fn pattern_snapshot(pattern: &RsPattern) -> Pattern {
    Pattern {
        id: pattern.id.to_string(),
        multipliers: pattern.multipliers.clone(),
    }
}

pub fn curve_snapshot(curve: &RsCurve) -> Curve {
    Curve {
        id: curve.id.to_string(),
        x: curve.x.clone(),
        y: curve.y.clone(),
    }
}

pub fn control_snapshot(
    control: &RsControl,
    network: &RsNetwork,
    options: &SimulationOptions,
) -> Control {
    let mut condition = control.condition.clone();
    condition.convert_from_standard(options);

    let py_condition = match condition {
        RsControlCondition::HighPressure { node_index, target } => ControlCondition::HighPressure {
            node_id: network.nodes[node_index].id.to_string(),
            target,
        },
        RsControlCondition::LowPressure { node_index, target } => ControlCondition::LowPressure {
            node_id: network.nodes[node_index].id.to_string(),
            target,
        },
        RsControlCondition::HighLevel { tank_index, target } => ControlCondition::HighLevel {
            tank_id: network.nodes[tank_index].id.to_string(),
            target,
        },
        RsControlCondition::LowLevel { tank_index, target } => ControlCondition::LowLevel {
            tank_id: network.nodes[tank_index].id.to_string(),
            target,
        },
        RsControlCondition::Time { seconds } => ControlCondition::Time { seconds },
        RsControlCondition::ClockTime { seconds } => ControlCondition::ClockTime { seconds },
    };

    let mut setting = None;
    if let Some(value) = control.setting {
        let mut converted = RsControl {
            condition: control.condition.clone(),
            link_id: control.link_id.clone(),
            setting: Some(value),
            status: control.status,
        };
        if let Some(&link_index) = network.link_map.get(&control.link_id) {
            converted.convert_setting_from_standard(&network.links[link_index], options);
        }
        setting = converted.setting;
    }

    Control {
        link_id: control.link_id.to_string(),
        condition: py_condition,
        setting,
        status: control.status.map(Into::into),
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Demand>()?;
    m.add_class::<Junction>()?;
    m.add_class::<Tank>()?;
    m.add_class::<Reservoir>()?;
    m.add_class::<Node>()?;
    m.add_class::<Pipe>()?;
    m.add_class::<Pump>()?;
    m.add_class::<Valve>()?;
    m.add_class::<Link>()?;
    m.add_class::<Pattern>()?;
    m.add_class::<Curve>()?;
    m.add_class::<ControlCondition>()?;
    m.add_class::<Control>()?;
    Ok(())
}
