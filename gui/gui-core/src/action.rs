//! The `Action` enum: the only way `AppState` may be mutated (see `state::AppState::apply`).

use std::path::PathBuf;

use epanet_rs::model::link::LinkStatus;
use epanet_rs::model::network::{
    JunctionUpdate, LinkUpdate, NodeUpdate, PipeData, PipeUpdate, PumpData, PumpUpdate,
    ReservoirData, ReservoirUpdate, TankData, TankUpdate, ValveData, ValveUpdate,
};
use epanet_rs::model::valve::ValveType;

use crate::geo::GeoTransform;

/// The active editing tool. Determines what a canvas click/drag does.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum Tool {
    #[default]
    Select,
    Pan,
    AddJunction,
    AddTank,
    AddReservoir,
    /// Connector tool for pipes/pumps/valves; `kind` picks what gets created
    /// when the drag completes on a second node.
    AddLink {
        kind: NewLinkKind,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NewLinkKind {
    Pipe,
    Pump,
    Valve(ValveType),
}

/// Every possible mutation of `AppState`. UI code only ever constructs these
/// and hands them to `AppState::apply`; it never mutates domain state
/// directly. (Not `Debug`: several variants embed epanet-rs `*Update`
/// structs that don't derive `Debug`.)
pub enum Action {
    // --- selection ---
    SelectNode(String),
    SelectLink(String),
    /// Add to the current multi-selection instead of replacing it.
    AddNodeToSelection(String),
    AddLinkToSelection(String),
    SelectInBox {
        nodes: Vec<String>,
        links: Vec<String>,
    },
    ClearSelection,

    // --- tool / camera (not undo-tracked) ---
    SetTool(Tool),
    PanBy {
        dx: f64,
        dy: f64,
    },
    ZoomAt {
        factor: f64,
        screen_x: f64,
        screen_y: f64,
    },
    SetViewport {
        width: f64,
        height: f64,
    },
    FitToContent,
    /// Fit the camera to the bounds of the current selection only.
    FitToSelection,
    /// Set the camera's zoom factor directly, keeping the current center
    /// fixed (used for "zoom to 100%").
    SetZoom(f64),
    /// Pan so `(x, y)` (world space) becomes the center of the viewport,
    /// without changing zoom (e.g. centering on a single selected node,
    /// where fitting-to-bounds would otherwise zoom all the way in).
    CenterOn {
        x: f64,
        y: f64,
    },

    // --- basemap georeferencing (view/session setting, not undo-tracked:
    // it describes the network's relationship to the earth, not the
    // network itself) ---
    /// Sets (or clears, via `None`) the affine transform used to place a
    /// loaded `.pmtiles` basemap under the network - see `crate::geo`.
    SetGeoTransform(Option<GeoTransform>),

    // --- topology edits (undo snapshot pushed before applying) ---
    AddJunctionAt {
        x: f64,
        y: f64,
    },
    AddTankAt {
        x: f64,
        y: f64,
    },
    AddReservoirAt {
        x: f64,
        y: f64,
    },
    ConnectLink {
        from: String,
        to: String,
        kind: NewLinkKind,
    },
    DeleteSelection,

    // --- move (no snapshot per-frame; caller pushes CommitDrag on mouse-up) ---
    MoveNodeTo {
        id: String,
        x: f64,
        y: f64,
    },
    /// Marks the end of a drag gesture: pushes a single undo snapshot for
    /// everything that happened since the drag started.
    CommitDrag,

    // --- property edits (typed update structs re-used from epanet-rs).
    // `ids` may hold a single id (ordinary edit) or several (bulk edit of a
    // multi-selection); `AppState::apply` applies the same `update` to every
    // id in one undo step either way.
    /// Generic node fields available on every node type (elevation,
    /// coordinates, disabled) — the only way to (de)activate or reposition
    /// a tank/reservoir, since their typed `*Update` structs don't carry
    /// `disabled`.
    UpdateNodes {
        ids: Vec<String>,
        update: NodeUpdate,
    },
    UpdateJunctions {
        ids: Vec<String>,
        update: JunctionUpdate,
    },
    UpdateTanks {
        ids: Vec<String>,
        update: TankUpdate,
    },
    UpdateReservoirs {
        ids: Vec<String>,
        update: ReservoirUpdate,
    },
    /// Generic link fields available on every link type (endpoints,
    /// vertices, initial status).
    UpdateLinks {
        ids: Vec<String>,
        update: LinkUpdate,
    },
    UpdatePipes {
        ids: Vec<String>,
        update: PipeUpdate,
    },
    UpdatePumps {
        ids: Vec<String>,
        update: PumpUpdate,
    },
    UpdateValves {
        ids: Vec<String>,
        update: ValveUpdate,
    },

    // --- simulation ---
    RunSimulation {
        parallel: bool,
    },
    SimulationProgress(f32),
    SimulationFinished(Result<epanet_rs::solver::result::SolverResult, String>),
    SetReportStep(usize),

    // --- undo/redo ---
    Undo,
    Redo,

    // --- import/export ---
    ImportInp(PathBuf),
    ExportInp(PathBuf),

    ClearStatus,
}

/// Data payloads for the add-node/add-link actions, kept out of `Action`
/// itself to avoid an unreadably large enum variant; these mirror the
/// epanet-rs `*Data` structs 1:1 for the fields the GUI currently exposes.
impl NewLinkKind {
    pub fn label(&self) -> &'static str {
        match self {
            NewLinkKind::Pipe => "Pipe",
            NewLinkKind::Pump => "Pump",
            NewLinkKind::Valve(_) => "Valve",
        }
    }
}

pub(crate) fn default_pipe_data(start: &str, end: &str) -> PipeData {
    PipeData {
        start_node: start.into(),
        end_node: end.into(),
        length: 100.0,
        diameter: 12.0,
        roughness: 100.0,
        minor_loss: 0.0,
        check_valve: false,
        initial_status: LinkStatus::Open,
        vertices: None,
    }
}

pub(crate) fn default_pump_data(start: &str, end: &str) -> PumpData {
    PumpData {
        start_node: start.into(),
        end_node: end.into(),
        speed: 1.0,
        head_curve_id: None,
        power: 0.0,
        initial_status: LinkStatus::Open,
        vertices: None,
    }
}

pub(crate) fn default_valve_data(start: &str, end: &str, valve_type: ValveType) -> ValveData {
    ValveData {
        start_node: start.into(),
        end_node: end.into(),
        diameter: 12.0,
        valve_type,
        setting: 0.0,
        curve_id: None,
        minor_loss: 0.0,
        initial_status: LinkStatus::Active,
        vertices: None,
    }
}

pub(crate) fn default_tank_data(x: f64, y: f64) -> TankData {
    TankData {
        elevation: 0.0,
        initial_level: 10.0,
        min_level: 0.0,
        max_level: 20.0,
        diameter: 50.0,
        min_volume: 0.0,
        volume_curve_id: None,
        overflow: false,
        coordinates: Some((x, y)),
    }
}

pub(crate) fn default_reservoir_data(x: f64, y: f64) -> ReservoirData {
    ReservoirData {
        elevation: 0.0,
        head_pattern: None,
        coordinates: Some((x, y)),
    }
}
