//! `gui-core`: UI-framework-agnostic application state for the EPANET
//! network editor. No `egui`/`eframe`/`wgpu` dependency lives here — only
//! `gui-app` may hold or dispatch against a rendering/UI framework; this
//! crate defines the `Action`s it may dispatch and the single `AppState`
//! reducer that applies them, funneling every domain mutation through
//! `epanet_rs::model::network::Network`'s existing add/update/remove API.

pub mod action;
pub mod camera;
pub mod geo;
pub mod selection;
pub mod sim_runner;
pub mod spatial;
pub mod state;

pub use action::{Action, NewLinkKind, Tool};
pub use camera::{Bounds, Camera};
pub use geo::{GeoTransform, ORIGIN_SHIFT, lonlat_to_mercator, tile_bounds_mercator};
pub use selection::Selection;
pub use sim_runner::{SimHandle, SimMessage};
pub use spatial::SpatialIndex;
pub use state::{AppState, LinkSimResult, NodeSimResult};

// Re-export the epanet-rs types the GUI works with most often, so gui-app
// generally only needs `use gui_core::*` plus the specific update-struct
// types it constructs for property edits.
pub use epanet_rs;
