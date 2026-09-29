//! `AppState`: the single owner of domain + UI state, and `apply`, the only
//! function allowed to mutate it. Every editing action funnels through the
//! existing `epanet_rs::model::network::Network` add/update/remove API — this
//! reducer never duplicates validation or unit-conversion logic.

use std::path::Path;

use epanet_rs::model::network::{
    JunctionUpdate, LinkUpdate, Network, NodeUpdate, ReservoirUpdate, TankUpdate,
};
use epanet_rs::simulation::Simulation;

use crate::action::{
    Action, NewLinkKind, Tool, default_pipe_data, default_pump_data, default_reservoir_data,
    default_tank_data, default_valve_data,
};
use crate::camera::{Bounds, Camera};
use crate::geo::GeoTransform;
use crate::selection::Selection;
use crate::sim_runner::{SimHandle, SimMessage};
use crate::spatial::SpatialIndex;

/// Maximum number of undo snapshots retained (each is a full `Network`
/// clone; capping bounds memory use on very long editing sessions).
const MAX_UNDO_DEPTH: usize = 100;

pub struct AppState {
    pub simulation: Simulation,
    pub selection: Selection,
    pub tool: Tool,
    pub camera: Camera,
    /// The network -> Web Mercator calibration for the basemap layer, or
    /// `None` if the network hasn't been georeferenced (in which case
    /// `gui-app`'s basemap layer has nowhere to place tile geometry and
    /// draws nothing). Not undo-tracked: it describes how this network
    /// session relates to the earth, not an edit to the network itself.
    pub geo_transform: Option<GeoTransform>,

    pub sim_results: Option<epanet_rs::solver::result::SolverResult>,
    pub sim_running: bool,
    pub sim_progress: Option<f32>,
    pub sim_error: Option<String>,
    pub report_step: usize,

    pub status_message: Option<String>,

    undo_stack: Vec<Network>,
    redo_stack: Vec<Network>,
    /// True while a drag gesture (node move) is in progress; ensures the
    /// whole drag only pushes a single undo snapshot, taken on the first
    /// `MoveNodeTo` of the gesture.
    drag_in_progress: bool,

    spatial: SpatialIndex,
    spatial_dirty: bool,

    sim_handle: Option<SimHandle>,
}

/// Whether `update` touches any field that feeds into the hydraulic solve,
/// as opposed to purely cosmetic/display fields (`coordinates`) that a
/// stale simulation result is still perfectly valid for. Used to decide
/// whether a property edit should invalidate `sim_results`.
fn node_update_affects_hydraulics(update: &NodeUpdate) -> bool {
    update.elevation.is_some() || update.disabled.is_some()
}

fn junction_update_affects_hydraulics(update: &JunctionUpdate) -> bool {
    update.elevation.is_some()
        || update.demands.is_some()
        || update.basedemand.is_some()
        || update.pattern.is_some()
        || update.emitter_coefficient.is_some()
        || update.disabled.is_some()
}

fn tank_update_affects_hydraulics(update: &TankUpdate) -> bool {
    update.elevation.is_some()
        || update.initial_level.is_some()
        || update.min_level.is_some()
        || update.max_level.is_some()
        || update.diameter.is_some()
        || update.min_volume.is_some()
        || update.overflow.is_some()
}

fn reservoir_update_affects_hydraulics(update: &ReservoirUpdate) -> bool {
    update.elevation.is_some() || update.head_pattern.is_some()
}

/// `vertices` are a display-only polyline and don't affect the solve;
/// `start_node`/`end_node`/`initial_status` all do.
fn link_update_affects_hydraulics(update: &LinkUpdate) -> bool {
    update.start_node.is_some() || update.end_node.is_some() || update.initial_status.is_some()
}

impl AppState {
    pub fn new(network: Network) -> Self {
        let simulation = Simulation::new(network);
        let spatial = SpatialIndex::build(&simulation.network, 100.0);
        Self {
            simulation,
            selection: Selection::default(),
            tool: Tool::default(),
            camera: Camera::default(),
            geo_transform: None,
            sim_results: None,
            sim_running: false,
            sim_progress: None,
            sim_error: None,
            report_step: 0,
            status_message: None,
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
            drag_in_progress: false,
            spatial,
            spatial_dirty: false,
            sim_handle: None,
        }
    }

    pub fn network(&self) -> &Network {
        &self.simulation.network
    }

    /// Rebuild the spatial index if anything since the last call may have
    /// invalidated it (topology change, coordinate change, undo/redo,
    /// import). Cheap no-op otherwise. Call once per frame before hit
    /// testing / culling.
    pub fn ensure_spatial_index(&mut self) -> &SpatialIndex {
        if self.spatial_dirty {
            self.spatial = SpatialIndex::build(&self.simulation.network, 100.0);
            self.spatial_dirty = false;
        }
        &self.spatial
    }

    pub fn spatial_index(&self) -> &SpatialIndex {
        &self.spatial
    }

    pub fn can_undo(&self) -> bool {
        !self.undo_stack.is_empty()
    }
    pub fn can_redo(&self) -> bool {
        !self.redo_stack.is_empty()
    }

    /// Poll the background simulation (if any) and turn its messages into
    /// actions. Call once per frame from the UI event loop; apply the
    /// returned actions the same way any other action is applied.
    pub fn drain_simulation_messages(&mut self) -> Vec<Action> {
        let Some(handle) = &self.sim_handle else {
            return Vec::new();
        };
        handle
            .poll()
            .into_iter()
            .map(|msg| match msg {
                SimMessage::Progress(p) => Action::SimulationProgress(p),
                SimMessage::Finished(r) => Action::SimulationFinished(r),
            })
            .collect()
    }

    fn push_undo_snapshot(&mut self) {
        if self.undo_stack.len() >= MAX_UNDO_DEPTH {
            self.undo_stack.remove(0);
        }
        self.undo_stack.push(self.simulation.network.clone());
        self.redo_stack.clear();
    }

    fn mark_spatial_dirty(&mut self) {
        self.spatial_dirty = true;
    }

    fn next_id(&self, prefix: &str) -> String {
        let mut n = self.simulation.network.nodes.len() + self.simulation.network.links.len() + 1;
        loop {
            let candidate = format!("{prefix}{n}");
            if !self
                .simulation
                .network
                .node_map
                .contains_key(candidate.as_str())
                && !self
                    .simulation
                    .network
                    .link_map
                    .contains_key(candidate.as_str())
            {
                return candidate;
            }
            n += 1;
        }
    }

    fn set_error(&mut self, message: impl std::fmt::Display) {
        self.status_message = Some(message.to_string());
    }

    /// Apply `f` (one of `Network::update_*`) to every id in `ids`, as part
    /// of a single undo step the caller has already pushed via
    /// `push_undo_snapshot`. Used for both ordinary edits (`ids.len() ==
    /// 1`) and bulk edits of a multi-selection. If every id fails (e.g. the
    /// whole selection turned out to be the wrong entity kind for this
    /// update), the just-pushed snapshot is discarded since nothing
    /// actually changed; partial failures in a mixed-type bulk selection are
    /// reported but don't roll back the ids that did succeed. Returns
    /// whether at least one id succeeded, so callers can decide whether the
    /// edit actually changed anything (e.g. to invalidate stale simulation
    /// results only when it did).
    fn apply_to_each<F, E>(&mut self, ids: &[String], mut f: F) -> bool
    where
        F: FnMut(&mut Network, &str) -> Result<(), E>,
        E: std::fmt::Display,
    {
        let mut succeeded = 0usize;
        let mut last_error: Option<String> = None;
        for id in ids {
            match f(&mut self.simulation.network, id) {
                Ok(()) => succeeded += 1,
                Err(e) => last_error = Some(e.to_string()),
            }
        }
        let any_succeeded = succeeded > 0;
        if !any_succeeded {
            self.undo_stack.pop();
        } else {
            self.mark_spatial_dirty();
        }
        if let Some(e) = last_error {
            self.set_error(e);
        }
        any_succeeded
    }

    /// Discard any simulation results (and related bookkeeping) once the
    /// network's topology or hydraulic properties change under it, so a
    /// stale result is never shown as if it still applied to the edited
    /// network. Not called for pure display edits (moving a node,
    /// editing link vertices) since those don't affect the hydraulic model.
    fn invalidate_sim_results(&mut self) {
        self.sim_results = None;
        self.sim_progress = None;
        self.sim_error = None;
        self.report_step = 0;
    }

    /// The single mutation point for `AppState`. UI code (or replayed
    /// background-thread messages) only ever constructs `Action`s and hands
    /// them here.
    pub fn apply(&mut self, action: Action) {
        match action {
            // --- selection ---
            Action::SelectNode(id) => self.selection = Selection::single_node(id),
            Action::SelectLink(id) => self.selection = Selection::single_link(id),
            Action::AddNodeToSelection(id) => {
                if !self.selection.contains_node(&id) {
                    self.selection.nodes.push(id);
                }
            }
            Action::AddLinkToSelection(id) => {
                if !self.selection.contains_link(&id) {
                    self.selection.links.push(id);
                }
            }
            Action::SelectInBox { nodes, links } => {
                self.selection = Selection { nodes, links };
            }
            Action::ClearSelection => self.selection.clear(),

            // --- tool / camera ---
            Action::SetTool(tool) => self.tool = tool,
            Action::PanBy { dx, dy } => self.camera.pan_screen(dx, dy),
            Action::ZoomAt {
                factor,
                screen_x,
                screen_y,
            } => {
                self.camera.zoom_at(factor, (screen_x, screen_y));
            }
            Action::SetViewport { width, height } => self.camera.viewport = (width, height),
            Action::FitToContent => {
                if let Some(bounds) = self.content_bounds() {
                    self.camera.fit_to_bounds(bounds);
                }
            }
            Action::FitToSelection => {
                if let Some(bounds) = self.selection_bounds() {
                    self.camera.fit_to_bounds(bounds);
                }
            }
            Action::SetZoom(zoom) => {
                self.camera.zoom = zoom.clamp(1e-6, 1e6);
            }
            Action::CenterOn { x, y } => {
                self.camera.center = (x, y);
            }

            Action::SetGeoTransform(transform) => self.geo_transform = transform,

            // --- topology edits ---
            Action::AddJunctionAt { x, y } => {
                let id = self.next_id("J");
                self.push_undo_snapshot();
                let data = epanet_rs::model::network::JunctionData {
                    elevation: 0.0,
                    coordinates: Some((x, y)),
                    ..Default::default()
                };
                match self.simulation.network.add_junction(&id, &data) {
                    Ok(()) => {
                        self.selection = Selection::single_node(id);
                        self.mark_spatial_dirty();
                        self.invalidate_sim_results();
                    }
                    Err(e) => {
                        self.undo_stack.pop();
                        self.set_error(e);
                    }
                }
            }
            Action::AddTankAt { x, y } => {
                let id = self.next_id("T");
                self.push_undo_snapshot();
                match self
                    .simulation
                    .network
                    .add_tank(&id, &default_tank_data(x, y))
                {
                    Ok(()) => {
                        self.selection = Selection::single_node(id);
                        self.mark_spatial_dirty();
                        self.invalidate_sim_results();
                    }
                    Err(e) => {
                        self.undo_stack.pop();
                        self.set_error(e);
                    }
                }
            }
            Action::AddReservoirAt { x, y } => {
                let id = self.next_id("R");
                self.push_undo_snapshot();
                match self
                    .simulation
                    .network
                    .add_reservoir(&id, &default_reservoir_data(x, y))
                {
                    Ok(()) => {
                        self.selection = Selection::single_node(id);
                        self.mark_spatial_dirty();
                        self.invalidate_sim_results();
                    }
                    Err(e) => {
                        self.undo_stack.pop();
                        self.set_error(e);
                    }
                }
            }
            Action::ConnectLink { from, to, kind } => {
                if from == to {
                    self.set_error("Cannot connect a node to itself");
                    return;
                }
                let prefix = match &kind {
                    NewLinkKind::Pipe => "P",
                    NewLinkKind::Pump => "PU",
                    NewLinkKind::Valve(_) => "V",
                };
                let id = self.next_id(prefix);
                self.push_undo_snapshot();
                let result = match &kind {
                    NewLinkKind::Pipe => self
                        .simulation
                        .network
                        .add_pipe(&id, &default_pipe_data(&from, &to)),
                    NewLinkKind::Pump => self
                        .simulation
                        .network
                        .add_pump(&id, &default_pump_data(&from, &to)),
                    NewLinkKind::Valve(vt) => self
                        .simulation
                        .network
                        .add_valve(&id, &default_valve_data(&from, &to, vt.clone())),
                };
                match result {
                    Ok(()) => {
                        self.selection = Selection::single_link(id);
                        self.mark_spatial_dirty();
                        self.invalidate_sim_results();
                    }
                    Err(e) => {
                        self.undo_stack.pop();
                        self.set_error(e);
                    }
                }
            }
            Action::DeleteSelection => {
                if self.selection.is_empty() {
                    return;
                }
                self.push_undo_snapshot();
                for link_id in self.selection.links.clone() {
                    let _ = self.simulation.network.remove_link(&link_id, true);
                }
                for node_id in self.selection.nodes.clone() {
                    let _ = self.simulation.network.remove_node(&node_id, true);
                }
                self.selection.clear();
                self.mark_spatial_dirty();
                self.invalidate_sim_results();
            }

            // --- move ---
            Action::MoveNodeTo { id, x, y } => {
                if !self.drag_in_progress {
                    self.push_undo_snapshot();
                    self.drag_in_progress = true;
                }
                let _ = self.simulation.network.update_node(
                    &id,
                    &NodeUpdate {
                        coordinates: Some((x, y)),
                        ..Default::default()
                    },
                );
                self.mark_spatial_dirty();
            }
            Action::CommitDrag => self.drag_in_progress = false,

            // --- property edits ---
            // Each of these applies `update` to every id in `ids` (one id for
            // an ordinary edit, several for a bulk edit of a multi-selection)
            // as a single undo step. Individual failures (e.g. a `Tank`-only
            // field sent to a `Junction` in a mixed bulk selection) are
            // reported but don't stop the rest of the batch. Any simulation
            // result is invalidated once something that actually feeds the
            // solve changed (not for a coordinates/vertices-only edit).
            Action::UpdateNodes { ids, update } => {
                self.push_undo_snapshot();
                let hydraulic = node_update_affects_hydraulics(&update);
                let changed = self.apply_to_each(&ids, |network, id| network.update_node(id, &update));
                if changed && hydraulic {
                    self.invalidate_sim_results();
                }
            }
            Action::UpdateJunctions { ids, update } => {
                self.push_undo_snapshot();
                let hydraulic = junction_update_affects_hydraulics(&update);
                let changed = self.apply_to_each(&ids, |network, id| network.update_junction(id, &update));
                if changed && hydraulic {
                    self.invalidate_sim_results();
                }
            }
            Action::UpdateTanks { ids, update } => {
                self.push_undo_snapshot();
                let hydraulic = tank_update_affects_hydraulics(&update);
                let changed = self.apply_to_each(&ids, |network, id| network.update_tank(id, &update));
                if changed && hydraulic {
                    self.invalidate_sim_results();
                }
            }
            Action::UpdateReservoirs { ids, update } => {
                self.push_undo_snapshot();
                let hydraulic = reservoir_update_affects_hydraulics(&update);
                let changed = self.apply_to_each(&ids, |network, id| network.update_reservoir(id, &update));
                if changed && hydraulic {
                    self.invalidate_sim_results();
                }
            }
            Action::UpdateLinks { ids, update } => {
                self.push_undo_snapshot();
                let hydraulic = link_update_affects_hydraulics(&update);
                let changed = self.apply_to_each(&ids, |network, id| network.update_link(id, &update));
                if changed && hydraulic {
                    self.invalidate_sim_results();
                }
            }
            Action::UpdatePipes { ids, update } => {
                self.push_undo_snapshot();
                // Every `PipeUpdate` field (length/diameter/roughness/minor
                // loss/check valve) feeds the solve, unlike the generic node
                // and link updates above.
                if self.apply_to_each(&ids, |network, id| network.update_pipe(id, &update)) {
                    self.invalidate_sim_results();
                }
            }
            Action::UpdatePumps { ids, update } => {
                self.push_undo_snapshot();
                if self.apply_to_each(&ids, |network, id| network.update_pump(id, &update)) {
                    self.invalidate_sim_results();
                }
            }
            Action::UpdateValves { ids, update } => {
                self.push_undo_snapshot();
                if self.apply_to_each(&ids, |network, id| network.update_valve(id, &update)) {
                    self.invalidate_sim_results();
                }
            }

            // --- simulation ---
            Action::RunSimulation { parallel } => {
                if self.sim_running {
                    return;
                }
                self.sim_running = true;
                self.sim_progress = Some(0.0);
                self.sim_error = None;
                let network_snapshot = self.simulation.network.clone();
                self.sim_handle = Some(SimHandle::spawn(network_snapshot, parallel));
            }
            Action::SimulationProgress(p) => self.sim_progress = Some(p),
            Action::SimulationFinished(result) => {
                self.sim_running = false;
                self.sim_handle = None;
                match result {
                    Ok(results) => {
                        self.sim_results = Some(results);
                        self.report_step = 0;
                        self.sim_progress = Some(1.0);
                    }
                    Err(e) => {
                        self.sim_error = Some(e);
                        self.sim_progress = None;
                    }
                }
            }
            Action::SetReportStep(step) => {
                let max_step = self
                    .sim_results
                    .as_ref()
                    .map(|r| r.heads.len().saturating_sub(1))
                    .unwrap_or(0);
                self.report_step = step.min(max_step);
            }

            // --- undo/redo ---
            Action::Undo => {
                if let Some(prev) = self.undo_stack.pop() {
                    let current = std::mem::replace(&mut self.simulation.network, prev);
                    self.redo_stack.push(current);
                    self.after_network_replaced();
                }
            }
            Action::Redo => {
                if let Some(next) = self.redo_stack.pop() {
                    let current = std::mem::replace(&mut self.simulation.network, next);
                    self.undo_stack.push(current);
                    self.after_network_replaced();
                }
            }

            // --- import/export ---
            Action::ImportInp(path) => match load_network(&path) {
                Ok(network) => {
                    self.simulation = Simulation::new(network);
                    self.undo_stack.clear();
                    self.redo_stack.clear();
                    self.selection.clear();
                    self.invalidate_sim_results();
                    // A freshly imported network's `[COORDINATES]` are
                    // generally a different local/national grid than
                    // whatever the previous network was calibrated
                    // against, so any existing basemap alignment is
                    // presumed stale rather than silently kept.
                    self.geo_transform = None;
                    self.mark_spatial_dirty();
                    self.ensure_spatial_index();
                    if let Some(bounds) = self.content_bounds() {
                        self.camera.fit_to_bounds(bounds);
                    }
                    self.status_message = Some(format!("Imported {}", path.display()));
                }
                Err(e) => self.set_error(e),
            },
            Action::ExportInp(path) => {
                let Some(path_str) = path.to_str() else {
                    self.set_error("Export path is not valid UTF-8");
                    return;
                };
                match self.simulation.network.save_network(path_str) {
                    Ok(()) => self.status_message = Some(format!("Exported {}", path.display())),
                    Err(e) => self.set_error(e),
                }
            }

            Action::ClearStatus => self.status_message = None,
        }
    }

    /// Called after `Undo`/`Redo` swap in a whole different `Network`
    /// snapshot. Any simulation result is invalidated unconditionally here:
    /// unlike a single tracked edit, we have no cheap way to tell whether
    /// the swapped-in network differs only cosmetically (e.g. undoing a
    /// pure node move) from the one the result was computed for.
    fn after_network_replaced(&mut self) {
        self.mark_spatial_dirty();
        self.ensure_spatial_index();
        self.selection
            .nodes
            .retain(|id| self.simulation.network.node_map.contains_key(id.as_str()));
        self.selection
            .links
            .retain(|id| self.simulation.network.link_map.contains_key(id.as_str()));
        self.drag_in_progress = false;
        self.invalidate_sim_results();
    }

    /// World-space bounds of every node with coordinates, used by
    /// `FitToContent` and on import.
    pub fn content_bounds(&self) -> Option<Bounds> {
        Bounds::from_points(
            self.simulation
                .network
                .nodes
                .iter()
                .filter_map(|n| n.coordinates),
        )
    }

    /// World-space bounds of the current selection (selected nodes' points,
    /// plus the endpoints and vertices of selected links), used by
    /// `FitToSelection`. `None` if nothing in the selection has coordinates.
    pub fn selection_bounds(&self) -> Option<Bounds> {
        let network = &self.simulation.network;
        let mut points: Vec<(f64, f64)> = Vec::new();

        for id in &self.selection.nodes {
            if let Some(&idx) = network.node_map.get(id.as_str())
                && let Some(p) = network.nodes[idx].coordinates
            {
                points.push(p);
            }
        }

        for id in &self.selection.links {
            if let Some(&idx) = network.link_map.get(id.as_str()) {
                let link = &network.links[idx];
                if let Some(&si) = network.node_map.get(link.start_node_id.as_ref())
                    && let Some(p) = network.nodes[si].coordinates
                {
                    points.push(p);
                }
                if let Some(&ei) = network.node_map.get(link.end_node_id.as_ref())
                    && let Some(p) = network.nodes[ei].coordinates
                {
                    points.push(p);
                }
                if let Some(vertices) = &link.vertices {
                    points.extend(vertices.iter().copied());
                }
            }
        }

        Bounds::from_points(points)
    }

    /// Simulation results for node `id` at the current report step, in the
    /// network's configured units (matches what `node.elevation` etc. are
    /// displayed in elsewhere). `None` if there's no simulation result, `id`
    /// isn't a node, or the result no longer matches the live network's node
    /// count (topology changed since the run).
    pub fn node_result(&self, id: &str) -> Option<NodeSimResult> {
        let results = self.sim_results.as_ref()?;
        let network = &self.simulation.network;
        let &index = network.node_map.get(id)?;
        let max_step = results.heads.len().checked_sub(1)?;
        let step = self.report_step.min(max_step);
        let heads = results.heads.get(step)?;
        let demands = results.demands.get(step)?;
        if heads.len() != network.nodes.len() || demands.len() != network.nodes.len() {
            return None;
        }
        let head = *heads.get(index)?;
        let demand = *demands.get(index)?;
        let per_feet = network.options.unit_system.per_feet();
        let pressure = head - network.nodes[index].elevation * per_feet;
        Some(NodeSimResult { head, pressure, demand })
    }

    /// Simulation results for link `id` at the current report step, in the
    /// network's configured flow units (positive = flowing start->end).
    /// `None` under the same conditions as `node_result`.
    pub fn link_result(&self, id: &str) -> Option<LinkSimResult> {
        let results = self.sim_results.as_ref()?;
        let network = &self.simulation.network;
        let &index = network.link_map.get(id)?;
        let max_step = results.heads.len().checked_sub(1)?;
        let step = self.report_step.min(max_step);
        let flows = results.flows.get(step)?;
        if flows.len() != network.links.len() {
            return None;
        }
        let flow = *flows.get(index)?;
        Some(LinkSimResult { flow })
    }
}

/// Simulation results for a single node at the app's current report step;
/// see `AppState::node_result`.
pub struct NodeSimResult {
    pub head: f64,
    pub pressure: f64,
    pub demand: f64,
}

/// Simulation results for a single link at the app's current report step;
/// see `AppState::link_result`.
pub struct LinkSimResult {
    pub flow: f64,
}

fn load_network(path: &Path) -> Result<Network, String> {
    let path_str = path.to_str().ok_or("Import path is not valid UTF-8")?;
    Network::from_file(path_str).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use epanet_rs::model::network::JunctionUpdate;
    use epanet_rs::model::options::HeadlossFormula;
    use epanet_rs::model::units::FlowUnits;
    use epanet_rs::solver::result::SolverResult;

    fn empty_state() -> AppState {
        AppState::new(Network::new(FlowUnits::CFS, HeadlossFormula::HazenWilliams))
    }

    /// A trivially-shaped "sim just ran" result, for tests that only care
    /// about whether `sim_results` gets invalidated, not its contents.
    fn fake_result(n_nodes: usize, n_links: usize) -> SolverResult {
        SolverResult::new(n_links, n_nodes, 1)
    }

    #[test]
    fn add_junction_creates_node_and_selects_it() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 1.0, y: 2.0 });
        assert_eq!(state.network().nodes.len(), 1);
        let id = state.selection.as_single_node().unwrap().to_string();
        assert_eq!(state.network().nodes[0].coordinates, Some((1.0, 2.0)));
        assert!(state.can_undo());
        assert_eq!(state.network().node_map[id.as_str()], 0);
    }

    #[test]
    fn connect_link_adds_pipe_between_nodes() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        let a = state.selection.as_single_node().unwrap().to_string();
        state.apply(Action::AddJunctionAt { x: 10.0, y: 0.0 });
        let b = state.selection.as_single_node().unwrap().to_string();

        state.apply(Action::ConnectLink {
            from: a.clone(),
            to: b.clone(),
            kind: NewLinkKind::Pipe,
        });

        assert_eq!(state.network().links.len(), 1);
        let link = &state.network().links[0];
        assert_eq!(link.start_node_id.as_ref(), a.as_str());
        assert_eq!(link.end_node_id.as_ref(), b.as_str());
    }

    #[test]
    fn delete_selection_removes_node_and_its_links() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        let a = state.selection.as_single_node().unwrap().to_string();
        state.apply(Action::AddJunctionAt { x: 10.0, y: 0.0 });
        let b = state.selection.as_single_node().unwrap().to_string();
        state.apply(Action::ConnectLink {
            from: a.clone(),
            to: b,
            kind: NewLinkKind::Pipe,
        });

        state.apply(Action::SelectNode(a));
        state.apply(Action::DeleteSelection);

        assert_eq!(state.network().nodes.len(), 1);
        assert_eq!(state.network().links.len(), 0);
        assert!(state.selection.is_empty());
    }

    #[test]
    fn update_junction_property_goes_through_network_api() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        let id = state.selection.as_single_node().unwrap().to_string();

        state.apply(Action::UpdateJunctions {
            ids: vec![id.clone()],
            update: JunctionUpdate {
                elevation: Some(42.0),
                ..Default::default()
            },
        });

        let idx = state.network().node_map[id.as_str()];
        assert!((state.network().nodes[idx].elevation - 42.0).abs() < 1e-9);
    }

    #[test]
    fn bulk_update_applies_to_every_selected_id_in_one_undo_step() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        let a = state.selection.as_single_node().unwrap().to_string();
        state.apply(Action::AddJunctionAt { x: 1.0, y: 0.0 });
        let b = state.selection.as_single_node().unwrap().to_string();
        let undo_depth_before = state.undo_stack.len();

        state.apply(Action::UpdateJunctions {
            ids: vec![a.clone(), b.clone()],
            update: JunctionUpdate {
                elevation: Some(7.5),
                ..Default::default()
            },
        });

        assert_eq!(state.undo_stack.len(), undo_depth_before + 1);
        for id in [&a, &b] {
            let idx = state.network().node_map[id.as_str()];
            assert!((state.network().nodes[idx].elevation - 7.5).abs() < 1e-9);
        }

        state.apply(Action::Undo);
        for id in [&a, &b] {
            let idx = state.network().node_map[id.as_str()];
            assert!((state.network().nodes[idx].elevation - 0.0).abs() < 1e-9);
        }
    }

    #[test]
    fn bulk_update_on_generic_node_fields_works_across_mixed_node_types() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        let junction = state.selection.as_single_node().unwrap().to_string();
        state.apply(Action::AddReservoirAt { x: 5.0, y: 0.0 });
        let reservoir = state.selection.as_single_node().unwrap().to_string();

        // `disabled` isn't on `JunctionUpdate`/`ReservoirUpdate`, only the
        // generic `NodeUpdate`, so this exercises the one bulk action that
        // must work across heterogeneous node types.
        state.apply(Action::UpdateNodes {
            ids: vec![junction.clone(), reservoir.clone()],
            update: epanet_rs::model::network::NodeUpdate {
                disabled: Some(true),
                ..Default::default()
            },
        });

        for id in [&junction, &reservoir] {
            let idx = state.network().node_map[id.as_str()];
            assert!(state.network().nodes[idx].disabled);
        }
    }

    #[test]
    fn undo_redo_restores_and_replays_network_state() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        assert_eq!(state.network().nodes.len(), 1);

        state.apply(Action::Undo);
        assert_eq!(state.network().nodes.len(), 0);
        assert!(state.can_redo());

        state.apply(Action::Redo);
        assert_eq!(state.network().nodes.len(), 1);
    }

    #[test]
    fn drag_gesture_pushes_a_single_undo_snapshot() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        let id = state.selection.as_single_node().unwrap().to_string();
        let undo_depth_before_drag = state.undo_stack.len();

        state.apply(Action::MoveNodeTo {
            id: id.clone(),
            x: 1.0,
            y: 1.0,
        });
        state.apply(Action::MoveNodeTo {
            id: id.clone(),
            x: 2.0,
            y: 2.0,
        });
        state.apply(Action::MoveNodeTo {
            id: id.clone(),
            x: 3.0,
            y: 3.0,
        });
        state.apply(Action::CommitDrag);

        assert_eq!(state.undo_stack.len(), undo_depth_before_drag + 1);
        let idx = state.network().node_map[id.as_str()];
        assert_eq!(state.network().nodes[idx].coordinates, Some((3.0, 3.0)));

        state.apply(Action::Undo);
        let idx = state.network().node_map[id.as_str()];
        assert_eq!(state.network().nodes[idx].coordinates, Some((0.0, 0.0)));
    }

    #[test]
    fn selection_helpers_add_and_clear() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        let a = state.selection.as_single_node().unwrap().to_string();
        state.apply(Action::AddJunctionAt { x: 1.0, y: 1.0 });
        let b = state.selection.as_single_node().unwrap().to_string();

        state.apply(Action::AddNodeToSelection(a.clone()));
        assert!(state.selection.contains_node(&a));
        assert!(state.selection.contains_node(&b));

        state.apply(Action::ClearSelection);
        assert!(state.selection.is_empty());
    }

    #[test]
    fn property_edit_invalidates_stale_sim_results() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        let id = state.selection.as_single_node().unwrap().to_string();
        state.sim_results = Some(fake_result(1, 0));
        state.report_step = 0;

        state.apply(Action::UpdateJunctions {
            ids: vec![id],
            update: JunctionUpdate {
                elevation: Some(5.0),
                ..Default::default()
            },
        });

        assert!(state.sim_results.is_none(), "changing elevation should invalidate results");
    }

    #[test]
    fn coordinates_only_edit_does_not_invalidate_sim_results() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        let id = state.selection.as_single_node().unwrap().to_string();
        state.sim_results = Some(fake_result(1, 0));

        // A pure position edit (dragging a node, or the position field in
        // the properties panel) is display-only and shouldn't touch a
        // still-valid simulation result.
        state.apply(Action::MoveNodeTo { id: id.clone(), x: 5.0, y: 5.0 });
        state.apply(Action::CommitDrag);
        assert!(state.sim_results.is_some(), "moving a node shouldn't invalidate results");

        state.apply(Action::UpdateNodes {
            ids: vec![id],
            update: NodeUpdate {
                coordinates: Some((9.0, 9.0)),
                ..Default::default()
            },
        });
        assert!(
            state.sim_results.is_some(),
            "a coordinates-only UpdateNodes shouldn't invalidate results"
        );
    }

    #[test]
    fn disabling_a_node_invalidates_sim_results() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        let id = state.selection.as_single_node().unwrap().to_string();
        state.sim_results = Some(fake_result(1, 0));

        state.apply(Action::UpdateNodes {
            ids: vec![id],
            update: NodeUpdate {
                disabled: Some(true),
                ..Default::default()
            },
        });

        assert!(state.sim_results.is_none(), "disabling a node changes the topology it solves over");
    }

    #[test]
    fn topology_edits_invalidate_sim_results() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        let a = state.selection.as_single_node().unwrap().to_string();
        state.apply(Action::AddJunctionAt { x: 10.0, y: 0.0 });
        let b = state.selection.as_single_node().unwrap().to_string();
        state.sim_results = Some(fake_result(2, 0));

        state.apply(Action::ConnectLink {
            from: a.clone(),
            to: b,
            kind: NewLinkKind::Pipe,
        });
        assert!(state.sim_results.is_none(), "adding a link should invalidate results");

        state.sim_results = Some(fake_result(2, 1));
        state.apply(Action::SelectNode(a));
        state.apply(Action::DeleteSelection);
        assert!(state.sim_results.is_none(), "deleting a node should invalidate results");
    }

    #[test]
    fn undo_invalidates_sim_results() {
        let mut state = empty_state();
        state.apply(Action::AddJunctionAt { x: 0.0, y: 0.0 });
        state.sim_results = Some(fake_result(1, 0));

        state.apply(Action::Undo);

        assert!(state.sim_results.is_none());
    }
}
