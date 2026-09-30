use eframe::egui::{self, Color32, Pos2, Rect, Sense, Stroke, Vec2};

use gui_core::AppState;
use gui_core::action::{Action, NewLinkKind, Tool};
use gui_core::camera::Bounds;

/// Ephemeral, UI-only interaction state that doesn't belong in `AppState`
/// (it's meaningless outside of "what is the mouse doing right now").
#[derive(Default)]
pub struct Interaction {
    dragging_node: Option<String>,
    connecting_from: Option<String>,
    box_select_start: Option<Pos2>,
    /// True while a universal pan gesture (spacebar or middle-mouse drag) is
    /// active, regardless of the currently selected `Tool`.
    panning: bool,
    /// `(camera.center.x, camera.center.y, camera.zoom)` as of the last
    /// frame - see `camera_is_settled` below.
    last_camera: Option<(f64, f64, f64)>,
}

impl Interaction {
    /// Cancel any in-progress add-link connection or box-select drag (e.g.
    /// when Escape backs out to the Select tool), so no stale rubber-band
    /// preview is left hanging around on the canvas. Doesn't touch
    /// `dragging_node`/`panning`: those end naturally when the mouse is
    /// released and shouldn't be interrupted mid-gesture.
    pub fn cancel_transient_gesture(&mut self) {
        self.connecting_from = None;
        self.box_select_start = None;
    }

    /// Returns whether `camera` is unchanged from the last call (i.e. the
    /// view has "settled"), updating the stored value either way.
    ///
    /// Used to skip the direction-arrow/value-label overlay while the
    /// camera is actively panning or zooming: that overlay is drawn with
    /// `egui::Painter` (unlike node/link shapes, which Galileo now owns -
    /// see `network_layers.rs`), so it CPU-tessellates from scratch on
    /// every call - up to `MAX_RENDERED_SHAPES` link arrows, every single
    /// frame. A camera-only change (pan/zoom) doesn't touch any feature
    /// data, so skipping this overlay for those frames costs nothing
    /// visually except the arrows/labels themselves blinking out for the
    /// gesture's duration and reappearing once it settles - a standard
    /// "defer expensive decoration during interaction" tradeoff, not a
    /// correctness issue.
    fn camera_is_settled(&mut self, camera: &gui_core::camera::Camera) -> bool {
        let key = (camera.center.0, camera.center.1, camera.zoom);
        let settled = self.last_camera == Some(key);
        self.last_camera = Some(key);
        settled
    }
}

const HIT_RADIUS_PX: f32 = 10.0;
const NODE_RADIUS_PX: f32 = 4.5;
/// Offset used to place a value label off to the side of a pipe's
/// direction arrow (drawn by Galileo - see `network_layers.rs`), so the
/// two don't sit directly on top of each other.
const ARROW_SPACING_PX: f32 = 46.0;
/// Simulation-result value labels ("the color alone is not clear at all")
/// only draw while at most this many nodes+links are on screen at once;
/// past that, text would overlap into an illegible smear and tank the frame
/// rate on very large networks (some real INP files have 100k+ junctions).
/// Zooming in (fewer entities visible) brings the labels back.
const MAX_LABELED_ENTITIES: usize = 300;
/// Hard ceiling on how many pipe polyline segments plus node symbols we ask
/// egui to tessellate in a single frame. Each `painter.line_segment` /
/// `circle` / `rect` / `convex_polygon` call becomes its own `Shape`,
/// independently tessellated by epaint into anti-aliased triangle geometry
/// with no batching; drawing every entity in a full, zoomed-out view of a
/// large (100k+ node/link) INP file can overflow wgpu's single-buffer size
/// limit and crash the whole app (`egui_vertex_buffer` exceeding the
/// device's max buffer size). Past this budget we draw only an evenly
/// spaced subset of the visible nodes/links instead of everything; zooming
/// in shrinks the visible set and brings back full detail automatically.
const MAX_RENDERED_SHAPES: usize = 60_000;

pub fn draw_canvas(
    ui: &mut egui::Ui,
    state: &mut AppState,
    interaction: &mut Interaction,
    instanced: &mut Option<crate::instanced_renderer::InstancedRenderer>,
    galileo_target: &mut Option<crate::galileo_target::GalileoTarget>,
    actions: &mut Vec<Action>,
) {
    let rect = ui.available_rect_before_wrap();
    actions.push(Action::SetViewport {
        width: rect.width() as f64,
        height: rect.height() as f64,
    });

    let response = ui.allocate_rect(rect, Sense::click_and_drag());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, Color32::from_gray(24));

    // Keep the spatial index current before using it for culling/hit-testing.
    state.ensure_spatial_index();

    let camera = state.camera;
    let to_screen = |world: (f64, f64)| -> Pos2 {
        let (x, y) = camera.world_to_screen(world);
        Pos2::new(rect.left() + x as f32, rect.top() + y as f32)
    };
    let to_world = |screen: Pos2| -> (f64, f64) {
        camera.screen_to_world((
            (screen.x - rect.left()) as f64,
            (screen.y - rect.top()) as f64,
        ))
    };

    // --- scroll to zoom, anchored at the pointer ---
    if response.hovered() {
        let scroll_y = ui.input(|i| i.smooth_scroll_delta.y);
        if scroll_y.abs() > 0.01 {
            if let Some(pos) = response.hover_pos() {
                let factor = 2.0_f64.powf(scroll_y as f64 / 200.0).clamp(0.1, 10.0);
                actions.push(Action::ZoomAt {
                    factor,
                    screen_x: (pos.x - rect.left()) as f64,
                    screen_y: (pos.y - rect.top()) as f64,
                });
            }
        }
    }

    let hit_radius_world = (HIT_RADIUS_PX / camera.zoom.max(1e-6) as f32) as f64;

    // --- universal pan (spacebar or middle-mouse drag), regardless of tool ---
    let panning = handle_universal_pan(ui, &response, interaction, actions);

    // --- gesture handling, dispatched by the active tool (suppressed while panning) ---
    if !panning {
        match &state.tool {
            Tool::Pan => handle_pan(&response, actions),
            Tool::Select => handle_select(
                &response,
                state,
                to_world,
                hit_radius_world,
                interaction,
                actions,
            ),
            Tool::AddJunction => handle_add_at(&response, to_world, actions, |x, y| {
                Action::AddJunctionAt { x, y }
            }),
            Tool::AddTank => {
                handle_add_at(&response, to_world, actions, |x, y| Action::AddTankAt { x, y })
            }
            Tool::AddReservoir => handle_add_at(&response, to_world, actions, |x, y| {
                Action::AddReservoirAt { x, y }
            }),
            Tool::AddLink { kind } => handle_add_link(
                &response,
                state,
                to_world,
                hit_radius_world,
                kind.clone(),
                interaction,
                actions,
            ),
        }
    }

    // --- apply this frame's camera-changing actions immediately ---
    // `app.rs`'s "collect then apply" pattern (see its own doc comment)
    // is correct for domain edits - it needs every action gesture-handled
    // above to be visible together, after layout, exactly once. But for a
    // *continuous* drag/scroll gesture like pan/zoom, deferring the
    // resulting `Action::PanBy`/`ZoomAt` to that end-of-frame batch means
    // this frame's rendering below (which already captured `camera`,
    // above) draws with *last* frame's camera - the on-screen network
    // would always trail the mouse by one frame, worse under any
    // additional GPU/present-queue latency on top. Pulled out and applied
    // right here instead, then `camera` is re-read fresh before anything
    // below uses it for drawing.
    //
    // Not undo-tracked either way (see `gui_core::state::AppState`), so
    // applying immediately - rather than through the batched path domain
    // edits go through - has no correctness/undo implications. The
    // gesture handling *above* this point (hit-testing for select/connect/
    // box-select) correctly still used the *pre*-pan camera: a click or
    // drag that started this frame should hit-test against what was
    // actually on screen at the top of the frame, not a camera that
    // hasn't been drawn yet.
    let mut remaining_actions = Vec::with_capacity(actions.len());
    for action in actions.drain(..) {
        match action {
            Action::PanBy { .. } | Action::ZoomAt { .. } | Action::CenterOn { .. } => {
                state.apply(action);
            }
            other => remaining_actions.push(other),
        }
    }
    *actions = remaining_actions;
    let camera = state.camera;

    // --- basemap backdrop (currently empty - see `galileo_target.rs`) ---
    // NOT painted right now, on purpose: `target.paint`'s render pass
    // clears its texture (to opaque white, as it happens) even with zero
    // layers, and blitting that over the background fill above replaced
    // it wholesale - the reported white canvas background (and the
    // "selected" white highlight becoming invisible against it) traced
    // back to exactly this, not to anything in `instanced_renderer.rs`.
    // Skipping the paint call - not just skipping *adding* a basemap
    // layer - avoids drawing an empty backdrop that actively overwrites
    // the real background with nothing. `target`'s `Map`/view setup
    // stays wired below (harmless, and camera-synced) so the only change
    // needed when a real `.pmtiles` layer exists is re-enabling this
    // block, not rebuilding the integration.
    if let Some(target) = galileo_target.as_mut() {
        let map_center = match state.geo_transform.as_ref() {
            Some(transform) => transform.to_mercator(camera.center),
            None => camera.center,
        };
        let resolution = 1.0 / camera.zoom.max(1e-9);
        let view = galileo::MapView::new_projected_with_crs(
            &galileo_types::cartesian::Point2::new(map_center.0, map_center.1),
            resolution,
            galileo_types::geo::Crs::EPSG3857,
        )
        .with_size(galileo_types::cartesian::Size::new(
            rect.width().max(1.0) as f64,
            rect.height().max(1.0) as f64,
        ));
        target.map().set_view(view);
    }

    // --- network topology: true GPU instancing (see `instanced_renderer.rs`) ---
    // Two draw calls total, regardless of node/link count - one instanced
    // draw for node symbols + direction arrows, one for links. `sync` is
    // cheap when nothing changed; only the changed instances' bytes get
    // re-uploaded (`queue.write_buffer`), not a CPU re-tessellation pass.
    if let Some(renderer) = instanced.as_mut() {
        renderer.sync(state);
        renderer.set_camera(camera.center, camera.zoom, (rect.width(), rect.height()));
        for shape in renderer.paint_shape(rect) {
            painter.add(shape);
        }
    }

    // --- render ---
    let visible_bounds = camera.visible_bounds();
    let network = state.network();
    let spatial = state.spatial_index();
    let results = ResultOverlay::new(state);
    let visible_links = spatial.visible_links(&visible_bounds);
    let visible_nodes = spatial.visible_nodes(&visible_bounds);
    // See `MAX_RENDERED_SHAPES`: count the actual shape-generating work this
    // frame would do (one shape per node, one per pipe polyline segment)
    // and, if it's over budget, only draw an evenly spaced subset of the
    // visible nodes/links rather than everything - keeps the per-frame
    // vertex buffer bounded no matter how large/zoomed-out the network is.
    //
    // This budget now only bounds the *overlay* pass below (direction
    // arrows + value labels via `egui::Painter`) - node/link shapes and
    // colors themselves are drawn by the Galileo layer above, which has no
    // decimation cutoff (see `network_layers.rs`).
    let total_link_segments: usize = visible_links
        .iter()
        .filter_map(|&i| spatial.link_line(i))
        .map(|line| line.len().saturating_sub(1).max(1))
        .sum();
    let total_shapes = total_link_segments + visible_nodes.len();
    let decimation_stride = if total_shapes > MAX_RENDERED_SHAPES {
        total_shapes / MAX_RENDERED_SHAPES + 1
    } else {
        1
    };
    // See `MAX_LABELED_ENTITIES`: only label results when the view isn't
    // already crowded, so labels stay legible instead of turning into an
    // unreadable smear on large/zoomed-out networks.
    let show_value_labels =
        results.enabled && visible_links.len() + visible_nodes.len() <= MAX_LABELED_ENTITIES;

    // Direction arrows + value labels are drawn with `egui::Painter`
    // (unlike node/link shapes, which Galileo owns - see
    // `network_layers.rs`), so - unlike Galileo's own redraw-skip in
    // `galileo_target.rs` - they CPU-tessellate from scratch on every call
    // regardless of whether anything changed. Skip this overlay entirely
    // while the camera is actively panning/zooming (see
    // `Interaction::camera_is_settled`'s doc comment): a pan/zoom touches
    // no feature data, so there's nothing for this overlay to redraw
    // *correctly* anyway, only redundantly - it reappears once the camera
    // settles.
    let overlay_enabled = interaction.camera_is_settled(&camera);

    // Value labels only - node/link shapes, colors, and direction arrows
    // are all drawn by the Galileo layer painted above (see
    // `network_layers.rs`; `LinkSymbol::render` draws the arrow as part of
    // each link's own cached geometry, not per-frame here).
    if overlay_enabled {
        for (pos, &link_index) in visible_links.iter().enumerate() {
            if pos % decimation_stride != 0 {
                continue;
            }
            if !show_value_labels {
                continue;
            }
            let Some(line) = spatial.link_line(link_index) else {
                continue;
            };
            let screen_points: Vec<Pos2> = line.iter().map(|&p| to_screen(p)).collect();
            if let Some(label) = results.link_label(link_index) {
                // Offset away from the pipe's exact midpoint (roughly
                // where its direction arrow sits) so the label doesn't
                // sit directly on top of it.
                let total = polyline_length(&screen_points);
                let label_dist = (total * 0.5 + ARROW_SPACING_PX * 0.5).min(total);
                let (label_anchor, tangent) = point_and_tangent_at(&screen_points, label_dist);
                let normal = Vec2::new(-tangent.y, tangent.x);
                let label_pos = label_anchor + normal * 15.0;
                draw_value_label(&painter, label_pos, egui::Align2::CENTER_CENTER, &label);
            }
        }
    }

    if overlay_enabled && show_value_labels {
        for (pos, &node_index) in visible_nodes.iter().enumerate() {
            if pos % decimation_stride != 0 {
                continue;
            }
            let Some(point) = spatial.node_point(node_index) else {
                continue;
            };
            let center = to_screen(point);
            if let Some(label) = results.node_label(node_index) {
                draw_value_label(
                    &painter,
                    center + Vec2::new(NODE_RADIUS_PX + 4.0, -NODE_RADIUS_PX - 4.0),
                    egui::Align2::LEFT_BOTTOM,
                    &label,
                );
            }
        }
    }


    // connect-tool rubber band preview
    if let Some(from_id) = &interaction.connecting_from {
        if let Some(&idx) = network.node_map.get(from_id.as_str()) {
            if let Some(point) = network.nodes[idx].coordinates {
                if let Some(pos) = response.hover_pos().or(response.interact_pointer_pos()) {
                    painter
                        .line_segment([to_screen(point), pos], Stroke::new(1.5, Color32::YELLOW));
                }
            }
        }
    }

    // box-select preview
    if let Some(start) = interaction.box_select_start {
        if let Some(pos) = response.interact_pointer_pos().or(response.hover_pos()) {
            let selection_rect = Rect::from_two_pos(start, pos);
            painter.rect_stroke(
                selection_rect,
                0.0,
                Stroke::new(1.0, Color32::LIGHT_BLUE),
                egui::StrokeKind::Middle,
            );
        }
    }

    results.draw_legend(&painter, rect, show_value_labels, visible_links.len() + visible_nodes.len());
    if decimation_stride > 1 {
        painter.text(
            Pos2::new(rect.left() + 8.0, rect.bottom() - 20.0),
            egui::Align2::LEFT_BOTTOM,
            format!(
                "Showing 1 in {decimation_stride} nodes/links ({total_shapes} visible) - zoom in for full detail"
            ),
            egui::FontId::monospace(12.0),
            Color32::from_gray(200),
        );
    }
}

/// Pan the camera while the spacebar or the middle mouse button is held,
/// regardless of the currently active `Tool` — so the user never has to
/// switch out of Select (or Add-*) mode just to pan around.
///
/// Returns `true` if a pan gesture consumed this frame's pointer input, in
/// which case the caller must skip its normal per-tool gesture handling
/// (otherwise a space+left-drag would simultaneously move a node/box-select
/// under the Select tool while also panning).
fn handle_universal_pan(
    ui: &egui::Ui,
    response: &egui::Response,
    interaction: &mut Interaction,
    actions: &mut Vec<Action>,
) -> bool {
    let (space_down, middle_down, primary_down, delta) = ui.input(|i| {
        (
            i.key_down(egui::Key::Space),
            i.pointer.button_down(egui::PointerButton::Middle),
            i.pointer.button_down(egui::PointerButton::Primary),
            i.pointer.delta(),
        )
    });

    let engage = middle_down || (space_down && primary_down);

    if engage && (interaction.panning || response.hovered()) {
        interaction.panning = true;
    }

    if interaction.panning {
        if engage {
            if delta != Vec2::ZERO {
                actions.push(Action::PanBy { dx: delta.x as f64, dy: delta.y as f64 });
            }
            ui.ctx().set_cursor_icon(egui::CursorIcon::Grabbing);
            return true;
        }
        interaction.panning = false;
    } else if (space_down || middle_down) && response.hovered() {
        ui.ctx().set_cursor_icon(egui::CursorIcon::Grab);
    }

    false
}

fn handle_pan(response: &egui::Response, actions: &mut Vec<Action>) {
    if response.dragged() {
        let delta = response.drag_delta();
        actions.push(Action::PanBy { dx: delta.x as f64, dy: delta.y as f64 });
    }
}

fn handle_add_at(
    response: &egui::Response,
    to_world: impl Fn(Pos2) -> (f64, f64),
    actions: &mut Vec<Action>,
    make_action: impl Fn(f64, f64) -> Action,
) {
    if response.clicked() {
        if let Some(pos) = response.interact_pointer_pos() {
            let (x, y) = to_world(pos);
            actions.push(make_action(x, y));
        }
    }
}

fn handle_select(
    response: &egui::Response,
    state: &AppState,
    to_world: impl Fn(Pos2) -> (f64, f64),
    hit_radius_world: f64,
    interaction: &mut Interaction,
    actions: &mut Vec<Action>,
) {
    let spatial = state.spatial_index();
    let extend = response_extends_selection(response);

    if response.drag_started() {
        if let Some(pos) = response.interact_pointer_pos() {
            let world = to_world(pos);
            if let Some(node_index) = spatial.nearest_node(world, hit_radius_world) {
                let id = state.network().nodes[node_index].id.to_string();
                if !state.selection.contains_node(&id) {
                    if extend {
                        actions.push(Action::AddNodeToSelection(id.clone()));
                    } else {
                        actions.push(Action::SelectNode(id.clone()));
                    }
                }
                interaction.dragging_node = Some(id);
            } else {
                interaction.box_select_start = Some(pos);
            }
        }
    } else if response.dragged() {
        if let Some(id) = interaction.dragging_node.clone() {
            if let Some(pos) = response.interact_pointer_pos() {
                let (x, y) = to_world(pos);
                actions.push(Action::MoveNodeTo { id, x, y });
            }
        }
        // box select rectangle is drawn from `box_select_start` to the
        // current pointer position each frame in `draw_canvas`.
    } else if response.drag_stopped() {
        if interaction.dragging_node.take().is_some() {
            actions.push(Action::CommitDrag);
        } else if let Some(start) = interaction.box_select_start.take() {
            if let Some(end) = response.interact_pointer_pos() {
                let world_start = to_world(start);
                let world_end = to_world(end);
                let bounds = Bounds {
                    min_x: world_start.0.min(world_end.0),
                    min_y: world_start.1.min(world_end.1),
                    max_x: world_start.0.max(world_end.0),
                    max_y: world_start.1.max(world_end.1),
                };
                let mut ids: Vec<String> = spatial
                    .nodes_in_bounds(bounds)
                    .into_iter()
                    .map(|i| state.network().nodes[i].id.to_string())
                    .collect();
                if extend {
                    for existing in &state.selection.nodes {
                        if !ids.contains(existing) {
                            ids.push(existing.clone());
                        }
                    }
                }
                actions.push(Action::SelectInBox { nodes: ids, links: vec![] });
            }
        }
    } else if response.clicked() {
        if let Some(pos) = response.interact_pointer_pos() {
            let world = to_world(pos);
            if let Some(node_index) = spatial.nearest_node(world, hit_radius_world) {
                let id = state.network().nodes[node_index].id.to_string();
                if extend {
                    actions.push(Action::AddNodeToSelection(id));
                } else {
                    actions.push(Action::SelectNode(id));
                }
            } else if let Some(link_index) = spatial.nearest_link(world, hit_radius_world) {
                let id = state.network().links[link_index].id.to_string();
                if extend {
                    actions.push(Action::AddLinkToSelection(id));
                } else {
                    actions.push(Action::SelectLink(id));
                }
            } else if !extend {
                actions.push(Action::ClearSelection);
            }
        }
    }
}

/// Shift or Ctrl/Cmd held: extend the selection instead of replacing it.
fn response_extends_selection(response: &egui::Response) -> bool {
    response.ctx.input(|i| i.modifiers.shift || i.modifiers.command)
}

fn handle_add_link(
    response: &egui::Response,
    state: &AppState,
    to_world: impl Fn(Pos2) -> (f64, f64),
    hit_radius_world: f64,
    kind: NewLinkKind,
    interaction: &mut Interaction,
    actions: &mut Vec<Action>,
) {
    if !response.clicked() {
        return;
    }
    let Some(pos) = response.interact_pointer_pos() else {
        return;
    };
    let world = to_world(pos);
    let hit = state
        .spatial_index()
        .nearest_node(world, hit_radius_world)
        .map(|i| state.network().nodes[i].id.to_string());

    match (interaction.connecting_from.take(), hit) {
        (None, Some(id)) => interaction.connecting_from = Some(id),
        (Some(from), Some(to)) if from != to => {
            actions.push(Action::ConnectLink { from, to, kind });
        }
        // clicked empty space, or clicked the same node twice: cancel
        _ => {}
    }
}

/// Total on-screen length of a polyline (sum of its segment lengths).
fn polyline_length(points: &[Pos2]) -> f32 {
    points.windows(2).map(|w| w[0].distance(w[1])).sum()
}

/// The point at arc-length `target_dist` along a screen-space polyline
/// (clamped to the polyline's actual length), and the (normalized) local
/// tangent direction there, pointing from `points[0]` toward `points[last]`.
/// Walks the actual vertex chain rather than assuming a straight line
/// between the two endpoints, so placement on a curved/multi-vertex pipe
/// lands on the drawn route instead of floating off to one side of it.
fn point_and_tangent_at(points: &[Pos2], target_dist: f32) -> (Pos2, Vec2) {
    let fallback_tangent = Vec2::new(1.0, 0.0);
    if points.len() < 2 {
        return (points.first().copied().unwrap_or(Pos2::ZERO), fallback_tangent);
    }
    let target_dist = target_dist.max(0.0);
    let mut walked = 0.0f32;
    for (i, w) in points.windows(2).enumerate() {
        let seg_len = w[0].distance(w[1]);
        let is_last_segment = i + 2 == points.len();
        if walked + seg_len >= target_dist || is_last_segment {
            let t = if seg_len > 0.0 {
                ((target_dist - walked) / seg_len).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let pos = w[0] + (w[1] - w[0]) * t;
            let dir = w[1] - w[0];
            let tangent = if dir.length_sq() < 1e-6 { fallback_tangent } else { dir.normalized() };
            return (pos, tangent);
        }
        walked += seg_len;
    }
    (points[0], fallback_tangent)
}

/// Draws `text` anchored at `pos` with a cheap dark outline (the same text
/// re-drawn 1px off in each direction, then once more on top in white) so
/// it stays legible over any node/link/background color without having to
/// measure text extents to paint a proper backdrop rect. Used for the
/// per-entity simulation-result value labels - "the color alone is not
/// clear at all" - since the heatmap color alone doesn't say a number.
fn draw_value_label(painter: &egui::Painter, pos: Pos2, anchor: egui::Align2, text: &str) {
    let font = egui::FontId::monospace(10.0);
    for offset in [
        Vec2::new(-1.0, 0.0),
        Vec2::new(1.0, 0.0),
        Vec2::new(0.0, -1.0),
        Vec2::new(0.0, 1.0),
    ] {
        painter.text(pos + offset, anchor, text, font.clone(), Color32::from_black_alpha(235));
    }
    painter.text(pos, anchor, text, font, Color32::WHITE);
}

/// Computes (and caches, for one frame) the value-mapped node/link colors
/// for the currently selected simulation report step, plus the legend text.
/// `None` whenever there's no result to show, or the result no longer
/// matches the live network's node/link count (e.g. topology was edited
/// after the run) — callers fall back to the static per-type colors.
/// Simulation-result color/label overlay for the current report step.
/// `pub(crate)`: also used by `network_layers.rs` to color Galileo
/// features, not just this module's own `egui::Painter` legend.
pub(crate) struct ResultOverlay {
    node_range: Option<(f64, f64)>,
    link_range: Option<(f64, f64)>,
    node_values: Vec<f64>,
    link_values: Vec<f64>,
    /// Signed flow per link at the current report step (`None` when there's
    /// no result), used to flip the direction arrow to match simulated flow
    /// rather than the model's static start->end convention.
    link_flow_signs: Vec<f64>,
    /// Head/pressure unit abbreviation ("ft"/"m") and flow/demand unit
    /// abbreviation (e.g. "GPM", "LPS") for the network's configured units,
    /// so every displayed value carries the right unit alongside it.
    head_unit: &'static str,
    flow_unit: String,
    enabled: bool,
}

impl ResultOverlay {
    pub(crate) fn new(state: &AppState) -> Self {
        let Some(results) = &state.sim_results else {
            return Self::disabled();
        };
        let step = state.report_step.min(results.heads.len().saturating_sub(1));
        let Some(heads) = results.heads.get(step) else {
            return Self::disabled();
        };
        let Some(flows) = results.flows.get(step) else {
            return Self::disabled();
        };
        let network = state.network();
        if heads.len() != network.nodes.len() || flows.len() != network.links.len() {
            // Network topology changed since this result was computed.
            return Self::disabled();
        }

        let per_feet = network.options.unit_system.per_feet();
        let node_values: Vec<f64> = heads
            .iter()
            .zip(network.nodes.iter())
            .map(|(head, node)| head - node.elevation * per_feet)
            .collect();
        let link_values: Vec<f64> = flows.iter().map(|f| f.abs()).collect();

        Self {
            node_range: value_range(&node_values),
            link_range: value_range(&link_values),
            node_values,
            link_values,
            link_flow_signs: flows.clone(),
            head_unit: crate::panels::head_unit_label(network),
            flow_unit: crate::panels::flow_unit_label(network),
            enabled: true,
        }
    }

    fn disabled() -> Self {
        Self {
            node_range: None,
            link_range: None,
            node_values: Vec::new(),
            link_values: Vec::new(),
            link_flow_signs: Vec::new(),
            head_unit: "",
            flow_unit: String::new(),
            enabled: false,
        }
    }

    /// The result-driven fill color for `node_index`'s heat value, if any
    /// results are active - `network_layers.rs` wraps this into a
    /// `galileo::Color` for the Galileo-rendered node fill (there's no
    /// `egui::Color32` dependency here since that module has none).
    pub(crate) fn node_color_rgb(&self, node_index: usize) -> Option<(u8, u8, u8)> {
        let (lo, hi) = self.node_range?;
        let value = *self.node_values.get(node_index)?;
        Some(heat_color_rgb(normalize(value, lo, hi)))
    }

    /// Same as `node_color_rgb`, for links.
    pub(crate) fn link_color_rgb(&self, link_index: usize) -> Option<(u8, u8, u8)> {
        let (lo, hi) = self.link_range?;
        let value = *self.link_values.get(link_index)?;
        Some(heat_color_rgb(normalize(value, lo, hi)))
    }

    /// Formatted pressure value (with unit) for the on-canvas label next to
    /// `node_index` ("the color alone is not clear at all").
    fn node_label(&self, node_index: usize) -> Option<String> {
        if !self.enabled {
            return None;
        }
        let value = *self.node_values.get(node_index)?;
        Some(format!("{value:.1} {}", self.head_unit))
    }

    /// Formatted signed-flow value (with unit) for the on-canvas label next
    /// to `link_index`.
    fn link_label(&self, link_index: usize) -> Option<String> {
        if !self.enabled {
            return None;
        }
        let value = *self.link_flow_signs.get(link_index)?;
        Some(format!("{value:.2} {}", self.flow_unit))
    }

    /// Whether the simulated flow on `link_index` runs from the link's end
    /// node to its start node (i.e. opposite the model's stored direction).
    /// Whether flow direction on `link_index` is the reverse of the
    /// model's start -> end convention (negative flow) - `pub(crate)`:
    /// also used by `network_layers.rs`'s GPU-rendered direction arrows,
    /// not just this module's own (removed) `egui::Painter` ones.
    pub(crate) fn flow_reversed(&self, link_index: usize) -> bool {
        self.enabled
            && self
                .link_flow_signs
                .get(link_index)
                .is_some_and(|&f| f < 0.0)
    }

    /// Draws one range's legend entry: a title line, a horizontal gradient
    /// bar showing the actual `heat_color` scale ("the color alone is not
    /// clear at all" - this spells out what each color means), and the
    /// low/high value (with unit) at the bar's two ends so the scale is
    /// visibly anchored to the current min/max rather than some fixed
    /// range. Returns the y position just below the drawn block, so
    /// callers can stack multiple legend entries.
    fn draw_range_legend(
        painter: &egui::Painter,
        top_left: Pos2,
        bar_width: f32,
        title: &str,
        lo: f64,
        hi: f64,
        unit: &str,
        decimals: usize,
    ) -> f32 {
        let title_font = egui::FontId::monospace(12.0);
        let value_font = egui::FontId::monospace(11.0);
        painter.text(top_left, egui::Align2::LEFT_TOP, title, title_font, Color32::WHITE);

        let bar_top = top_left.y + 15.0;
        let bar_rect = Rect::from_min_size(Pos2::new(top_left.x, bar_top), egui::vec2(bar_width, 10.0));
        draw_gradient_bar(painter, bar_rect);
        painter.rect_stroke(
            bar_rect,
            0.0,
            Stroke::new(1.0, Color32::from_gray(80)),
            egui::StrokeKind::Outside,
        );

        let value_y = bar_rect.bottom() + 2.0;
        painter.text(
            Pos2::new(bar_rect.left(), value_y),
            egui::Align2::LEFT_TOP,
            format!("{lo:.decimals$} {unit}"),
            value_font.clone(),
            Color32::WHITE,
        );
        painter.text(
            Pos2::new(bar_rect.right(), value_y),
            egui::Align2::RIGHT_TOP,
            format!("{hi:.decimals$} {unit}"),
            value_font,
            Color32::WHITE,
        );

        value_y + 14.0
    }

    fn draw_legend(&self, painter: &egui::Painter, rect: Rect, showing_value_labels: bool, visible_count: usize) {
        if !self.enabled {
            return;
        }
        let x = rect.left() + 8.0;
        let bar_width = 160.0_f32.min(rect.width() - 16.0).max(0.0);
        let mut y = rect.top() + 8.0;
        if let Some((lo, hi)) = self.node_range {
            y = Self::draw_range_legend(
                painter,
                Pos2::new(x, y),
                bar_width,
                "Pressure (head - elevation)",
                lo,
                hi,
                self.head_unit,
                1,
            );
        }
        if let Some((lo, hi)) = self.link_range {
            y = Self::draw_range_legend(
                painter,
                Pos2::new(x, y),
                bar_width,
                "|Flow|  (arrows show simulated direction)",
                lo,
                hi,
                &self.flow_unit,
                2,
            );
        }
        if !showing_value_labels {
            painter.text(
                Pos2::new(x, y),
                egui::Align2::LEFT_TOP,
                format!(
                    "Zoom in to see value labels ({visible_count} items visible, limit {MAX_LABELED_ENTITIES})"
                ),
                egui::FontId::monospace(12.0),
                Color32::from_gray(200),
            );
        }
    }
}

fn value_range(values: &[f64]) -> Option<(f64, f64)> {
    if values.is_empty() {
        return None;
    }
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in values {
        if v.is_finite() {
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    if !lo.is_finite() || !hi.is_finite() {
        return None;
    }
    Some((lo, hi))
}

fn normalize(value: f64, lo: f64, hi: f64) -> f64 {
    if (hi - lo).abs() < 1e-12 {
        0.5
    } else {
        ((value - lo) / (hi - lo)).clamp(0.0, 1.0)
    }
}

/// Blue (low) -> yellow (mid) -> red (high) heatmap, the same 3-stop
/// gradient used by the legend text.
fn heat_color(t: f64) -> Color32 {
    let (r, g, b) = heat_color_rgb(t);
    Color32::from_rgb(r, g, b)
}

/// The color math behind `heat_color`, without the `egui::Color32`
/// dependency - shared with `network_layers.rs`'s `galileo::Color`-based
/// feature fill, which has no reason to depend on `egui` at all.
fn heat_color_rgb(t: f64) -> (u8, u8, u8) {
    let t = t.clamp(0.0, 1.0) as f32;
    const LOW: (u8, u8, u8) = (60, 90, 220);
    const MID: (u8, u8, u8) = (250, 210, 60);
    const HIGH: (u8, u8, u8) = (220, 60, 60);
    fn lerp(a: (u8, u8, u8), b: (u8, u8, u8), t: f32) -> (u8, u8, u8) {
        let l = |x: u8, y: u8| (x as f32 + (y as f32 - x as f32) * t).round() as u8;
        (l(a.0, b.0), l(a.1, b.1), l(a.2, b.2))
    }
    if t < 0.5 { lerp(LOW, MID, t / 0.5) } else { lerp(MID, HIGH, (t - 0.5) / 0.5) }
}

/// Paints `rect` as a left-to-right swatch of the `heat_color` gradient
/// (low value at the left edge, high value at the right edge) - the actual
/// color scale for the legend, not just min/max text, since "the color
/// alone is not clear at all" cuts both ways: a number needs the scale to
/// read the color, and the color needs the scale to read the number.
fn draw_gradient_bar(painter: &egui::Painter, rect: Rect) {
    const STEPS: usize = 48;
    if rect.width() <= 0.0 {
        return;
    }
    let step_width = rect.width() / STEPS as f32;
    for i in 0..STEPS {
        let t = i as f64 / (STEPS - 1) as f64;
        let x0 = rect.left() + step_width * i as f32;
        // Slightly overlap adjacent strips so there's no antialiasing seam.
        let strip = Rect::from_min_max(
            Pos2::new(x0, rect.top()),
            Pos2::new((x0 + step_width + 0.5).min(rect.right()), rect.bottom()),
        );
        painter.rect_filled(strip, 0.0, heat_color(t));
    }
}
