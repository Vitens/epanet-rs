//! The `eframe::App` implementation: collects `Action`s from every panel and
//! the canvas each frame, then applies them to `AppState` exactly once,
//! after layout — the "collect then apply" pattern from the plan.

use eframe::egui;

use gui_core::AppState;
use gui_core::action::Action;

use crate::canvas::{self, Interaction};
use crate::command_palette::{self, PaletteState};
use crate::galileo_target::GalileoTarget;
use crate::instanced_renderer::InstancedRenderer;
use crate::panels;
use crate::shortcuts;

pub struct EditorApp {
    state: AppState,
    interaction: Interaction,
    palette: PaletteState,
    instanced: Option<InstancedRenderer>,
    /// Kept around empty (no layers) for now - see `galileo_target.rs`'s
    /// doc comment: node/link topology rendering moved to
    /// `instanced_renderer.rs`'s true GPU instancing, but a future
    /// `.pmtiles`/vector-tile basemap layer still belongs here, composited
    /// as the backdrop underneath the instanced topology.
    galileo_target: Option<GalileoTarget>,
    /// Currently unused: `GalileoTarget::paint` isn't called right now
    /// (see `canvas.rs`'s basemap-backdrop comment - painting an empty
    /// map actively overwrote the real background). Kept for when that's
    /// re-enabled alongside a real basemap layer, rather than removed and
    /// re-added.
    #[allow(dead_code)]
    render_state: Option<egui_wgpu::RenderState>,
}

impl EditorApp {
    pub fn new(cc: &eframe::CreationContext<'_>, state: AppState) -> Self {
        let render_state = cc.wgpu_render_state.clone();

        // `render_state` is only `None` if eframe somehow initialized a
        // non-wgpu backend despite `NativeOptions::renderer =
        // Renderer::Wgpu` (see `main.rs`) - shouldn't happen in practice,
        // but degrade to "nothing GPU-drawn" rather than panicking if it
        // ever does.
        let instanced = render_state
            .as_ref()
            .map(|rs| InstancedRenderer::new(rs, &state));
        let galileo_target = render_state.as_ref().map(|rs| {
            let view = galileo::MapView::new_projected_with_crs(
                &galileo_types::cartesian::Point2::new(0.0, 0.0),
                1.0,
                galileo_types::geo::Crs::EPSG3857,
            );
            let map = galileo::Map::new(view, Vec::new(), None);
            GalileoTarget::new(rs, map)
        });

        Self {
            state,
            interaction: Interaction::default(),
            palette: PaletteState::default(),
            instanced,
            galileo_target,
            render_state,
        }
    }
}

impl eframe::App for EditorApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();

        // Drain any progress/finished messages from a background simulation
        // thread first, so they're applied in the same batch as this
        // frame's UI-driven actions.
        let mut actions: Vec<Action> = self.state.drain_simulation_messages();

        shortcuts::handle_shortcuts(
            &ctx,
            &self.state,
            &mut self.palette,
            &mut self.interaction,
            &mut actions,
        );

        panels::top_menu(ui, &self.state, &mut actions);
        panels::left_panel(ui, &self.state, &mut actions);
        panels::right_panel(ui, &self.state, &mut actions);
        panels::bottom_toolbar(&ctx, &self.state, &mut actions);

        egui::CentralPanel::default().show_inside(ui, |ui| {
            canvas::draw_canvas(
                ui,
                &mut self.state,
                &mut self.interaction,
                &mut self.instanced,
                &mut self.galileo_target,
                &mut actions,
            );
        });

        command_palette::draw(&ctx, &self.state, &mut self.palette, &mut actions);

        for action in actions {
            self.state.apply(action);
        }

        // Keep polling the background simulation thread even if the mouse
        // isn't moving, so the progress bar / result actually update.
        if self.state.sim_running {
            ctx.request_repaint();
        }
    }
}
