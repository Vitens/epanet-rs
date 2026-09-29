//! The `eframe::App` implementation: collects `Action`s from every panel and
//! the canvas each frame, then applies them to `AppState` exactly once,
//! after layout — the "collect then apply" pattern from the plan.

use eframe::egui;

use gui_core::AppState;
use gui_core::action::Action;

use crate::canvas::{self, Interaction};
use crate::command_palette::{self, PaletteState};
use crate::panels;
use crate::shortcuts;

pub struct EditorApp {
    state: AppState,
    interaction: Interaction,
    palette: PaletteState,
}

impl EditorApp {
    pub fn new(state: AppState) -> Self {
        Self {
            state,
            interaction: Interaction::default(),
            palette: PaletteState::default(),
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

        egui::CentralPanel::default().show(ui, |ui| {
            canvas::draw_canvas(ui, &mut self.state, &mut self.interaction, &mut actions);
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
