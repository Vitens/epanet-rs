//! Global keyboard shortcuts, turned into `Action`s (or a palette toggle)
//! once per frame. Plain letter shortcuts and Delete/Backspace are
//! suppressed entirely while any widget has keyboard focus (e.g. editing a
//! numeric field or the command palette's search box), so typing never gets
//! hijacked into a tool switch or a deletion.

use eframe::egui::{self, Key, Modifiers};

use epanet_rs::model::valve::ValveType;
use gui_core::AppState;
use gui_core::action::{Action, NewLinkKind, Tool};

use crate::canvas::Interaction;
use crate::command_palette::PaletteState;

pub fn handle_shortcuts(
    ctx: &egui::Context,
    state: &AppState,
    palette: &mut PaletteState,
    interaction: &mut Interaction,
    actions: &mut Vec<Action>,
) {
    // Ctrl/Cmd+P always toggles the quick-actions palette, regardless of
    // focus, and swallows the key so nothing else reacts to it.
    let ctrl_p = ctx.input_mut(|i| i.consume_key(Modifiers::COMMAND, Key::P));
    if ctrl_p {
        palette.toggle();
    }

    // While the palette is open it owns the keyboard: its search box is
    // focused, and its own `draw` reads Escape/Enter/arrow keys directly.
    if palette.open {
        return;
    }

    if ctx.memory(|m| m.focused().is_some()) {
        return;
    }

    ctx.input_mut(|i| {
        // Escape backs out of whatever tool/gesture is active and returns to
        // the plain Select tool, cancelling any in-progress add-link
        // connection or box-select so no stale rubber-band is left behind.
        if i.consume_key(Modifiers::NONE, Key::Escape) {
            interaction.cancel_transient_gesture();
            if state.tool != Tool::Select {
                actions.push(Action::SetTool(Tool::Select));
            }
        }

        // `consume_key`'s modifier match is "at least these modifiers, extra
        // Shift/Alt ignored" (see `Modifiers::matches_logically`), so a more
        // specific shortcut that is a superset of a plainer one (Ctrl+Shift+Z
        // vs Ctrl+Z, Shift+P vs P) must be checked - and have its event
        // consumed - first, or the specific one can never fire.
        if i.consume_key(Modifiers::COMMAND.plus(Modifiers::SHIFT), Key::Z)
            || i.consume_key(Modifiers::COMMAND, Key::Y)
        {
            actions.push(Action::Redo);
        } else if i.consume_key(Modifiers::COMMAND, Key::Z) {
            actions.push(Action::Undo);
        }

        if i.consume_key(Modifiers::COMMAND, Key::Space) && !state.sim_running {
            actions.push(Action::RunSimulation { parallel: false });
        }

        if i.consume_key(Modifiers::SHIFT, Key::Num0) {
            actions.push(Action::SetZoom(1.0));
        }
        if i.consume_key(Modifiers::SHIFT, Key::Num1) {
            actions.push(Action::FitToContent);
        }
        if i.consume_key(Modifiers::SHIFT, Key::Num2) && !state.selection.is_empty() {
            actions.push(Action::FitToSelection);
        }

        if !state.selection.is_empty()
            && (i.consume_key(Modifiers::NONE, Key::Delete)
                || i.consume_key(Modifiers::NONE, Key::Backspace))
        {
            actions.push(Action::DeleteSelection);
        }

        // Shift+P (add pump) before plain P (add pipe): see the ordering
        // note above.
        if i.consume_key(Modifiers::SHIFT, Key::P) {
            actions.push(Action::SetTool(Tool::AddLink {
                kind: NewLinkKind::Pump,
            }));
        } else if i.consume_key(Modifiers::NONE, Key::P) {
            actions.push(Action::SetTool(Tool::AddLink {
                kind: NewLinkKind::Pipe,
            }));
        }
        if i.consume_key(Modifiers::NONE, Key::J) {
            actions.push(Action::SetTool(Tool::AddJunction));
        }
        if i.consume_key(Modifiers::NONE, Key::T) {
            actions.push(Action::SetTool(Tool::AddTank));
        }
        if i.consume_key(Modifiers::NONE, Key::R) {
            actions.push(Action::SetTool(Tool::AddReservoir));
        }
        if i.consume_key(Modifiers::NONE, Key::V) {
            actions.push(Action::SetTool(Tool::AddLink {
                kind: NewLinkKind::Valve(ValveType::TCV),
            }));
        }
    });
}
