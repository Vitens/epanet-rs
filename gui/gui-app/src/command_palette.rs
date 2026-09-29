//! Ctrl+P quick-actions palette: a fuzzy-filterable list of every action
//! that makes sense to run right now (tool switches, undo/redo, view
//! commands, simulation, import/export, ...), mirroring the menu/toolbar
//! `Action`s from a single place. It doubles as a search bar: typing any
//! text also matches node/link IDs in the network. Picking a match always
//! centers the view on it; a link additionally zooms to fit it (it has
//! genuine extent to frame), while a node only re-centers, leaving the zoom
//! alone (a point has no "size" to fit to, and jumping the zoom around
//! while browsing single nodes would be disorienting).

use eframe::egui;

use epanet_rs::model::network::LinkUpdate;
use epanet_rs::model::valve::ValveType;
use gui_core::AppState;
use gui_core::action::{Action, NewLinkKind, Tool};

use crate::panels::{link_type_label, node_type_label};

/// Cap on how many node/link ID matches are shown at once, so a huge
/// network with a short/common query string doesn't flood the list.
const MAX_ENTITY_MATCHES: usize = 40;

/// Ephemeral UI-only state for the palette; not part of `AppState` since it
/// has no meaning outside "what is this modal currently showing".
#[derive(Default)]
pub struct PaletteState {
    pub open: bool,
    query: String,
    selected: usize,
    just_opened: bool,
}

impl PaletteState {
    /// Open the palette (resetting its query/selection) if closed, or close
    /// it if already open.
    pub fn toggle(&mut self) {
        if self.open {
            self.close();
        } else {
            self.open = true;
            self.query.clear();
            self.selected = 0;
            self.just_opened = true;
        }
    }

    pub fn close(&mut self) {
        self.open = false;
    }
}

/// What happens when a palette entry is chosen. Most entries just push a
/// single `Action`; import/export need a native file dialog first, and the
/// node/link search matches need to chain a select with a camera move, so
/// they get their own variants instead of trying to squeeze everything
/// through one `Action`.
enum Effect {
    Action(Action),
    Actions(Vec<Action>),
    ImportInp,
    ExportInp,
}

struct Entry {
    label: String,
    shortcut: Option<&'static str>,
    enabled: bool,
    effect: Effect,
}

impl Entry {
    fn action(label: &str, shortcut: Option<&'static str>, enabled: bool, action: Action) -> Self {
        Entry {
            label: label.to_string(),
            shortcut,
            enabled,
            effect: Effect::Action(action),
        }
    }
}

/// Draw the palette (a no-op if it isn't open) and apply the chosen entry's
/// effect directly into `actions` (or via a file dialog for import/export).
pub fn draw(ctx: &egui::Context, state: &AppState, palette: &mut PaletteState, actions: &mut Vec<Action>) {
    if !palette.open {
        return;
    }

    let mut entries = build_entries(state, &palette.query);
    let query_lower = palette.query.to_lowercase();
    let filtered_indices: Vec<usize> = entries
        .iter()
        .enumerate()
        .filter(|(_, e)| {
            e.enabled && (query_lower.is_empty() || e.label.to_lowercase().contains(&query_lower))
        })
        .map(|(i, _)| i)
        .collect();

    if filtered_indices.is_empty() {
        palette.selected = 0;
    } else if palette.selected >= filtered_indices.len() {
        palette.selected = filtered_indices.len() - 1;
    }

    let mut close = false;
    let mut chosen_index: Option<usize> = None;

    egui::Window::new("Quick actions")
        .id(egui::Id::new("quick_actions_palette"))
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_TOP, egui::vec2(0.0, 80.0))
        .order(egui::Order::Foreground)
        .fixed_size(egui::vec2(440.0, 380.0))
        .show(ctx, |ui| {
            let text_response = ui.add(
                egui::TextEdit::singleline(&mut palette.query)
                    .hint_text("Type a command, or a node/link ID to select it...")
                    .desired_width(f32::INFINITY),
            );
            if palette.just_opened {
                text_response.request_focus();
                palette.just_opened = false;
            }

            ui.input(|i| {
                if i.key_pressed(egui::Key::Escape) {
                    close = true;
                }
                if !filtered_indices.is_empty() {
                    if i.key_pressed(egui::Key::ArrowDown) {
                        palette.selected = (palette.selected + 1).min(filtered_indices.len() - 1);
                    }
                    if i.key_pressed(egui::Key::ArrowUp) {
                        palette.selected = palette.selected.saturating_sub(1);
                    }
                    if i.key_pressed(egui::Key::Enter) {
                        chosen_index = Some(filtered_indices[palette.selected]);
                    }
                }
            });

            ui.separator();
            egui::ScrollArea::vertical()
                .max_height(320.0)
                .show(ui, |ui| {
                    if filtered_indices.is_empty() {
                        ui.weak("No matching commands or IDs.");
                    }
                    for (row, &entry_idx) in filtered_indices.iter().enumerate() {
                        let entry = &entries[entry_idx];
                        let row_selected = row == palette.selected;
                        let text = match entry.shortcut {
                            Some(sc) => format!("{}   [{sc}]", entry.label),
                            None => entry.label.clone(),
                        };
                        if ui.selectable_label(row_selected, text).clicked() {
                            chosen_index = Some(entry_idx);
                        }
                    }
                });
        });

    if let Some(idx) = chosen_index {
        let entry = entries.swap_remove(idx);
        apply_effect(entry.effect, actions);
        close = true;
    }

    if close {
        palette.close();
    }
}

fn apply_effect(effect: Effect, actions: &mut Vec<Action>) {
    match effect {
        Effect::Action(action) => actions.push(action),
        Effect::Actions(actions_to_push) => actions.extend(actions_to_push),
        Effect::ImportInp => {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("EPANET INP", &["inp"])
                .pick_file()
            {
                actions.push(Action::ImportInp(path));
            }
        }
        Effect::ExportInp => {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("EPANET INP", &["inp"])
                .save_file()
            {
                actions.push(Action::ExportInp(path));
            }
        }
    }
}

fn build_entries(state: &AppState, query: &str) -> Vec<Entry> {
    let mut entries = Vec::new();

    entries.push(Entry::action("Undo", Some("Ctrl+Z"), state.can_undo(), Action::Undo));
    entries.push(Entry::action(
        "Redo",
        Some("Ctrl+Shift+Z"),
        state.can_redo(),
        Action::Redo,
    ));
    entries.push(Entry::action(
        "Delete selection",
        Some("Delete"),
        !state.selection.is_empty(),
        Action::DeleteSelection,
    ));
    entries.push(Entry::action(
        "Clear selection",
        None,
        !state.selection.is_empty(),
        Action::ClearSelection,
    ));

    entries.push(Entry::action("Select tool", None, true, Action::SetTool(Tool::Select)));
    entries.push(Entry::action(
        "Pan tool",
        None,
        true,
        Action::SetTool(Tool::Pan),
    ));
    entries.push(Entry::action(
        "Add Junction",
        Some("J"),
        true,
        Action::SetTool(Tool::AddJunction),
    ));
    entries.push(Entry::action(
        "Add Tank",
        Some("T"),
        true,
        Action::SetTool(Tool::AddTank),
    ));
    entries.push(Entry::action(
        "Add Reservoir",
        Some("R"),
        true,
        Action::SetTool(Tool::AddReservoir),
    ));
    entries.push(Entry::action(
        "Add Pipe",
        Some("P"),
        true,
        Action::SetTool(Tool::AddLink { kind: NewLinkKind::Pipe }),
    ));
    entries.push(Entry::action(
        "Add Pump",
        Some("Shift+P"),
        true,
        Action::SetTool(Tool::AddLink { kind: NewLinkKind::Pump }),
    ));
    for vt in [
        ValveType::PRV,
        ValveType::PSV,
        ValveType::PBV,
        ValveType::FCV,
        ValveType::TCV,
        ValveType::PCV,
        ValveType::GPV,
    ] {
        let shortcut = if vt == ValveType::TCV { Some("V") } else { None };
        entries.push(Entry::action(
            &format!("Add Valve ({vt})"),
            shortcut,
            true,
            Action::SetTool(Tool::AddLink {
                kind: NewLinkKind::Valve(vt),
            }),
        ));
    }

    entries.push(Entry::action(
        "Run simulation",
        Some("Ctrl+Space"),
        !state.sim_running,
        Action::RunSimulation { parallel: false },
    ));

    entries.push(Entry::action("Fit all", Some("Shift+1"), true, Action::FitToContent));
    entries.push(Entry::action(
        "Fit to selection",
        Some("Shift+2"),
        !state.selection.is_empty(),
        Action::FitToSelection,
    ));
    entries.push(Entry::action(
        "Zoom to 100%",
        Some("Shift+0"),
        true,
        Action::SetZoom(1.0),
    ));

    if let Some(id) = state.selection.as_single_link()
        && let Some(&idx) = state.network().link_map.get(id)
    {
        let link = &state.network().links[idx];
        let start = link.start_node_id.to_string();
        let end = link.end_node_id.to_string();
        // Reverse the display vertices along with the endpoints, or the
        // polyline ends up crossing itself (start/end swapped but the
        // vertex order still reads in the old direction).
        let vertices = link.vertices.as_ref().map(|v| v.iter().rev().copied().collect());
        entries.push(Entry::action(
            "Swap pipe direction",
            None,
            true,
            Action::UpdateLinks {
                ids: vec![id.to_string()],
                update: LinkUpdate {
                    start_node: Some(end.into()),
                    end_node: Some(start.into()),
                    vertices,
                    ..Default::default()
                },
            },
        ));
    }

    entries.push(Entry {
        label: "Import INP...".to_string(),
        shortcut: None,
        enabled: true,
        effect: Effect::ImportInp,
    });
    entries.push(Entry {
        label: "Export INP...".to_string(),
        shortcut: None,
        enabled: true,
        effect: Effect::ExportInp,
    });

    entries.extend(entity_matches(state, query));

    entries
}

/// Rank how well `haystack` (already lowercased) matches `needle` (already
/// lowercased): exact match ranks best, then prefix, then substring
/// anywhere. `None` if it doesn't match at all.
fn match_rank(haystack: &str, needle: &str) -> Option<u8> {
    if haystack == needle {
        Some(0)
    } else if haystack.starts_with(needle) {
        Some(1)
    } else if haystack.contains(needle) {
        Some(2)
    } else {
        None
    }
}

/// Node/link IDs whose text matches `query`, turned into "select" entries.
/// Picking a node selects it and centers the camera on its coordinates
/// without touching zoom; picking a link selects it and fits/zooms the
/// camera to it (a link has real extent to frame; a node is just a point).
/// Empty (rather than listing every node/link) when there's no query, since
/// that could be thousands of entries on a large network.
fn entity_matches(state: &AppState, query: &str) -> Vec<Entry> {
    if query.trim().is_empty() {
        return Vec::new();
    }
    let query_lower = query.to_lowercase();
    let network = state.network();

    let mut ranked: Vec<(u8, Entry)> = Vec::new();

    for node in &network.nodes {
        let id_lower = node.id.to_lowercase();
        if let Some(rank) = match_rank(&id_lower, &query_lower) {
            let label = format!("Select node {}  [{}]", node.id, node_type_label(&node.node_type));
            let mut acts = vec![Action::SelectNode(node.id.to_string())];
            if let Some((x, y)) = node.coordinates {
                acts.push(Action::CenterOn { x, y });
            }
            ranked.push((
                rank,
                Entry {
                    label,
                    shortcut: None,
                    enabled: true,
                    effect: Effect::Actions(acts),
                },
            ));
        }
    }
    for link in &network.links {
        let id_lower = link.id.to_lowercase();
        if let Some(rank) = match_rank(&id_lower, &query_lower) {
            let label = format!("Select link {}  [{}]", link.id, link_type_label(&link.link_type));
            let acts = vec![Action::SelectLink(link.id.to_string()), Action::FitToSelection];
            ranked.push((
                rank,
                Entry {
                    label,
                    shortcut: None,
                    enabled: true,
                    effect: Effect::Actions(acts),
                },
            ));
        }
    }

    ranked.sort_by_key(|(rank, entry)| (*rank, entry.label.clone()));
    ranked
        .into_iter()
        .take(MAX_ENTITY_MATCHES)
        .map(|(_, entry)| entry)
        .collect()
}
