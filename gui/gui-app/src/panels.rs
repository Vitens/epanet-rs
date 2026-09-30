//! Top menu, left element list, right property inspector, and the floating
//! bottom tool palette. Every widget here only ever pushes `Action`s into
//! the shared `actions` buffer; `AppState` is read-only from this module.

use eframe::egui;

use epanet_rs::model::link::{Link, LinkStatus, LinkType};
use epanet_rs::model::network::{
    JunctionUpdate, LinkUpdate, Network, NodeUpdate, PipeUpdate, PumpUpdate, ReservoirUpdate,
    TankUpdate, ValveUpdate,
};
use epanet_rs::model::node::{Node, NodeType};
use epanet_rs::model::options::HeadlossFormula;
use epanet_rs::model::units::{UnitConversion, UnitSystem};
use epanet_rs::model::valve::ValveType;
use gui_core::AppState;
use gui_core::action::{Action, NewLinkKind, Tool};

/// Every `ValveType`, in the order offered by pickers.
const VALVE_TYPES: [ValveType; 7] = [
    ValveType::PRV,
    ValveType::PSV,
    ValveType::PBV,
    ValveType::FCV,
    ValveType::TCV,
    ValveType::PCV,
    ValveType::GPV,
];

/// The three initial statuses a user may explicitly choose for a link (the
/// rest of `LinkStatus` are solver-computed outcomes, not valid inputs).
const LINK_STATUSES: [LinkStatus; 3] = [LinkStatus::Open, LinkStatus::Closed, LinkStatus::Active];

pub fn top_menu(ui: &mut egui::Ui, state: &AppState, actions: &mut Vec<Action>) {
    egui::Panel::top("top_menu").show_inside(ui, |ui| {
        egui::MenuBar::new().ui(ui, |ui| {
            ui.menu_button("File", |ui| {
                if ui.button("Import INP...").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("EPANET INP", &["inp"])
                        .pick_file()
                    {
                        actions.push(Action::ImportInp(path));
                    }
                    ui.close_kind(egui::UiKind::Menu);
                }
                if ui.button("Export INP...").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("EPANET INP", &["inp"])
                        .save_file()
                    {
                        actions.push(Action::ExportInp(path));
                    }
                    ui.close_kind(egui::UiKind::Menu);
                }
            });
            ui.menu_button("Edit", |ui| {
                if ui
                    .add_enabled(state.can_undo(), egui::Button::new("Undo (Ctrl+Z)"))
                    .clicked()
                {
                    actions.push(Action::Undo);
                    ui.close_kind(egui::UiKind::Menu);
                }
                if ui
                    .add_enabled(state.can_redo(), egui::Button::new("Redo (Ctrl+Shift+Z)"))
                    .clicked()
                {
                    actions.push(Action::Redo);
                    ui.close_kind(egui::UiKind::Menu);
                }
                ui.separator();
                if ui
                    .add_enabled(
                        !state.selection.is_empty(),
                        egui::Button::new("Delete selection (Delete)"),
                    )
                    .clicked()
                {
                    actions.push(Action::DeleteSelection);
                    ui.close_kind(egui::UiKind::Menu);
                }
            });
            ui.menu_button("View", |ui| {
                if ui.button("Fit all (Shift+1)").clicked() {
                    actions.push(Action::FitToContent);
                    ui.close_kind(egui::UiKind::Menu);
                }
                if ui
                    .add_enabled(
                        !state.selection.is_empty(),
                        egui::Button::new("Fit to selection (Shift+2)"),
                    )
                    .clicked()
                {
                    actions.push(Action::FitToSelection);
                    ui.close_kind(egui::UiKind::Menu);
                }
                if ui.button("Zoom to 100% (Shift+0)").clicked() {
                    actions.push(Action::SetZoom(1.0));
                    ui.close_kind(egui::UiKind::Menu);
                }
                ui.separator();
                ui.label("Quick actions & search... (Ctrl+P)");
            });

            ui.separator();
            ui.label(format!(
                "{} nodes, {} links",
                state.network().nodes.len(),
                state.network().links.len()
            ));

            if let Some(msg) = &state.status_message {
                ui.separator();
                ui.label(msg);
                if ui.small_button("x").clicked() {
                    actions.push(Action::ClearStatus);
                }
            }
        });
    });
}

pub fn left_panel(ui: &mut egui::Ui, state: &AppState, actions: &mut Vec<Action>) {
    egui::Panel::left("left_panel")
        .resizable(true)
        .default_size(260.0)
        .min_size(180.0)
        .show_inside(ui, |ui| {
            ui.heading("Elements");
            ui.add_space(4.0);

            egui::CollapsingHeader::new(format!("Nodes ({})", state.network().nodes.len()))
                .default_open(true)
                .show(ui, |ui| node_list(ui, state, actions));

            egui::CollapsingHeader::new(format!("Links ({})", state.network().links.len()))
                .default_open(true)
                .show(ui, |ui| link_list(ui, state, actions));
        });
}

/// A selectable row that always spans the full width of its container
/// (`ScrollArea` reports the scrollbar to whichever content is widest, so
/// without this the row buttons - and the scrollbar right after them - end
/// up only as wide as the longest id/label, floating short of the actual
/// panel edge instead of filling it).
fn full_width_selectable_row(ui: &mut egui::Ui, selected: bool, label: String) -> egui::Response {
    let full_width = ui.available_width();
    ui.add(egui::Button::selectable(selected, label).min_size(egui::vec2(full_width, 0.0)))
}

fn node_list(ui: &mut egui::Ui, state: &AppState, actions: &mut Vec<Action>) {
    let nodes = &state.network().nodes;
    let row_height = ui.text_style_height(&egui::TextStyle::Body);
    egui::ScrollArea::vertical()
        .id_salt("nodes_scroll")
        .max_height(240.0)
        .auto_shrink([false, true])
        .show_rows(ui, row_height, nodes.len(), |ui, range| {
            for i in range {
                let node = &nodes[i];
                let selected = state.selection.contains_node(&node.id);
                let label = format!("{}  [{}]", node.id, node_type_label(&node.node_type));
                let extend = ui.input(|i| i.modifiers.shift || i.modifiers.command);
                if full_width_selectable_row(ui, selected, label).clicked() {
                    if extend {
                        actions.push(Action::AddNodeToSelection(node.id.to_string()));
                    } else {
                        actions.push(Action::SelectNode(node.id.to_string()));
                    }
                }
            }
        });
}

fn link_list(ui: &mut egui::Ui, state: &AppState, actions: &mut Vec<Action>) {
    let links = &state.network().links;
    let row_height = ui.text_style_height(&egui::TextStyle::Body);
    egui::ScrollArea::vertical()
        .id_salt("links_scroll")
        .max_height(240.0)
        .auto_shrink([false, true])
        .show_rows(ui, row_height, links.len(), |ui, range| {
            for i in range {
                let link = &links[i];
                let selected = state.selection.contains_link(&link.id);
                let label = format!("{}  [{}]", link.id, link_type_label(&link.link_type));
                let extend = ui.input(|i| i.modifiers.shift || i.modifiers.command);
                if full_width_selectable_row(ui, selected, label).clicked() {
                    if extend {
                        actions.push(Action::AddLinkToSelection(link.id.to_string()));
                    } else {
                        actions.push(Action::SelectLink(link.id.to_string()));
                    }
                }
            }
        });
}

pub(crate) fn node_type_label(node_type: &NodeType) -> &'static str {
    match node_type {
        NodeType::Junction(_) => "Junction",
        NodeType::Tank(_) => "Tank",
        NodeType::Reservoir(_) => "Reservoir",
    }
}

pub(crate) fn link_type_label(link_type: &LinkType) -> &'static str {
    match link_type {
        LinkType::Pipe(_) => "Pipe",
        LinkType::Pump(_) => "Pump",
        LinkType::Valve(_) => "Valve",
    }
}

/// Head/pressure unit abbreviation for the network's configured unit system
/// ("ft" for US customary, "m" for SI) - used everywhere a head or pressure
/// value is displayed.
pub(crate) fn head_unit_label(network: &Network) -> &'static str {
    if network.options.unit_system.per_feet() == 1.0 {
        "ft"
    } else {
        "m"
    }
}

/// Plain length unit abbreviation ("ft"/"m") - same scale as
/// `head_unit_label` (elevation, tank levels/diameter, pipe length,
/// Darcy-Weisbach roughness); named separately so call sites read by
/// physical quantity rather than borrowing the "head" name.
pub(crate) fn length_unit_label(network: &Network) -> &'static str {
    head_unit_label(network)
}

/// Pipe/valve *diameter* unit ("in" for US, "mm" for SI) - diameters are
/// converted with their own extra /12 (in) or /1000 (mm) step on top of the
/// plain length conversion (see `Pipe`/`Valve::convert_from_standard`), so
/// they need a different label than a plain length even within the same
/// unit system. Note `Tank.diameter` is *not* one of these - it uses the
/// plain length conversion, so tanks use `length_unit_label` instead.
pub(crate) fn diameter_unit_label(network: &Network) -> &'static str {
    if network.options.unit_system == UnitSystem::US {
        "in"
    } else {
        "mm"
    }
}

/// Tank volume unit ("ft³"/"m³").
pub(crate) fn volume_unit_label(network: &Network) -> &'static str {
    if network.options.unit_system == UnitSystem::US {
        "ft³"
    } else {
        "m³"
    }
}

/// Pump power unit ("hp" for US, "kW" for SI).
pub(crate) fn power_unit_label(network: &Network) -> &'static str {
    if network.options.unit_system == UnitSystem::US {
        "hp"
    } else {
        "kW"
    }
}

/// Pipe roughness unit: a length, only for Darcy-Weisbach (`Pipe` treats it
/// as a length there, see `Pipe::convert_from_standard`); the
/// Hazen-Williams C-factor and Chezy-Manning coefficient are dimensionless
/// and never unit-converted, so they get no unit suffix.
pub(crate) fn roughness_unit_label(network: &Network) -> &'static str {
    if network.options.headloss_formula == HeadlossFormula::DarcyWeisbach {
        length_unit_label(network)
    } else {
        ""
    }
}

/// PRV/PSV/PBV valve setting unit. `Valve::convert_*` always converts
/// through PSI for a US network - hardcoded, independent of
/// `options.pressure_units` - and through meters for SI; this matches that
/// actual behavior exactly rather than a theoretically "nicer" but wrong
/// label derived from `pressure_units`.
pub(crate) fn valve_pressure_setting_unit_label(network: &Network) -> &'static str {
    if network.options.unit_system == UnitSystem::US {
        "psi"
    } else {
        "m"
    }
}

/// The unit for a valve's `setting` field, which depends on its
/// `ValveType`: a pressure for PRV/PSV/PBV, a flow for FCV, and
/// dimensionless (never unit-converted) for TCV/PCV/GPV.
pub(crate) fn valve_setting_unit_label(network: &Network, valve_type: &ValveType) -> String {
    match valve_type {
        ValveType::PRV | ValveType::PSV | ValveType::PBV => {
            valve_pressure_setting_unit_label(network).to_string()
        }
        ValveType::FCV => flow_unit_label(network),
        ValveType::TCV | ValveType::PCV | ValveType::GPV => String::new(),
    }
}

/// Flow/demand unit abbreviation (e.g. "GPM", "LPS") for the network's
/// configured flow units - used everywhere a flow or demand value is
/// displayed.
pub(crate) fn flow_unit_label(network: &Network) -> String {
    network.options.flow_units.to_string()
}

/// Emitter coefficient unit: flow per pressure to the `emitter_exponent`
/// power (`Q = K * P^n`), e.g. "GPM/psi^0.5" - built dynamically since the
/// exponent is itself a per-network option. Uses the same (hardcoded PSI
/// for US / meters for SI) pressure basis as `Junction::convert_*` actually
/// uses, not `options.pressure_units`.
pub(crate) fn emitter_unit_label(network: &Network) -> String {
    let flow = flow_unit_label(network);
    let pressure = valve_pressure_setting_unit_label(network);
    let exponent = network.options.emitter_exponent;
    format!("{flow}/{pressure}^{exponent}")
}

pub fn right_panel(ui: &mut egui::Ui, state: &AppState, actions: &mut Vec<Action>) {
    egui::Panel::right("right_panel")
        .resizable(true)
        .default_size(320.0)
        .show_inside(ui, |ui| {
            ui.heading("Properties");
            ui.add_space(4.0);

            egui::ScrollArea::vertical()
                .id_salt("properties_scroll")
                .show(ui, |ui| {
                    if let Some(id) = state.selection.as_single_node() {
                        node_inspector(ui, state, id, actions);
                    } else if let Some(id) = state.selection.as_single_link() {
                        link_inspector(ui, state, id, actions);
                    } else if state.selection.is_empty() {
                        ui.weak("Nothing selected.");
                        ui.label("Click a node or link, or drag a box to select several.");
                        ui.label("Hold Shift/Ctrl while clicking to extend the selection.");
                    } else {
                        bulk_inspector(ui, state, actions);
                    }
                });
        });
}

/// A labeled numeric drag field, showing `unit` alongside the label (e.g.
/// "Elevation (ft)") when non-empty - "show correct units for properties".
fn drag_field(ui: &mut egui::Ui, label: &str, unit: &str, value: &mut f64) -> bool {
    ui.horizontal(|ui| {
        if unit.is_empty() {
            ui.label(label);
        } else {
            ui.label(format!("{label} ({unit})"));
        }
        ui.add(egui::DragValue::new(value).speed(0.1))
    })
    .inner
    .changed()
}

/// Shows a compact "Simulation results (step N/M):" block for the currently
/// selected entity, right under its id/type header - so results sit next to
/// the properties they explain instead of only being visible as a canvas
/// color. `compute` returns the rows to show (label, formatted value), or
/// `None` when there's no result for this entity yet (no simulation has
/// run, or the network topology changed since the last run); a no-op in
/// that case.
fn simulation_result_section(
    ui: &mut egui::Ui,
    state: &AppState,
    compute: impl FnOnce(&AppState) -> Option<Vec<(&'static str, String)>>,
) {
    let Some(rows) = compute(state) else {
        return;
    };
    ui.separator();
    let step = state.report_step + 1;
    let total_steps = state.sim_results.as_ref().map(|r| r.heads.len()).unwrap_or(0);
    ui.label(format!("Simulation results (step {step}/{total_steps})"));
    egui::Grid::new("sim_result_grid").num_columns(2).show(ui, |ui| {
        for (label, value) in rows {
            ui.weak(label);
            ui.monospace(value);
            ui.end_row();
        }
    });
}

/// Combo box following the `Option<Option<Box<str>>>` tri-state convention
/// used by patterns/curves: returns `Some(new_value)` (itself an
/// `Option<Box<str>>`, `None` meaning "clear") when the user picks a
/// different entry, or `None` if nothing changed this frame.
fn pattern_or_curve_picker(
    ui: &mut egui::Ui,
    label: &str,
    id_salt: &str,
    current: Option<&str>,
    options: impl Iterator<Item = String>,
) -> Option<Option<Box<str>>> {
    let mut result = None;
    let selected_text = current.unwrap_or("(none)").to_string();
    egui::ComboBox::from_id_salt(id_salt)
        .selected_text(selected_text)
        .show_ui(ui, |ui| {
            if ui.selectable_label(current.is_none(), "(none)").clicked() && current.is_some() {
                result = Some(None);
            }
            for candidate in options {
                let selected = current == Some(candidate.as_str());
                if ui.selectable_label(selected, &candidate).clicked() && !selected {
                    result = Some(Some(candidate.as_str().into()));
                }
            }
        });
    ui.horizontal(|ui| {
        ui.weak(label);
    });
    result
}

fn node_inspector(ui: &mut egui::Ui, state: &AppState, id: &str, actions: &mut Vec<Action>) {
    let Some(&index) = state.network().node_map.get(id) else {
        return;
    };
    // A display-unit copy: everything read from `node` below (other than
    // `coordinates`, which is never unit-converted at all) is in the
    // network's configured display units, matching exactly what
    // `Action::UpdateNodes`/`UpdateJunctions`/etc. expect back - reusing
    // the model's own `UnitConversion` impls instead of re-deriving the
    // same conversion factors here.
    let mut node: Node = state.network().nodes[index].clone();
    node.convert_from_standard(&state.network().options);
    let head_unit = head_unit_label(state.network());

    ui.label(format!("Node: {}", node.id));
    ui.label(node_type_label(&node.node_type));

    simulation_result_section(ui, state, |state| {
        let result = state.node_result(id)?;
        let head_unit = head_unit_label(state.network());
        let flow_unit = flow_unit_label(state.network());
        Some(vec![
            ("Head", format!("{:.2} {head_unit}", result.head)),
            ("Pressure", format!("{:.2} {head_unit}", result.pressure)),
            ("Demand", format!("{:.3} {flow_unit}", result.demand)),
        ])
    });
    ui.separator();

    // --- generic fields, available on every node type ---
    let mut elevation = node.elevation;
    if drag_field(ui, "Elevation", head_unit, &mut elevation) {
        actions.push(Action::UpdateNodes {
            ids: vec![id.to_string()],
            update: NodeUpdate {
                elevation: Some(elevation),
                ..Default::default()
            },
        });
    }
    let (mut x, mut y) = node.coordinates.unwrap_or((0.0, 0.0));
    let (cx, cy) = ui
        .horizontal(|ui| {
            ui.label("Position");
            let cx = ui.add(egui::DragValue::new(&mut x).prefix("x: ").speed(0.5));
            let cy = ui.add(egui::DragValue::new(&mut y).prefix("y: ").speed(0.5));
            (cx.changed(), cy.changed())
        })
        .inner;
    if cx || cy {
        actions.push(Action::UpdateNodes {
            ids: vec![id.to_string()],
            update: NodeUpdate {
                coordinates: Some((x, y)),
                ..Default::default()
            },
        });
    }
    let mut disabled = node.disabled;
    if ui.checkbox(&mut disabled, "Disabled").changed() {
        actions.push(Action::UpdateNodes {
            ids: vec![id.to_string()],
            update: NodeUpdate {
                disabled: Some(disabled),
                ..Default::default()
            },
        });
    }
    ui.separator();

    match &node.node_type {
        NodeType::Junction(junction) => {
            let flow_unit = flow_unit_label(state.network());
            let mut base_demand = junction
                .demands
                .first()
                .map(|d| d.basedemand)
                .unwrap_or(0.0);
            if drag_field(ui, "Base demand", &flow_unit, &mut base_demand) {
                actions.push(Action::UpdateJunctions {
                    ids: vec![id.to_string()],
                    update: JunctionUpdate {
                        basedemand: Some(base_demand),
                        ..Default::default()
                    },
                });
            }
            let current_pattern = junction.demands.first().and_then(|d| d.pattern.as_deref());
            if let Some(pattern) = pattern_or_curve_picker(
                ui,
                "Demand pattern",
                "junction_pattern",
                current_pattern,
                state.network().patterns.iter().map(|p| p.id.to_string()),
            ) {
                actions.push(Action::UpdateJunctions {
                    ids: vec![id.to_string()],
                    update: JunctionUpdate {
                        pattern: Some(pattern),
                        ..Default::default()
                    },
                });
            }
            let emitter_unit = emitter_unit_label(state.network());
            let mut emitter = junction.emitter_coefficient;
            if drag_field(ui, "Emitter coeff.", &emitter_unit, &mut emitter) {
                actions.push(Action::UpdateJunctions {
                    ids: vec![id.to_string()],
                    update: JunctionUpdate {
                        emitter_coefficient: Some(emitter),
                        ..Default::default()
                    },
                });
            }
        }
        NodeType::Tank(tank) => {
            let volume_unit = volume_unit_label(state.network());
            let mut initial_level = tank.initial_level;
            if drag_field(ui, "Initial level", head_unit, &mut initial_level) {
                actions.push(Action::UpdateTanks {
                    ids: vec![id.to_string()],
                    update: TankUpdate {
                        initial_level: Some(initial_level),
                        ..Default::default()
                    },
                });
            }
            let mut min_level = tank.min_level;
            if drag_field(ui, "Min level", head_unit, &mut min_level) {
                actions.push(Action::UpdateTanks {
                    ids: vec![id.to_string()],
                    update: TankUpdate {
                        min_level: Some(min_level),
                        ..Default::default()
                    },
                });
            }
            let mut max_level = tank.max_level;
            if drag_field(ui, "Max level", head_unit, &mut max_level) {
                actions.push(Action::UpdateTanks {
                    ids: vec![id.to_string()],
                    update: TankUpdate {
                        max_level: Some(max_level),
                        ..Default::default()
                    },
                });
            }
            // Unlike Pipe/Valve diameter, `Tank.diameter` is converted with
            // the plain length factor (ft/m), not the extra in/mm step -
            // see `length_unit_label`'s doc comment.
            let mut diameter = tank.diameter;
            if drag_field(ui, "Diameter", head_unit, &mut diameter) {
                actions.push(Action::UpdateTanks {
                    ids: vec![id.to_string()],
                    update: TankUpdate {
                        diameter: Some(diameter),
                        ..Default::default()
                    },
                });
            }
            let mut min_volume = tank.min_volume;
            if drag_field(ui, "Min volume", volume_unit, &mut min_volume) {
                actions.push(Action::UpdateTanks {
                    ids: vec![id.to_string()],
                    update: TankUpdate {
                        min_volume: Some(min_volume),
                        ..Default::default()
                    },
                });
            }
            let mut overflow = tank.overflow;
            if ui.checkbox(&mut overflow, "Overflow").changed() {
                actions.push(Action::UpdateTanks {
                    ids: vec![id.to_string()],
                    update: TankUpdate {
                        overflow: Some(overflow),
                        ..Default::default()
                    },
                });
            }
        }
        NodeType::Reservoir(reservoir) => {
            if let Some(pattern) = pattern_or_curve_picker(
                ui,
                "Head pattern",
                "reservoir_head_pattern",
                reservoir.head_pattern.as_deref(),
                state.network().patterns.iter().map(|p| p.id.to_string()),
            ) {
                actions.push(Action::UpdateReservoirs {
                    ids: vec![id.to_string()],
                    update: ReservoirUpdate {
                        head_pattern: Some(pattern),
                        ..Default::default()
                    },
                });
            }
        }
    }

    ui.separator();
    if ui.button("Delete").clicked() {
        actions.push(Action::DeleteSelection);
    }
}

fn link_inspector(ui: &mut egui::Ui, state: &AppState, id: &str, actions: &mut Vec<Action>) {
    let Some(&index) = state.network().link_map.get(id) else {
        return;
    };
    // The original, standard-unit link - kept around only for the pipe
    // minor-loss un-normalization below, which needs the actual standard
    // (feet) diameter the stored `minor_loss` was normalized against, not
    // the display-unit one.
    let standard_link: &Link = &state.network().links[index];
    // A display-unit copy: see the equivalent comment in `node_inspector`.
    let mut link: Link = standard_link.clone();
    link.convert_from_standard(&state.network().options);

    ui.label(format!("Link: {}", link.id));
    ui.label(link_type_label(&link.link_type));

    simulation_result_section(ui, state, |state| {
        let result = state.link_result(id)?;
        let flow_unit = flow_unit_label(state.network());
        let direction = if result.flow < 0.0 { "end → start" } else { "start → end" };
        Some(vec![
            ("Flow", format!("{:.3} {flow_unit}", result.flow)),
            ("Direction", direction.to_string()),
        ])
    });
    ui.separator();

    // --- direction: which way is start, which way is end ---
    ui.label("Direction");
    ui.horizontal(|ui| {
        ui.monospace(link.start_node_id.as_ref());
        ui.label("→ (start to end) →");
        ui.monospace(link.end_node_id.as_ref());
    });
    if ui
        .button("⇄ Swap direction")
        .on_hover_text("Swap start/end nodes; flips the arrow drawn on the canvas.")
        .clicked()
    {
        let start = link.start_node_id.to_string();
        let end = link.end_node_id.to_string();
        // Reverse the display vertices too, or the polyline would still
        // read start(old A)->v1->..->vn->end(old B) with A/B swapped,
        // jumping straight across from the new start to whichever vertex
        // used to be nearest the old start - a crossed, wrong-looking path.
        let vertices = link.vertices.as_ref().map(|v| v.iter().rev().copied().collect());
        actions.push(Action::UpdateLinks {
            ids: vec![id.to_string()],
            update: LinkUpdate {
                start_node: Some(end.into()),
                end_node: Some(start.into()),
                vertices,
                ..Default::default()
            },
        });
    }
    ui.separator();

    let mut status = link.initial_status;
    egui::ComboBox::from_label("Status")
        .selected_text(status.to_string())
        .show_ui(ui, |ui| {
            for candidate in LINK_STATUSES {
                if ui
                    .selectable_value(&mut status, candidate, candidate.to_string())
                    .changed()
                {
                    actions.push(Action::UpdateLinks {
                        ids: vec![id.to_string()],
                        update: LinkUpdate {
                            initial_status: Some(status),
                            ..Default::default()
                        },
                    });
                }
            }
        });
    ui.separator();

    match &link.link_type {
        LinkType::Pipe(pipe) => {
            let length_unit = length_unit_label(state.network());
            let diameter_unit = diameter_unit_label(state.network());
            let roughness_unit = roughness_unit_label(state.network());
            let mut length = pipe.length;
            if drag_field(ui, "Length", length_unit, &mut length) {
                actions.push(Action::UpdatePipes {
                    ids: vec![id.to_string()],
                    update: PipeUpdate {
                        length: Some(length),
                        ..Default::default()
                    },
                });
            }
            let mut diameter = pipe.diameter;
            if drag_field(ui, "Diameter", diameter_unit, &mut diameter) {
                actions.push(Action::UpdatePipes {
                    ids: vec![id.to_string()],
                    update: PipeUpdate {
                        diameter: Some(diameter),
                        ..Default::default()
                    },
                });
            }
            let mut roughness = pipe.roughness;
            if drag_field(ui, "Roughness", roughness_unit, &mut roughness) {
                actions.push(Action::UpdatePipes {
                    ids: vec![id.to_string()],
                    update: PipeUpdate {
                        roughness: Some(roughness),
                        ..Default::default()
                    },
                });
            }
            // `pipe.minor_loss` is stored pre-normalized against diameter^4
            // (see `Network::update_pipe`); un-normalize it back to the
            // user-facing K coefficient for display, matching what
            // `PipeUpdate::minor_loss` expects on the way back in. This
            // needs the *standard*-unit (feet) diameter the stored value
            // was normalized against, not the display-unit one above.
            let standard_diameter = match &standard_link.link_type {
                LinkType::Pipe(standard_pipe) => standard_pipe.diameter,
                _ => 0.0,
            };
            let mut minor_loss = if standard_diameter.abs() > 1e-12 {
                pipe.minor_loss * standard_diameter.powi(4) / 0.02517
            } else {
                0.0
            };
            if drag_field(ui, "Minor loss (K)", "", &mut minor_loss) {
                actions.push(Action::UpdatePipes {
                    ids: vec![id.to_string()],
                    update: PipeUpdate {
                        minor_loss: Some(minor_loss),
                        ..Default::default()
                    },
                });
            }
            let mut check_valve = pipe.check_valve;
            if ui.checkbox(&mut check_valve, "Check valve").changed() {
                actions.push(Action::UpdatePipes {
                    ids: vec![id.to_string()],
                    update: PipeUpdate {
                        check_valve: Some(check_valve),
                        ..Default::default()
                    },
                });
            }
        }
        LinkType::Pump(pump) => {
            let power_unit = power_unit_label(state.network());
            let mut speed = pump.speed;
            if drag_field(ui, "Speed", "", &mut speed) {
                actions.push(Action::UpdatePumps {
                    ids: vec![id.to_string()],
                    update: PumpUpdate {
                        speed: Some(speed),
                        ..Default::default()
                    },
                });
            }
            let mut power = pump.power;
            if drag_field(ui, "Power", power_unit, &mut power) {
                actions.push(Action::UpdatePumps {
                    ids: vec![id.to_string()],
                    update: PumpUpdate {
                        power: Some(power),
                        ..Default::default()
                    },
                });
            }
            if let Some(curve) = pattern_or_curve_picker(
                ui,
                "Head curve",
                "pump_head_curve",
                pump.head_curve_id.as_deref(),
                state.network().curves.iter().map(|c| c.id.to_string()),
            ) {
                actions.push(Action::UpdatePumps {
                    ids: vec![id.to_string()],
                    update: PumpUpdate {
                        head_curve_id: Some(curve),
                        ..Default::default()
                    },
                });
            }
        }
        LinkType::Valve(valve) => {
            let diameter_unit = diameter_unit_label(state.network());
            let setting_unit = valve_setting_unit_label(state.network(), &valve.valve_type);
            let mut diameter = valve.diameter;
            if drag_field(ui, "Diameter", diameter_unit, &mut diameter) {
                actions.push(Action::UpdateValves {
                    ids: vec![id.to_string()],
                    update: ValveUpdate {
                        diameter: Some(diameter),
                        ..Default::default()
                    },
                });
            }
            let mut setting = valve.setting;
            if drag_field(ui, "Setting", &setting_unit, &mut setting) {
                actions.push(Action::UpdateValves {
                    ids: vec![id.to_string()],
                    update: ValveUpdate {
                        setting: Some(setting),
                        ..Default::default()
                    },
                });
            }
            let mut minor_loss = valve.minor_loss;
            if drag_field(ui, "Minor loss (K)", "", &mut minor_loss) {
                actions.push(Action::UpdateValves {
                    ids: vec![id.to_string()],
                    update: ValveUpdate {
                        minor_loss: Some(minor_loss),
                        ..Default::default()
                    },
                });
            }

            let mut valve_type = valve.valve_type.clone();
            egui::ComboBox::from_label("Valve type")
                .selected_text(valve_type.to_string())
                .show_ui(ui, |ui| {
                    for candidate in VALVE_TYPES {
                        if ui
                            .selectable_value(&mut valve_type, candidate.clone(), candidate.to_string())
                            .changed()
                        {
                            // Changing valve_type requires a fresh setting;
                            // reuse the current one as a sane default.
                            actions.push(Action::UpdateValves {
                                ids: vec![id.to_string()],
                                update: ValveUpdate {
                                    valve_type: Some(valve_type.clone()),
                                    setting: Some(valve.setting),
                                    ..Default::default()
                                },
                            });
                        }
                    }
                });

            if matches!(valve.valve_type, ValveType::GPV | ValveType::PCV) {
                if let Some(curve) = pattern_or_curve_picker(
                    ui,
                    "Valve curve",
                    "valve_curve",
                    valve.curve_id.as_deref(),
                    state.network().curves.iter().map(|c| c.id.to_string()),
                ) {
                    actions.push(Action::UpdateValves {
                        ids: vec![id.to_string()],
                        update: ValveUpdate {
                            curve_id: Some(curve),
                            ..Default::default()
                        },
                    });
                }
            }
        }
    }

    ui.separator();
    if ui.button("Delete").clicked() {
        actions.push(Action::DeleteSelection);
    }
}

/// Shown when the selection contains more than one node/link (of any mix of
/// types). Generic fields apply to every id in the selection at once;
/// type-specific fields only appear when every selected node/link shares the
/// same concrete type. Numeric fields here are "set all selected to..."
/// actions, not a reflection of a (possibly non-uniform) current value.
fn bulk_inspector(ui: &mut egui::Ui, state: &AppState, actions: &mut Vec<Action>) {
    let selection = state.selection.clone();
    ui.label(format!(
        "{} node(s), {} link(s) selected",
        selection.nodes.len(),
        selection.links.len()
    ));
    ui.weak("Bulk edit: fields below apply to every selected entity.");
    ui.separator();

    if !selection.nodes.is_empty() {
        ui.heading("Nodes");
        bulk_node_fields(ui, state, &selection.nodes, actions);
        ui.separator();
    }
    if !selection.links.is_empty() {
        ui.heading("Links");
        bulk_link_fields(ui, state, &selection.links, actions);
        ui.separator();
    }

    if ui.button("Delete selection").clicked() {
        actions.push(Action::DeleteSelection);
    }
}

fn bulk_node_fields(ui: &mut egui::Ui, state: &AppState, ids: &[String], actions: &mut Vec<Action>) {
    let head_unit = head_unit_label(state.network());
    let mut elevation = 0.0;
    if drag_field(ui, "Set elevation (all)", head_unit, &mut elevation) {
        actions.push(Action::UpdateNodes {
            ids: ids.to_vec(),
            update: NodeUpdate {
                elevation: Some(elevation),
                ..Default::default()
            },
        });
    }
    ui.horizontal(|ui| {
        if ui.button("Enable all").clicked() {
            actions.push(Action::UpdateNodes {
                ids: ids.to_vec(),
                update: NodeUpdate {
                    disabled: Some(false),
                    ..Default::default()
                },
            });
        }
        if ui.button("Disable all").clicked() {
            actions.push(Action::UpdateNodes {
                ids: ids.to_vec(),
                update: NodeUpdate {
                    disabled: Some(true),
                    ..Default::default()
                },
            });
        }
    });

    let network = state.network();
    let mut types = ids
        .iter()
        .filter_map(|id| network.node_map.get(id.as_str()))
        .map(|&i| &network.nodes[i].node_type);
    let Some(first) = types.next() else { return };
    let same_type = types.all(|t| std::mem::discriminant(t) == std::mem::discriminant(first));
    if !same_type {
        ui.weak("Mixed node types: only the generic fields above apply.");
        return;
    }

    match first {
        NodeType::Junction(_) => {
            let flow_unit = flow_unit_label(state.network());
            let emitter_unit = emitter_unit_label(state.network());
            let mut base_demand = 0.0;
            if drag_field(ui, "Set base demand (all)", &flow_unit, &mut base_demand) {
                actions.push(Action::UpdateJunctions {
                    ids: ids.to_vec(),
                    update: JunctionUpdate {
                        basedemand: Some(base_demand),
                        ..Default::default()
                    },
                });
            }
            let mut emitter = 0.0;
            if drag_field(ui, "Set emitter coeff. (all)", &emitter_unit, &mut emitter) {
                actions.push(Action::UpdateJunctions {
                    ids: ids.to_vec(),
                    update: JunctionUpdate {
                        emitter_coefficient: Some(emitter),
                        ..Default::default()
                    },
                });
            }
        }
        NodeType::Tank(_) => {
            let mut initial_level = 0.0;
            if drag_field(ui, "Set initial level (all)", head_unit, &mut initial_level) {
                actions.push(Action::UpdateTanks {
                    ids: ids.to_vec(),
                    update: TankUpdate {
                        initial_level: Some(initial_level),
                        ..Default::default()
                    },
                });
            }
            let mut min_level = 0.0;
            if drag_field(ui, "Set min level (all)", head_unit, &mut min_level) {
                actions.push(Action::UpdateTanks {
                    ids: ids.to_vec(),
                    update: TankUpdate {
                        min_level: Some(min_level),
                        ..Default::default()
                    },
                });
            }
            let mut max_level = 0.0;
            if drag_field(ui, "Set max level (all)", head_unit, &mut max_level) {
                actions.push(Action::UpdateTanks {
                    ids: ids.to_vec(),
                    update: TankUpdate {
                        max_level: Some(max_level),
                        ..Default::default()
                    },
                });
            }
            // Plain length unit, not `diameter_unit_label` - see
            // `length_unit_label`'s doc comment (Tank is not Pipe/Valve).
            let mut diameter = 0.0;
            if drag_field(ui, "Set diameter (all)", head_unit, &mut diameter) {
                actions.push(Action::UpdateTanks {
                    ids: ids.to_vec(),
                    update: TankUpdate {
                        diameter: Some(diameter),
                        ..Default::default()
                    },
                });
            }
        }
        NodeType::Reservoir(_) => {
            let mut head = 0.0;
            if drag_field(ui, "Set elevation/head (all)", head_unit, &mut head) {
                actions.push(Action::UpdateReservoirs {
                    ids: ids.to_vec(),
                    update: ReservoirUpdate {
                        elevation: Some(head),
                        ..Default::default()
                    },
                });
            }
        }
    }
}

fn bulk_link_fields(ui: &mut egui::Ui, state: &AppState, ids: &[String], actions: &mut Vec<Action>) {
    let mut status = LinkStatus::Open;
    egui::ComboBox::from_label("Set status (all)")
        .selected_text(status.to_string())
        .show_ui(ui, |ui| {
            for candidate in LINK_STATUSES {
                if ui
                    .selectable_value(&mut status, candidate, candidate.to_string())
                    .changed()
                {
                    actions.push(Action::UpdateLinks {
                        ids: ids.to_vec(),
                        update: LinkUpdate {
                            initial_status: Some(candidate),
                            ..Default::default()
                        },
                    });
                }
            }
        });

    let network = state.network();
    let mut types = ids
        .iter()
        .filter_map(|id| network.link_map.get(id.as_str()))
        .map(|&i| &network.links[i].link_type);
    let Some(first) = types.next() else { return };
    let same_type = types.all(|t| std::mem::discriminant(t) == std::mem::discriminant(first));
    if !same_type {
        ui.weak("Mixed link types: only status above applies.");
        return;
    }

    match first {
        LinkType::Pipe(_) => {
            let diameter_unit = diameter_unit_label(state.network());
            let roughness_unit = roughness_unit_label(state.network());
            let length_unit = length_unit_label(state.network());
            let mut diameter = 0.0;
            if drag_field(ui, "Set diameter (all)", diameter_unit, &mut diameter) {
                actions.push(Action::UpdatePipes {
                    ids: ids.to_vec(),
                    update: PipeUpdate {
                        diameter: Some(diameter),
                        ..Default::default()
                    },
                });
            }
            let mut roughness = 0.0;
            if drag_field(ui, "Set roughness (all)", roughness_unit, &mut roughness) {
                actions.push(Action::UpdatePipes {
                    ids: ids.to_vec(),
                    update: PipeUpdate {
                        roughness: Some(roughness),
                        ..Default::default()
                    },
                });
            }
            let mut length = 0.0;
            if drag_field(ui, "Set length (all)", length_unit, &mut length) {
                actions.push(Action::UpdatePipes {
                    ids: ids.to_vec(),
                    update: PipeUpdate {
                        length: Some(length),
                        ..Default::default()
                    },
                });
            }
        }
        LinkType::Pump(_) => {
            let mut speed = 0.0;
            if drag_field(ui, "Set speed (all)", "", &mut speed) {
                actions.push(Action::UpdatePumps {
                    ids: ids.to_vec(),
                    update: PumpUpdate {
                        speed: Some(speed),
                        ..Default::default()
                    },
                });
            }
        }
        LinkType::Valve(valve) => {
            let setting_unit = valve_setting_unit_label(state.network(), &valve.valve_type);
            let mut setting = 0.0;
            if drag_field(ui, "Set setting (all)", &setting_unit, &mut setting) {
                actions.push(Action::UpdateValves {
                    ids: ids.to_vec(),
                    update: ValveUpdate {
                        setting: Some(setting),
                        ..Default::default()
                    },
                });
            }
        }
    }
}

pub fn bottom_toolbar(ctx: &egui::Context, state: &AppState, actions: &mut Vec<Action>) {
    egui::Area::new(egui::Id::new("bottom_toolbar"))
        .anchor(egui::Align2::CENTER_BOTTOM, egui::vec2(0.0, -12.0))
        .order(egui::Order::Foreground)
        .show(ctx, |ui| {
            egui::Frame::popup(ui.style()).show(ui, |ui| {
                ui.horizontal(|ui| {
                    tool_button(ui, state, actions, Tool::Select, "Select", Some("Hold Space/middle mouse button to pan"));
                    ui.separator();
                    tool_button(ui, state, actions, Tool::AddJunction, "+ Junction", Some("J"));
                    tool_button(ui, state, actions, Tool::AddTank, "+ Tank", Some("T"));
                    tool_button(ui, state, actions, Tool::AddReservoir, "+ Reservoir", Some("R"));
                    tool_button(
                        ui,
                        state,
                        actions,
                        Tool::AddLink { kind: NewLinkKind::Pipe },
                        "+ Pipe",
                        Some("P"),
                    );
                    tool_button(
                        ui,
                        state,
                        actions,
                        Tool::AddLink { kind: NewLinkKind::Pump },
                        "+ Pump",
                        Some("Shift+P"),
                    );
                    valve_tool_control(ui, state, actions);
                    ui.separator();

                    if ui
                        .add_enabled(!state.sim_running, egui::Button::new("Run"))
                        .on_hover_text("Ctrl+Space")
                        .clicked()
                    {
                        actions.push(Action::RunSimulation { parallel: false });
                    }
                    if state.sim_running {
                        let p = state.sim_progress.unwrap_or(0.0);
                        ui.add(
                            egui::ProgressBar::new(p)
                                .desired_width(120.0)
                                .show_percentage(),
                        );
                    }
                    if let Some(err) = &state.sim_error {
                        ui.colored_label(egui::Color32::from_rgb(230, 90, 90), err);
                    }
                    if let Some(results) = &state.sim_results {
                        let max_step = results.heads.len().saturating_sub(1);
                        let mut step = state.report_step;
                        ui.separator();
                        if ui
                            .add(egui::Slider::new(&mut step, 0..=max_step).text("Report step"))
                            .changed()
                        {
                            actions.push(Action::SetReportStep(step));
                        }
                    }
                    ui.separator();
                    ui.weak("Ctrl+P: quick actions & search");
                });
            });
        });
}

/// The "+ Valve" tool button plus an inline combo to pick which `ValveType`
/// the connector tool will create (the toolbar previously hard-coded TCV).
fn valve_tool_control(ui: &mut egui::Ui, state: &AppState, actions: &mut Vec<Action>) {
    let active_valve_type = match &state.tool {
        Tool::AddLink { kind: NewLinkKind::Valve(vt) } => Some(vt.clone()),
        _ => None,
    };
    let selected = active_valve_type.is_some();
    let label = match &active_valve_type {
        Some(vt) => format!("+ Valve ({vt})"),
        None => "+ Valve".to_string(),
    };
    if ui.selectable_label(selected, label).on_hover_text("V").clicked() {
        let vt = active_valve_type.clone().unwrap_or(ValveType::TCV);
        actions.push(Action::SetTool(Tool::AddLink {
            kind: NewLinkKind::Valve(vt),
        }));
    }
    egui::ComboBox::from_id_salt("valve_type_picker")
        .selected_text(
            active_valve_type
                .map(|v| v.to_string())
                .unwrap_or_else(|| "type".to_string()),
        )
        .width(60.0)
        .show_ui(ui, |ui| {
            for candidate in VALVE_TYPES {
                if ui.selectable_label(false, candidate.to_string()).clicked() {
                    actions.push(Action::SetTool(Tool::AddLink {
                        kind: NewLinkKind::Valve(candidate),
                    }));
                }
            }
        });
}

fn tool_button(
    ui: &mut egui::Ui,
    state: &AppState,
    actions: &mut Vec<Action>,
    tool: Tool,
    label: &str,
    shortcut: Option<&str>,
) {
    let selected = state.tool == tool;
    let response = ui.selectable_label(selected, label);
    let response = if let Some(sc) = shortcut {
        response.on_hover_text(sc)
    } else {
        response
    };
    if response.clicked() {
        actions.push(Action::SetTool(tool));
    }
}
