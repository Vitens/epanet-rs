//! FFI accessors for simple controls: `EN_addcontrol`, `EN_getcontrol`, etc.

use crate::ffi::enums::ControlType;
use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::ffi::util::write_out;

use crate::model::control::{Control, ControlCondition};
use crate::model::link::{LinkStatus, LinkType};
use crate::model::network::Network;
use crate::model::node::NodeType;
use crate::model::units::UnitConversion;
use crate::simulation::Simulation;

use std::os::raw::{c_double, c_int};

/// The EPANET control type that matches a stored condition.
fn control_type_of(condition: &ControlCondition) -> ControlType {
    match condition {
        ControlCondition::LowLevel { .. } | ControlCondition::LowPressure { .. } => {
            ControlType::LowLevel
        }
        ControlCondition::HighLevel { .. } | ControlCondition::HighPressure { .. } => {
            ControlType::HiLevel
        }
        ControlCondition::Time { .. } => ControlType::Timer,
        ControlCondition::ClockTime { .. } => ControlType::TimeOfDay,
    }
}

/// Index of the node a level/pressure control watches, if any.
fn control_node_index(condition: &ControlCondition) -> Option<usize> {
    match condition {
        ControlCondition::LowLevel { tank_index, .. }
        | ControlCondition::HighLevel { tank_index, .. } => Some(*tank_index),
        ControlCondition::LowPressure { node_index, .. }
        | ControlCondition::HighPressure { node_index, .. } => Some(*node_index),
        _ => None,
    }
}

/// Builds a control in *user* units and converts it to internal units.
///
/// `level` is a tank level, a nodal pressure or a time in seconds depending on
/// `control_type` and the kind of node referenced by `node_index`.
fn build_control(
    network: &Network,
    control_type: ControlType,
    link_index: usize,
    setting: f64,
    node_index: usize,
    level: f64,
) -> Result<Control, ErrorCode> {
    let link = network
        .links
        .get(link_index)
        .ok_or(ErrorCode::UndefinedLink)?;

    let condition = match control_type {
        ControlType::Timer => ControlCondition::Time {
            seconds: level.max(0.0) as usize,
        },
        ControlType::TimeOfDay => ControlCondition::ClockTime {
            seconds: level.max(0.0) as usize,
        },
        ControlType::LowLevel | ControlType::HiLevel => {
            let node = network
                .nodes
                .get(node_index)
                .ok_or(ErrorCode::UndefinedNode)?;
            let is_tank = matches!(node.node_type, NodeType::Tank(_));
            let high = control_type == ControlType::HiLevel;
            match (is_tank, high) {
                (true, false) => ControlCondition::LowLevel {
                    tank_index: node_index,
                    target: level,
                },
                (true, true) => ControlCondition::HighLevel {
                    tank_index: node_index,
                    target: level,
                },
                (false, false) => ControlCondition::LowPressure {
                    node_index,
                    target: level,
                },
                (false, true) => ControlCondition::HighPressure {
                    node_index,
                    target: level,
                },
            }
        }
    };

    // a pipe can only be switched on or off, pumps and valves take a setting
    let (status, setting) = match &link.link_type {
        LinkType::Pipe(_) => {
            if setting == 0.0 {
                (Some(LinkStatus::Closed), None)
            } else {
                (Some(LinkStatus::Open), None)
            }
        }
        LinkType::Pump(_) => {
            if setting == 0.0 {
                (Some(LinkStatus::Closed), None)
            } else {
                (None, Some(setting))
            }
        }
        LinkType::Valve(_) => (None, Some(setting)),
    };

    let mut control = Control {
        condition,
        link_id: link.id.clone(),
        setting,
        status,
    };

    control.convert_to_standard(&network.options);
    control.convert_setting_to_standard(link, &network.options);

    Ok(control)
}

/// Reads back a control in user units.
fn read_control(
    simulation: &Simulation,
    index: usize,
) -> Result<(ControlType, usize, f64, usize, f64), ErrorCode> {
    let network = &simulation.network;
    let control = network
        .controls
        .get(index)
        .ok_or(ErrorCode::NonexistentControl)?;

    let link_index = *network
        .link_map
        .get(&control.link_id)
        .ok_or(ErrorCode::UndefinedLink)?;

    // work on a copy so the stored control keeps its internal units
    let mut control = control.clone();
    control.convert_setting_from_standard(&network.links[link_index], &network.options);
    control.convert_from_standard(&network.options);

    let setting = match (control.setting, control.status) {
        (Some(setting), _) => setting,
        (None, Some(LinkStatus::Closed) | None | Some(LinkStatus::FixedClosed)) => 0.0,
        (None, Some(_)) => 1.0,
    };

    let node_index = control_node_index(&control.condition);

    let level = match &control.condition {
        ControlCondition::LowLevel { target, .. }
        | ControlCondition::HighLevel { target, .. }
        | ControlCondition::LowPressure { target, .. }
        | ControlCondition::HighPressure { target, .. } => *target,
        ControlCondition::Time { seconds } | ControlCondition::ClockTime { seconds } => {
            *seconds as f64
        }
    };

    Ok((
        control_type_of(&control.condition),
        link_index + 1,
        setting,
        node_index.map(|i| i + 1).unwrap_or(0),
        level,
    ))
}

/// Adds a simple control to the project.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_index` must either be null or point to a writable `c_int`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_addcontrol(
    ph: *mut Project,
    control_type: c_int,
    link_index: c_int,
    setting: c_double,
    node_index: c_int,
    level: c_double,
    out_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let Some(control_type) = ControlType::from_repr(control_type) else {
        return ErrorCode::InvalidParameterCode;
    };

    if link_index < 1 {
        return ErrorCode::UndefinedLink;
    }

    let control = match build_control(
        &simulation.network,
        control_type,
        (link_index - 1) as usize,
        setting,
        (node_index - 1).max(0) as usize,
        level,
    ) {
        Ok(control) => control,
        Err(code) => return code,
    };

    simulation.network.controls.push(control);
    simulation.network.properties_version += 1;

    unsafe { write_out(out_index, simulation.network.controls.len() as c_int) };

    ErrorCode::Ok
}

/// Deletes a simple control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_deletecontrol(ph: *mut Project, index: c_int) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let index = (index - 1) as usize;
    if index >= simulation.network.controls.len() {
        return ErrorCode::NonexistentControl;
    }

    simulation.network.controls.remove(index);
    simulation.network.properties_version += 1;

    ErrorCode::Ok
}

/// Retrieves the properties of a simple control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// Every out-parameter must either be null or point to writable memory.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getcontrol(
    ph: *mut Project,
    index: c_int,
    out_type: *mut c_int,
    out_link_index: *mut c_int,
    out_setting: *mut c_double,
    out_node_index: *mut c_int,
    out_level: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    if index < 1 {
        return ErrorCode::NonexistentControl;
    }

    match read_control(simulation, (index - 1) as usize) {
        Ok((control_type, link_index, setting, node_index, level)) => unsafe {
            write_out(out_type, control_type as c_int);
            write_out(out_link_index, link_index as c_int);
            write_out(out_setting, setting as c_double);
            write_out(out_node_index, node_index as c_int);
            write_out(out_level, level as c_double);
            ErrorCode::Ok
        },
        Err(code) => code,
    }
}

/// Sets the properties of an existing simple control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setcontrol(
    ph: *mut Project,
    index: c_int,
    control_type: c_int,
    link_index: c_int,
    setting: c_double,
    node_index: c_int,
    level: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let index = (index - 1) as usize;
    if index >= simulation.network.controls.len() {
        return ErrorCode::NonexistentControl;
    }

    let Some(control_type) = ControlType::from_repr(control_type) else {
        return ErrorCode::InvalidParameterCode;
    };

    if link_index < 1 {
        return ErrorCode::UndefinedLink;
    }

    let control = match build_control(
        &simulation.network,
        control_type,
        (link_index - 1) as usize,
        setting,
        (node_index - 1).max(0) as usize,
        level,
    ) {
        Ok(control) => control,
        Err(code) => return code,
    };

    simulation.network.controls[index] = control;
    simulation.network.properties_version += 1;

    ErrorCode::Ok
}

/// Gets the enabled/disabled status of a simple control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getcontrolenabled(
    _ph: *mut Project,
    _index: c_int,
    _out_enabled: *mut c_int,
) -> ErrorCode {
    // TODO: controls in epanet-rs are always active; there is no enabled flag
    // to report.
    ErrorCode::NotImplemented
}

/// Enables or disables a simple control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setcontrolenabled(
    _ph: *mut Project,
    _index: c_int,
    _enabled: c_int,
) -> ErrorCode {
    // TODO: controls in epanet-rs are always active; use EN_deletecontrol to
    // remove one instead.
    ErrorCode::NotImplemented
}
