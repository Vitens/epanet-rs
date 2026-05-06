//! FFI link accessors: `EN_addlink`, `EN_getlinkvalue`, `EN_setlinkvalue`, etc.

use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::model::link::LinkType;
use crate::model::valve::ValveType;

use crate::constants::MperFT;
use crate::ffi::enums::LinkProperty;
use crate::ffi::enums::LinkType as ENLinkType;
use crate::ffi::enums::MISSING_VALUE;

use crate::model::network::modify::{LinkUpdate, PipeUpdate, PumpUpdate, ValveUpdate};
use crate::model::network::modify::{PipeData, PumpData, ValveData};

use crate::model::link::LinkStatus;
use crate::model::options::HeadlossFormula;
use crate::model::units::UnitSystem;

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_double, c_int};

/// Delete the link from the network
/// TODO: Implement action code
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_deletelink(
    ph: *mut Project,
    index: c_int,
    action_code: c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let link_id = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link.id.clone(),
        None => return ErrorCode::UndefinedLink,
    };

    let unconditional = match action_code {
        0 => true,
        1 => false,
        _ => return ErrorCode::InvalidParameterCode,
    };

    match simulation.network.remove_link(&link_id, unconditional) {
        Ok(_) => ErrorCode::Ok,
        Err(_) => ErrorCode::DeleteNodeOrLinkInControl,
    }
}

/// Gets the index of a node given its ID name.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
/// `out_index` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getlinkindex(
    ph: *mut Project,
    id: *const c_char,
    out_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let c_str = unsafe { CStr::from_ptr(id) };
    let link_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    // get the link index from the network
    let link_index = match simulation.network.link_map.get(link_id) {
        Some(&index) => index,
        None => return ErrorCode::UndefinedLink,
    };

    // EPANET indexes from 1, so we need to add 1 to the index
    unsafe { *out_index = (link_index + 1) as c_int };

    ErrorCode::Ok
}

// Gets the ID name of a node given its index.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_id` must point to a buffer large enough for the result string including NUL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getlinkid(
    ph: *mut Project,
    index: c_int,
    out_id: *mut c_char,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let link_id = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link.id.as_ref(),
        None => return ErrorCode::UndefinedLink,
    };

    let c_str = CString::new(link_id).unwrap();

    unsafe {
        std::ptr::copy_nonoverlapping(c_str.as_ptr(), out_id, c_str.as_bytes_with_nul().len());
    }

    ErrorCode::Ok
}

// Sets the ID name of a node given its index.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setlinkid(
    ph: *mut Project,
    index: c_int,
    id: *const c_char,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    // EPANET indexes from 1, so we need to subtract 1 from the index
    let link = match simulation.network.links.get_mut((index - 1) as usize) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };

    let c_str = unsafe { CStr::from_ptr(id) };

    let new_link_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    // check if the new node id is already in use
    if simulation.network.link_map.contains_key(new_link_id) {
        return ErrorCode::DuplicateId;
    }

    // remove the old link id from the link map
    simulation.network.link_map.remove(&link.id);

    // update the link id
    link.id = new_link_id.into();

    // update the link map
    simulation
        .network
        .link_map
        .insert(new_link_id.into(), index as usize);

    ErrorCode::Ok
}

// Get the link type given its index.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_type` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getlinktype(
    ph: *mut Project,
    index: c_int,
    out_type: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let link = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };

    let link_type_int = match &link.link_type {
        LinkType::Pipe(pipe) => match pipe.check_valve {
            true => ENLinkType::CVPipe as i32,
            false => ENLinkType::Pipe as i32,
        },
        LinkType::Pump(_) => ENLinkType::Pump as i32,
        LinkType::Valve(valve) => match valve.valve_type {
            ValveType::PRV => ENLinkType::PRV as i32,
            ValveType::PSV => ENLinkType::PSV as i32,
            ValveType::PBV => ENLinkType::PBV as i32,
            ValveType::FCV => ENLinkType::FCV as i32,
            ValveType::GPV => ENLinkType::GPV as i32,
            ValveType::TCV => ENLinkType::TCV as i32,
            ValveType::PCV => ENLinkType::PCV as i32,
        },
    };

    unsafe { *out_type = link_type_int };

    ErrorCode::Ok
}

// Get node index of the start and end nodes of a link given its index.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_start_node` must be a valid non-null writable pointer.
/// `out_end_node` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getlinknodes(
    ph: *mut Project,
    index: c_int,
    out_start_node: *mut c_int,
    out_end_node: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let link = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };

    unsafe { *out_start_node = (link.start_node + 1) as c_int };
    unsafe { *out_end_node = (link.end_node + 1) as c_int };

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setlinknodes(
    ph: *mut Project,
    index: c_int,
    start_node: c_int,
    end_node: c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    // check if the node ids are not the same
    if start_node == end_node {
        return ErrorCode::LinkSameStartEnd;
    }

    let link_id = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link.id.clone(),
        None => return ErrorCode::UndefinedLink,
    };

    // lookup start/end node ids
    let start_node_id = match simulation.network.nodes.get((start_node - 1) as usize) {
        Some(node) => node.id.clone(),
        None => return ErrorCode::UndefinedNode,
    };

    let end_node_id = match simulation.network.nodes.get((end_node - 1) as usize) {
        Some(node) => node.id.clone(),
        None => return ErrorCode::UndefinedNode,
    };

    match simulation.network.update_link(
        &link_id,
        &LinkUpdate {
            start_node: Some(start_node_id),
            end_node: Some(end_node_id),
            ..Default::default()
        },
    ) {
        Ok(_) => ErrorCode::Ok,
        Err(_) => ErrorCode::InvalidParameterCode,
    };

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_value` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getlinkvalue(
    ph: *mut Project,
    index: c_int,
    property: c_int,
    out_value: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let index = (index - 1) as usize;

    let link = match simulation.network.links.get(index) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };

    let options = &simulation.network.options;

    let property = match LinkProperty::from_repr(property) {
        Some(property) => property,
        None => return ErrorCode::InvalidParameterCode,
    };

    let value = match property {
        LinkProperty::Diameter => {
            // convert from standard units (Ft) to in (US) or mm (SI)
            let conversion = match options.unit_system {
                UnitSystem::US => 12.0,         // ft to in
                UnitSystem::SI => MperFT * 1e3, // ft to mm
            };

            match &link.link_type {
                LinkType::Pipe(pipe) => pipe.diameter * conversion,
                LinkType::Valve(valve) => valve.diameter * conversion,
                LinkType::Pump(_) => 0.0,
            }
        }
        LinkProperty::Length => match &link.link_type {
            LinkType::Pipe(pipe) => pipe.length * options.unit_system.per_feet(),
            _ => 0.0,
        },
        LinkProperty::Roughness => match &link.link_type {
            LinkType::Pipe(pipe) => {
                if pipe.headloss_formula == HeadlossFormula::DarcyWeisbach {
                    pipe.roughness * options.unit_system.per_feet()
                } else {
                    pipe.roughness
                }
            }
            _ => 0.0,
        },
        LinkProperty::MinorLoss => match &link.link_type {
            LinkType::Pipe(pipe) => pipe.minor_loss * pipe.diameter.powi(4) / 0.02517,
            LinkType::Valve(valve) => valve.minor_loss * valve.diameter.powi(4) / 0.02517,
            _ => 0.0,
        },
        LinkProperty::HeadLoss => {
            if let Some(state) = simulation.solved_state() {
                let headloss = state.heads[link.start_node] - state.heads[link.end_node];
                match &link.link_type {
                    LinkType::Pipe(_) | LinkType::Valve(_) => {
                        headloss.abs() * options.unit_system.per_feet()
                    }
                    LinkType::Pump(_) => headloss,
                }
            } else {
                0.0
            }
        }
        LinkProperty::Status => {
            if let Some(state) = simulation.solved_state() {
                match state.statuses[index] {
                    LinkStatus::Xhead => 2.0, // pump closed
                    LinkStatus::Closed => 0.0,
                    LinkStatus::FixedClosed => 0.0,
                    _ => 1.0,
                }
            } else {
                1.0
            }
        }
        LinkProperty::InitStatus => match &link.initial_status {
            LinkStatus::Closed => 0.0,
            _ => 1.0,
        },
        LinkProperty::Setting | LinkProperty::InitSetting => match &link.link_type {
            LinkType::Pipe(_) => 0.0,
            LinkType::Valve(valve) => match valve.valve_type {
                ValveType::FCV => valve.setting * options.flow_units.per_cfs(),
                ValveType::TCV => valve.setting,
                ValveType::PCV => valve.setting,
                ValveType::PRV => {
                    let downstream_node = simulation.network.nodes.get(link.end_node).unwrap();
                    (valve.setting + downstream_node.elevation) * options.pressure_units.per_feet()
                }
                ValveType::PSV => {
                    let upstream_node = simulation.network.nodes.get(link.start_node).unwrap();
                    (valve.setting + upstream_node.elevation) * options.pressure_units.per_feet()
                }
                _ => 0.0,
            },
            LinkType::Pump(pump) => pump.speed,
        },
        // TODO: Implement KBulk and KWall
        LinkProperty::KBulk => 0.0,
        LinkProperty::KWall => 0.0,

        // Flow, Velocity, Headloss (dynamic properties) are only available if the state is available and solved, otherwise return MISSING_VALUE
        LinkProperty::Flow => simulation.solved_state().map_or(MISSING_VALUE, |state| {
            state.flows[index] * options.flow_units.per_cfs()
        }),
        LinkProperty::Velocity => simulation.solved_state().map_or(MISSING_VALUE, |state| {
            let flow = state.flows[index].abs();
            match &link.link_type {
                LinkType::Pipe(pipe) => {
                    flow / (pipe.diameter.powi(2) * std::f64::consts::PI / 4.0)
                        * options.unit_system.per_feet()
                }
                LinkType::Valve(valve) => {
                    flow / (valve.diameter.powi(2) * std::f64::consts::PI / 4.0)
                        * options.unit_system.per_feet()
                }
                LinkType::Pump(_) => 0.0,
            }
        }),

        // TODO: expose other link properties
        _ => 0.0,
    };

    unsafe { *out_value = value };

    ErrorCode::Ok
}

// Set the property value of a link
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setlinkvalue(
    ph: *mut Project,
    index: c_int,
    property: c_int,
    value: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let index = (index - 1) as usize;

    let link = match simulation.network.links.get(index) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };
    let link_id = link.id.clone();

    let property = match LinkProperty::from_repr(property) {
        Some(property) => property,
        None => return ErrorCode::InvalidParameterCode,
    };

    let result = match property {
        LinkProperty::Diameter => match &link.link_type {
            LinkType::Pipe(_) => simulation.network.update_pipe(
                &link_id,
                &PipeUpdate {
                    diameter: Some(value),
                    ..Default::default()
                },
            ),
            LinkType::Valve(_) => simulation.network.update_valve(
                &link_id,
                &ValveUpdate {
                    diameter: Some(value),
                    ..Default::default()
                },
            ),
            LinkType::Pump(_) => {
                return ErrorCode::InvalidParameterCode;
            }
        },
        LinkProperty::Length => simulation.network.update_pipe(
            &link_id,
            &PipeUpdate {
                length: Some(value),
                ..Default::default()
            },
        ),
        LinkProperty::Roughness => simulation.network.update_pipe(
            &link_id,
            &PipeUpdate {
                roughness: Some(value),
                ..Default::default()
            },
        ),
        LinkProperty::MinorLoss => match &link.link_type {
            LinkType::Pipe(_) => simulation.network.update_pipe(
                &link_id,
                &PipeUpdate {
                    minor_loss: Some(value),
                    ..Default::default()
                },
            ),
            LinkType::Valve(_) => simulation.network.update_valve(
                &link_id,
                &ValveUpdate {
                    minor_loss: Some(value),
                    ..Default::default()
                },
            ),
            LinkType::Pump(_) => {
                return ErrorCode::InvalidParameterCode;
            }
        },
        LinkProperty::InitStatus | LinkProperty::Status => {
            let link_status = match value {
                0.0 => LinkStatus::Closed,
                1.0 => LinkStatus::Open,
                _ => return ErrorCode::InvalidParameterCode,
            };

            // if the link is a pipe with check_valve, return an error
            if let LinkType::Pipe(pipe) = &link.link_type
                && pipe.check_valve
            {
                return ErrorCode::IllegalValveControl;
            }
            simulation.network.update_link(
                &link_id,
                &LinkUpdate {
                    initial_status: Some(link_status),
                    ..Default::default()
                },
            )
        }
        LinkProperty::InitSetting | LinkProperty::Setting => {
            match &link.link_type {
                // for pipes, the setting is the roughness
                LinkType::Pipe(_) => simulation.network.update_pipe(
                    &link_id,
                    &PipeUpdate {
                        roughness: Some(value),
                        ..Default::default()
                    },
                ),
                // for valves, the setting is the valve setting
                LinkType::Valve(_) => simulation.network.update_valve(
                    &link_id,
                    &ValveUpdate {
                        setting: Some(value),
                        ..Default::default()
                    },
                ),
                // for pumps, the setting is the pump speed
                LinkType::Pump(_) => simulation.network.update_pump(
                    &link_id,
                    &PumpUpdate {
                        speed: Some(value),
                        ..Default::default()
                    },
                ),
            }
        }

        LinkProperty::GpvCurve | LinkProperty::PcvCurve => {
            let curve_id = match simulation.network.curves.get(value as usize - 1) {
                Some(curve) => curve.id.clone(),
                None => return ErrorCode::UndefinedCurve,
            };

            simulation.network.update_valve(
                &link_id,
                &ValveUpdate {
                    curve_id: Some(Some(curve_id)),
                    ..Default::default()
                },
            )
        }

        _ => return ErrorCode::InvalidParameterCode,
    };

    if result.is_err() {
        return ErrorCode::IllegalLinkProperty;
    }

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
/// `start_node` must be a valid non-null pointer to a NUL-terminated C string.
/// `end_node` must be a valid non-null pointer to a NUL-terminated C string.
/// `out_index` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_addlink(
    ph: *mut Project,
    id: *const c_char,
    link_type: c_int,
    start_node: *const c_char,
    end_node: *const c_char,
    out_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let c_str = unsafe { CStr::from_ptr(id) };
    let link_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    let link_type = match ENLinkType::from_repr(link_type) {
        Some(link_type) => link_type,
        None => return ErrorCode::InvalidParameterCode,
    };

    let c_str = unsafe { CStr::from_ptr(start_node) };
    let start_node_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };
    let c_str = unsafe { CStr::from_ptr(end_node) };
    let end_node_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    let mut new_pipe = PipeData {
        start_node: start_node_id.into(),
        end_node: end_node_id.into(),
        length: 330.0,
        diameter: 10.0 / 12.0 * simulation.network.options.unit_system.per_feet(), // 10 inches to
        roughness: match simulation.network.options.headloss_formula {
            HeadlossFormula::DarcyWeisbach => 0.0015,
            HeadlossFormula::HazenWilliams => 130.0,
            HeadlossFormula::ChezyManning => 0.01,
        },
        minor_loss: 0.0,
        check_valve: false,
        vertices: None,
        initial_status: LinkStatus::Open,
    };

    let mut new_valve = ValveData {
        start_node: start_node_id.into(),
        end_node: end_node_id.into(),
        diameter: 10.0 / 12.0 * simulation.network.options.unit_system.per_feet(), // 10 inches to
        valve_type: ValveType::PRV,
        setting: 0.0,
        minor_loss: 0.0,
        initial_status: LinkStatus::Active,
        vertices: None,
        curve_id: None,
    };

    let response = match link_type {
        ENLinkType::CVPipe => {
            new_pipe.check_valve = true;
            simulation.network.add_pipe(link_id, &new_pipe)
        }
        ENLinkType::Pipe => simulation.network.add_pipe(link_id, &new_pipe),
        ENLinkType::PRV => {
            new_valve.valve_type = ValveType::PRV;
            simulation.network.add_valve(link_id, &new_valve)
        }
        ENLinkType::PSV => {
            new_valve.valve_type = ValveType::PSV;
            simulation.network.add_valve(link_id, &new_valve)
        }
        ENLinkType::PBV => {
            new_valve.valve_type = ValveType::PBV;
            simulation.network.add_valve(link_id, &new_valve)
        }
        ENLinkType::FCV => {
            new_valve.valve_type = ValveType::FCV;
            simulation.network.add_valve(link_id, &new_valve)
        }
        ENLinkType::TCV => {
            new_valve.valve_type = ValveType::TCV;
            simulation.network.add_valve(link_id, &new_valve)
        }
        ENLinkType::PCV => {
            new_valve.valve_type = ValveType::PCV;
            simulation.network.add_valve(link_id, &new_valve)
        }
        ENLinkType::GPV => {
            new_valve.valve_type = ValveType::GPV;
            simulation.network.add_valve(link_id, &new_valve)
        }
        ENLinkType::Pump => {
            let new_pump = PumpData {
                start_node: start_node_id.into(),
                end_node: end_node_id.into(),
                speed: 1.0,
                head_curve_id: None,
                power: 0.0,
                initial_status: LinkStatus::Open,
                vertices: None,
            };
            simulation.network.add_pump(link_id, &new_pump)
        }
    };

    match response {
        Ok(_) => {
            unsafe {
                *out_index = (*simulation.network.link_map.get(link_id).unwrap() + 1) as c_int
            };
            ErrorCode::Ok
        }
        Err(_) => ErrorCode::InvalidParameterCode,
    }
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_index` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getheadcurveindex(
    ph: *mut Project,
    index: c_int,
    out_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let link = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };

    // no head curve, return 0
    unsafe { *out_index = 0 };

    match &link.link_type {
        LinkType::Pump(pump) => {
            if let Some(head_curve_id) = &pump.head_curve_id
                && let Some(index) = simulation.network.curve_map.get(head_curve_id)
            {
                unsafe { *out_index = (*index + 1) as c_int };
            }
        }
        _ => return ErrorCode::UndefinedPump,
    }

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setheadcurveindex(
    ph: *mut Project,
    index: c_int,
    head_curve_index: c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let (link_id, link_type) = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => (link.id.clone(), link.link_type.clone()),
        None => return ErrorCode::UndefinedLink,
    };

    let result = match link_type {
        LinkType::Pump(_) => {
            if head_curve_index == 0 {
                // clear the head curve
                simulation.network.update_pump(
                    &link_id,
                    &PumpUpdate {
                        head_curve_id: Some(None),
                        ..Default::default()
                    },
                )
            } else {
                // set the head curve
                let head_curve_id = match simulation
                    .network
                    .curves
                    .get((head_curve_index - 1) as usize)
                {
                    Some(curve) => curve.id.clone(),
                    None => return ErrorCode::UndefinedCurve,
                };

                simulation.network.update_pump(
                    &link_id,
                    &PumpUpdate {
                        head_curve_id: Some(Some(head_curve_id)),
                        ..Default::default()
                    },
                )
            }
        }
        _ => return ErrorCode::UndefinedPump,
    };

    if result.is_err() {
        return ErrorCode::IllegalLinkProperty;
    }

    ErrorCode::Ok
}
