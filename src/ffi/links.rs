//! FFI link accessors: `EN_addlink`, `EN_getlinkvalue`, `EN_setlinkvalue`, etc.

use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut, input_error_code};
use crate::ffi::util::{MAX_ID, write_out, write_str};
use crate::model::link::LinkType;
use crate::model::valve::ValveType;
use crate::simulation::Simulation;

use crate::constants::MperFT;
use crate::ffi::enums::LinkProperty;
use crate::ffi::enums::LinkType as ENLinkType;
use crate::ffi::enums::MISSING_VALUE;
use crate::ffi::enums::PumpStateType;
use crate::ffi::enums::PumpType;

use crate::model::curve::HeadCurveType;
use crate::model::network::modify::{LinkUpdate, PipeUpdate, PumpUpdate, ValveUpdate};
use crate::model::network::modify::{PipeData, PumpData, ValveData};
use crate::model::options::SimulationOptions;

use crate::model::link::LinkStatus;
use crate::model::options::HeadlossFormula;
use crate::model::units::UnitSystem;

use std::ffi::CStr;
use std::os::raw::{c_char, c_double, c_int};

/// Diameter of a newly created pipe or valve, in the project's diameter units
/// (inches or millimetres), matching legacy EPANET's `EN_addlink`.
const DEFAULT_DIAMETER: f64 = 10.0;
/// Length of a newly created pipe, in the project's length units.
const DEFAULT_LENGTH: f64 = 330.0;

/// Builds the default data of a newly created pipe, in the project's units.
fn default_pipe_data(
    start_node: &str,
    end_node: &str,
    options: &SimulationOptions,
    check_valve: bool,
) -> PipeData {
    PipeData {
        start_node: start_node.into(),
        end_node: end_node.into(),
        length: DEFAULT_LENGTH,
        diameter: DEFAULT_DIAMETER,
        roughness: match options.headloss_formula {
            HeadlossFormula::DarcyWeisbach => 0.0015,
            HeadlossFormula::HazenWilliams => 130.0,
            HeadlossFormula::ChezyManning => 0.01,
        },
        minor_loss: 0.0,
        check_valve,
        vertices: None,
        initial_status: LinkStatus::Open,
    }
}

/// Builds the default data of a newly created valve, in the project's units.
fn default_valve_data(start_node: &str, end_node: &str, valve_type: ValveType) -> ValveData {
    ValveData {
        start_node: start_node.into(),
        end_node: end_node.into(),
        diameter: DEFAULT_DIAMETER,
        valve_type,
        setting: 0.0,
        minor_loss: 0.0,
        initial_status: LinkStatus::Active,
        vertices: None,
        curve_id: None,
    }
}

/// Builds the default data of a newly created pump.
fn default_pump_data(start_node: &str, end_node: &str) -> PumpData {
    PumpData {
        start_node: start_node.into(),
        end_node: end_node.into(),
        speed: 1.0,
        head_curve_id: None,
        power: 0.0,
        initial_status: LinkStatus::Open,
        vertices: None,
    }
}

/// Maps a toolkit link-type code onto the valve type it represents.
fn valve_type_of(link_type: ENLinkType) -> Option<ValveType> {
    match link_type {
        ENLinkType::PRV => Some(ValveType::PRV),
        ENLinkType::PSV => Some(ValveType::PSV),
        ENLinkType::PBV => Some(ValveType::PBV),
        ENLinkType::FCV => Some(ValveType::FCV),
        ENLinkType::TCV => Some(ValveType::TCV),
        ENLinkType::PCV => Some(ValveType::PCV),
        ENLinkType::GPV => Some(ValveType::GPV),
        _ => None,
    }
}

/// Maps a link onto its toolkit link-type code.
fn link_type_code(link_type: &LinkType) -> ENLinkType {
    match link_type {
        LinkType::Pipe(pipe) => {
            if pipe.check_valve {
                ENLinkType::CVPipe
            } else {
                ENLinkType::Pipe
            }
        }
        LinkType::Pump(_) => ENLinkType::Pump,
        LinkType::Valve(valve) => match valve.valve_type {
            ValveType::PRV => ENLinkType::PRV,
            ValveType::PSV => ENLinkType::PSV,
            ValveType::PBV => ENLinkType::PBV,
            ValveType::FCV => ENLinkType::FCV,
            ValveType::GPV => ENLinkType::GPV,
            ValveType::TCV => ENLinkType::TCV,
            ValveType::PCV => ENLinkType::PCV,
        },
    }
}

/// Converts a diameter from feet to the project's diameter units.
fn diameter_conversion(options: &SimulationOptions) -> f64 {
    match options.unit_system {
        UnitSystem::US => 12.0,         // ft to in
        UnitSystem::SI => MperFT * 1e3, // ft to mm
    }
}

/// Delete the link from the network
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
/// `out_id` must point to a buffer of at least `EN_MAXID + 1` bytes.
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

    unsafe { write_str(out_id, link_id, MAX_ID) };

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
    let index = (index - 1) as usize;
    let link = match simulation.network.links.get_mut(index) {
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
        .insert(new_link_id.into(), index);

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

    unsafe { write_out(out_type, link_type_code(&link.link_type) as c_int) };

    ErrorCode::Ok
}

/// Changes a link's type.
///
/// Converting between a pipe and a pump or valve replaces the link, so its
/// index may change; the new index is written back to `inout_index`. All
/// type-specific properties are reset to the defaults used by [`EN_addlink`].
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `inout_index` must be a valid non-null readable and writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setlinktype(
    ph: *mut Project,
    inout_index: *mut c_int,
    link_type: c_int,
    action_code: c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    if inout_index.is_null() {
        return ErrorCode::InvalidFormat;
    }
    let index = unsafe { *inout_index };

    let new_type = match ENLinkType::from_repr(link_type) {
        Some(link_type) => link_type,
        None => return ErrorCode::InvalidParameterCode,
    };

    let unconditional = match action_code {
        0 => true,
        1 => false,
        _ => return ErrorCode::InvalidParameterCode,
    };

    let link = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };

    let old_type = link_type_code(&link.link_type);
    if old_type == new_type {
        return ErrorCode::Ok;
    }

    let link_id = link.id.clone();
    let start_node_id = link.start_node_id.clone();
    let end_node_id = link.end_node_id.clone();
    let initial_status = link.initial_status;
    let vertices = link.vertices.clone();

    // a pipe only differs from a check valve pipe by a flag
    if matches!(old_type, ENLinkType::Pipe | ENLinkType::CVPipe)
        && matches!(new_type, ENLinkType::Pipe | ENLinkType::CVPipe)
    {
        return match simulation.network.update_pipe(
            &link_id,
            &PipeUpdate {
                check_valve: Some(new_type == ENLinkType::CVPipe),
                ..Default::default()
            },
        ) {
            Ok(_) => ErrorCode::Ok,
            Err(_) => ErrorCode::IllegalLinkProperty,
        };
    }

    // valves keep their identity (and index) when only their kind changes
    if let (Some(_), Some(new_valve_type)) = (valve_type_of(old_type), valve_type_of(new_type)) {
        if new_valve_type == ValveType::GPV {
            // TODO: a GPV needs a head loss curve, which this call cannot supply.
            return ErrorCode::NotImplemented;
        }
        return match simulation.network.update_valve(
            &link_id,
            &ValveUpdate {
                valve_type: Some(new_valve_type),
                setting: Some(0.0),
                curve_id: Some(None),
                ..Default::default()
            },
        ) {
            Ok(_) => ErrorCode::Ok,
            Err(_) => ErrorCode::IllegalLinkProperty,
        };
    }

    if new_type == ENLinkType::GPV {
        // TODO: a GPV needs a head loss curve, which this call cannot supply.
        return ErrorCode::NotImplemented;
    }

    // any other conversion recreates the link, which moves it to the end of
    // the link list
    if let Err(error) = simulation.network.remove_link(&link_id, unconditional) {
        return match input_error_code(&error) {
            ErrorCode::UndefinedLink => ErrorCode::UndefinedLink,
            _ => ErrorCode::DeleteNodeOrLinkInControl,
        };
    }

    let result = match new_type {
        ENLinkType::Pipe | ENLinkType::CVPipe => {
            let mut data = default_pipe_data(
                &start_node_id,
                &end_node_id,
                &simulation.network.options,
                new_type == ENLinkType::CVPipe,
            );
            data.initial_status = initial_status;
            data.vertices = vertices;
            simulation.network.add_pipe(&link_id, &data)
        }
        ENLinkType::Pump => {
            let mut data = default_pump_data(&start_node_id, &end_node_id);
            data.initial_status = initial_status;
            data.vertices = vertices;
            simulation.network.add_pump(&link_id, &data)
        }
        valve => {
            let valve_type = match valve_type_of(valve) {
                Some(valve_type) => valve_type,
                None => return ErrorCode::InvalidParameterCode,
            };
            let mut data = default_valve_data(&start_node_id, &end_node_id, valve_type);
            data.vertices = vertices;
            simulation.network.add_valve(&link_id, &data)
        }
    };

    match result {
        Ok(_) => {
            let new_index = simulation
                .network
                .link_map
                .get(&link_id)
                .copied()
                .unwrap_or(0);
            unsafe { write_out(inout_index, (new_index + 1) as c_int) };
            ErrorCode::Ok
        }
        Err(error) => input_error_code(&error),
    }
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
    }
}

/// Reads a single link property, converted to the project's units.
///
/// `index` is zero-based. Properties that only apply to another kind of link
/// read as zero; properties whose underlying feature is unsupported return
/// [`ErrorCode::NotImplemented`].
fn link_value(
    simulation: &Simulation,
    index: usize,
    property: LinkProperty,
) -> Result<f64, ErrorCode> {
    let link = match simulation.network.links.get(index) {
        Some(link) => link,
        None => return Err(ErrorCode::UndefinedLink),
    };

    let options = &simulation.network.options;

    let value = match property {
        LinkProperty::Diameter => {
            let conversion = diameter_conversion(options);
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
            let status = match simulation.solved_state() {
                Some(state) => state.statuses[index],
                None => link.initial_status,
            };
            match status {
                LinkStatus::Xhead
                | LinkStatus::TempClosed
                | LinkStatus::Closed
                | LinkStatus::FixedClosed => 0.0,
                _ => 1.0,
            }
        }
        LinkProperty::InitStatus => match &link.initial_status {
            LinkStatus::Closed | LinkStatus::FixedClosed => 0.0,
            _ => 1.0,
        },
        LinkProperty::Setting | LinkProperty::InitSetting => match &link.link_type {
            LinkType::Pipe(_) => 0.0,
            LinkType::Valve(valve) => match valve.valve_type {
                ValveType::FCV => valve.setting * options.flow_units.per_cfs(),
                ValveType::TCV => valve.setting,
                ValveType::PCV => valve.setting,
                // PRV/PSV/PBV settings are stored as pressure (in feet of
                // head); just convert to the user's pressure units.
                ValveType::PRV | ValveType::PSV | ValveType::PBV => {
                    valve.setting * options.pressure_units.per_feet()
                }
                ValveType::GPV => 0.0,
            },
            LinkType::Pump(pump) => pump.speed,
        },

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

        LinkProperty::PumpState => match &link.link_type {
            LinkType::Pump(_) => {
                let status = match simulation.solved_state() {
                    Some(state) => state.statuses[index],
                    None => link.initial_status,
                };
                match status {
                    LinkStatus::Xhead => PumpStateType::XHead as i32 as f64,
                    LinkStatus::Xflow => PumpStateType::XFlow as i32 as f64,
                    LinkStatus::Closed | LinkStatus::TempClosed | LinkStatus::FixedClosed => {
                        PumpStateType::Closed as i32 as f64
                    }
                    _ => PumpStateType::Open as i32 as f64,
                }
            }
            _ => 0.0,
        },
        LinkProperty::PumpPower => match &link.link_type {
            LinkType::Pump(pump) => pump.power * options.unit_system.per_horsepower(),
            _ => 0.0,
        },
        LinkProperty::PumpHCurve => match &link.link_type {
            LinkType::Pump(pump) => pump
                .head_curve_id
                .as_ref()
                .and_then(|id| simulation.network.curve_map.get(id))
                .map(|index| (index + 1) as f64)
                .unwrap_or(0.0),
            _ => 0.0,
        },
        LinkProperty::GpvCurve | LinkProperty::PcvCurve => match &link.link_type {
            LinkType::Valve(valve) => valve
                .curve_id
                .as_ref()
                .and_then(|id| simulation.network.curve_map.get(id))
                .map(|index| (index + 1) as f64)
                .unwrap_or(0.0),
            _ => 0.0,
        },
        LinkProperty::ValveType => match &link.link_type {
            LinkType::Valve(_) => link_type_code(&link.link_type) as i32 as f64,
            _ => return Err(ErrorCode::NotAValve),
        },
        LinkProperty::LinkInControl => {
            let in_control = simulation
                .network
                .controls
                .iter()
                .any(|control| control.link_id == link.id);
            if in_control { 1.0 } else { 0.0 }
        }

        // TODO: pump energy use, efficiency and cost are not computed.
        LinkProperty::Energy
        | LinkProperty::PumpEffic
        | LinkProperty::PumpECurve
        | LinkProperty::PumpECost
        | LinkProperty::PumpEPat
        | LinkProperty::LinkPattern => return Err(ErrorCode::NotImplemented),
        // TODO: pipe leakage is not modelled.
        LinkProperty::LeakArea | LinkProperty::LeakExpan | LinkProperty::LinkLeakage => {
            return Err(ErrorCode::NotImplemented);
        }
        // TODO: water quality analysis is not implemented.
        LinkProperty::KBulk | LinkProperty::KWall | LinkProperty::LinkQual => {
            return Err(ErrorCode::NotImplemented);
        }
    };

    Ok(value)
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

    let property = match LinkProperty::from_repr(property) {
        Some(property) => property,
        None => return ErrorCode::InvalidParameterCode,
    };

    if index < 1 {
        return ErrorCode::UndefinedLink;
    }

    match link_value(simulation, (index - 1) as usize, property) {
        Ok(value) => {
            unsafe { write_out(out_value, value) };
            ErrorCode::Ok
        }
        Err(code) => code,
    }
}

/// Retrieves a property value for all links, in index order.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_values` must point to writable memory for `EN_getcount(EN_LINKCOUNT)`
/// `c_double` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getlinkvalues(
    ph: *mut Project,
    property: c_int,
    out_values: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let property = match LinkProperty::from_repr(property) {
        Some(property) => property,
        None => return ErrorCode::InvalidParameterCode,
    };

    if out_values.is_null() {
        return ErrorCode::InvalidFormat;
    }

    for index in 0..simulation.network.links.len() {
        match link_value(simulation, index, property) {
            Ok(value) => unsafe { *out_values.add(index) = value },
            Err(code) => return code,
        }
    }

    ErrorCode::Ok
}

/// Writes a single link property, interpreting `value` in the project's units.
fn set_link_value(
    simulation: &mut Simulation,
    index: usize,
    property: LinkProperty,
    value: f64,
) -> ErrorCode {
    let link = match simulation.network.links.get(index) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };
    let link_id = link.id.clone();
    let link_type = link.link_type.clone();
    let is_check_valve = matches!(&link.link_type, LinkType::Pipe(pipe) if pipe.check_valve);

    let result = match property {
        LinkProperty::Diameter => {
            if value <= 0.0 {
                return ErrorCode::IllegalLinkProperty;
            }
            match &link_type {
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
                LinkType::Pump(_) => return ErrorCode::InvalidParameterCode,
            }
        }
        LinkProperty::Length => {
            if value <= 0.0 {
                return ErrorCode::IllegalLinkProperty;
            }
            simulation.network.update_pipe(
                &link_id,
                &PipeUpdate {
                    length: Some(value),
                    ..Default::default()
                },
            )
        }
        LinkProperty::Roughness => {
            if value <= 0.0 {
                return ErrorCode::IllegalLinkProperty;
            }
            simulation.network.update_pipe(
                &link_id,
                &PipeUpdate {
                    roughness: Some(value),
                    ..Default::default()
                },
            )
        }
        LinkProperty::MinorLoss => {
            if value < 0.0 {
                return ErrorCode::IllegalLinkProperty;
            }
            match &link_type {
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
                LinkType::Pump(_) => return ErrorCode::InvalidParameterCode,
            }
        }
        LinkProperty::InitStatus | LinkProperty::Status => {
            let is_valve = matches!(&link_type, LinkType::Valve(_));
            let link_status = match (value, is_valve) {
                (0.0, false) => LinkStatus::Closed,
                (1.0, false) => LinkStatus::Open,
                (0.0, true) => LinkStatus::FixedClosed,
                (1.0, true) => LinkStatus::FixedOpen,
                _ => return ErrorCode::InvalidParameterCode,
            };

            // a check valve pipe controls its own status
            if is_check_valve {
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
            match &link_type {
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

        LinkProperty::PumpPower => {
            if value <= 0.0 {
                return ErrorCode::IllegalLinkProperty;
            }
            simulation.network.update_pump(
                &link_id,
                &PumpUpdate {
                    power: Some(value),
                    ..Default::default()
                },
            )
        }

        LinkProperty::PumpHCurve => {
            let curve_id = if value < 1.0 {
                None
            } else {
                match simulation.network.curves.get(value as usize - 1) {
                    Some(curve) => Some(curve.id.clone()),
                    None => return ErrorCode::UndefinedCurve,
                }
            };
            simulation.network.update_pump(
                &link_id,
                &PumpUpdate {
                    head_curve_id: Some(curve_id),
                    ..Default::default()
                },
            )
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

        // TODO: pump energy data, pipe leakage and water quality reaction
        // coefficients are not part of the network model.
        LinkProperty::LinkPattern
        | LinkProperty::PumpEffic
        | LinkProperty::PumpECurve
        | LinkProperty::PumpECost
        | LinkProperty::PumpEPat
        | LinkProperty::LeakArea
        | LinkProperty::LeakExpan
        | LinkProperty::KBulk
        | LinkProperty::KWall
        | LinkProperty::LinkQual => return ErrorCode::NotImplemented,

        // computed results cannot be assigned
        _ => return ErrorCode::InvalidParameterCode,
    };

    match result {
        Ok(_) => ErrorCode::Ok,
        Err(error) => match input_error_code(&error) {
            ErrorCode::UndefinedLink => ErrorCode::UndefinedLink,
            ErrorCode::UndefinedCurve => ErrorCode::UndefinedCurve,
            _ => ErrorCode::IllegalLinkProperty,
        },
    }
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

    let property = match LinkProperty::from_repr(property) {
        Some(property) => property,
        None => return ErrorCode::InvalidParameterCode,
    };

    if index < 1 {
        return ErrorCode::UndefinedLink;
    }

    set_link_value(simulation, (index - 1) as usize, property, value)
}

/// Sets a property value for all links, in index order.
///
/// On failure the 1-based index of the offending link is written to
/// `out_bad_index`.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `values` must point to `EN_getcount(EN_LINKCOUNT)` readable `c_double` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setlinkvalues(
    ph: *mut Project,
    property: c_int,
    values: *const c_double,
    out_bad_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let property = match LinkProperty::from_repr(property) {
        Some(property) => property,
        None => return ErrorCode::InvalidParameterCode,
    };

    if values.is_null() {
        return ErrorCode::InvalidFormat;
    }

    for index in 0..simulation.network.links.len() {
        let value = unsafe { *values.add(index) };
        let code = set_link_value(simulation, index, property, value);
        if code != ErrorCode::Ok {
            unsafe { write_out(out_bad_index, (index + 1) as c_int) };
            return code;
        }
    }

    ErrorCode::Ok
}

/// Sets a group of properties for a pipe link.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setpipedata(
    ph: *mut Project,
    index: c_int,
    length: c_double,
    diameter: c_double,
    roughness: c_double,
    minor_loss: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let link_id = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link.id.clone(),
        None => return ErrorCode::UndefinedLink,
    };

    if length <= 0.0 || diameter <= 0.0 || roughness <= 0.0 || minor_loss < 0.0 {
        return ErrorCode::IllegalLinkProperty;
    }

    match simulation.network.update_pipe(
        &link_id,
        &PipeUpdate {
            length: Some(length),
            diameter: Some(diameter),
            roughness: Some(roughness),
            minor_loss: Some(minor_loss),
            ..Default::default()
        },
    ) {
        Ok(_) => ErrorCode::Ok,
        Err(_) => ErrorCode::IllegalLinkProperty,
    }
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

    let response = match link_type {
        ENLinkType::Pipe | ENLinkType::CVPipe => {
            let data = default_pipe_data(
                start_node_id,
                end_node_id,
                &simulation.network.options,
                link_type == ENLinkType::CVPipe,
            );
            simulation.network.add_pipe(link_id, &data)
        }
        ENLinkType::Pump => {
            let data = default_pump_data(start_node_id, end_node_id);
            simulation.network.add_pump(link_id, &data)
        }
        valve => {
            let valve_type = match valve_type_of(valve) {
                Some(valve_type) => valve_type,
                None => return ErrorCode::InvalidParameterCode,
            };
            let data = default_valve_data(start_node_id, end_node_id, valve_type);
            simulation.network.add_valve(link_id, &data)
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

/// Retrieves the type of head curve used by a pump.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_type` must be null or a valid writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getpumptype(
    ph: *mut Project,
    index: c_int,
    out_type: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let link = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };

    let LinkType::Pump(pump) = &link.link_type else {
        return ErrorCode::UndefinedPump;
    };

    let pump_type = match &pump.head_curve {
        Some(curve) => match curve.curve_type {
            HeadCurveType::Custom => PumpType::Custom,
            _ => PumpType::PowerFunc,
        },
        None if pump.power > 0.0 => PumpType::ConstHp,
        None => PumpType::NoCurve,
    };

    unsafe { write_out(out_type, pump_type as c_int) };

    ErrorCode::Ok
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

/// Retrieves the number of internal vertex points assigned to a link.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_count` must be null or a valid writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getvertexcount(
    ph: *mut Project,
    index: c_int,
    out_count: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let link = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };

    let count = link.vertices.as_ref().map_or(0, |v| v.len());
    unsafe { write_out(out_count, count as c_int) };

    ErrorCode::Ok
}

/// Retrieves the coordinates of a link's vertex point.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_x` and `out_y` must be null or valid writable pointers.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getvertex(
    ph: *mut Project,
    index: c_int,
    vertex: c_int,
    out_x: *mut c_double,
    out_y: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let link = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };

    let vertices = match &link.vertices {
        Some(vertices) => vertices,
        None => return ErrorCode::LinkNoVertices,
    };

    if vertex < 1 || vertex as usize > vertices.len() {
        return ErrorCode::LinkNoVertices;
    }

    let (x, y) = vertices[(vertex - 1) as usize];
    unsafe { write_out(out_x, x) };
    unsafe { write_out(out_y, y) };

    ErrorCode::Ok
}

/// Sets the coordinates of a link's vertex point.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setvertex(
    ph: *mut Project,
    index: c_int,
    vertex: c_int,
    x: c_double,
    y: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let link = match simulation.network.links.get_mut((index - 1) as usize) {
        Some(link) => link,
        None => return ErrorCode::UndefinedLink,
    };

    let vertices = match &mut link.vertices {
        Some(vertices) => vertices,
        None => return ErrorCode::LinkNoVertices,
    };

    if vertex < 1 || vertex as usize > vertices.len() {
        return ErrorCode::LinkNoVertices;
    }

    vertices[(vertex - 1) as usize] = (x, y);

    ErrorCode::Ok
}

/// Assigns a set of internal vertex points to a link.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `x` and `y` must each point to `count` readable `c_double` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setvertices(
    ph: *mut Project,
    index: c_int,
    x: *const c_double,
    y: *const c_double,
    count: c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let link_id = match simulation.network.links.get((index - 1) as usize) {
        Some(link) => link.id.clone(),
        None => return ErrorCode::UndefinedLink,
    };

    if count < 0 || (count > 0 && (x.is_null() || y.is_null())) {
        return ErrorCode::InvalidFormat;
    }

    let count = count as usize;
    let mut vertices = Vec::with_capacity(count);
    for i in 0..count {
        vertices.push(unsafe { (*x.add(i), *y.add(i)) });
    }

    match simulation.network.update_link(
        &link_id,
        &LinkUpdate {
            vertices: Some(vertices),
            ..Default::default()
        },
    ) {
        Ok(_) => ErrorCode::Ok,
        Err(_) => ErrorCode::UndefinedLink,
    }
}
