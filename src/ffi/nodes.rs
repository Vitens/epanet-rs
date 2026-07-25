//! FFI node accessors: `EN_addnode`, `EN_getnodevalue`, `EN_setnodevalue`, etc.

use crate::ffi::enums::NodeProperty;
use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut, input_error_code};
use crate::ffi::util::{MAX_ID, read_str, write_out, write_str};
use crate::model::control::ControlCondition;
use crate::model::demand::Demand;
use crate::model::network::modify::{
    JunctionData, JunctionUpdate, NodeUpdate, ReservoirData, ReservoirUpdate, TankData, TankUpdate,
};
use crate::model::node::NodeType;
use crate::simulation::Simulation;

use crate::ffi::enums::NodeType as ENNodeType;

use std::ffi::CStr;
use std::os::raw::{c_char, c_double, c_int};

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
/// `out_index` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_addnode(
    ph: *mut Project,
    id: *const c_char,
    node_type: c_int,
    out_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let c_str = unsafe { CStr::from_ptr(id) };
    let node_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    let node_type = match ENNodeType::from_repr(node_type) {
        Some(node_type) => node_type,
        None => return ErrorCode::InvalidParameterCode,
    };

    let result = match node_type {
        ENNodeType::Junction => simulation.network.add_junction(
            node_id,
            &JunctionData {
                elevation: 0.0,
                demands: vec![Demand {
                    basedemand: 0.0,
                    pattern: None,
                    pattern_index: None,
                    name: None,
                }],
                emitter_coefficient: 0.0,
                coordinates: None,
            },
        ),
        ENNodeType::Reservoir => simulation.network.add_reservoir(
            node_id,
            &ReservoirData {
                elevation: 0.0,
                head_pattern: None,
                coordinates: None,
            },
        ),
        ENNodeType::Tank => simulation.network.add_tank(
            node_id,
            &TankData {
                elevation: 0.0,
                initial_level: 0.0,
                min_level: 0.0,
                max_level: 0.0,
                diameter: 0.0,
                min_volume: 0.0,
                volume_curve_id: None,
                overflow: false,
                coordinates: None,
            },
        ),
    };

    if result.is_err() {
        return ErrorCode::IllegalNodeProperty;
    }
    unsafe { *out_index = (*simulation.network.node_map.get(node_id).unwrap() + 1) as c_int };
    ErrorCode::Ok
}

/// Gets the index of a node given its ID name.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
/// `out_index` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getnodeindex(
    ph: *mut Project,
    id: *const c_char,
    out_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let c_str = unsafe { CStr::from_ptr(id) };
    let node_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    // get the node index from the network
    let node_index = match simulation.network.node_map.get(node_id) {
        Some(&index) => index,
        None => return ErrorCode::UndefinedNode,
    };

    // EPANET indexes from 1, so we need to add 1 to the index
    unsafe { *out_index = (node_index + 1) as c_int };

    ErrorCode::Ok
}

// Gets the ID name of a node given its index.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_id` must point to a buffer of at least `EN_MAXID + 1` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getnodeid(
    ph: *mut Project,
    index: c_int,
    out_id: *mut c_char,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    // EPANET indexes from 1, so we need to subtract 1 from the index
    let index = index - 1;

    let node_id = match simulation.network.nodes.get(index as usize) {
        Some(node) => node.id.as_ref(),
        None => return ErrorCode::UndefinedNode,
    };

    unsafe { write_str(out_id, node_id, MAX_ID) };

    ErrorCode::Ok
}

// Sets the ID name of a node given its index.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setnodeid(
    ph: *mut Project,
    index: c_int,
    id: *const c_char,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let c_str = unsafe { CStr::from_ptr(id) };

    let new_node_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    // EPANET indexes from 1, so we need to subtract 1 from the index
    let index = index - 1;

    let node = match simulation.network.nodes.get_mut(index as usize) {
        Some(node) => node,
        None => return ErrorCode::UndefinedNode,
    };

    // check if the new node id is already in use
    if simulation.network.node_map.contains_key(new_node_id) {
        return ErrorCode::DuplicateId;
    }

    let old_node_id = node.id.clone();

    // remove the old node id from the node map
    simulation.network.node_map.remove(&node.id);

    // update the node id
    node.id = new_node_id.into();

    // update the node map
    simulation
        .network
        .node_map
        .insert(new_node_id.into(), index as usize);

    // update all links that point to the old node id to point to the new node id
    for link in simulation.network.links.iter_mut() {
        if link.start_node_id == old_node_id {
            link.start_node_id = new_node_id.into();
        }
        if link.end_node_id == old_node_id {
            link.end_node_id = new_node_id.into();
        }
    }

    ErrorCode::Ok
}

// Get the node type given its index.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_type` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getnodetype(
    ph: *mut Project,
    index: c_int,
    out_type: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    // EPANET indexes from 1, so we need to add 1 to the index
    let index = index - 1;

    let node_type = match simulation.network.nodes.get(index as usize) {
        Some(node) => &node.node_type,
        None => return ErrorCode::UndefinedNode,
    };

    let node_type_int = match node_type {
        NodeType::Junction(_) => ENNodeType::Junction as i32,
        NodeType::Reservoir(_) => ENNodeType::Reservoir as i32,
        NodeType::Tank(_) => ENNodeType::Tank as i32,
    };

    unsafe { *out_type = node_type_int };

    ErrorCode::Ok
}

/// Reads a single node property, converted to the project's units.
///
/// `index` is zero-based. Properties whose underlying feature is unsupported
/// return [`ErrorCode::NotImplemented`].
fn node_value(
    simulation: &Simulation,
    index: usize,
    property: NodeProperty,
) -> Result<f64, ErrorCode> {
    let node = match simulation.network.nodes.get(index) {
        Some(node) => node,
        None => return Err(ErrorCode::UndefinedNode),
    };

    let options = &simulation.network.options;
    let unit_system = &options.unit_system;
    let flow_units = &options.flow_units;
    let pressure_units = &options.pressure_units;

    let value = match property {
        NodeProperty::Elevation => node.elevation * unit_system.per_feet(),
        NodeProperty::BaseDemand => match &node.node_type {
            NodeType::Junction(junction) => junction
                .demands
                .first()
                .map(|d| d.basedemand * flow_units.per_cfs())
                .unwrap_or(0.0),
            _ => 0.0,
        },
        NodeProperty::Pattern => match &node.node_type {
            NodeType::Junction(junction) => junction
                .demands
                .first()
                .and_then(|d| d.pattern_index)
                .map(|index| (index + 1) as f64)
                .unwrap_or(0.0),
            NodeType::Reservoir(reservoir) => reservoir
                .head_pattern_index
                .map(|index| (index + 1) as f64)
                .unwrap_or(0.0),
            _ => 0.0,
        },
        NodeProperty::Emitter => match &node.node_type {
            NodeType::Junction(junction) => {
                let e = junction.emitter_coefficient;
                if e > 0.0 {
                    flow_units.per_cfs()
                        / (pressure_units.per_feet() * e).powf(1.0 / options.emitter_exponent)
                } else {
                    0.0
                }
            }
            _ => 0.0,
        },

        NodeProperty::TankLevel => match &node.node_type {
            NodeType::Tank(tank) => match simulation.solved_state() {
                Some(state) => (state.heads[index] - node.elevation) * unit_system.per_feet(),
                None => tank.initial_level * unit_system.per_feet(),
            },
            _ => 0.0,
        },

        NodeProperty::InitVolume => match &node.node_type {
            NodeType::Tank(tank) => {
                tank.volume_at_level(tank.initial_level) * unit_system.per_cubic_feet()
            }
            _ => 0.0,
        },

        NodeProperty::Demand => simulation
            .solved_state()
            .map_or(0.0, |state| state.demands[index] * flow_units.per_cfs()),
        NodeProperty::Head => simulation
            .solved_state()
            .map_or(0.0, |state| state.heads[index] * unit_system.per_feet()),
        NodeProperty::Pressure => simulation.solved_state().map_or(0.0, |state| {
            (state.heads[index] - node.elevation) * pressure_units.per_feet()
        }),
        NodeProperty::EmitterFlow => simulation.solved_state().map_or(0.0, |state| {
            state.emitter_flows[index] * flow_units.per_cfs()
        }),
        // the converged demand also carries the emitter outflow
        NodeProperty::DemandFlow => simulation.solved_state().map_or(0.0, |state| {
            (state.demands[index] - state.emitter_flows[index]) * flow_units.per_cfs()
        }),

        NodeProperty::TankDiam => match &node.node_type {
            NodeType::Tank(tank) => tank.diameter * unit_system.per_feet(),
            _ => 0.0,
        },
        NodeProperty::MinVolume => match &node.node_type {
            NodeType::Tank(tank) => tank.min_volume() * unit_system.per_cubic_feet(),
            _ => 0.0,
        },
        NodeProperty::MaxVolume => match &node.node_type {
            NodeType::Tank(tank) => tank.max_volume() * unit_system.per_cubic_feet(),
            _ => 0.0,
        },
        NodeProperty::MinLevel => match &node.node_type {
            NodeType::Tank(tank) => tank.min_level * unit_system.per_feet(),
            _ => 0.0,
        },
        NodeProperty::MaxLevel => match &node.node_type {
            NodeType::Tank(tank) => tank.max_level * unit_system.per_feet(),
            _ => 0.0,
        },

        NodeProperty::TankVolume => match &node.node_type {
            NodeType::Tank(tank) => {
                if let Some(state) = simulation.solved_state() {
                    tank.volume_at_head(state.heads[index]) * unit_system.per_cubic_feet()
                } else {
                    tank.volume_at_level(tank.initial_level) * unit_system.per_cubic_feet()
                }
            }
            _ => 0.0,
        },

        NodeProperty::VolCurve => match &node.node_type {
            NodeType::Tank(tank) => tank
                .volume_curve_id
                .as_ref()
                .and_then(|id| simulation.network.curve_map.get(id))
                .map(|index| (index + 1) as f64)
                .unwrap_or(0.0),
            _ => 0.0,
        },

        NodeProperty::CanOverflow => match &node.node_type {
            NodeType::Tank(tank) if tank.overflow => 1.0,
            _ => 0.0,
        },

        NodeProperty::NodeInControl => {
            let in_control =
                simulation
                    .network
                    .controls
                    .iter()
                    .any(|control| match control.condition {
                        ControlCondition::HighLevel { tank_index, .. }
                        | ControlCondition::LowLevel { tank_index, .. } => tank_index == index,
                        ControlCondition::HighPressure { node_index, .. }
                        | ControlCondition::LowPressure { node_index, .. } => node_index == index,
                        _ => false,
                    });
            if in_control { 1.0 } else { 0.0 }
        }

        // TODO: pressure-driven demand deficits are not retained after the
        // solver converges.
        NodeProperty::DemandDeficit | NodeProperty::FullDemand => {
            return Err(ErrorCode::NotImplemented);
        }
        // TODO: pipe/node leakage is not modelled.
        NodeProperty::LeakageFlow => return Err(ErrorCode::NotImplemented),
        // TODO: water quality analysis (including tank mixing) is not implemented.
        NodeProperty::InitQual
        | NodeProperty::SourceQual
        | NodeProperty::SourcePat
        | NodeProperty::SourceType
        | NodeProperty::SourceMass
        | NodeProperty::Quality
        | NodeProperty::MixModel
        | NodeProperty::MixZoneVol
        | NodeProperty::MixFraction
        | NodeProperty::TankKBulk => return Err(ErrorCode::NotImplemented),
    };

    Ok(value)
}

/// Retrieves the property value of a node
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_value` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getnodevalue(
    ph: *mut Project,
    index: c_int,
    property: c_int,
    out_value: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let property = match NodeProperty::from_repr(property) {
        Some(property) => property,
        None => return ErrorCode::InvalidParameterCode,
    };

    // EPANET indexes from 1, so we need to subtract 1 from the index
    if index < 1 {
        return ErrorCode::UndefinedNode;
    }
    match node_value(simulation, (index - 1) as usize, property) {
        Ok(value) => {
            unsafe { write_out(out_value, value as c_double) };
            ErrorCode::Ok
        }
        Err(code) => code,
    }
}

/// Retrieves a property value for all nodes, in index order.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_values` must point to writable memory for `EN_getcount(EN_NODECOUNT)`
/// `c_double` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getnodevalues(
    ph: *mut Project,
    property: c_int,
    out_values: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let property = match NodeProperty::from_repr(property) {
        Some(property) => property,
        None => return ErrorCode::InvalidParameterCode,
    };

    if out_values.is_null() {
        return ErrorCode::InvalidFormat;
    }

    for index in 0..simulation.network.nodes.len() {
        match node_value(simulation, index, property) {
            Ok(value) => unsafe { *out_values.add(index) = value as c_double },
            Err(code) => return code,
        }
    }

    ErrorCode::Ok
}

/// Writes a single node property, interpreting `value` in the project's units.
fn set_node_value(
    simulation: &mut Simulation,
    index: usize,
    property: NodeProperty,
    value: f64,
) -> ErrorCode {
    let network = &mut simulation.network;

    let node_id = match network.nodes.get(index) {
        Some(node) => node.id.clone(),
        None => return ErrorCode::UndefinedNode,
    };
    let is_reservoir = matches!(network.nodes[index].node_type, NodeType::Reservoir(_));

    let result = match property {
        NodeProperty::Elevation => network.update_node(
            &node_id,
            &NodeUpdate {
                elevation: Some(value),
                ..Default::default()
            },
        ),
        NodeProperty::BaseDemand => network.update_junction(
            &node_id,
            &JunctionUpdate {
                basedemand: Some(value),
                ..Default::default()
            },
        ),
        NodeProperty::Emitter => network.update_junction(
            &node_id,
            &JunctionUpdate {
                emitter_coefficient: Some(value),
                ..Default::default()
            },
        ),
        NodeProperty::TankDiam => {
            if value <= 0.0 {
                return ErrorCode::IllegalNodeProperty;
            }
            network.update_tank(
                &node_id,
                &TankUpdate {
                    diameter: Some(value),
                    ..Default::default()
                },
            )
        }
        NodeProperty::MinLevel => network.update_tank(
            &node_id,
            &TankUpdate {
                min_level: Some(value),
                ..Default::default()
            },
        ),
        NodeProperty::MaxLevel => network.update_tank(
            &node_id,
            &TankUpdate {
                max_level: Some(value),
                ..Default::default()
            },
        ),
        NodeProperty::TankLevel => network.update_tank(
            &node_id,
            &TankUpdate {
                initial_level: Some(value),
                ..Default::default()
            },
        ),
        NodeProperty::MinVolume => {
            if value < 0.0 {
                return ErrorCode::IllegalNodeProperty;
            }
            network.update_tank(
                &node_id,
                &TankUpdate {
                    min_volume: Some(value),
                    ..Default::default()
                },
            )
        }
        NodeProperty::CanOverflow => network.update_tank(
            &node_id,
            &TankUpdate {
                overflow: Some(value != 0.0),
                ..Default::default()
            },
        ),
        NodeProperty::Pattern => {
            // a zero index clears the node's pattern
            let pattern_id = if value < 1.0 {
                None
            } else {
                match network.patterns.get(value as usize - 1) {
                    Some(pattern) => Some(pattern.id.clone()),
                    None => return ErrorCode::UndefinedPattern,
                }
            };
            if is_reservoir {
                network.update_reservoir(
                    &node_id,
                    &ReservoirUpdate {
                        head_pattern: Some(pattern_id),
                        ..Default::default()
                    },
                )
            } else {
                network.update_junction(
                    &node_id,
                    &JunctionUpdate {
                        pattern: Some(pattern_id),
                        ..Default::default()
                    },
                )
            }
        }
        // TODO: tank volume curves, water quality sources and tank mixing are
        // not implemented.
        NodeProperty::VolCurve
        | NodeProperty::InitQual
        | NodeProperty::SourceQual
        | NodeProperty::SourcePat
        | NodeProperty::SourceType
        | NodeProperty::SourceMass
        | NodeProperty::MixModel
        | NodeProperty::MixFraction
        | NodeProperty::TankKBulk => return ErrorCode::NotImplemented,
        // computed results cannot be assigned
        _ => return ErrorCode::InvalidParameterCode,
    };

    match result {
        Ok(_) => ErrorCode::Ok,
        Err(error) => match input_error_code(&error) {
            ErrorCode::UndefinedNode => ErrorCode::UndefinedNode,
            ErrorCode::NotATank => ErrorCode::NotATank,
            _ => ErrorCode::IllegalNodeProperty,
        },
    }
}

// Set the property value of a node
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setnodevalue(
    ph: *mut Project,
    index: c_int,
    property: c_int,
    value: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let property = match NodeProperty::from_repr(property) {
        Some(property) => property,
        None => return ErrorCode::InvalidParameterCode,
    };

    if index < 1 {
        return ErrorCode::UndefinedNode;
    }

    set_node_value(simulation, (index - 1) as usize, property, value)
}

/// Sets a property value for all nodes, in index order.
///
/// On failure the 1-based index of the offending node is written to
/// `out_bad_index`.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `values` must point to `EN_getcount(EN_NODECOUNT)` readable `c_double` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setnodevalues(
    ph: *mut Project,
    property: c_int,
    values: *const c_double,
    out_bad_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let property = match NodeProperty::from_repr(property) {
        Some(property) => property,
        None => return ErrorCode::InvalidParameterCode,
    };

    if values.is_null() {
        return ErrorCode::InvalidFormat;
    }

    for index in 0..simulation.network.nodes.len() {
        let value = unsafe { *values.add(index) };
        let code = set_node_value(simulation, index, property, value);
        if code != ErrorCode::Ok {
            unsafe { write_out(out_bad_index, (index + 1) as c_int) };
            return code;
        }
    }

    ErrorCode::Ok
}

/// Sets a junction's elevation, primary base demand and demand pattern in one call.
///
/// An empty or null `demand_pattern` clears the junction's demand pattern.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `demand_pattern` must be null or a valid pointer to a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setjuncdata(
    ph: *mut Project,
    index: c_int,
    elevation: c_double,
    demand: c_double,
    demand_pattern: *const c_char,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let node_id = match simulation.network.nodes.get((index - 1) as usize) {
        Some(node) => node.id.clone(),
        None => return ErrorCode::UndefinedNode,
    };

    let pattern = match unsafe { read_str(demand_pattern) } {
        Some(pattern) if !pattern.is_empty() => {
            if !simulation.network.pattern_map.contains_key(pattern) {
                return ErrorCode::UndefinedPattern;
            }
            Some(pattern.into())
        }
        _ => None,
    };

    match simulation.network.update_junction(
        &node_id,
        &JunctionUpdate {
            elevation: Some(elevation),
            basedemand: Some(demand),
            pattern: Some(pattern),
            ..Default::default()
        },
    ) {
        Ok(_) => ErrorCode::Ok,
        Err(_) => ErrorCode::IllegalNodeProperty,
    }
}

/// Sets a group of properties for a tank node.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `volume_curve` must be null or a valid pointer to a NUL-terminated C string.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn EN_settankdata(
    ph: *mut Project,
    index: c_int,
    elevation: c_double,
    initial_level: c_double,
    min_level: c_double,
    max_level: c_double,
    diameter: c_double,
    min_volume: c_double,
    volume_curve: *const c_char,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let node_id = match simulation.network.nodes.get((index - 1) as usize) {
        Some(node) => node.id.clone(),
        None => return ErrorCode::UndefinedNode,
    };

    if let Some(curve) = unsafe { read_str(volume_curve) }
        && !curve.is_empty()
    {
        // TODO: tanks with a volume curve are not supported by the solver.
        return ErrorCode::NotImplemented;
    }

    if min_level < 0.0
        || min_level > initial_level
        || max_level < initial_level
        || diameter < 0.0
        || min_volume < 0.0
    {
        return ErrorCode::IllegalNodeProperty;
    }

    match simulation.network.update_tank(
        &node_id,
        &TankUpdate {
            elevation: Some(elevation),
            initial_level: Some(initial_level),
            min_level: Some(min_level),
            max_level: Some(max_level),
            diameter: Some(diameter),
            min_volume: Some(min_volume),
            ..Default::default()
        },
    ) {
        Ok(_) => ErrorCode::Ok,
        Err(error) => match input_error_code(&error) {
            ErrorCode::NotATank => ErrorCode::NotATank,
            _ => ErrorCode::IllegalNodeProperty,
        },
    }
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_x` must be a valid non-null writable pointer.
/// `out_y` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getcoord(
    ph: *mut Project,
    index: c_int,
    out_x: *mut c_double,
    out_y: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    // EPANET indexes from 1, so we need to subtract 1 from the index
    let index = (index - 1) as usize;

    let node = match simulation.network.nodes.get(index) {
        Some(node) => node,
        None => return ErrorCode::UndefinedNode,
    };

    if let Some(coordinates) = node.coordinates {
        unsafe { *out_x = coordinates.0 as c_double };
        unsafe { *out_y = coordinates.1 as c_double };
    } else {
        return ErrorCode::NodeNoCoordinates;
    }

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setcoord(
    ph: *mut Project,
    index: c_int,
    x: c_double,
    y: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    // EPANET indexes from 1, so we need to subtract 1 from the index
    let index = (index - 1) as usize;

    let node = match simulation.network.nodes.get_mut(index) {
        Some(node) => node,
        None => return ErrorCode::UndefinedNode,
    };
    node.coordinates = Some((x, y));
    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_deletenode(
    ph: *mut Project,
    index: c_int,
    action_code: c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    // EPANET indexes from 1, so we need to subtract 1 from the index
    let index = (index - 1) as usize;

    let node_id = match simulation.network.nodes.get(index) {
        Some(node) => node.id.clone(),
        None => return ErrorCode::UndefinedNode,
    };

    let result = simulation.network.remove_node(&node_id, action_code == 0);
    if result.is_err() {
        return ErrorCode::DeleteNodeOrLinkInControl;
    }

    ErrorCode::Ok
}
