//! FFI accessors for a junction's demand categories.

use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::ffi::util::{MAX_ID, read_str, write_out, write_str};

use crate::model::demand::Demand;
use crate::model::network::Network;
use crate::model::node::NodeType;

use std::os::raw::{c_char, c_double, c_int};

/// Returns the demand list of the junction at `index` (0-based).
fn demands_of(network: &Network, index: usize) -> Result<&Vec<Demand>, ErrorCode> {
    let node = network.nodes.get(index).ok_or(ErrorCode::UndefinedNode)?;
    match &node.node_type {
        NodeType::Junction(junction) => Ok(&junction.demands),
        // only junctions carry demand categories
        _ => Err(ErrorCode::IllegalNodeProperty),
    }
}

fn demands_of_mut(network: &mut Network, index: usize) -> Result<&mut Vec<Demand>, ErrorCode> {
    let node = network
        .nodes
        .get_mut(index)
        .ok_or(ErrorCode::UndefinedNode)?;
    match &mut node.node_type {
        NodeType::Junction(junction) => Ok(&mut junction.demands),
        _ => Err(ErrorCode::IllegalNodeProperty),
    }
}

/// Retrieves the number of demand categories for a junction.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getnumdemands(
    ph: *mut Project,
    node_index: c_int,
    out_num_demands: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    match demands_of(&simulation.network, (node_index - 1) as usize) {
        Ok(demands) => {
            unsafe { write_out(out_num_demands, demands.len() as c_int) };
            ErrorCode::Ok
        }
        Err(code) => code,
    }
}

/// Retrieves the base demand of one of a junction's demand categories.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getbasedemand(
    ph: *mut Project,
    node_index: c_int,
    demand_index: c_int,
    out_base_demand: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let flow_factor = simulation.network.options.flow_units.per_cfs();

    let demands = match demands_of(&simulation.network, (node_index - 1) as usize) {
        Ok(demands) => demands,
        Err(code) => return code,
    };

    let Some(demand) = demands.get((demand_index - 1) as usize) else {
        return ErrorCode::NonexistentDemandCategory;
    };

    unsafe { write_out(out_base_demand, demand.basedemand * flow_factor) };

    ErrorCode::Ok
}

/// Sets the base demand of one of a junction's demand categories.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setbasedemand(
    ph: *mut Project,
    node_index: c_int,
    demand_index: c_int,
    base_demand: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let flow_factor = simulation.network.options.flow_units.per_cfs();
    let node_index = (node_index - 1) as usize;

    let demands = match demands_of_mut(&mut simulation.network, node_index) {
        Ok(demands) => demands,
        Err(code) => return code,
    };

    let Some(demand) = demands.get_mut((demand_index - 1) as usize) else {
        return ErrorCode::NonexistentDemandCategory;
    };

    demand.basedemand = base_demand / flow_factor;

    simulation.network.updated_nodes.insert(node_index);
    simulation.network.properties_version += 1;

    ErrorCode::Ok
}

/// Retrieves the index of the time pattern used by a demand category.
/// A value of 0 means the category has no pattern assigned.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getdemandpattern(
    ph: *mut Project,
    node_index: c_int,
    demand_index: c_int,
    out_pattern_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let demands = match demands_of(&simulation.network, (node_index - 1) as usize) {
        Ok(demands) => demands,
        Err(code) => return code,
    };

    let Some(demand) = demands.get((demand_index - 1) as usize) else {
        return ErrorCode::NonexistentDemandCategory;
    };

    let pattern_index = demand
        .pattern
        .as_deref()
        .and_then(|id| simulation.network.pattern_map.get(id))
        .map(|&i| i + 1)
        .unwrap_or(0);

    unsafe { write_out(out_pattern_index, pattern_index as c_int) };

    ErrorCode::Ok
}

/// Assigns a time pattern to a demand category. Pattern index 0 clears the
/// pattern.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setdemandpattern(
    ph: *mut Project,
    node_index: c_int,
    demand_index: c_int,
    pattern_index: c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let pattern = if pattern_index == 0 {
        None
    } else {
        match simulation
            .network
            .patterns
            .get((pattern_index - 1) as usize)
        {
            Some(pattern) => Some((pattern.id.clone(), (pattern_index - 1) as usize)),
            None => return ErrorCode::UndefinedPattern,
        }
    };

    let node_index = (node_index - 1) as usize;

    let demands = match demands_of_mut(&mut simulation.network, node_index) {
        Ok(demands) => demands,
        Err(code) => return code,
    };

    let Some(demand) = demands.get_mut((demand_index - 1) as usize) else {
        return ErrorCode::NonexistentDemandCategory;
    };

    match pattern {
        Some((id, index)) => {
            demand.pattern = Some(id);
            demand.pattern_index = Some(index);
        }
        None => {
            demand.pattern = None;
            demand.pattern_index = None;
        }
    }

    simulation.network.updated_nodes.insert(node_index);
    simulation.network.properties_version += 1;

    ErrorCode::Ok
}

/// Retrieves the name assigned to a demand category.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_name` must point to a buffer of at least `EN_MAXID + 1` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getdemandname(
    ph: *mut Project,
    node_index: c_int,
    demand_index: c_int,
    out_name: *mut c_char,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let demands = match demands_of(&simulation.network, (node_index - 1) as usize) {
        Ok(demands) => demands,
        Err(code) => return code,
    };

    let Some(demand) = demands.get((demand_index - 1) as usize) else {
        return ErrorCode::NonexistentDemandCategory;
    };

    unsafe { write_str(out_name, demand.name.as_deref().unwrap_or(""), MAX_ID) };

    ErrorCode::Ok
}

/// Assigns a name to a demand category.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `name` must be a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setdemandname(
    ph: *mut Project,
    node_index: c_int,
    demand_index: c_int,
    name: *const c_char,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let Some(name) = (unsafe { read_str(name) }) else {
        return ErrorCode::InvalidIdName;
    };
    let name: Box<str> = name.into();

    let demands = match demands_of_mut(&mut simulation.network, (node_index - 1) as usize) {
        Ok(demands) => demands,
        Err(code) => return code,
    };

    let Some(demand) = demands.get_mut((demand_index - 1) as usize) else {
        return ErrorCode::NonexistentDemandCategory;
    };

    demand.name = if name.is_empty() { None } else { Some(name) };

    ErrorCode::Ok
}

/// Retrieves the index of a demand category given its name.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `demand_name` must be a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getdemandindex(
    ph: *mut Project,
    node_index: c_int,
    demand_name: *const c_char,
    out_demand_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let Some(name) = (unsafe { read_str(demand_name) }) else {
        return ErrorCode::InvalidIdName;
    };

    let demands = match demands_of(&simulation.network, (node_index - 1) as usize) {
        Ok(demands) => demands,
        Err(code) => return code,
    };

    let index = demands
        .iter()
        .position(|demand| demand.name.as_deref().unwrap_or("") == name);

    match index {
        Some(index) => {
            unsafe { write_out(out_demand_index, (index + 1) as c_int) };
            ErrorCode::Ok
        }
        None => ErrorCode::NonexistentDemandCategory,
    }
}

/// Appends a new demand category to a junction.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `demand_name` must either be null or a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_adddemand(
    ph: *mut Project,
    node_index: c_int,
    base_demand: c_double,
    demand_pattern: *const c_char,
    demand_name: *const c_char,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let flow_factor = simulation.network.options.flow_units.per_cfs();

    // resolve the pattern before borrowing the node mutably
    let pattern = match unsafe { read_str(demand_pattern) } {
        None | Some("") => None,
        Some(pattern_id) => match simulation.network.pattern_map.get(pattern_id) {
            Some(&index) => Some((Box::<str>::from(pattern_id), index)),
            None => return ErrorCode::UndefinedPattern,
        },
    };

    let name = unsafe { read_str(demand_name) }
        .filter(|name| !name.is_empty())
        .map(Box::<str>::from);

    let node_index = (node_index - 1) as usize;

    let demands = match demands_of_mut(&mut simulation.network, node_index) {
        Ok(demands) => demands,
        Err(code) => return code,
    };

    demands.push(Demand {
        basedemand: base_demand / flow_factor,
        pattern: pattern.as_ref().map(|(id, _)| id.clone()),
        pattern_index: pattern.as_ref().map(|(_, index)| *index),
        name,
    });

    simulation.network.updated_nodes.insert(node_index);
    simulation.network.properties_version += 1;

    ErrorCode::Ok
}

/// Deletes a demand category from a junction.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_deletedemand(
    ph: *mut Project,
    node_index: c_int,
    demand_index: c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let node_index = (node_index - 1) as usize;

    let demands = match demands_of_mut(&mut simulation.network, node_index) {
        Ok(demands) => demands,
        Err(code) => return code,
    };

    let demand_index = (demand_index - 1) as usize;
    if demand_index >= demands.len() {
        return ErrorCode::NonexistentDemandCategory;
    }

    demands.remove(demand_index);

    simulation.network.updated_nodes.insert(node_index);
    simulation.network.properties_version += 1;

    ErrorCode::Ok
}
