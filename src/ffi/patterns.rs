//! FFI pattern accessors: `EN_addpattern`, `EN_getpatternvalue`, `EN_setpatternvalue`, etc.

use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::ffi::util::{MAX_ID, write_out, write_str};

use std::ffi::CStr;
use std::os::raw::{c_char, c_double, c_int};

use crate::model::network::modify::PatternData;
use crate::model::node::NodeType;

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_addpattern(ph: *mut Project, id: *const c_char) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let c_str = unsafe { CStr::from_ptr(id) };
    let pattern_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    // add the pattern to the network
    let result = simulation.network.add_pattern(
        pattern_id,
        &PatternData {
            multipliers: vec![1.0],
        },
    );

    // return error if the pattern already exists
    if result.is_err() {
        return ErrorCode::DuplicateId;
    }

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_deletepattern(ph: *mut Project, index: c_int) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let pattern_id = match simulation.network.patterns.get((index - 1) as usize) {
        Some(pattern) => pattern.id.clone(),
        None => return ErrorCode::UndefinedPattern,
    };

    // delete the pattern from the network
    let result = simulation.network.remove_pattern(&pattern_id);

    // return error if the pattern is not found
    if result.is_err() {
        return ErrorCode::UndefinedPattern;
    }

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_value` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getaveragepatternvalue(
    ph: *mut Project,
    index: c_int,
    out_value: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let index = (index - 1) as usize;

    let pattern = match simulation.network.patterns.get(index) {
        Some(pattern) => pattern,
        None => return ErrorCode::UndefinedPattern,
    };

    let value = pattern.multipliers.iter().sum::<f64>() / pattern.multipliers.len() as f64;

    unsafe { *out_value = value as c_double };

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_id` must point to a buffer large enough for the result string including NUL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getpatternid(
    ph: *mut Project,
    index: c_int,
    out_id: *mut c_char,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let index = (index - 1) as usize;

    let pattern_id = match simulation.network.patterns.get(index) {
        Some(pattern) => pattern.id.clone(),
        None => return ErrorCode::UndefinedPattern,
    };

    unsafe { write_str(out_id, &pattern_id, MAX_ID) };

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setpatternid(
    ph: *mut Project,
    index: c_int,
    id: *const c_char,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let index = (index - 1) as usize;

    let pattern = match simulation.network.patterns.get_mut(index) {
        Some(pattern) => pattern,
        None => return ErrorCode::UndefinedPattern,
    };

    let old_pattern_id = pattern.id.clone();

    let c_str = unsafe { CStr::from_ptr(id) };
    let new_pattern_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    // check if the new pattern id is already in use
    if simulation.network.pattern_map.contains_key(new_pattern_id) {
        return ErrorCode::DuplicateId;
    }
    // remove the old pattern id from the pattern map
    simulation.network.pattern_map.remove(&pattern.id);

    // update the pattern id
    pattern.id = new_pattern_id.into();

    // update the pattern map
    simulation
        .network
        .pattern_map
        .insert(new_pattern_id.into(), index);

    // update all nodes that point to the old pattern id to point to the new pattern id
    for node in simulation.network.nodes.iter_mut() {
        match &mut node.node_type {
            NodeType::Junction(junction) => {
                for demand in junction.demands.iter_mut() {
                    if demand.pattern.as_deref() == Some(&old_pattern_id) {
                        demand.pattern = Some(new_pattern_id.into());
                    }
                }
            }
            NodeType::Reservoir(reservoir)
                if reservoir.head_pattern.as_deref() == Some(&old_pattern_id) =>
            {
                reservoir.head_pattern = Some(new_pattern_id.into());
            }
            _ => {}
        }
    }

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
/// `out_index` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getpatternindex(
    ph: *mut Project,
    id: *const c_char,
    out_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let c_str = unsafe { CStr::from_ptr(id) };
    let pattern_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    let index = match simulation.network.pattern_map.get(pattern_id) {
        Some(&index) => index,
        None => return ErrorCode::UndefinedPattern,
    };

    unsafe { *out_index = (index + 1) as c_int };

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_count` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getpatternlen(
    ph: *mut Project,
    index: c_int,
    out_count: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let index = (index - 1) as usize;

    let pattern = match simulation.network.patterns.get(index) {
        Some(pattern) => pattern,
        None => return ErrorCode::UndefinedPattern,
    };

    unsafe { *out_count = pattern.multipliers.len() as c_int };
    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_value` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getpatternvalue(
    ph: *mut Project,
    index: c_int,
    period: c_int,
    out_value: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let index = (index - 1) as usize;

    let pattern = match simulation.network.patterns.get(index) {
        Some(pattern) => pattern,
        None => return ErrorCode::UndefinedPattern,
    };

    if period < 1 || period as usize > pattern.multipliers.len() {
        return ErrorCode::InvalidParameterCode;
    }

    let value = pattern.multipliers[(period - 1) as usize];

    unsafe { write_out(out_value, value as c_double) };
    ErrorCode::Ok
}

/// Sets a time pattern's multiplier for a given period.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setpatternvalue(
    ph: *mut Project,
    index: c_int,
    period: c_int,
    value: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let index = (index - 1) as usize;

    let pattern = match simulation.network.patterns.get_mut(index) {
        Some(pattern) => pattern,
        None => return ErrorCode::UndefinedPattern,
    };

    if period < 1 || period as usize > pattern.multipliers.len() {
        return ErrorCode::InvalidParameterCode;
    }

    pattern.multipliers[(period - 1) as usize] = value;

    ErrorCode::Ok
}

/// Loads a set of multipliers from a file into a time pattern.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_loadpatternfile(
    _ph: *mut Project,
    _filename: *const c_char,
    _id: *const c_char,
) -> ErrorCode {
    // TODO: epanet-rs does not read EPANET's standalone pattern file format.
    ErrorCode::NotImplemented
}

///
/// # Safety
///
/// `ph` must be a valid project handle. `multipliers` must point to `count`
/// readable `c_double` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setpattern(
    ph: *mut Project,
    index: c_int,
    multipliers: *const c_double,
    count: c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let index = (index - 1) as usize;

    let pattern = match simulation.network.patterns.get_mut(index) {
        Some(pattern) => pattern,
        None => return ErrorCode::UndefinedPattern,
    };

    if count < 1 || multipliers.is_null() {
        return ErrorCode::InvalidParameterCode;
    }

    pattern.multipliers =
        unsafe { std::slice::from_raw_parts(multipliers, count as usize).to_vec() };
    ErrorCode::Ok
}
