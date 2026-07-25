//! FFI curve accessors: `EN_addcurve`, `EN_getcurvevalue`, `EN_setcurvevalue`, etc.

use crate::ffi::enums::CurveType;
use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::ffi::util::{MAX_ID, write_out, write_str};

use crate::model::link::LinkType;
use crate::model::network::modify::{CurveData, CurveUpdate};
use crate::model::node::NodeType;
use crate::model::valve::ValveType;

use std::ffi::CStr;
use std::os::raw::{c_char, c_double, c_int};

// Adds a new curve to the network.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_addcurve(ph: *mut Project, id: *const c_char) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let c_str = unsafe { CStr::from_ptr(id) };
    let curve_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    let result = simulation.network.add_curve(
        curve_id,
        &CurveData {
            x: vec![1.0],
            y: vec![1.0],
        },
    );

    if result.is_err() {
        return ErrorCode::DuplicateId;
    }

    ErrorCode::Ok
}

/// Gets the index of a curve given its ID name.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
/// `out_index` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getcurveindex(
    ph: *mut Project,
    id: *const c_char,
    out_index: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let c_str = unsafe { CStr::from_ptr(id) };
    let curve_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    let curve_index = match simulation.network.curve_map.get(curve_id) {
        Some(&index) => index,
        None => return ErrorCode::UndefinedCurve,
    };

    // EPANET indexes from 1
    unsafe { *out_index = (curve_index + 1) as c_int };

    ErrorCode::Ok
}

/// Gets the ID name of a curve given its index.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_id` must point to a buffer large enough for the result string including NUL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getcurveid(
    ph: *mut Project,
    index: c_int,
    out_id: *mut c_char,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    // EPANET indexes from 1, so we need to subtract 1 from the index
    let index = index - 1;

    let curve_id = match simulation.network.curves.get(index as usize) {
        Some(curve) => curve.id.as_ref(),
        None => return ErrorCode::UndefinedCurve,
    };

    unsafe { write_str(out_id, curve_id, MAX_ID) };

    ErrorCode::Ok
}

/// Deletes a curve from the network. Fails while the curve is still referenced
/// by a pump, valve or tank.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_deletecurve(ph: *mut Project, index: c_int) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let curve_id = match simulation.network.curves.get((index - 1) as usize) {
        Some(curve) => curve.id.clone(),
        None => return ErrorCode::UndefinedCurve,
    };

    match simulation.network.remove_curve(&curve_id) {
        Ok(()) => ErrorCode::Ok,
        Err(_) => ErrorCode::IllegalLinkProperty,
    }
}

/// Returns how a curve is used within the network.
///
/// `epanet-rs` does not store an explicit curve type, so the type is derived
/// from the objects that reference the curve; unreferenced curves report
/// [`CurveType::Generic`].
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getcurvetype(
    ph: *mut Project,
    index: c_int,
    out_type: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let curve_id = match simulation.network.curves.get((index - 1) as usize) {
        Some(curve) => curve.id.clone(),
        None => return ErrorCode::UndefinedCurve,
    };

    let mut curve_type = CurveType::Generic;

    for link in simulation.network.links.iter() {
        match &link.link_type {
            LinkType::Pump(pump) if pump.head_curve_id.as_deref() == Some(&*curve_id) => {
                curve_type = CurveType::Pump;
            }
            LinkType::Valve(valve) if valve.curve_id.as_deref() == Some(&*curve_id) => {
                curve_type = match valve.valve_type {
                    ValveType::GPV => CurveType::HLoss,
                    _ => CurveType::Valve,
                };
            }
            _ => continue,
        }
        break;
    }

    if curve_type == CurveType::Generic {
        let used_by_tank = simulation
            .network
            .nodes
            .iter()
            .any(|node| match &node.node_type {
                NodeType::Tank(tank) => tank.volume_curve_id.as_deref() == Some(&*curve_id),
                _ => false,
            });
        if used_by_tank {
            curve_type = CurveType::Volume;
        }
    }

    unsafe { write_out(out_type, curve_type as c_int) };

    ErrorCode::Ok
}

/// Assigns a type to a curve.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setcurvetype(
    _ph: *mut Project,
    _index: c_int,
    _curve_type: c_int,
) -> ErrorCode {
    // TODO: curves in epanet-rs have no explicit type; their meaning follows
    // from the pump/valve/tank that references them.
    ErrorCode::NotImplemented
}

/// Sets the ID name of a curve given its index.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `id` must be a valid non-null pointer to a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setcurveid(
    ph: *mut Project,
    index: c_int,
    id: *const c_char,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let c_str = unsafe { CStr::from_ptr(id) };
    let new_curve_id = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::InvalidIdName,
    };

    // EPANET indexes from 1, so we need to subtract 1 from the index
    let index = index - 1;

    let curve = match simulation.network.curves.get_mut(index as usize) {
        Some(curve) => curve,
        None => return ErrorCode::UndefinedCurve,
    };

    // check if the new curve id is already in use
    if simulation.network.curve_map.contains_key(new_curve_id) {
        return ErrorCode::DuplicateId;
    }

    // remove the old curve id from the curve map
    simulation.network.curve_map.remove(&curve.id);

    // update the curve id
    curve.id = new_curve_id.into();

    // update the curve map
    simulation
        .network
        .curve_map
        .insert(new_curve_id.into(), index as usize);

    ErrorCode::Ok
}

// Get the length of a curve given its index.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_count` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getcurvelen(
    ph: *mut Project,
    index: c_int,
    out_count: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    // EPANET indexes from 1, so we need to subtract 1 from the index
    let index = index - 1;

    let curve = match simulation.network.curves.get(index as usize) {
        Some(curve) => curve,
        None => return ErrorCode::UndefinedCurve,
    };

    let count = curve.x.len();

    unsafe { *out_count = count as c_int };

    ErrorCode::Ok
}

// Retrieve the value of a single point on a curve
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_x` must be a valid non-null writable pointer.
/// `out_y` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getcurvevalue(
    ph: *mut Project,
    index: c_int,
    point_index: c_int,
    out_x: *mut c_double,
    out_y: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    // EPANET indexes from 1, so we need to subtract 1 from the index
    let index = index - 1;

    let curve = match simulation.network.curves.get(index as usize) {
        Some(curve) => curve,
        None => return ErrorCode::UndefinedCurve,
    };

    let point_index = (point_index - 1) as usize;

    if point_index >= curve.x.len() || curve.x.is_empty() {
        return ErrorCode::UndefinedCurve;
    }

    let x = curve.x[point_index];
    let y = curve.y[point_index];

    unsafe { write_out(out_x, x as c_double) };
    unsafe { write_out(out_y, y as c_double) };

    ErrorCode::Ok
}

/// Sets the value of a single data point on a curve.
///
/// Any pump or valve using the curve has its derived coefficients rebuilt.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setcurvevalue(
    ph: *mut Project,
    index: c_int,
    point_index: c_int,
    x: c_double,
    y: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let curve = match simulation.network.curves.get((index - 1) as usize) {
        Some(curve) => curve,
        None => return ErrorCode::UndefinedCurve,
    };

    let point_index = (point_index - 1) as usize;
    if point_index >= curve.x.len() {
        return ErrorCode::InvalidParameterCode;
    }

    let curve_id = curve.id.clone();
    let mut new_x = curve.x.clone();
    let mut new_y = curve.y.clone();
    new_x[point_index] = x;
    new_y[point_index] = y;

    // go through update_curve so any derived pump/valve curve is rebuilt
    match simulation.network.update_curve(
        &curve_id,
        &CurveUpdate {
            x: Some(new_x),
            y: Some(new_y),
        },
    ) {
        Ok(()) => ErrorCode::Ok,
        Err(_) => ErrorCode::InvalidParameterCode,
    }
}

/// Retrieves all of a curve's data.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_id` must point to a buffer of at least `EN_MAXID + 1` bytes.
/// `out_x` and `out_y` must each point to a buffer able to hold the curve's
/// points (use [`EN_getcurvelen`] to size them).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getcurve(
    ph: *mut Project,
    index: c_int,
    out_id: *mut c_char,
    out_n_points: *mut c_int,
    out_x: *mut c_double,
    out_y: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    // EPANET indexes from 1, so we need to subtract 1 from the index
    let index = index - 1;

    let curve = match simulation.network.curves.get(index as usize) {
        Some(curve) => curve,
        None => return ErrorCode::UndefinedCurve,
    };

    unsafe { write_str(out_id, &curve.id, MAX_ID) };
    unsafe { write_out(out_n_points, curve.x.len() as c_int) };

    if !out_x.is_null() {
        unsafe { std::ptr::copy_nonoverlapping(curve.x.as_ptr(), out_x, curve.x.len()) };
    }
    if !out_y.is_null() {
        unsafe { std::ptr::copy_nonoverlapping(curve.y.as_ptr(), out_y, curve.y.len()) };
    }

    ErrorCode::Ok
}

// Sets the values of a curve given its index.
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `x` must point to at least `count` readable `c_double` values.
/// `y` must point to at least `count` readable `c_double` values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setcurve(
    ph: *mut Project,
    index: c_int,
    x: *const c_double,
    y: *const c_double,
    count: c_int,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let index = index - 1;

    let curve_id = match simulation.network.curves.get_mut(index as usize) {
        Some(curve) => curve.id.clone(),
        None => return ErrorCode::UndefinedCurve,
    };

    if count < 1 || x.is_null() || y.is_null() {
        return ErrorCode::InvalidParameterCode;
    }

    let x = unsafe { std::slice::from_raw_parts(x, count as usize).to_vec() };
    let y = unsafe { std::slice::from_raw_parts(y, count as usize).to_vec() };

    let result = simulation.network.update_curve(
        &curve_id,
        &CurveUpdate {
            x: Some(x),
            y: Some(y),
        },
    );

    if result.is_err() {
        return ErrorCode::InvalidParameterCode;
    }

    ErrorCode::Ok
}
