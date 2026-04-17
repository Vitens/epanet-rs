use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_double};

/// Gets the index of a curve given its ID name.
#[unsafe(no_mangle)]
pub extern "C" fn EN_getcurveindex(ph: *mut Project, id: *const c_char, out_index: *mut c_int) -> ErrorCode {
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
#[unsafe(no_mangle)]
pub extern "C" fn EN_getcurveid(ph: *mut Project, index: c_int, out_id: *mut c_char) -> ErrorCode {
  let simulation = get_simulation!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = index-1;

  let curve_id = match simulation.network.curves.get(index as usize) {
    Some(curve) => curve.id.as_ref(),
    None => return ErrorCode::UndefinedCurve,
  };

  let c_str = CString::new(curve_id).unwrap();

  unsafe {
    std::ptr::copy_nonoverlapping(c_str.as_ptr(), out_id, c_str.as_bytes_with_nul().len());
  }

  ErrorCode::Ok
}

/// Sets the ID name of a curve given its index.
#[unsafe(no_mangle)]
pub extern "C" fn EN_setcurveid(ph: *mut Project, index: c_int, id: *const c_char) -> ErrorCode {
  let simulation = get_simulation_mut!(ph);

  let c_str = unsafe { CStr::from_ptr(id) };
  let new_curve_id = match c_str.to_str() {
    Ok(s) => s,
    Err(_) => return ErrorCode::InvalidIdName,
  };

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = index-1;

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
  simulation.network.curve_map.insert(new_curve_id.into(), index as usize);
  
  ErrorCode::Ok
}

// Get the length of a curve given its index.
#[unsafe(no_mangle)]
pub extern "C" fn EN_getcurvelen(ph: *mut Project, index: c_int, out_count: *mut c_int) -> ErrorCode {
  let simulation = get_simulation!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = index-1;

  let curve = match simulation.network.curves.get(index as usize) {
    Some(curve) => curve,
    None => return ErrorCode::UndefinedCurve,
  };

  let count = curve.x.len();

  unsafe { *out_count = count as c_int };

  ErrorCode::Ok
}

// Retrieve the value of a single point on a curve
#[unsafe(no_mangle)]
pub extern "C" fn EN_getcurvevalue(ph: *mut Project, index: c_int, point_index: c_int, out_x: *mut c_double, out_y: *mut c_double) -> ErrorCode {

  let simulation = get_simulation!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = index-1;

  let curve = match simulation.network.curves.get(index as usize) {
    Some(curve) => curve,
    None => return ErrorCode::UndefinedCurve,
  };

  let point_index = (point_index-1) as usize;

  if point_index >= curve.x.len() || curve.x.is_empty() {
    return ErrorCode::UndefinedCurve;
  }

  let x = curve.x[point_index as usize];
  let y = curve.y[point_index as usize];

  unsafe { *out_x = x as c_double };
  unsafe { *out_y = y as c_double };

  ErrorCode::Ok

}

// Retrieve all of a curve's points
#[unsafe(no_mangle)]
pub extern "C" fn EN_getcurve(ph: *mut Project, index: c_int, out_n_points: *mut c_int, out_x: *mut *mut c_double, out_y: *mut *mut c_double) -> ErrorCode {

  let simulation = get_simulation!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = index-1;

  let curve = match simulation.network.curves.get(index as usize) {
    Some(curve) => curve,
    None => return ErrorCode::UndefinedCurve,
  };

  let n_points = curve.x.len();

  unsafe { *out_n_points = n_points as c_int };

  unsafe { *out_x = curve.x.as_ptr() as *mut c_double };
  unsafe { *out_y = curve.y.as_ptr() as *mut c_double };

  ErrorCode::Ok
}
