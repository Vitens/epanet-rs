use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::model::node::NodeType;

use crate::ffi::enums::NodeType as ENNodeType;

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};

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

