use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation};
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};

/// Gets the index of a node given its ID name.
#[unsafe(no_mangle)]
pub extern "C" fn EN_getnodeindex(ph: *mut Project, id: *const c_char, out_index: *mut c_int) -> ErrorCode {

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

  unsafe { *out_index = node_index as c_int };

  ErrorCode::Ok
}

// EN_getnodeid — Gets the ID name of a node given its index.

#[unsafe(no_mangle)]
pub extern "C" fn EN_getnodeid(ph: *mut Project, index: c_int, out_id: *mut c_char) -> ErrorCode {
  let simulation = get_simulation!(ph);

  let node_id = match simulation.network.nodes.get(index as usize) {
    Some(node) => node.id.as_ref(),
    None => return ErrorCode::UndefinedNode,
  };

  let c_str = CString::new(node_id).unwrap();

  unsafe {
    std::ptr::copy_nonoverlapping(c_str.as_ptr(), out_id, c_str.as_bytes_with_nul().len());
  }

  ErrorCode::Ok
}
