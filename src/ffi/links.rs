use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::model::link::LinkType;
use crate::model::valve::ValveType;

use crate::ffi::enums::LinkType as ENLinkType;

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};

/// Gets the index of a node given its ID name.
#[unsafe(no_mangle)]
pub extern "C" fn EN_getlinkindex(ph: *mut Project, id: *const c_char, out_index: *mut c_int) -> ErrorCode {

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
  unsafe { *out_index = (link_index+1) as c_int };

  ErrorCode::Ok
}

// Gets the ID name of a node given its index.
#[unsafe(no_mangle)]
pub extern "C" fn EN_getlinkid(ph: *mut Project, index: c_int, out_id: *mut c_char) -> ErrorCode {
  let simulation = get_simulation!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = index-1;

  let link_id = match simulation.network.links.get(index as usize) {
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
#[unsafe(no_mangle)]
pub extern "C" fn EN_setlinkid(ph: *mut Project, index: c_int, id: *const c_char) -> ErrorCode {
  let simulation = get_simulation_mut!(ph);

  let c_str = unsafe { CStr::from_ptr(id) };

  let new_link_id = match c_str.to_str() {
    Ok(s) => s,
    Err(_) => return ErrorCode::InvalidIdName,
  };

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = index-1;

  let link = match simulation.network.links.get_mut(index as usize) {
    Some(link) => link,
    None => return ErrorCode::UndefinedLink,
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
  simulation.network.link_map.insert(new_link_id.into(), index as usize);

  ErrorCode::Ok
}

// Get the link type given its index.
#[unsafe(no_mangle)]
pub extern "C" fn EN_getlinktype(ph: *mut Project, index: c_int, out_type: *mut c_int) -> ErrorCode {
  let simulation = get_simulation!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = index-1;

  let link_type = match simulation.network.links.get(index as usize) {
    Some(link) => &link.link_type,
    None => return ErrorCode::UndefinedNode,
  };


  let link_type_int = match link_type {
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
#[unsafe(no_mangle)]
pub extern "C" fn EN_getlinknodes(ph: *mut Project, index: c_int, out_start_node: *mut c_int, out_end_node: *mut c_int) -> ErrorCode {
  let simulation = get_simulation!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = index-1;

  let link = match simulation.network.links.get(index as usize) {
    Some(link) => link,
    None => return ErrorCode::UndefinedLink,
  };

  unsafe { *out_start_node = (link.start_node+1) as c_int };
  unsafe { *out_end_node = (link.end_node+1) as c_int };

  ErrorCode::Ok
}