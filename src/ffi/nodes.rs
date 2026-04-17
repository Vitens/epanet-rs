use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::model::node::NodeType;

use crate::ffi::enums::NodeType as ENNodeType;

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

  // EPANET indexes from 1, so we need to add 1 to the index
  unsafe { *out_index = (node_index+1) as c_int };

  ErrorCode::Ok
}

// Gets the ID name of a node given its index.
#[unsafe(no_mangle)]
pub extern "C" fn EN_getnodeid(ph: *mut Project, index: c_int, out_id: *mut c_char) -> ErrorCode {
  let simulation = get_simulation!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = index-1;

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

// Sets the ID name of a node given its index.
#[unsafe(no_mangle)]
pub extern "C" fn EN_setnodeid(ph: *mut Project, index: c_int, id: *const c_char) -> ErrorCode {
  let simulation = get_simulation_mut!(ph);

  let c_str = unsafe { CStr::from_ptr(id) };

  let new_node_id = match c_str.to_str() {
    Ok(s) => s,
    Err(_) => return ErrorCode::InvalidIdName,
  };

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = index-1;

  let node = match simulation.network.nodes.get_mut(index as usize) {
    Some(node) => node,
    None => return ErrorCode::UndefinedNode,
  };

  // check if the new node id is already in use
  if simulation.network.node_map.contains_key(new_node_id) {
    return ErrorCode::DuplicateId;
  }

  // remove the old node id from the node map
  simulation.network.node_map.remove(&node.id);

  // update the node id
  node.id = new_node_id.into();

  // update the node map
  simulation.network.node_map.insert(new_node_id.into(), index as usize);

  // update all links that point to the old node id to point to the new node id
  for link in simulation.network.links.iter_mut() {
    if link.start_node_id == node.id {
      link.start_node_id = new_node_id.into();
    }
    if link.end_node_id == node.id {
      link.end_node_id = new_node_id.into();
    }
  }

  ErrorCode::Ok
}

// Get the node type given its index.
#[unsafe(no_mangle)]
pub extern "C" fn EN_getnodetype(ph: *mut Project, index: c_int, out_type: *mut c_int) -> ErrorCode {
  let simulation = get_simulation!(ph);

  // EPANET indexes from 1, so we need to add 1 to the index
  let index = index-1;

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