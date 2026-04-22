use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::ffi::enums::NodeProperty;
use crate::model::node::NodeType;
use crate::model::network::modify::{NodeUpdate, JunctionUpdate, TankUpdate, ReservoirUpdate};

use crate::ffi::enums::NodeType as ENNodeType;

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_double};

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

/// Retrieves the property value of a node
#[unsafe(no_mangle)]
pub extern "C" fn EN_getnodevalue(ph: *mut Project, index: c_int, property: c_int, out_value: *mut c_double) -> ErrorCode {
  let simulation = get_simulation!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = (index-1) as usize;

  let node = match simulation.network.nodes.get(index) {
    Some(node) => node,
    None => return ErrorCode::UndefinedNode,
  };

  let options = &simulation.network.options;
  let unit_system = &options.unit_system;
  let flow_units = &options.flow_units;
  let pressure_units = &options.pressure_units;

  let property = match NodeProperty::from_repr(property) {
    Some(property) => property,
    None => return ErrorCode::InvalidParameterCode,
  };

  let value = match property {
    NodeProperty::Elevation => node.elevation * unit_system.per_feet(),
    NodeProperty::BaseDemand => match &node.node_type {
      NodeType::Junction(junction) => junction.basedemand * flow_units.per_cfs(),
      _ => 0.0,
    },
    NodeProperty::Pattern => match &node.node_type {
      NodeType::Junction(junction) => junction.pattern_index.map(|index| (index+1) as f64).unwrap_or(0.0),
      NodeType::Reservoir(reservoir) => reservoir.head_pattern_index.map(|index| (index+1) as f64).unwrap_or(0.0),
      _ => 0.0,
    },
    NodeProperty::Emitter => match &node.node_type {
      NodeType::Junction(junction) => {
        let e = junction.emitter_coefficient;
        if e > 0.0 {
          flow_units.per_cfs() / (pressure_units.per_feet() * e).powf(1.0 / options.emitter_exponent)
        } else {
          0.0
        }
      }
      _ => 0.0,
    },
    NodeProperty::InitQual => 0.0,  // TODO: quality not implemented yet
    NodeProperty::SourceQual => 0.0,  // TODO: quality not implemented yet
    NodeProperty::SourcePat => 0.0,  // TODO: quality not implemented yet
    NodeProperty::SourceType => 0.0,  // TODO: quality not implemented yet

    NodeProperty::TankLevel => match &node.node_type {
      NodeType::Tank(_) => simulation.solved_state().map_or(0.0, |state| state.heads[index] - node.elevation) * unit_system.per_feet(),
      _ => 0.0,
    },

    NodeProperty::Demand => simulation.solved_state().map_or(0.0, |state| state.demands[index] * flow_units.per_cfs()),
    NodeProperty::Head => simulation.solved_state().map_or(0.0, |state| state.heads[index] * unit_system.per_feet()),
    NodeProperty::Pressure => simulation.solved_state().map_or(0.0, |state| (state.heads[index] - node.elevation) * pressure_units.per_feet()),
    NodeProperty::Quality => 0.0,  // TODO: quality not implemented yet
    NodeProperty::SourceMass => 0.0,  // TODO: mass not implemented yet
    _ => 0.0,
  };

  unsafe { *out_value = value as c_double };

  ErrorCode::Ok
}

// Set the property value of a node
#[unsafe(no_mangle)]
pub extern "C" fn EN_setnodevalue(ph: *mut Project, index: c_int, property: c_int, value: c_double) -> ErrorCode {
  let simulation = get_simulation_mut!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = (index-1) as usize;

  let network = &mut simulation.network;

  let node_id = match network.nodes.get(index) {
    Some(node) => node.id.clone(),
    None => return ErrorCode::UndefinedNode,
  };

  let property = match NodeProperty::from_repr(property) {
    Some(property) => property,
    None => return ErrorCode::InvalidParameterCode,
  };


  let result = match property {
    NodeProperty::Elevation => {
      network.update_node(&node_id, &NodeUpdate { elevation: Some(value as f64), ..Default::default() })
    }
    NodeProperty::BaseDemand => {
      network.update_junction(&node_id, &JunctionUpdate { basedemand: Some(value as f64), ..Default::default() })
    }
    NodeProperty::Emitter => {
      network.update_junction(&node_id, &JunctionUpdate { emitter_coefficient: Some(value as f64), ..Default::default() })
    }
    NodeProperty::TankDiam => {
      network.update_tank(&node_id, &TankUpdate { diameter: Some(value as f64), ..Default::default() })
    }
    NodeProperty::MinLevel => {
      network.update_tank(&node_id, &TankUpdate { min_level: Some(value as f64), ..Default::default() })
    }
    NodeProperty::MaxLevel => {
      network.update_tank(&node_id, &TankUpdate { max_level: Some(value as f64), ..Default::default() })
    }
    // TODO: implement missing properties
    _ => return ErrorCode::InvalidParameterCode,
  };

  if result.is_err() {
    return ErrorCode::IllegalNodeProperty;
  }

  ErrorCode::Ok
}