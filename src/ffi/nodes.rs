use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::ffi::enums::NodeProperty;
use crate::model::node::NodeType;
use crate::model::network::modify::{NodeUpdate, JunctionUpdate, TankUpdate, ReservoirData, JunctionData, TankData};

use crate::ffi::enums::NodeType as ENNodeType;

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_double};

#[unsafe(no_mangle)]
pub extern "C" fn EN_addnode(ph: *mut Project, id: *const c_char, node_type: c_int, out_index: *mut c_int) -> ErrorCode {
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
    ENNodeType::Junction => simulation.network.add_junction(node_id, &JunctionData {
      elevation: 0.0,
      basedemand: 0.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    }),
    ENNodeType::Reservoir => simulation.network.add_reservoir(node_id, &ReservoirData {
      elevation: 0.0,
      head_pattern: None,
      coordinates: None,
    }),
    ENNodeType::Tank => simulation.network.add_tank(node_id, &TankData {
      elevation: 0.0,
      initial_level: 0.0,
      min_level: 0.0,
      max_level: 0.0,
      diameter: 0.0,
      min_volume: 0.0,
      volume_curve_id: None,
      overflow: false,
      coordinates: None,
    })
  };

  if result.is_err() {
    return ErrorCode::IllegalNodeProperty;
  }
  unsafe { *out_index = (*simulation.network.node_map.get(node_id).unwrap()+1) as c_int };
  ErrorCode::Ok
}

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
      NodeType::Junction(junction) => junction.pattern_index.map(|index| (index+1) as f64).unwrap_or(123.0),
      NodeType::Reservoir(reservoir) => reservoir.head_pattern_index.map(|index| (index+1) as f64).unwrap_or(123.0),
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

    NodeProperty::InitVolume => match &node.node_type {
      NodeType::Tank(tank) => tank.volume_at_level(tank.initial_level) * options.unit_system.per_cubic_feet(),
      _ => 0.0,
    },

    NodeProperty::Demand => simulation.solved_state().map_or(0.0, |state| state.demands[index] * flow_units.per_cfs()),
    NodeProperty::Head => simulation.solved_state().map_or(0.0, |state| state.heads[index] * unit_system.per_feet()),
    NodeProperty::Pressure => simulation.solved_state().map_or(0.0, |state| (state.heads[index] - node.elevation) * pressure_units.per_feet()),
    NodeProperty::Quality => 0.0,  // TODO: quality not implemented yet
    NodeProperty::TankDiam => match &node.node_type {
      NodeType::Tank(tank) => tank.diameter * options.unit_system.per_feet(),
      _ => 0.0,
    },
    NodeProperty::MinVolume => match &node.node_type {
      NodeType::Tank(tank) => tank.min_volume() * options.unit_system.per_cubic_feet(),
      _ => 0.0,
    },
    NodeProperty::MaxVolume => match &node.node_type {
      NodeType::Tank(tank) => tank.max_volume() * options.unit_system.per_cubic_feet(),
      _ => 0.0,
    },
    NodeProperty::MinLevel => match &node.node_type {
      NodeType::Tank(tank) => tank.min_level * options.unit_system.per_feet(),
      _ => 0.0,
    },
    NodeProperty::MaxLevel => match &node.node_type {
      NodeType::Tank(tank) => tank.max_level * options.unit_system.per_feet(),
      _ => 0.0,
    },

    NodeProperty::TankVolume => match &node.node_type {
      NodeType::Tank(tank) => {
        if let Some(state) = simulation.solved_state() {
          tank.volume_at_head(state.heads[index]) * options.unit_system.per_cubic_feet()
        } else {
          tank.volume_at_level(tank.initial_level) * options.unit_system.per_cubic_feet()
        }
      }
      _ => 0.0,
    },

    NodeProperty::CanOverflow => match &node.node_type {
      NodeType::Tank(tank) => if tank.overflow { 1.0 } else { 0.0 },
      _ => 0.0,
    },
    NodeProperty::DemandDeficit => 0.0,
    NodeProperty::NodeInControl => 0.0,
    NodeProperty::EmitterFlow => 0.0,
    NodeProperty::LeakageFlow => 0.0,
    NodeProperty::DemandFlow => 0.0,
    NodeProperty::FullDemand => 0.0,
    NodeProperty::SourceMass => 0.0,  // TODO: mass not implemented yet
    _ => -123.0,
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
    NodeProperty::TankLevel => {
      network.update_tank(&node_id, &TankUpdate { initial_level: Some(value as f64), ..Default::default() })
    },
    NodeProperty::Pattern => {
      let pattern_id = match network.patterns.get(value as usize - 1 ) {
        Some(pattern) => pattern.id.clone(),
        None => return ErrorCode::UndefinedPattern,
      };
      network.update_junction(&node_id, &JunctionUpdate { pattern: Some(Some(pattern_id)), ..Default::default() })
    }
    // TODO: implement missing properties
    _ => return ErrorCode::InvalidParameterCode,
  };

  if result.is_err() {
    return ErrorCode::IllegalNodeProperty;
  }

  ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_getcoord(ph: *mut Project, index: c_int, out_x: *mut c_double, out_y: *mut c_double) -> ErrorCode {
  let simulation = get_simulation!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = (index-1) as usize;

  let node = match simulation.network.nodes.get(index) {
    Some(node) => node,
    None => return ErrorCode::UndefinedNode,
  };

  if let Some(coordinates) = node.coordinates {
    unsafe { *out_x = coordinates.0 as c_double };
    unsafe { *out_y = coordinates.1 as c_double };
  } else {
    return ErrorCode::NodeNoCoordinates
  }

  ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_setcoord(ph: *mut Project, index: c_int, x: c_double, y: c_double) -> ErrorCode {
  let simulation = get_simulation_mut!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = (index-1) as usize;

  let node = match simulation.network.nodes.get_mut(index) {
    Some(node) => node,
    None => return ErrorCode::UndefinedNode,
  };
  node.coordinates = Some((x as f64, y as f64));
  ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_deletenode(ph: *mut Project, index: c_int, action_code: c_int) -> ErrorCode {
  let simulation = get_simulation_mut!(ph);

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let index = (index-1) as usize;

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