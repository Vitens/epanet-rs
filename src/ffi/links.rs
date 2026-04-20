use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::model::link::LinkType;
use crate::model::valve::ValveType;

use crate::ffi::enums::LinkType as ENLinkType;
use crate::ffi::enums::LinkProperty;
use crate::ffi::enums::MISSING_VALUE;
use crate::constants::MperFT;

use crate::model::units::UnitSystem;
use crate::model::options::HeadlossFormula;
use crate::model::link::LinkStatus;


use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_double};

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

  let link_id = match simulation.network.links.get((index - 1) as usize) {
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

  // EPANET indexes from 1, so we need to subtract 1 from the index
  let link = match simulation.network.links.get_mut((index - 1) as usize) {
    Some(link) => link,
    None => return ErrorCode::UndefinedLink,
  };

  let c_str = unsafe { CStr::from_ptr(id) };

  let new_link_id = match c_str.to_str() {
    Ok(s) => s,
    Err(_) => return ErrorCode::InvalidIdName,
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

  let link = match simulation.network.links.get((index - 1) as usize) {
    Some(link) => link,
    None => return ErrorCode::UndefinedLink,
  };

  let link_type_int = match &link.link_type {
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

  let link = match simulation.network.links.get((index - 1) as usize) {
    Some(link) => link,
    None => return ErrorCode::UndefinedLink,
  };

  unsafe { *out_start_node = (link.start_node+1) as c_int };
  unsafe { *out_end_node = (link.end_node+1) as c_int };

  ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_getlinkvalue(ph: *mut Project, index: c_int, property: c_int, out_value: *mut c_double) -> ErrorCode {
  let simulation = get_simulation!(ph);

  let index = (index - 1) as usize;

  let link = match simulation.network.links.get(index) {
    Some(link) => link,
    None => return ErrorCode::UndefinedLink,
  };

  let options = &simulation.network.options;

  let property = match LinkProperty::from_repr(property) {
    Some(property) => property,
    None => return ErrorCode::InvalidParameterCode,
  };

  let value = match property {
    LinkProperty::Diameter => {
      // convert from standard units (Ft) to in (US) or mm (SI)
      let conversion = match options.unit_system {
        UnitSystem::US => 12.0, // ft to in
        UnitSystem::SI => MperFT * 1e3, // ft to mm
      };

      match &link.link_type {
        LinkType::Pipe(pipe) => pipe.diameter * conversion,
        LinkType::Valve(valve) => valve.diameter * conversion,
        LinkType::Pump(_) => 0.0,
      }
    },
    LinkProperty::Length => {
      match &link.link_type {
        LinkType::Pipe(pipe) => pipe.length * options.unit_system.per_feet(),
        _ => 0.0,
      }
    },
    LinkProperty::Roughness => {
      match &link.link_type {
        LinkType::Pipe(pipe) => {
          if pipe.headloss_formula == HeadlossFormula::DarcyWeisbach {
            pipe.roughness * options.unit_system.per_feet()
          } else {
            pipe.roughness
          }
        }
        _ => 0.0,
      }
    },
    LinkProperty::MinorLoss => {
      match &link.link_type {
        LinkType::Pipe(pipe) => pipe.minor_loss * pipe.diameter.powi(4) / 0.02517,
        LinkType::Valve(valve) => valve.minor_loss * valve.diameter.powi(4) / 0.02517,
        _ => 0.0,
      }
    },
    LinkProperty::InitStatus => {
      match &link.initial_status {
        LinkStatus::Closed => 0.0,
        _ => 1.0,
      }
    }
    // TODO: Implement init setting
    LinkProperty::InitSetting => 0.0,
    // TODO: Implement KBulk and KWall
    LinkProperty::KBulk => 0.0,
    LinkProperty::KWall => 0.0,

    // Flow, Velocity, Headloss (dynamic properties) are only available if the state is available and solved, otherwise return MISSING_VALUE
    LinkProperty::Flow => simulation.solved_state().map_or(MISSING_VALUE, |state| state.flows[index] * options.flow_units.per_cfs()),
    LinkProperty::Velocity => simulation.solved_state().map_or(MISSING_VALUE, |state| {
      let flow = state.flows[index].abs();
      match &link.link_type {
        LinkType::Pipe(pipe) => { flow / (pipe.diameter.powi(2) * std::f64::consts::PI / 4.0) * options.unit_system.per_feet() },
        LinkType::Valve(valve) => { flow / (valve.diameter.powi(2) * std::f64::consts::PI / 4.0) * options.unit_system.per_feet() },
        LinkType::Pump(_) => 0.0,
      }
    }),

    // TODO: expose other link properties
    _ => 0.0,
  };

  unsafe { *out_value = value };

  ErrorCode::Ok
}