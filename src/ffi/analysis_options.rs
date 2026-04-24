//! FFI analysis-options accessors: `EN_gettimeparam`, `EN_setoption`, `EN_setdemandmodel`, etc.

use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation_mut};

use crate::ffi::enums::{TimeParameter, DemandModel as ENDemandModel, SimOption};
use crate::model::options::DemandModel;

use std::os::raw::{c_int, c_long};


#[unsafe(no_mangle)]
pub extern "C" fn EN_settimeparam(ph: *mut Project, param: c_int, value: c_long) -> ErrorCode {
  let simulation = get_simulation_mut!(ph);

  let time_options = &mut simulation.network.options.time_options;

  let value = value as usize;

  let param = match TimeParameter::from_repr(param) {
    Some(param) => param,
    None => return ErrorCode::InvalidParameterCode,
  };

  match param {
    TimeParameter::Duration => {
      time_options.duration = value;
      if time_options.duration > time_options.pattern_start {
        time_options.pattern_start = time_options.duration;
      }
    },
    TimeParameter::HydStep => {
      if value == 0 { return ErrorCode::InvalidOptionValue; }
      time_options.hydraulic_timestep = value;
      time_options.hydraulic_timestep = time_options.hydraulic_timestep.min(time_options.pattern_timestep);
      time_options.hydraulic_timestep = time_options.hydraulic_timestep.min(time_options.report_timestep);
    },
    TimeParameter::PatternStep => {
      if value == 0 { return ErrorCode::InvalidOptionValue; }
      time_options.pattern_timestep = value;
      if time_options.hydraulic_timestep > time_options.pattern_timestep {
        time_options.hydraulic_timestep = time_options.pattern_timestep;
      }
    },
    TimeParameter::PatternStart => {
      println!("Pattern start: {}", value);
      time_options.pattern_start = value;
    }
    TimeParameter::ReportStep => {
      if value == 0 { return ErrorCode::InvalidOptionValue; }
      time_options.report_timestep = value;
      if time_options.hydraulic_timestep > time_options.report_timestep {
        time_options.hydraulic_timestep = time_options.report_timestep;
      }
    },
    _ => return ErrorCode::InvalidParameterCode,
  }
  
  ErrorCode::Ok

}

#[unsafe(no_mangle)]
pub extern "C" fn EN_setdemandmodel(ph: *mut Project, demand_model: c_int) -> ErrorCode {
  let simulation = get_simulation_mut!(ph);

  let demand_model = match ENDemandModel::from_repr(demand_model) {
    Some(demand_model) => demand_model,
    None => return ErrorCode::InvalidParameterCode,
  };

  match demand_model {
    ENDemandModel::Pda => { simulation.network.options.demand_model = DemandModel::PDA; }
    ENDemandModel::Dda => { simulation.network.options.demand_model = DemandModel::DDA; }
  }

  ErrorCode::Ok

}

#[unsafe(no_mangle)]
pub extern "C" fn EN_setoption(ph: *mut Project, option: c_int, _value: c_long) -> ErrorCode {
  let _simulation = get_simulation_mut!(ph);

  let _option = match SimOption::from_repr(option) {
    Some(option) => option,
    None => return ErrorCode::InvalidParameterCode,
  };
  // TODO: implement setoption
  ErrorCode::Ok
}