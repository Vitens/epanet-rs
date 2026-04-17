use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation_mut};

use crate::ffi::enums::TimeParameter;

use std::os::raw::{c_int, c_double};


#[unsafe(no_mangle)]
pub extern "C" fn EN_settimeparam(ph: *mut Project, param: c_int, value: c_double) -> ErrorCode {
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