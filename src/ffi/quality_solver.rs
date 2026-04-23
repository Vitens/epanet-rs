use crate::ffi::error_codes::ErrorCode;
use std::os::raw::{c_int, c_long};

use crate::ffi::project::Project;

#[unsafe(no_mangle)]
pub extern "C" fn EN_openQ(ph: *mut Project) -> ErrorCode {
  ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_initQ(ph: *mut Project, _initflag: c_int) -> ErrorCode {
  ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_runQ(ph: *mut Project, _time: c_long) -> ErrorCode {
  ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_nextQ(ph: *mut Project, _time: c_long) -> ErrorCode {
  ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_stepQ(ph: *mut Project, _time: c_long) -> ErrorCode {
  ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_closeQ(ph: *mut Project) -> ErrorCode {
  ErrorCode::Ok
}