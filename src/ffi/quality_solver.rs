//! FFI water-quality-solver stubs: `EN_openQ`, `EN_initQ`, `EN_runQ`, `EN_nextQ`, `EN_closeQ`.

use crate::ffi::error_codes::ErrorCode;
use std::os::raw::{c_int, c_long};

use crate::ffi::project::Project;

#[unsafe(no_mangle)]
pub extern "C" fn EN_openQ(_ph: *mut Project) -> ErrorCode {
    ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_initQ(_ph: *mut Project, _initflag: c_int) -> ErrorCode {
    ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_runQ(_ph: *mut Project, _time: c_long) -> ErrorCode {
    ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_nextQ(_ph: *mut Project, _time: c_long) -> ErrorCode {
    ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_stepQ(_ph: *mut Project, _time: c_long) -> ErrorCode {
    ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_closeQ(_ph: *mut Project) -> ErrorCode {
    ErrorCode::Ok
}
