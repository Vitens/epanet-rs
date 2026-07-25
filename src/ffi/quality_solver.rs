//! FFI water-quality solver and quality-option accessors.
//!
//! `epanet-rs` is a hydraulic-only solver: every function that would run or
//! configure a water quality analysis reports
//! [`ErrorCode::NotImplemented`]. The two read-only accessors below report
//! that no quality analysis is configured, which is always accurate.

use crate::ffi::enums::QualityType;
use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation};
use crate::ffi::util::{MAX_ID, write_out, write_str};

use std::os::raw::{c_char, c_int, c_long};

/// Runs a complete water quality simulation.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_solveQ(_ph: *mut Project) -> ErrorCode {
    // TODO: water quality analysis is not implemented.
    ErrorCode::NotImplemented
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_openQ(_ph: *mut Project) -> ErrorCode {
    // TODO: water quality analysis is not implemented.
    ErrorCode::NotImplemented
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_initQ(_ph: *mut Project, _save_flag: c_int) -> ErrorCode {
    // TODO: water quality analysis is not implemented.
    ErrorCode::NotImplemented
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_runQ(_ph: *mut Project, _out_current_time: *mut c_long) -> ErrorCode {
    // TODO: water quality analysis is not implemented.
    ErrorCode::NotImplemented
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_nextQ(_ph: *mut Project, _out_time_step: *mut c_long) -> ErrorCode {
    // TODO: water quality analysis is not implemented.
    ErrorCode::NotImplemented
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_stepQ(_ph: *mut Project, _out_time_left: *mut c_long) -> ErrorCode {
    // TODO: water quality analysis is not implemented.
    ErrorCode::NotImplemented
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_closeQ(_ph: *mut Project) -> ErrorCode {
    // TODO: water quality analysis is not implemented.
    ErrorCode::NotImplemented
}

/// Retrieves the type of water quality analysis to be run.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// The output pointers must be null or valid writable pointers.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getqualtype(
    ph: *mut Project,
    out_qual_type: *mut c_int,
    out_trace_node: *mut c_int,
) -> ErrorCode {
    let _ = get_simulation!(ph);

    unsafe { write_out(out_qual_type, QualityType::None as c_int) };
    unsafe { write_out(out_trace_node, 0) };

    ErrorCode::Ok
}

/// Retrieves information about the type of water quality analysis to be run.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// The string buffers must be null or point to at least `EN_MAXID + 1` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getqualinfo(
    ph: *mut Project,
    out_qual_type: *mut c_int,
    out_chem_name: *mut c_char,
    out_chem_units: *mut c_char,
    out_trace_node: *mut c_int,
) -> ErrorCode {
    let _ = get_simulation!(ph);

    unsafe { write_out(out_qual_type, QualityType::None as c_int) };
    unsafe { write_str(out_chem_name, "", MAX_ID) };
    unsafe { write_str(out_chem_units, "", MAX_ID) };
    unsafe { write_out(out_trace_node, 0) };

    ErrorCode::Ok
}

/// Sets the type of water quality analysis to run.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setqualtype(
    _ph: *mut Project,
    _qual_type: c_int,
    _chem_name: *const c_char,
    _chem_units: *const c_char,
    _trace_node: *const c_char,
) -> ErrorCode {
    // TODO: water quality analysis is not implemented.
    ErrorCode::NotImplemented
}
