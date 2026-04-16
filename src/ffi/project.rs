use std::ffi::CStr;
use std::os::raw::c_char;

use crate::model::network::Network;
use crate::simulation::Simulation;

use super::ErrorCode;


// Alias EPANET type Project to EPANET-RS style Simulation
pub type Project = Simulation;

/// Creates an EPANET project.
///
/// The caller receives an opaque handle via `*ph` that must be passed to all
/// subsequent API calls and eventually freed with [`EN_deleteproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_createproject(ph: *mut *mut Project) -> ErrorCode {
    if ph.is_null() {
        return ErrorCode::InvalidHandle;
    }
    unsafe { *ph = std::ptr::null_mut() };
    ErrorCode::Ok
}

/// Deletes a currently opened EPANET project and frees all of its memory.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_deleteproject(ph: *mut Project) -> ErrorCode {
    if ph.is_null() {
        return ErrorCode::Ok;
    }
    unsafe { drop(Box::from_raw(ph)) };
    ErrorCode::Ok
}

/// Opens an EPANET input file, reads in the network data, and initializes
/// the hydraulic solver and state so the project is ready for simulation.
///
/// `rptFile` and `outFile` are accepted for API compatibility but currently ignored.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_open(
    ph: *mut *mut Project,
    inp_file: *const c_char,
    _rpt_file: *const c_char,
    _out_file: *const c_char,
) -> ErrorCode {
    if ph.is_null() {
        return ErrorCode::InvalidHandle;
    }
    if inp_file.is_null() {
        return ErrorCode::FileNotFound;
    }

    let c_str = unsafe { CStr::from_ptr(inp_file) };
    let path = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::FileNotFound,
    };

    let mut network = Network::default();
    if let Err(_) = network.read_file(path) {
        return ErrorCode::InputError;
    }

    let simulation = match Simulation::new(network) {
        Ok(s) => s,
        Err(_) => return ErrorCode::InputError,
    };

    unsafe { *ph = Box::into_raw(Box::new(simulation)) };
    ErrorCode::Ok
}
