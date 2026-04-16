use std::ffi::CStr;
use std::os::raw::{c_char, c_int};

use crate::error::InputError;
use crate::model::network::Network;
use crate::simulation::Simulation;

use super::ErrorCode;

/// Opaque EPANET project handle wrapping an optional [`Simulation`].
///
/// - After [`EN_createproject`] the handle exists but is empty (`None`).
/// - After [`EN_open`] it contains a ready-to-run `Simulation`.
pub struct Project{
  simulation: Option<Simulation>,
}

/// Creates an EPANET project.
///
/// Allocates a new project handle via `*ph` that must be passed to all
/// subsequent API calls and eventually freed with [`EN_deleteproject`].
#[unsafe(no_mangle)]
pub extern "C" fn EN_createproject(ph: *mut *mut Project) -> ErrorCode {
    unsafe { *ph = Box::into_raw(Box::new(Project { simulation: None })) };
    ErrorCode::Ok
}

/// Deletes a project and frees all of its memory.
#[unsafe(no_mangle)]
pub extern "C" fn EN_deleteproject(ph: *mut Project) -> ErrorCode {
    unsafe { drop(Box::from_raw(ph)) };
    ErrorCode::Ok
}

/// Opens an EPANET input file, reads in the network data, and initializes
/// the hydraulic solver and state so the project is ready for simulation.
///
/// `rptFile` and `outFile` are accepted for API compatibility but currently ignored.
#[unsafe(no_mangle)]
pub extern "C" fn EN_open(
    ph: *mut Project,
    inp_file: *const c_char,
    _rpt_file: *const c_char,
    _out_file: *const c_char,
) -> ErrorCode {

    if ph.is_null() {
      return ErrorCode::InvalidHandle;
    }

    if inp_file.is_null() {
        return ErrorCode::CannotOpenInputFile;
    }

    let c_str = unsafe { CStr::from_ptr(inp_file) };

    let path = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::CannotOpenInputFile,
    };

    let mut network = Network::default();
    if let Err(e) = network.read_file(path) {
        return match e {
            InputError::FileOpen { .. } | InputError::FileRead(_) => ErrorCode::CannotOpenInputFile,
            InputError::Parse { .. } => ErrorCode::InputError,
        };
    }

    let simulation = match Simulation::new(network) {
        Ok(s) => s,
        Err(_) => return ErrorCode::CannotSolveHydraulics,
    };

    let project = unsafe { &mut *ph };
    project.simulation = Some(simulation);
    ErrorCode::Ok
}

/// Retrieves the text of an error message given its error code.
///
/// Writes up to `maxLen` bytes (including the null terminator) into `errmsg`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_geterror(errcode: c_int, errmsg: *mut c_char, max_len: c_int) -> ErrorCode {
    if errmsg.is_null() || max_len <= 0 {
        return ErrorCode::InvalidFormat;
    }

    let msg = match ErrorCode::from_repr(errcode) {
        Some(code) => code.to_string(),
        None => format!("Unknown error code: {}", errcode),
    };

    let buf = unsafe { std::slice::from_raw_parts_mut(errmsg as *mut u8, max_len as usize) };
    let bytes = msg.as_bytes();
    let copy_len = bytes.len().min(buf.len() - 1);
    buf[..copy_len].copy_from_slice(&bytes[..copy_len]);
    buf[copy_len] = 0;

    ErrorCode::Ok
}
