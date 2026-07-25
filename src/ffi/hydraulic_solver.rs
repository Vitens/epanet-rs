//! FFI hydraulic-solver workflow: `EN_openH`, `EN_initH`, `EN_runH`, `EN_nextH`, `EN_solveH`, `EN_closeH`.

use crate::ffi::enums::InitHydOption;
use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation_mut};
use crate::ffi::util::write_out;
use crate::solver::state::SolverState;

use std::os::raw::{c_char, c_int, c_long};

/// Runs a complete hydraulic simulation for every time step.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_solveH(ph: *mut Project) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);
    match simulation.solve_hydraulics(false) {
        Ok(_) => ErrorCode::Ok,
        Err(_) => ErrorCode::CannotSolveHydraulics,
    }
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_openH(ph: *mut Project) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);
    let result = simulation.initialize_hydraulics();
    match result {
        Ok(_) => ErrorCode::Ok,
        Err(_) => ErrorCode::CannotSolveHydraulics,
    }
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_closeH(ph: *mut Project) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);
    simulation.solver = None;
    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_initH(ph: *mut Project, init_flag: c_int) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);
    if simulation.solver.is_none() {
        return ErrorCode::HydraulicSolverNotOpened;
    }
    if InitHydOption::from_repr(init_flag).is_none() {
        return ErrorCode::InvalidParameterCode;
    }
    // Saving results to a binary hydraulics file is not supported, so both the
    // save and the re-initialize flags behave the same: reset the simulation to
    // its initial conditions.
    simulation.state = Some(SolverState::new_with_initial_values(&simulation.network));
    simulation.time = 0;
    simulation.solved = false;
    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_current_time` must be null or a valid writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_runH(ph: *mut Project, out_current_time: *mut c_long) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);
    if simulation.solver.is_none() {
        return ErrorCode::HydraulicSolverNotOpened;
    }
    let result = simulation.run_hydraulics();

    match result {
        Ok(time) => {
            unsafe { write_out(out_current_time, time as c_long) };
            ErrorCode::Ok
        }
        Err(_) => ErrorCode::CannotSolveHydraulics,
    }
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_time` must be null or a valid writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_nextH(ph: *mut Project, out_time: *mut c_long) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);
    if simulation.solver.is_none() {
        return ErrorCode::HydraulicSolverNotOpened;
    }
    let result = simulation.next_hydraulic_timestep();
    unsafe { write_out(out_time, result as c_long) };
    ErrorCode::Ok
}

/// Saves the hydraulic results of the current time step to the binary output file.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_saveH(_ph: *mut Project) -> ErrorCode {
    // TODO: epanet-rs does not write the EPANET binary output file.
    ErrorCode::NotImplemented
}

/// Saves the current contents of the binary hydraulics file to a named file.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_savehydfile(_ph: *mut Project, _filename: *const c_char) -> ErrorCode {
    // TODO: epanet-rs does not write the EPANET binary hydraulics file.
    ErrorCode::NotImplemented
}

/// Uses a previously saved binary hydraulics file to supply a project's hydraulics.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_usehydfile(_ph: *mut Project, _filename: *const c_char) -> ErrorCode {
    // TODO: epanet-rs does not read the EPANET binary hydraulics file.
    ErrorCode::NotImplemented
}
