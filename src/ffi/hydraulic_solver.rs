//! FFI hydraulic-solver workflow: `EN_openH`, `EN_initH`, `EN_runH`, `EN_nextH`, `EN_solveH`, `EN_closeH`.

use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation_mut};
use crate::solver::state::SolverState;

use std::os::raw::{c_int, c_long};

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
pub unsafe extern "C" fn EN_initH(ph: *mut Project, _initflag: c_int) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);
    if simulation.solver.is_none() {
        return ErrorCode::HydraulicSolverNotOpened;
    }
    // reset the simulation to initial conditions
    simulation.state = Some(SolverState::new_with_initial_values(&simulation.network));
    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_runH(ph: *mut Project, _time: c_int) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);
    if simulation.solver.is_none() {
        return ErrorCode::HydraulicSolverNotOpened;
    }
    let result = simulation.run_hydraulics();

    match result {
        Ok(_) => ErrorCode::Ok,
        Err(_) => ErrorCode::CannotSolveHydraulics,
    }
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_time` must be a valid non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_nextH(ph: *mut Project, out_time: *mut c_long) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);
    if simulation.solver.is_none() {
        return ErrorCode::HydraulicSolverNotOpened;
    }
    let result = simulation.next_hydraulic_timestep();
    if !out_time.is_null() {
        unsafe {
            *out_time = result as c_long;
        }
    }
    ErrorCode::Ok
}
