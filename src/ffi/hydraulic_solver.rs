use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation_mut};
use crate::solver::state::SolverState;

use std::os::raw::{c_int};

#[unsafe(no_mangle)]
pub extern "C" fn EN_openH(ph: *mut Project) -> ErrorCode {
  let simulation = get_simulation_mut!(ph);
  let result = simulation.initialize_hydraulics();
  match result {
    Ok(_) => ErrorCode::Ok,
    Err(_) => ErrorCode::CannotSolveHydraulics,
  }
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_closeH(ph: *mut Project) -> ErrorCode {
  let simulation = get_simulation_mut!(ph);
  simulation.solver = None;
  ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_initH(ph: *mut Project, _initflag: c_int) -> ErrorCode {
  let simulation = get_simulation_mut!(ph);
  if simulation.solver.is_none() {
    return ErrorCode::HydraulicSolverNotOpened;
  }
  // reset the simulation to initial conditions
  simulation.state = Some(SolverState::new_with_initial_values(&simulation.network));
  ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_runH(ph: *mut Project, _time: c_int) -> ErrorCode {
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