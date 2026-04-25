//! FFI project lifecycle: `EN_createproject`, `EN_deleteproject`, `EN_open`, `EN_close`, etc.

use std::ffi::CStr;
use std::os::raw::{c_char, c_int};

use crate::error::InputError;

use crate::ffi::enums::{FlowUnits as ENFlowUnits, HeadLossType as ENHeadLossType};
use crate::ffi::error_codes::ErrorCode;
use crate::simulation::Simulation;

use crate::model::options::HeadlossFormula;
use crate::model::units::FlowUnits;

/// Opaque EPANET project handle wrapping an optional [`Simulation`].
///
/// - After [`EN_createproject`] the handle exists but is empty (`None`).
/// - After [`EN_open`] it contains a ready-to-run `Simulation`.
pub struct Project {
    pub(crate) simulation: Option<Simulation>,
}

/// Macro to get the simulation from a project handle.
macro_rules! get_simulation {
    ($ph:expr) => {{
        if $ph.is_null() {
            return ErrorCode::InvalidHandle;
        }
        let project = unsafe { &*$ph };
        match project.simulation.as_ref() {
            Some(s) => s,
            None => return ErrorCode::NoNetworkData,
        }
    }};
}

/// Macro to get mutable simulation from a project handle.
macro_rules! get_simulation_mut {
    ($ph:expr) => {{
        if $ph.is_null() {
            return ErrorCode::InvalidHandle;
        }
        let project = unsafe { &mut *$ph };
        match project.simulation.as_mut() {
            Some(s) => s,
            None => return ErrorCode::NoNetworkData,
        }
    }};
}

pub(crate) use get_simulation;
#[allow(unused_imports)]
pub(crate) use get_simulation_mut;

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
/// Initializes a project with a new network.
#[unsafe(no_mangle)]
pub extern "C" fn EN_init(
    ph: *mut Project,
    _rpt_file: *const c_char,
    _out_file: *const c_char,
    flow_units: c_int,
    headloss_type: c_int,
) -> ErrorCode {
    if ph.is_null() {
        return ErrorCode::InvalidHandle;
    }

    let flow_units = match ENFlowUnits::from_repr(flow_units) {
        Some(flow_units) => flow_units,
        None => return ErrorCode::InvalidParameterCode,
    };
    let headloss_type = match ENHeadLossType::from_repr(headloss_type) {
        Some(headloss_formula) => headloss_formula,
        None => return ErrorCode::InvalidParameterCode,
    };

    // map flow_units to FlowUnits (rust enum)
    let flow_units = match flow_units {
        ENFlowUnits::Cfs => FlowUnits::CFS,
        ENFlowUnits::Gpm => FlowUnits::GPM,
        ENFlowUnits::Mgd => FlowUnits::MGD,
        ENFlowUnits::Imgd => FlowUnits::IMGD,
        ENFlowUnits::Afd => FlowUnits::AFD,
        ENFlowUnits::Lps => FlowUnits::LPS,
        ENFlowUnits::Lpm => FlowUnits::LPM,
        ENFlowUnits::Mld => FlowUnits::MLD,
        ENFlowUnits::Cmh => FlowUnits::CMH,
        ENFlowUnits::Cmd => FlowUnits::CMD,
        ENFlowUnits::Cms => FlowUnits::CMS,
    };

    let headloss_formula = match headloss_type {
        ENHeadLossType::HW => HeadlossFormula::HazenWilliams,
        ENHeadLossType::DW => HeadlossFormula::DarcyWeisbach,
        ENHeadLossType::CM => HeadlossFormula::ChezyManning,
    };

    let simulation = Simulation::init(flow_units, headloss_formula);

    let project = unsafe { &mut *ph };
    project.simulation = Some(simulation);

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

    let simulation = match Simulation::from_file(path) {
        Ok(s) => s,
        Err(e) => {
            return match e {
                InputError::FileOpen { .. } | InputError::FileRead(_) => {
                    ErrorCode::CannotOpenInputFile
                }
                InputError::Parse { .. } => ErrorCode::InputError,
                InputError::NodeExists { .. } => ErrorCode::DuplicateId,
                InputError::LinkExists { .. } => ErrorCode::DuplicateId,
                InputError::PatternNotFound { .. } => ErrorCode::UndefinedPattern,
                _ => ErrorCode::InputError,
            };
        }
    };

    let project = unsafe { &mut *ph };
    project.simulation = Some(simulation);
    ErrorCode::Ok
}

#[unsafe(no_mangle)]
pub extern "C" fn EN_saveinpfile(ph: *mut Project, inp_file: *const c_char) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let c_str = unsafe { CStr::from_ptr(inp_file) };
    let path = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return ErrorCode::CannotOpenInputFile,
    };

    let result = simulation.network.save_network(path);
    if result.is_err() {
        return ErrorCode::CannotOpenInputFile;
    }

    ErrorCode::Ok
}
