//! FFI project lifecycle: `EN_createproject`, `EN_deleteproject`, `EN_open`, `EN_close`, etc.

use std::ffi::CStr;
use std::os::raw::{c_char, c_int, c_void};

use crate::error::InputError;

use crate::ffi::enums::{FlowUnits as ENFlowUnits, HeadLossType as ENHeadLossType};
use crate::ffi::error_codes::ErrorCode;
use crate::ffi::util::{MAX_MSG, read_str, write_out, write_str};
use crate::simulation::Simulation;

use crate::model::options::HeadlossFormula;
use crate::model::units::FlowUnits;

/// Toolkit version reported by [`EN_getversion`]: EPANET 2.3.0.
const TOOLKIT_VERSION: c_int = 20300;

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
///
/// # Safety
///
/// `ph` must be a valid, non-null pointer to writable memory for one
/// `*mut Project`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_createproject(ph: *mut *mut Project) -> ErrorCode {
    if ph.is_null() {
        return ErrorCode::InvalidHandle;
    }

    unsafe { *ph = Box::into_raw(Box::new(Project { simulation: None })) };
    ErrorCode::Ok
}

/// Deletes a project and frees all of its memory.
///
/// # Safety
///
/// `ph` must be a valid project handle returned by [`EN_createproject`], and it
/// must not be used again after this call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_deleteproject(ph: *mut Project) -> ErrorCode {
    if ph.is_null() {
        return ErrorCode::InvalidHandle;
    }

    unsafe { drop(Box::from_raw(ph)) };
    ErrorCode::Ok
}
/// Initializes a project with a new network.
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle created by
/// [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_init(
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
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle. `inp_file` must be a
/// valid, non-null pointer to a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_open(
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
        Err(e) => return input_error_code(&e),
    };

    let project = unsafe { &mut *ph };
    project.simulation = Some(simulation);
    ErrorCode::Ok
}

/// Maps a parser error onto the closest EPANET error code.
pub(crate) fn input_error_code(error: &InputError) -> ErrorCode {
    match error {
        InputError::FileOpen { .. } | InputError::FileRead(_) => ErrorCode::CannotOpenInputFile,
        InputError::Parse { .. } => ErrorCode::InputError,
        InputError::NodeExists { .. } | InputError::LinkExists { .. } => ErrorCode::DuplicateId,
        InputError::NodeNotFound { .. } => ErrorCode::UndefinedNode,
        InputError::LinkNotFound { .. } => ErrorCode::UndefinedLink,
        InputError::PatternNotFound { .. } => ErrorCode::UndefinedPattern,
        InputError::CurveNotFound { .. } => ErrorCode::UndefinedCurve,
        InputError::NodeNotATank { .. } => ErrorCode::NotATank,
        InputError::LinkNotAValve { .. } => ErrorCode::NotAValve,
        InputError::LinkNotAPump { .. } => ErrorCode::UndefinedPump,
        InputError::TankLevelsInvalid { .. } => ErrorCode::InvalidTankLevels,
        InputError::NodeNotAJunction { .. } | InputError::NodeNotAReservoir { .. } => {
            ErrorCode::IllegalNodeProperty
        }
        InputError::LinkNotAPipe { .. } => ErrorCode::IllegalLinkProperty,
    }
}

/// Opens an EPANET input file, reporting the first error encountered instead of
/// aborting on it.
///
/// `epanet-rs` stops at the first parse error, so this behaves exactly like
/// [`EN_open`].
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle. `inp_file` must be a
/// valid, non-null pointer to a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_openX(
    ph: *mut Project,
    inp_file: *const c_char,
    rpt_file: *const c_char,
    out_file: *const c_char,
) -> ErrorCode {
    unsafe { EN_open(ph, inp_file, rpt_file, out_file) }
}

/// Closes a project's network data, leaving the handle itself allocated.
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle created by
/// [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_close(ph: *mut Project) -> ErrorCode {
    if ph.is_null() {
        return ErrorCode::InvalidHandle;
    }
    let project = unsafe { &mut *ph };
    project.simulation = None;
    ErrorCode::Ok
}

/// Runs a complete EPANET simulation from an input file.
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_runproject(
    _ph: *mut Project,
    _inp_file: *const c_char,
    _rpt_file: *const c_char,
    _out_file: *const c_char,
    _callback: Option<extern "C" fn(*mut c_char)>,
) -> ErrorCode {
    // TODO: requires report and binary output file writing, which epanet-rs does not support.
    ErrorCode::NotImplemented
}

/// Retrieves the title lines of a project.
///
/// `epanet-rs` stores the title as a single string; embedded newlines are
/// split over the three output lines.
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle. Each output pointer must be
/// null or point to a buffer of at least `EN_MAXMSG + 1` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_gettitle(
    ph: *mut Project,
    out_line1: *mut c_char,
    out_line2: *mut c_char,
    out_line3: *mut c_char,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let title = simulation.network.title.as_deref().unwrap_or("");
    let mut lines = title.lines();

    for out in [out_line1, out_line2, out_line3] {
        unsafe { write_str(out, lines.next().unwrap_or(""), MAX_MSG) };
    }

    ErrorCode::Ok
}

/// Sets the title lines of a project.
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle. Each line must be null or a
/// valid pointer to a NUL-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_settitle(
    ph: *mut Project,
    line1: *const c_char,
    line2: *const c_char,
    line3: *const c_char,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let mut lines: Vec<&str> = Vec::with_capacity(3);
    for line in [line1, line2, line3] {
        match unsafe { read_str(line) } {
            Some(text) if !text.is_empty() => lines.push(text),
            Some(_) => {}
            None if line.is_null() => {}
            None => return ErrorCode::InvalidFormat,
        }
    }

    simulation.network.title = if lines.is_empty() {
        None
    } else {
        Some(lines.join("\n").into())
    };

    ErrorCode::Ok
}

/// Retrieves the descriptive comment assigned to an object.
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getcomment(
    _ph: *mut Project,
    _object: c_int,
    _index: c_int,
    _out_comment: *mut c_char,
) -> ErrorCode {
    // TODO: the network model does not store per-object comments.
    ErrorCode::NotImplemented
}

/// Assigns a descriptive comment to an object.
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setcomment(
    _ph: *mut Project,
    _object: c_int,
    _index: c_int,
    _comment: *const c_char,
) -> ErrorCode {
    // TODO: the network model does not store per-object comments.
    ErrorCode::NotImplemented
}

/// Retrieves the tag assigned to an object.
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_gettag(
    _ph: *mut Project,
    _object: c_int,
    _index: c_int,
    _out_tag: *mut c_char,
) -> ErrorCode {
    // TODO: the network model does not store per-object tags.
    ErrorCode::NotImplemented
}

/// Assigns a tag to an object.
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_settag(
    _ph: *mut Project,
    _object: c_int,
    _index: c_int,
    _tag: *const c_char,
) -> ErrorCode {
    // TODO: the network model does not store per-object tags.
    ErrorCode::NotImplemented
}

/// Retrieves the toolkit version number.
///
/// # Safety
///
/// `out_version` must be a valid, non-null writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getversion(out_version: *mut c_int) -> ErrorCode {
    if out_version.is_null() {
        return ErrorCode::InvalidFormat;
    }
    unsafe { write_out(out_version, TOOLKIT_VERSION) };
    ErrorCode::Ok
}

/// Registers a callback that receives the project's report lines.
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setreportcallback(
    _ph: *mut Project,
    _callback: Option<extern "C" fn(*mut c_void, *mut c_void, *const c_char)>,
) -> ErrorCode {
    // TODO: epanet-rs does not generate report output.
    ErrorCode::NotImplemented
}

/// Supplies the user data passed to the report callback.
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setreportcallbackuserdata(
    _ph: *mut Project,
    _user_data: *mut c_void,
) -> ErrorCode {
    // TODO: epanet-rs does not generate report output.
    ErrorCode::NotImplemented
}

#[unsafe(no_mangle)]
///
/// # Safety
///
/// `ph` must be a valid, non-null project handle. `inp_file` must be a
/// valid, non-null pointer to a NUL-terminated C string.
pub unsafe extern "C" fn EN_saveinpfile(ph: *mut Project, inp_file: *const c_char) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    if inp_file.is_null() {
        return ErrorCode::CannotOpenInputFile;
    }

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
