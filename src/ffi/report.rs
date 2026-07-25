//! FFI reporting, counting and result-index functions.
//!
//! `epanet-rs` does not produce EPANET's textual report file, so the functions
//! that write, format or reset that report are marked as not implemented.

use crate::ffi::enums::{CountType, ObjectType};
use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation};
use crate::ffi::util::write_out;

use crate::model::node::NodeType;

use std::os::raw::{c_char, c_double, c_int, c_long};

/// Retrieves the number of objects of a given type in a project.
///
/// # Safety
///
/// `out_count` must be a valid, non-null pointer to writable memory for one `c_int`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getcount(
    ph: *mut Project,
    object: c_int,
    out_count: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let count_type = match CountType::from_repr(object) {
        Some(ct) => ct,
        None => return ErrorCode::InvalidParameterCode,
    };

    let net = &simulation.network;
    let count = match count_type {
        CountType::NodeCount => net.nodes.len(),
        CountType::TankCount => net
            .nodes
            .iter()
            .filter(|n| matches!(n.node_type, NodeType::Tank(_) | NodeType::Reservoir(_)))
            .count(),
        CountType::LinkCount => net.links.len(),
        CountType::PatCount => net.patterns.len(),
        CountType::CurveCount => net.curves.len(),
        CountType::ControlCount => net.controls.len(),
        CountType::RuleCount => net.rules.len(),
    };

    unsafe { write_out(out_count, count as c_int) };
    ErrorCode::Ok
}

/// Retrieves the text of an error message given its error code.
///
/// Writes up to `max_len` bytes (including the null terminator) into `errmsg`.
///
/// # Safety
///
/// `errmsg` must be a valid pointer to a buffer of at least `max_len` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_geterror(
    errcode: c_int,
    errmsg: *mut c_char,
    max_len: c_int,
) -> ErrorCode {
    if errmsg.is_null() || max_len <= 0 {
        return ErrorCode::InvalidFormat;
    }

    let msg = match ErrorCode::from_repr(errcode) {
        Some(code) => code.to_string(),
        None => format!("Unknown error code: {}", errcode),
    };

    unsafe { crate::ffi::util::write_str(errmsg, &msg, max_len as usize - 1) };

    ErrorCode::Ok
}

/// Retrieves the order in which a node or link appears in the binary output file.
///
/// `epanet-rs` keeps results for every node and link, so an object's result
/// index is identical to its object index.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_value` must be null or a valid writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getresultindex(
    ph: *mut Project,
    object_type: c_int,
    index: c_int,
    out_value: *mut c_int,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let count = match ObjectType::from_repr(object_type) {
        Some(ObjectType::Node) => simulation.network.nodes.len(),
        Some(ObjectType::Link) => simulation.network.links.len(),
        _ => return ErrorCode::InvalidParameterCode,
    };

    if index < 1 || index as usize > count {
        return match ObjectType::from_repr(object_type) {
            Some(ObjectType::Node) => ErrorCode::UndefinedNode,
            _ => ErrorCode::UndefinedLink,
        };
    }

    unsafe { write_out(out_value, index) };
    ErrorCode::Ok
}

/// Retrieves a particular simulation statistic.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getstatistic(
    _ph: *mut Project,
    _statistic_type: c_int,
    _out_value: *mut c_double,
) -> ErrorCode {
    // TODO: the solver discards its per-step convergence statistics.
    ErrorCode::NotImplemented
}

/// Retrieves the time until the next hydraulic event occurs.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_timetonextevent(
    _ph: *mut Project,
    _out_event_type: *mut c_int,
    _out_duration: *mut c_long,
    _out_element_index: *mut c_int,
) -> ErrorCode {
    // TODO: the simulation driver does not expose which event ends a time step.
    ErrorCode::NotImplemented
}

/// Writes a line of text to a project's report file.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_writeline(_ph: *mut Project, _line: *const c_char) -> ErrorCode {
    // TODO: epanet-rs does not write a report file.
    ErrorCode::NotImplemented
}

/// Writes simulation results in a tabular format to a project's report file.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_report(_ph: *mut Project) -> ErrorCode {
    // TODO: epanet-rs does not write a report file.
    ErrorCode::NotImplemented
}

/// Copies the current contents of a project's report file to another file.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_copyreport(_ph: *mut Project, _filename: *const c_char) -> ErrorCode {
    // TODO: epanet-rs does not write a report file.
    ErrorCode::NotImplemented
}

/// Clears the contents of a project's report file.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_clearreport(_ph: *mut Project) -> ErrorCode {
    // TODO: epanet-rs does not write a report file.
    ErrorCode::NotImplemented
}

/// Resets a project's report options to their default values.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_resetreport(_ph: *mut Project) -> ErrorCode {
    // TODO: epanet-rs does not write a report file.
    ErrorCode::NotImplemented
}

/// Processes a reporting format command.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setreport(_ph: *mut Project, _format: *const c_char) -> ErrorCode {
    // TODO: epanet-rs does not write a report file.
    ErrorCode::NotImplemented
}

/// Sets the level of hydraulic status reporting.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setstatusreport(_ph: *mut Project, _level: c_int) -> ErrorCode {
    // TODO: epanet-rs does not write a report file.
    ErrorCode::NotImplemented
}
