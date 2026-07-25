//! FFI accessors for rule-based controls.
//!
//! `epanet-rs` stores rules as structured conditions/actions rather than as the
//! flat premise/action lists of the original toolkit, and its rule parser is
//! internal to the input reader. Functions that would need to build or rewrite
//! individual clauses therefore return
//! [`ErrorCode::NotImplemented`](crate::ffi::error_codes::ErrorCode::NotImplemented).

use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::ffi::util::{MAX_ID, write_out, write_str};

use crate::model::link::LinkStatus;

use std::os::raw::{c_char, c_double, c_int};

/// Adds a new rule-based control, given its text representation.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_addrule(_ph: *mut Project, _rule: *const c_char) -> ErrorCode {
    // TODO: the rule text parser lives in the (private) input reader, so rules
    // can currently only be added by loading an .inp file.
    ErrorCode::NotImplemented
}

/// Deletes a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_deleterule(ph: *mut Project, index: c_int) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let index = (index - 1) as usize;
    if index >= simulation.network.rules.len() {
        return ErrorCode::NonexistentRule;
    }

    let rule_id = simulation.network.rules[index].id.clone();
    simulation.network.rules.remove(index);
    simulation.network.rule_map.remove(&rule_id);

    // the map stores positions, so everything after the removed rule shifts down
    for position in simulation.network.rule_map.values_mut() {
        if *position > index {
            *position -= 1;
        }
    }

    simulation.network.properties_version += 1;

    ErrorCode::Ok
}

/// Retrieves summary information about a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// Every out-parameter must either be null or point to writable memory.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getrule(
    ph: *mut Project,
    index: c_int,
    out_n_premises: *mut c_int,
    out_n_then_actions: *mut c_int,
    out_n_else_actions: *mut c_int,
    out_priority: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let Some(rule) = simulation.network.rules.get((index - 1) as usize) else {
        return ErrorCode::NonexistentRule;
    };

    let then_actions = rule
        .actions
        .iter()
        .filter(|action| !action.default_active)
        .count();

    unsafe {
        write_out(out_n_premises, rule.conditions.len() as c_int);
        write_out(out_n_then_actions, then_actions as c_int);
        write_out(
            out_n_else_actions,
            (rule.actions.len() - then_actions) as c_int,
        );
        write_out(out_priority, rule.priority.unwrap_or(0) as c_double);
    }

    ErrorCode::Ok
}

/// Retrieves the ID name of a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_id` must point to a buffer of at least `EN_MAXID + 1` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getruleID(
    ph: *mut Project,
    index: c_int,
    out_id: *mut c_char,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let Some(rule) = simulation.network.rules.get((index - 1) as usize) else {
        return ErrorCode::NonexistentRule;
    };

    unsafe { write_str(out_id, &rule.id, MAX_ID) };

    ErrorCode::Ok
}

/// Sets the priority of a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setrulepriority(
    ph: *mut Project,
    index: c_int,
    priority: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let Some(rule) = simulation.network.rules.get_mut((index - 1) as usize) else {
        return ErrorCode::NonexistentRule;
    };

    if priority < 0.0 {
        return ErrorCode::IllegalNumericValue;
    }

    rule.priority = if priority == 0.0 {
        None
    } else {
        Some(priority as usize)
    };

    simulation.network.properties_version += 1;

    ErrorCode::Ok
}

/// Retrieves the properties of a premise in a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn EN_getpremise(
    _ph: *mut Project,
    _rule_index: c_int,
    _premise_index: c_int,
    _out_logop: *mut c_int,
    _out_object: *mut c_int,
    _out_obj_index: *mut c_int,
    _out_variable: *mut c_int,
    _out_relop: *mut c_int,
    _out_status: *mut c_int,
    _out_value: *mut c_double,
) -> ErrorCode {
    // TODO: premises are stored as typed targets/attributes rather than the
    // toolkit's numeric object/variable codes.
    ErrorCode::NotImplemented
}

/// Sets the properties of a premise in a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn EN_setpremise(
    _ph: *mut Project,
    _rule_index: c_int,
    _premise_index: c_int,
    _logop: c_int,
    _object: c_int,
    _obj_index: c_int,
    _variable: c_int,
    _relop: c_int,
    _status: c_int,
    _value: c_double,
) -> ErrorCode {
    // TODO: see EN_getpremise.
    ErrorCode::NotImplemented
}

/// Sets the object index of a premise in a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setpremiseindex(
    _ph: *mut Project,
    _rule_index: c_int,
    _premise_index: c_int,
    _obj_index: c_int,
) -> ErrorCode {
    // TODO: see EN_getpremise.
    ErrorCode::NotImplemented
}

/// Sets the status a premise in a rule-based control compares against.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setpremisestatus(
    _ph: *mut Project,
    _rule_index: c_int,
    _premise_index: c_int,
    _status: c_int,
) -> ErrorCode {
    // TODO: see EN_getpremise.
    ErrorCode::NotImplemented
}

/// Sets the value a premise in a rule-based control compares against.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setpremisevalue(
    _ph: *mut Project,
    _rule_index: c_int,
    _premise_index: c_int,
    _value: c_double,
) -> ErrorCode {
    // TODO: see EN_getpremise.
    ErrorCode::NotImplemented
}

/// Shared implementation of `EN_getthenaction` / `EN_getelseaction`.
///
/// # Safety
///
/// See the callers.
unsafe fn get_action(
    ph: *mut Project,
    rule_index: c_int,
    action_index: c_int,
    out_link_index: *mut c_int,
    out_status: *mut c_int,
    out_setting: *mut c_double,
    else_action: bool,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let Some(rule) = simulation.network.rules.get((rule_index - 1) as usize) else {
        return ErrorCode::NonexistentRule;
    };

    let action = rule
        .actions
        .iter()
        .filter(|action| action.default_active == else_action)
        .nth((action_index - 1).max(0) as usize);

    let Some(action) = action else {
        return ErrorCode::NonexistentRuleClause;
    };

    let Some(&link_index) = simulation.network.link_map.get(&action.link_id) else {
        return ErrorCode::UndefinedLink;
    };

    // rule statuses use the same 1/2/3 numbering as EN_R_IS_OPEN/CLOSED/ACTIVE
    let status = match action.status {
        Some(LinkStatus::Open) | Some(LinkStatus::FixedOpen) => 1,
        Some(LinkStatus::Closed) | Some(LinkStatus::FixedClosed) => 2,
        Some(LinkStatus::Active) => 3,
        _ => 0,
    };

    unsafe {
        write_out(out_link_index, (link_index + 1) as c_int);
        write_out(out_status, status);
        write_out(out_setting, action.setting.unwrap_or(0.0));
    }

    ErrorCode::Ok
}

/// Retrieves a THEN action of a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getthenaction(
    ph: *mut Project,
    rule_index: c_int,
    action_index: c_int,
    out_link_index: *mut c_int,
    out_status: *mut c_int,
    out_setting: *mut c_double,
) -> ErrorCode {
    unsafe {
        get_action(
            ph,
            rule_index,
            action_index,
            out_link_index,
            out_status,
            out_setting,
            false,
        )
    }
}

/// Retrieves an ELSE action of a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getelseaction(
    ph: *mut Project,
    rule_index: c_int,
    action_index: c_int,
    out_link_index: *mut c_int,
    out_status: *mut c_int,
    out_setting: *mut c_double,
) -> ErrorCode {
    unsafe {
        get_action(
            ph,
            rule_index,
            action_index,
            out_link_index,
            out_status,
            out_setting,
            true,
        )
    }
}

/// Sets a THEN action of a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setthenaction(
    _ph: *mut Project,
    _rule_index: c_int,
    _action_index: c_int,
    _link_index: c_int,
    _status: c_int,
    _setting: c_double,
) -> ErrorCode {
    // TODO: rule actions store link ids and settings in internal units; safe
    // rewriting requires validation that epanet-rs does not expose yet.
    ErrorCode::NotImplemented
}

/// Sets an ELSE action of a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setelseaction(
    _ph: *mut Project,
    _rule_index: c_int,
    _action_index: c_int,
    _link_index: c_int,
    _status: c_int,
    _setting: c_double,
) -> ErrorCode {
    // TODO: see EN_setthenaction.
    ErrorCode::NotImplemented
}

/// Gets the enabled/disabled status of a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getruleenabled(
    _ph: *mut Project,
    _index: c_int,
    _out_enabled: *mut c_int,
) -> ErrorCode {
    // TODO: rules in epanet-rs are always active; there is no enabled flag.
    ErrorCode::NotImplemented
}

/// Enables or disables a rule-based control.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setruleenabled(
    _ph: *mut Project,
    _index: c_int,
    _enabled: c_int,
) -> ErrorCode {
    // TODO: see EN_getruleenabled.
    ErrorCode::NotImplemented
}
