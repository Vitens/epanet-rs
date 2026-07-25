//! FFI analysis-options accessors: `EN_gettimeparam`, `EN_setoption`, `EN_setdemandmodel`, etc.

use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};
use crate::ffi::util::write_out;

use crate::ffi::enums::{
    DemandModel as ENDemandModel, FlowUnits as ENFlowUnits, HeadLossType,
    PressUnits as ENPressUnits, SimOption, StatisticType, StatusReport, TimeParameter,
};
use crate::model::options::{DemandModel, HeadlossFormula};
use crate::model::units::{FlowUnits, PressureUnits, UnitSystem};

use std::os::raw::{c_double, c_int, c_long};

/// Minimum allowed difference between the required and minimum pressure of the
/// pressure-driven demand model (EPANET's `MINPDIFF`).
const MIN_PRESSURE_DIFFERENCE: f64 = 0.1;

/// Maps the internal flow units onto their toolkit code.
fn flow_units_code(units: &FlowUnits) -> ENFlowUnits {
    match units {
        FlowUnits::CFS => ENFlowUnits::Cfs,
        FlowUnits::GPM => ENFlowUnits::Gpm,
        FlowUnits::MGD => ENFlowUnits::Mgd,
        FlowUnits::IMGD => ENFlowUnits::Imgd,
        FlowUnits::AFD => ENFlowUnits::Afd,
        FlowUnits::LPS => ENFlowUnits::Lps,
        FlowUnits::LPM => ENFlowUnits::Lpm,
        FlowUnits::MLD => ENFlowUnits::Mld,
        FlowUnits::CMH => ENFlowUnits::Cmh,
        FlowUnits::CMD => ENFlowUnits::Cmd,
        FlowUnits::CMS => ENFlowUnits::Cms,
    }
}

/// Maps a toolkit flow-units code onto the internal flow units.
pub(crate) fn flow_units_from_code(code: ENFlowUnits) -> FlowUnits {
    match code {
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
    }
}

/// Maps the internal pressure units onto their toolkit code.
fn pressure_units_code(units: &PressureUnits) -> ENPressUnits {
    match units {
        PressureUnits::PSI => ENPressUnits::Psi,
        PressureUnits::KPA => ENPressUnits::Kpa,
        PressureUnits::METERS => ENPressUnits::Meters,
        PressureUnits::BAR => ENPressUnits::Bar,
        PressureUnits::FEET => ENPressUnits::Feet,
    }
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_value` must be null or a valid writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_gettimeparam(
    ph: *mut Project,
    param: c_int,
    out_value: *mut c_long,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let times = &simulation.network.options.time_options;

    let param = match TimeParameter::from_repr(param) {
        Some(param) => param,
        None => return ErrorCode::InvalidParameterCode,
    };

    let value = match param {
        TimeParameter::Duration => times.duration,
        TimeParameter::HydStep => times.hydraulic_timestep,
        TimeParameter::PatternStep => times.pattern_timestep,
        TimeParameter::PatternStart => times.pattern_start,
        TimeParameter::ReportStep => times.report_timestep,
        TimeParameter::StartTime => times.start_clocktime,
        TimeParameter::Periods => times.duration / times.report_timestep + 1,
        TimeParameter::HTime => simulation.time,
        // TODO: report start, rule step and the water quality clock are not
        // modelled; a statistic other than the full time series is never used.
        TimeParameter::ReportStart | TimeParameter::RuleStep => 0,
        TimeParameter::QualStep | TimeParameter::QTime => 0,
        TimeParameter::Statistic => StatisticType::Series as usize,
        TimeParameter::HaltFlag => 0,
        TimeParameter::NextEvent | TimeParameter::NextEventTank => {
            // TODO: the simulation driver does not expose the next event time.
            return ErrorCode::NotImplemented;
        }
    };

    unsafe { write_out(out_value, value as c_long) };
    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_settimeparam(
    ph: *mut Project,
    param: c_int,
    value: c_long,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let time_options = &mut simulation.network.options.time_options;

    if value < 0 {
        return ErrorCode::IllegalNumericValue;
    }
    let value = value as usize;

    let param = match TimeParameter::from_repr(param) {
        Some(param) => param,
        None => return ErrorCode::InvalidParameterCode,
    };

    match param {
        TimeParameter::Duration => {
            time_options.duration = value;
            if time_options.duration > time_options.pattern_start {
                time_options.pattern_start = time_options.duration;
            }
        }
        TimeParameter::HydStep => {
            if value == 0 {
                return ErrorCode::InvalidOptionValue;
            }
            time_options.hydraulic_timestep = value;
            time_options.hydraulic_timestep = time_options
                .hydraulic_timestep
                .min(time_options.pattern_timestep);
            time_options.hydraulic_timestep = time_options
                .hydraulic_timestep
                .min(time_options.report_timestep);
        }
        TimeParameter::PatternStep => {
            if value == 0 {
                return ErrorCode::InvalidOptionValue;
            }
            time_options.pattern_timestep = value;
            if time_options.hydraulic_timestep > time_options.pattern_timestep {
                time_options.hydraulic_timestep = time_options.pattern_timestep;
            }
        }
        TimeParameter::PatternStart => {
            time_options.pattern_start = value;
        }
        TimeParameter::ReportStep => {
            if value == 0 {
                return ErrorCode::InvalidOptionValue;
            }
            time_options.report_timestep = value;
            if time_options.hydraulic_timestep > time_options.report_timestep {
                time_options.hydraulic_timestep = time_options.report_timestep;
            }
        }
        TimeParameter::StartTime => {
            time_options.start_clocktime = value;
        }
        // TODO: report start, rule step, the reporting statistic and the water
        // quality time step are not part of the simulation options.
        TimeParameter::ReportStart
        | TimeParameter::RuleStep
        | TimeParameter::Statistic
        | TimeParameter::QualStep => return ErrorCode::NotImplemented,
        _ => return ErrorCode::InvalidParameterCode,
    }

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_units` must be null or a valid writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getflowunits(ph: *mut Project, out_units: *mut c_int) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let units = flow_units_code(&simulation.network.options.flow_units);
    unsafe { write_out(out_units, units as c_int) };

    ErrorCode::Ok
}

/// Sets a project's flow units.
///
/// Network data is stored internally in US customary units, so only the units
/// used to report and accept values change.
///
/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setflowunits(ph: *mut Project, units: c_int) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let units = match ENFlowUnits::from_repr(units) {
        Some(units) => flow_units_from_code(units),
        None => return ErrorCode::InvalidParameterCode,
    };

    let options = &mut simulation.network.options;
    let unit_system = match units {
        FlowUnits::CFS | FlowUnits::GPM | FlowUnits::MGD | FlowUnits::IMGD | FlowUnits::AFD => {
            UnitSystem::US
        }
        FlowUnits::LPS
        | FlowUnits::LPM
        | FlowUnits::MLD
        | FlowUnits::CMS
        | FlowUnits::CMH
        | FlowUnits::CMD => UnitSystem::SI,
    };

    // switching between unit systems also switches the default pressure units
    if unit_system != options.unit_system {
        options.pressure_units = match unit_system {
            UnitSystem::US => PressureUnits::PSI,
            UnitSystem::SI => PressureUnits::METERS,
        };
    }
    options.flow_units = units;
    options.unit_system = unit_system;

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// The output pointers must be null or valid writable pointers.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getdemandmodel(
    ph: *mut Project,
    out_type: *mut c_int,
    out_pmin: *mut c_double,
    out_preq: *mut c_double,
    out_pexp: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let options = &simulation.network.options;
    let demand_model = match options.demand_model {
        DemandModel::DDA => ENDemandModel::Dda,
        DemandModel::PDA => ENDemandModel::Pda,
    };
    let per_feet = options.pressure_units.per_feet();

    unsafe { write_out(out_type, demand_model as c_int) };
    unsafe { write_out(out_pmin, options.minimum_pressure * per_feet) };
    unsafe { write_out(out_preq, options.required_pressure * per_feet) };
    unsafe { write_out(out_pexp, options.pressure_exponent) };

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setdemandmodel(
    ph: *mut Project,
    demand_model: c_int,
    pmin: c_double,
    preq: c_double,
    pexp: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let demand_model = match ENDemandModel::from_repr(demand_model) {
        Some(demand_model) => demand_model,
        None => return ErrorCode::InvalidParameterCode,
    };

    let options = &mut simulation.network.options;

    if demand_model == ENDemandModel::Pda {
        if pmin < 0.0 || preq - pmin < MIN_PRESSURE_DIFFERENCE {
            return ErrorCode::IllegalPdaPressureLimits;
        }
        if pexp <= 0.0 {
            return ErrorCode::IllegalNumericValue;
        }

        let per_feet = options.pressure_units.per_feet();
        options.minimum_pressure = pmin / per_feet;
        options.required_pressure = preq / per_feet;
        options.pressure_exponent = pexp;
    }

    options.demand_model = match demand_model {
        ENDemandModel::Pda => DemandModel::PDA,
        ENDemandModel::Dda => DemandModel::DDA,
    };

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
/// `out_value` must be null or a valid writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getoption(
    ph: *mut Project,
    option: c_int,
    out_value: *mut c_double,
) -> ErrorCode {
    let simulation = get_simulation!(ph);

    let option = match SimOption::from_repr(option) {
        Some(option) => option,
        None => return ErrorCode::InvalidParameterCode,
    };

    let options = &simulation.network.options;

    let value = match option {
        SimOption::Trials => options.max_trials as f64,
        SimOption::Accuracy => options.accuracy,
        SimOption::EmitExpon => {
            // stored as the reciprocal of the user-facing exponent
            if options.emitter_exponent > 0.0 {
                1.0 / options.emitter_exponent
            } else {
                0.0
            }
        }
        SimOption::DemandMult => options.demand_multiplier,
        SimOption::FlowChange => options.max_flow_change.unwrap_or(0.0),
        SimOption::HeadLossForm => match options.headloss_formula {
            HeadlossFormula::HazenWilliams => HeadLossType::HW as i32 as f64,
            HeadlossFormula::DarcyWeisbach => HeadLossType::DW as i32 as f64,
            HeadlossFormula::ChezyManning => HeadLossType::CM as i32 as f64,
        },
        SimOption::SpGravity => options.specific_gravity,
        SimOption::SpViscos => options.viscosity,
        SimOption::CheckFreq => options.check_frequency as f64,
        SimOption::MaxCheck => options.max_check as f64,
        SimOption::DampLimit => options.damping_limit,
        SimOption::DemandPattern => options
            .pattern
            .as_ref()
            .and_then(|id| simulation.network.pattern_map.get(id))
            .map(|index| (index + 1) as f64)
            .unwrap_or(0.0),
        SimOption::PressUnits => pressure_units_code(&options.pressure_units) as i32 as f64,
        // an unbalanced solve always aborts and no status report is produced
        SimOption::Unbalanced => 0.0,
        SimOption::StatusReport => StatusReport::NoReport as i32 as f64,
        // emitters are always allowed to admit backflow
        SimOption::EmitBackflow => 1.0,
        // TODO: head-error convergence, energy costs and water quality
        // reaction options are not part of the simulation options.
        SimOption::HeadError
        | SimOption::GlobalEffic
        | SimOption::GlobalPrice
        | SimOption::GlobalPattern
        | SimOption::DemandCharge
        | SimOption::Tolerance
        | SimOption::SpDiffus
        | SimOption::BulkOrder
        | SimOption::WallOrder
        | SimOption::TankOrder
        | SimOption::ConcenLimit => return ErrorCode::NotImplemented,
    };

    unsafe { write_out(out_value, value) };
    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setoption(
    ph: *mut Project,
    option: c_int,
    value: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let option = match SimOption::from_repr(option) {
        Some(option) => option,
        None => return ErrorCode::InvalidParameterCode,
    };

    // resolve the demand pattern id before taking a mutable borrow of the options
    let demand_pattern = if option == SimOption::DemandPattern && value >= 1.0 {
        match simulation.network.patterns.get(value as usize - 1) {
            Some(pattern) => Some(pattern.id.clone()),
            None => return ErrorCode::UndefinedPattern,
        }
    } else {
        None
    };

    let options = &mut simulation.network.options;

    match option {
        SimOption::Trials => {
            if value < 1.0 {
                return ErrorCode::IllegalNumericValue;
            }
            options.max_trials = value as usize;
        }
        SimOption::Accuracy => {
            if value <= 0.0 {
                return ErrorCode::IllegalNumericValue;
            }
            options.accuracy = value;
        }
        SimOption::EmitExpon => {
            if value <= 0.0 {
                return ErrorCode::IllegalNumericValue;
            }
            // TODO: EPANET also rescales every junction's emitter coefficient
            // when the exponent changes.
            options.emitter_exponent = 1.0 / value;
        }
        SimOption::DemandMult => {
            if value <= 0.0 {
                return ErrorCode::IllegalNumericValue;
            }
            options.demand_multiplier = value;
        }
        SimOption::FlowChange => {
            if value < 0.0 {
                return ErrorCode::IllegalNumericValue;
            }
            options.max_flow_change = if value == 0.0 { None } else { Some(value) };
        }
        SimOption::SpGravity => {
            if value <= 0.0 {
                return ErrorCode::IllegalNumericValue;
            }
            options.specific_gravity = value;
        }
        SimOption::SpViscos => {
            if value <= 0.0 {
                return ErrorCode::IllegalNumericValue;
            }
            options.viscosity = value;
        }
        SimOption::CheckFreq => {
            if value < 0.0 {
                return ErrorCode::IllegalNumericValue;
            }
            options.check_frequency = value as usize;
        }
        SimOption::MaxCheck => {
            if value < 0.0 {
                return ErrorCode::IllegalNumericValue;
            }
            options.max_check = value as usize;
        }
        SimOption::DampLimit => {
            if value < 0.0 {
                return ErrorCode::IllegalNumericValue;
            }
            options.damping_limit = value;
        }
        SimOption::DemandPattern => {
            if value < 0.0 {
                return ErrorCode::IllegalNumericValue;
            }
            options.pattern = demand_pattern;
        }
        SimOption::PressUnits => {
            let units = match ENPressUnits::from_repr(value as c_int) {
                Some(units) => units,
                None => return ErrorCode::InvalidOptionValue,
            };
            options.pressure_units = match units {
                ENPressUnits::Psi => PressureUnits::PSI,
                ENPressUnits::Kpa => PressureUnits::KPA,
                ENPressUnits::Meters => PressureUnits::METERS,
                ENPressUnits::Bar => PressureUnits::BAR,
                ENPressUnits::Feet => PressureUnits::FEET,
            };
        }
        // TODO: the head loss formula is baked into every pipe when it is
        // created, and the remaining options (head-error convergence,
        // unbalanced handling, status reporting, energy costs and water
        // quality reactions) are not part of the simulation options.
        _ => return ErrorCode::NotImplemented,
    }

    ErrorCode::Ok
}
