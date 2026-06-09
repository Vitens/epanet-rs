//! FFI analysis-options accessors: `EN_gettimeparam`, `EN_getoption`, `EN_setoption`, `EN_setdemandmodel`, etc.

use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation, get_simulation_mut};

use crate::ffi::enums::{
    DemandModel as ENDemandModel, HeadLossType, PressUnits as ENPressUnits, SimOption,
    TimeParameter,
};
use crate::model::options::{DemandModel, HeadlossFormula};
use crate::model::units::PressureUnits;

use std::os::raw::{c_double, c_int, c_long};

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
        _ => return ErrorCode::InvalidParameterCode,
    }

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_setdemandmodel(
    ph: *mut Project,
    demand_model: c_int,
    minimum_pressure: c_double,
    required_pressure: c_double,
    pressure_exponent: c_double,
) -> ErrorCode {
    let simulation = get_simulation_mut!(ph);

    let demand_model = match ENDemandModel::from_repr(demand_model) {
        Some(demand_model) => demand_model,
        None => return ErrorCode::InvalidParameterCode,
    };

    if !minimum_pressure.is_finite()
        || !required_pressure.is_finite()
        || !pressure_exponent.is_finite()
    {
        return ErrorCode::IllegalNumericValue;
    }

    if minimum_pressure < 0.0 || required_pressure <= minimum_pressure || pressure_exponent <= 0.0 {
        return ErrorCode::IllegalPdaPressureLimits;
    }

    match demand_model {
        ENDemandModel::Pda => {
            simulation.network.options.demand_model = DemandModel::PDA;
        }
        ENDemandModel::Dda => {
            simulation.network.options.demand_model = DemandModel::DDA;
        }
    }

    simulation.network.options.minimum_pressure = minimum_pressure;
    simulation.network.options.required_pressure = required_pressure;
    simulation.network.options.pressure_exponent = pressure_exponent;

    ErrorCode::Ok
}

/// # Safety
///
/// `ph` must be a valid non-null project handle returned by [`EN_createproject`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn EN_getoption(
    ph: *mut Project,
    option: c_int,
    out_value: *mut c_double,
) -> ErrorCode {
    if !out_value.is_null() {
        unsafe { *out_value = 0.0 };
    }
    let simulation = get_simulation!(ph);

    let option = match SimOption::from_repr(option) {
        Some(o) => o,
        None => return ErrorCode::InvalidParameterCode,
    };

    let v: f64 = match option {
        SimOption::Trials => simulation.network.options.max_trials as f64,
        SimOption::Accuracy => simulation.network.options.accuracy,
        SimOption::EmitExpon => {
            // Stored internally as n = 1/exponent; return the user-facing exponent.
            if simulation.network.options.emitter_exponent > 0.0 {
                1.0 / simulation.network.options.emitter_exponent
            } else {
                0.0
            }
        }
        SimOption::DemandMult => simulation.network.options.demand_multiplier,
        SimOption::FlowChange => simulation.network.options.max_flow_change.unwrap_or(0.0),
        SimOption::HeadLossForm => match simulation.network.options.headloss_formula {
            HeadlossFormula::HazenWilliams => HeadLossType::HW as i32 as f64,
            HeadlossFormula::DarcyWeisbach => HeadLossType::DW as i32 as f64,
            HeadlossFormula::ChezyManning => HeadLossType::CM as i32 as f64,
        },
        SimOption::SpGravity => simulation.network.options.specific_gravity,
        SimOption::SpViscos => simulation.network.options.viscosity,
        SimOption::CheckFreq => simulation.network.options.check_frequency as f64,
        SimOption::MaxCheck => simulation.network.options.max_check as f64,
        SimOption::DemandPattern => match &simulation.network.options.pattern {
            None => 0.0,
            Some(id) => simulation
                .network
                .patterns
                .iter()
                .position(|p| p.id.as_ref() == id.as_ref())
                .map(|i| (i + 1) as f64)
                .unwrap_or(0.0),
        },
        SimOption::PressUnits => match simulation.network.options.pressure_units {
            PressureUnits::PSI => ENPressUnits::Psi as i32 as f64,
            PressureUnits::KPA => ENPressUnits::Kpa as i32 as f64,
            PressureUnits::METERS => ENPressUnits::Meters as i32 as f64,
            PressureUnits::BAR => ENPressUnits::Bar as i32 as f64,
            PressureUnits::FEET => ENPressUnits::Feet as i32 as f64,
        },
        // Not available in this implementation — return error, not a silent zero.
        _ => return ErrorCode::InvalidParameterCode,
    };

    if !out_value.is_null() {
        unsafe { *out_value = v };
    }
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

    if !value.is_finite() {
        return ErrorCode::IllegalNumericValue;
    }

    // EN_UNBALANCED is handled before the non-negative check in EPANET. It is not
    // available in this implementation (no ExtraIter field).
    if option == SimOption::Unbalanced as c_int {
        return ErrorCode::InvalidParameterCode;
    }

    // All other option values must be non-negative.
    if value < 0.0 {
        return ErrorCode::InvalidOptionValue;
    }

    let option = match SimOption::from_repr(option) {
        Some(option) => option,
        None => return ErrorCode::InvalidParameterCode,
    };

    match option {
        SimOption::Trials => {
            if value < 1.0 {
                return ErrorCode::InvalidOptionValue;
            }
            simulation.network.options.max_trials = value as usize;
        }

        SimOption::Accuracy => {
            if value < 1.0e-8 || value > 1.0e-1 {
                return ErrorCode::InvalidOptionValue;
            }
            simulation.network.options.accuracy = value;
        }

        SimOption::EmitExpon => {
            if value <= 0.0 {
                return ErrorCode::InvalidOptionValue;
            }
            // Store as n = 1/value to match EPANET internal convention.
            simulation.network.options.emitter_exponent = 1.0 / value;
        }

        SimOption::DemandMult => {
            simulation.network.options.demand_multiplier = value;
        }

        SimOption::FlowChange => {
            simulation.network.options.max_flow_change = Some(value);
        }

        SimOption::HeadLossForm => {
            // Cannot change while the hydraulic solver is open.
            if simulation.solver.is_some() {
                return ErrorCode::ModifyWhileSolverOpen;
            }

            let i = value.round() as i32;
            let head_loss_type = match HeadLossType::from_repr(i) {
                Some(f) => f,
                None => return ErrorCode::InvalidOptionValue,
            };
            simulation.network.options.headloss_formula = match head_loss_type {
                HeadLossType::HW => HeadlossFormula::HazenWilliams,
                HeadLossType::DW => HeadlossFormula::DarcyWeisbach,
                HeadLossType::CM => HeadlossFormula::ChezyManning,
            };
        }

        SimOption::SpGravity => {
            if value <= 0.0 {
                return ErrorCode::InvalidOptionValue;
            }
            simulation.network.options.specific_gravity = value;
        }

        SimOption::SpViscos => {
            if value <= 0.0 {
                return ErrorCode::InvalidOptionValue;
            }
            simulation.network.options.viscosity = value;
        }

        SimOption::CheckFreq => {
            simulation.network.options.check_frequency = value as usize;
        }

        SimOption::MaxCheck => {
            simulation.network.options.max_check = value as usize;
        }

        SimOption::DemandPattern => {
            let pat = value.round() as i32;
            if pat < 0 || pat > simulation.network.patterns.len() as i32 {
                return ErrorCode::UndefinedPattern;
            }
            simulation.network.options.pattern = if pat == 0 {
                None
            } else {
                Some(simulation.network.patterns[(pat - 1) as usize].id.clone())
            };
        }

        SimOption::PressUnits => {
            let unit = value.round() as i32;
            // C returns 205 (UndefinedPattern) for an out-of-range pressure unit.
            let press_units = match ENPressUnits::from_repr(unit) {
                Some(u) => u,
                None => return ErrorCode::UndefinedPattern,
            };
            simulation.network.options.pressure_units = match press_units {
                ENPressUnits::Psi => PressureUnits::PSI,
                ENPressUnits::Kpa => PressureUnits::KPA,
                ENPressUnits::Meters => PressureUnits::METERS,
                ENPressUnits::Bar => PressureUnits::BAR,
                ENPressUnits::Feet => PressureUnits::FEET,
            };
        }

        // Not available in this implementation.
        _ => return ErrorCode::InvalidParameterCode,
    }

    ErrorCode::Ok
}
