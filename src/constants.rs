//! Physical constants, unit-conversion factors and solver tolerances.

#![allow(non_upper_case_globals)]

/// Tolerance for zero flow
pub const Q_ZERO: f64 = 1e-6;
/// Small value for small numbers
pub const SMALL_VALUE: f64 = 1e-6;
/// Big value for large numbers
pub const BIG_VALUE: f64 = 1e8;
/// Tiny value for very small numbers
pub const TINY: f64 = 1e-6;
/// Minimum difference for Pressure Driven Analysis (PDA)
pub const PDA_MIN_DIFF : f64 = 0.1;

/// Viscosity of water
pub const VISCOSITY: f64 = 1.1e-5;
// PI (3.14159)
pub const PI : f64 = std::f64::consts::PI;

// Tolerance for head loss gradient
pub const RQ_TOL: f64 = 1e-7;
// Tolerance for flow difference
pub const Q_TOL : f64 = 0.0001;
/// Tolerance for head difference
pub const H_TOL : f64 = 0.0005;

/// Gallons per minute per cubic feet per second
pub const GPMperCFS: f64 =  448.831;
/// Million gallons per day per cubic feet per second
pub const MGDperCFS: f64 =  0.64632;
/// Imperial million gallons per day per cubic feet per second
pub const IMGDperCFS: f64 = 0.5382;
/// Acre-feet per day per cubic feet per second
pub const AFDperCFS: f64 =  1.9837;

/// Liters per second per cubic feet per second
pub const LPSperCFS: f64 =  28.317;
/// Liters per minute per cubic feet per second
pub const LPMperCFS: f64 =  1699.0;
/// Cubic meters per second per cubic feet per second
pub const CMSperCFS: f64 =  0.028317;
/// Cubic meters per hour per cubic feet per second
pub const CMHperCFS: f64 =  101.94;
/// Cubic meters per day per cubic feet per second
pub const CMDperCFS: f64 =  2446.6;
/// Million liters per day per cubic feet per second
pub const MLDperCFS: f64 =  2.4466;

// Conversion factors for SI quantities
/// Cubic meters per cubic feet per second
pub const M3perFT3 : f64  =  0.028317;
/// Meters per feet
pub const MperFT : f64    =  0.3048;
/// Kilowatts per horsepower
pub const KWperHP : f64   =  0.7457;

// Conversion factors for pressure units
/// Bar per feet
pub const BARperFT : f64 =  0.02953;
/// Kilopascals per feet
pub const KPAperFT : f64 =  2.988;
/// Pounds per square inch per feet
pub const PSIperFT : f64 =  0.4333;
