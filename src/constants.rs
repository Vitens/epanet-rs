#![allow(non_upper_case_globals)]

pub const SMALL_VALUE: f64 = 1e-6;
pub const BIG_VALUE: f64 = 1e8;
pub const TINY: f64 = 1e-6;

pub const VISCOSITY: f64 = 1.1e-5;
pub const PI : f64 = std::f64::consts::PI;

// Tolerance for head loss gradient
pub const RQ_TOL: f64 = 1e-7;
// Tolerance for flow difference
pub const Q_TOL : f64 = 0.0001;
pub const H_TOL : f64 = 0.0005;

// allow non uppercase consts
// Conversion factors for US flow units
pub const GPMperCFS: f64 =  448.831;
pub const MGDperCFS: f64 =  0.64632;
pub const IMGDperCFS: f64 = 0.5382;
pub const AFDperCFS: f64 =  1.9837;

// Conversion factors for SI flow units
pub const LPSperCFS: f64 =  28.317;
pub const LPMperCFS: f64 =  1699.0;
pub const CMSperCFS: f64 =  0.028317;
pub const CMHperCFS: f64 =  101.94;
pub const CMDperCFS: f64 =  2446.6;
pub const MLDperCFS: f64 =  2.4466;

// Conversion factors for SI quantities
// pub const LperFT3 : f64   =  28.317;
pub const M3perFT3 : f64  =  0.028317;
pub const MperFT : f64    =  0.3048;
// pub const KPAperPSI : f64 =  6.895;
// pub const KWperHP : f64   =  0.7457;

// Conversion factors for pressure units
pub const BARperFT : f64 =  0.02953;
pub const KPAperFT : f64 =  2.988;
pub const PSIperFT : f64 =  0.4333;
