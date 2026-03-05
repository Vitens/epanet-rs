use std::str::FromStr;

use crate::constants::*;
use serde::{Deserialize, Serialize};


// type aliases for internal units
// no enforcement (for now), only for additional clarity
pub type Ft = f64;
pub type Cfs = f64;
pub type Psi = f64;
pub type Ft3 = f64;


#[derive(Debug, PartialEq, Eq, Clone, Copy, Deserialize, Serialize)]
pub enum FlowUnits {
  // Imperial units
  CFS,  // Cubic feet per second
  GPM,  // Gallons per minute
  MGD,  // Million gallons per day
  IMGD, // Imperial million gallons per day
  AFD,  // Acre-feet per day
  // Metric units
  LPS,  // Liters per second
  LPM,  // Liters per minute
  MLD,  // Million liters per day
  CMS,  // Cubic meters per second
  CMH,  // Cubic meters per hour
  CMD   // Cubic meters per day
}

/// FromStr implementation for FlowUnits
impl FromStr for FlowUnits {
  type Err = String;
  fn from_str(s: &str) -> Result<Self, Self::Err> {
    match s.to_uppercase().as_str() {
      "CFS" => Ok(FlowUnits::CFS),
      "GPM" => Ok(FlowUnits::GPM),
      "MGD" => Ok(FlowUnits::MGD),
      "IMGD" => Ok(FlowUnits::IMGD),
      "AFD" => Ok(FlowUnits::AFD),
      "LPS" => Ok(FlowUnits::LPS),
      "LPM" => Ok(FlowUnits::LPM),
      "MLD" => Ok(FlowUnits::MLD),
      "CMS" => Ok(FlowUnits::CMS),
      "CMH" => Ok(FlowUnits::CMH),
      "CMD" => Ok(FlowUnits::CMD),
      _ => Err(format!("Invalid flow unit: {}", s)),
    }
  }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Deserialize, Serialize)]
pub enum PressureUnits {
  PSI,    // Pounds per square inch
  KPA,    // Kilopascals
  METERS, // Meters
  FEET,   // Feet
  BAR,    // Bar
}
impl FromStr for PressureUnits {
  type Err = String;
  fn from_str(s: &str) -> Result<Self, Self::Err> {
    match s.to_uppercase().as_str() {
      "PSI" => Ok(PressureUnits::PSI),
      "KPA" => Ok(PressureUnits::KPA),
      "METERS" => Ok(PressureUnits::METERS),
      "FEET" => Ok(PressureUnits::FEET),
      "BAR" => Ok(PressureUnits::BAR),
      _ => Err(format!("Invalid pressure unit: {}", s)),
    }
  }
}

#[derive(Debug, PartialEq, Eq, Deserialize, Serialize)]
pub enum UnitSystem {
  US, // US Customary units
  SI, // International System of Units (metric)
}

impl UnitSystem {
  pub fn per_feet(&self) -> f64 {
    match self {
      UnitSystem::US => 1.0,
      UnitSystem::SI => MperFT,
    }
  }
}

impl FlowUnits {
  // convert the flow units to standard units (US standard) and CFS
  pub fn per_cfs(&self) -> f64 {
    match self {
      FlowUnits::CFS => 1.0,
      FlowUnits::GPM => GPMperCFS,
      FlowUnits::MGD => MGDperCFS,
      FlowUnits::IMGD => IMGDperCFS,
      FlowUnits::AFD => AFDperCFS,
      FlowUnits::LPS => LPSperCFS,
      FlowUnits::LPM => LPMperCFS,
      FlowUnits::MLD => MLDperCFS,
      FlowUnits::CMS => CMSperCFS,
      FlowUnits::CMH => CMHperCFS,
      FlowUnits::CMD => CMDperCFS,
    }
  }
}

impl PressureUnits {
  pub fn per_feet(&self) -> f64 {
    match self {
      PressureUnits::PSI => PSIperFT,
      PressureUnits::KPA => KPAperFT,
      PressureUnits::METERS => MperFT,
      PressureUnits::FEET => 1.0,
      PressureUnits::BAR => BARperFT,
    }
  }
}

pub trait UnitConversion {
  fn convert_units(&mut self, flow: &FlowUnits, system: &UnitSystem, reverse: bool);
}
