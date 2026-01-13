use std::str::FromStr;

#[derive(Debug)]
pub enum FlowUnits {
  // Imperial units
  CFS,  // Cubic feet per second
  GMP,  // Gallons per minute
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
      "GMP" => Ok(FlowUnits::GMP),
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

#[derive(Debug)]
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

#[derive(Debug)]
pub enum UnitSystem {
  US, // US Customary units
  SI, // International System of Units (metric)
}

#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub enum HeadlossFormula {
  HazenWilliams, // H-W
  DarcyWeisbach, // D-W
  ChezyManning,  // C-M
}


#[derive(Debug)]
pub struct SimulationOptions {
  //
  pub flow_units: FlowUnits,
  pub pressure_units: PressureUnits,
  pub unit_system: UnitSystem,
  pub headloss_formula: HeadlossFormula,

  pub max_trials: usize,
  pub accuracy: f64,
  pub check_frequency: usize,
  pub max_check: usize,

  pub pattern: Option<Box<str>>,

}

/// Default implementation for SimulationOptions
impl Default for SimulationOptions {
  fn default() -> Self {
    Self {
      flow_units: FlowUnits::CFS,
      pressure_units: PressureUnits::PSI,
      unit_system: UnitSystem::US,
      headloss_formula: HeadlossFormula::HazenWilliams,
      max_trials: 40,
      accuracy: 0.001,
      check_frequency: 2,
      max_check: 10,
      pattern: None,
    }
  }
}