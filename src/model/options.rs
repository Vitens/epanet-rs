use crate::model::units::{FlowUnits, PressureUnits, UnitSystem};
use serde::{Deserialize, Serialize};

#[derive(Debug, Eq, PartialEq, Clone, Copy, Deserialize, Serialize)]
pub enum HeadlossFormula {
  HazenWilliams, // H-W
  DarcyWeisbach, // D-W
  ChezyManning,  // C-M
}

#[derive(Debug, Deserialize, Serialize)]
pub struct TimeOptions {
  pub duration: usize,      // duration of the simulation in hours
  pub hydraulic_timestep: usize, // hydraulic timestep in hours
  pub pattern_timestep: usize, // pattern timestep in hours
  pub pattern_start: usize, // pattern start in hours
  pub start_clocktime: usize, // start clocktime in hours
}

impl Default for TimeOptions {
  fn default() -> Self {
    Self {
      duration: 0,
      hydraulic_timestep: 3600,
      pattern_timestep: 3600,
      pattern_start: 0,
      start_clocktime: 0,
    }
  }
}


#[derive(Debug, Deserialize, Serialize)]
pub struct SimulationOptions {
  //
  pub flow_units: FlowUnits,
  pub pressure_units: PressureUnits,
  pub unit_system: UnitSystem,
  pub headloss_formula: HeadlossFormula,

  pub demand_multiplier: f64,

  pub max_trials: usize,
  pub accuracy: f64,
  pub max_flow_change: Option<f64>,
  pub check_frequency: usize,
  pub max_check: usize,

  pub pattern: Option<Box<str>>,

  pub time_options: TimeOptions,
}

/// Default implementation for SimulationOptions
impl Default for SimulationOptions {
  fn default() -> Self {
    Self {
      flow_units: FlowUnits::CFS,
      pressure_units: PressureUnits::FEET,
      unit_system: UnitSystem::US,
      headloss_formula: HeadlossFormula::HazenWilliams,
      demand_multiplier: 1.0,
      max_trials: 40,
      accuracy: 0.001,
      max_flow_change: None,
      check_frequency: 2,
      max_check: 10,
      pattern: None,
      time_options: TimeOptions::default(),
    }
  }
}