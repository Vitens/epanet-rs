use crate::model::units::{FlowUnits, PressureUnits, UnitSystem};
use serde::{Deserialize, Serialize};

#[derive(Debug, Eq, PartialEq, Clone, Copy, Deserialize, Serialize, strum::Display)]
pub enum HeadlossFormula {
  HazenWilliams, // H-W
  DarcyWeisbach, // D-W
  ChezyManning,  // C-M
}

#[derive(Debug, Eq, PartialEq, Clone, Copy, Deserialize, Serialize, strum::Display)]
pub enum DemandModel {
  // Pressure Driven Analysis
  PDA,
  // Demand Driven Analysis (default)
  DDA,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct TimeOptions {
  pub duration: usize,      // duration of the simulation in hours
  pub hydraulic_timestep: usize, // hydraulic timestep in hours
  pub report_timestep: usize, // report timestep in hours
  pub pattern_timestep: usize, // pattern timestep in hours
  pub pattern_start: usize, // pattern start in hours
  pub start_clocktime: usize, // start clocktime in hours
}

impl Default for TimeOptions {
  fn default() -> Self {
    Self {
      duration: 0,
      hydraulic_timestep: 3600,
      report_timestep: 3600,
      pattern_timestep: 3600,
      pattern_start: 0,
      start_clocktime: 0,
    }
  }
}


#[derive(Debug, Deserialize, Serialize, Clone)]
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

  pub specific_gravity: f64,
  pub viscosity: f64,

  pub emitter_exponent: f64,

  pub pattern: Option<Box<str>>,

  pub time_options: TimeOptions,

  pub demand_model: DemandModel,

  pub minimum_pressure: f64,
  pub required_pressure: f64,
  pub pressure_exponent: f64,
}

/// Default implementation for SimulationOptions
impl Default for SimulationOptions {
  fn default() -> Self {
    Self {
      flow_units: FlowUnits::CFS,
      pressure_units: PressureUnits::PSI,
      unit_system: UnitSystem::US,
      headloss_formula: HeadlossFormula::HazenWilliams,
      demand_multiplier: 1.0,
      max_trials: 40,
      accuracy: 0.001,
      max_flow_change: None,
      check_frequency: 2,
      emitter_exponent: 2.0,
      specific_gravity: 1.0,
      viscosity: 1.0,
      max_check: 10,
      pattern: None,
      time_options: TimeOptions::default(),
      demand_model: DemandModel::DDA,
      minimum_pressure: 0.0,
      required_pressure: 10.0,
      pressure_exponent: 0.5,
    }
  }
}

// unit conversion methods for simulationoptions
impl SimulationOptions {
  pub fn convert_to_standard(&mut self) {

    self.minimum_pressure = self.minimum_pressure / self.pressure_units.per_feet();
    self.required_pressure = self.required_pressure / self.pressure_units.per_feet();
  }
  pub fn convert_from_standard(&mut self) {
    self.minimum_pressure = self.minimum_pressure * self.pressure_units.per_feet();
    self.required_pressure = self.required_pressure * self.pressure_units.per_feet();
  }
}