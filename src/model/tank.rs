use crate::model::units::UnitConversion;
use crate::model::units::{FlowUnits, UnitSystem};
use crate::constants::*;
use crate::model::curve::Curve;
use std::sync::Arc;
use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize)]
pub struct Tank {
  pub elevation: f64,        // elevation of the tank (ft)
  pub initial_level: f64,    // initial level of the tank (ft)
  pub min_level: f64,        // minimum level of the tank (ft)
  pub max_level: f64,        // maximum level of the tank (ft)
  pub diameter: f64,        // nominal diameter of the tank (ft)
  pub min_volume: f64,      // minimum volume of the tank (ft^3)
  pub volume_curve_id: Option<Box<str>>, // id of the volume curve
  pub overflow: bool,                   // whether the tank has an overflow
  #[serde(skip)]
  pub volume_curve: Option<Arc<Curve>>, // volume curve
  #[serde(skip)]
  pub links_to: Vec<usize>, // indices of the links connected from N -> Tank
  #[serde(skip)]
  pub links_from: Vec<usize>, // indices of the links connected from Tank -> N
}

impl Tank {
  // calculate the time to fill or drain the tank given the head and flow
  pub fn time_to_fill_or_drain(&self, head: f64, flow: f64) -> usize {

    let current_volume = self.volume_at_head(head);

    let max_volume = self.max_volume();
    let min_volume = self.min_volume();

    // if flow is positive, the tank is filling
    if flow > 0.0 {
      let time_to_fill = (max_volume - current_volume) / flow;
      return time_to_fill as usize;
    }
    // if flow is negative, the tank is draining
    else {
      let time_to_drain = (current_volume - min_volume) / flow;
      return time_to_drain as usize;
    }

  }

  pub fn volume_at_head(&self, head: f64) -> f64 {
    let level = head - self.elevation;
    if self.volume_curve.is_some() {
      panic!("Tank volume curves not supported yet!");
    }
    else {
      return level * PI * self.diameter * self.diameter / 4.0;
    }
  }
  pub fn min_volume(&self) -> f64 {
    if self.volume_curve.is_some() {
      panic!("Tank volume curves not supported yet!");
    }
    else {
      return self.min_level * PI * self.diameter * self.diameter / 4.0;
    }
  }
  pub fn max_volume(&self) -> f64 {
    if self.volume_curve.is_some() {
      panic!("Tank volume curves not supported yet!");
    }
    else {
      return self.max_level * PI * self.diameter * self.diameter / 4.0;
    }
  }

  pub fn new_head(&self, delta_volume: f64, current_head: f64) -> f64 {

    let mut level = current_head - self.elevation;

    if self.volume_curve.is_some() {
      panic!("Tank volume curves not supported yet!");
    }
    else {
      // linear volume curve
      let surface_area = PI * self.diameter * self.diameter / 4.0; // in ft^2
      let new_level = level + delta_volume / surface_area; // in ft

      level = new_level.clamp(self.min_level, self.max_level);
    }

    return self.elevation + level;
  }
}

impl UnitConversion for Tank {
  fn convert_units(&mut self, _flow: &FlowUnits, system: &UnitSystem, _reverse: bool) {
    // convert the initial level, min level, max level, diameter, min volume
    if system == &UnitSystem::SI {
      self.initial_level = self.initial_level * MperFT;
      self.min_level = self.min_level * MperFT;
      self.max_level = self.max_level * MperFT;
      self.diameter = self.diameter * MperFT;
      self.min_volume = self.min_volume * M3perFT3;
    }
  }
}