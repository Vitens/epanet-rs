use crate::model::units::{Cfs, FlowUnits, UnitSystem, UnitConversion};
use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize)]
pub struct Junction {
  pub basedemand: Cfs,
  pub pattern: Option<Box<str>>,
}

impl UnitConversion for Junction {
  fn convert_units(&mut self, flow: &FlowUnits, _system: &UnitSystem, reverse: bool) {
    // convert from CFS to the given units
    if reverse {
      self.basedemand = self.basedemand * flow.per_cfs();
    } 
    // convert from the given units to CFS
    else {
      self.basedemand = self.basedemand / flow.per_cfs();
    }
  }
}