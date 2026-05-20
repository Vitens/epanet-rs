use serde::{Deserialize, Serialize};

use crate::model::units::{Cfs, UnitConversion};

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Demand {
    pub basedemand: Cfs,
    pub pattern: Option<Box<str>>,
    #[serde(skip)]
    pub pattern_index: Option<usize>,
    pub name: Option<Box<str>>,
}

impl UnitConversion for Demand {
    fn convert_to_standard(&mut self, options: &super::options::SimulationOptions) {
        self.basedemand /= options.flow_units.per_cfs();
    }
    fn convert_from_standard(&mut self, options: &super::options::SimulationOptions) {
        self.basedemand *= options.flow_units.per_cfs();
    }
}
