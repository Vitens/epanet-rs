//! Generic `Node` type with a `NodeType` variant for junctions, tanks and reservoirs.

use crate::model::junction::Junction;
use crate::model::options::SimulationOptions;
use crate::model::reservoir::Reservoir;
use crate::model::tank::Tank;
use crate::model::units::UnitConversion;

use serde::{Deserialize, Serialize};

/// Node struct
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Node {
    pub id: Box<str>,
    pub node_type: NodeType,
    pub elevation: f64,
    pub coordinates: Option<(f64, f64)>,
}

/// Node types
#[derive(Debug, Deserialize, Serialize, Clone)]
pub enum NodeType {
    Reservoir(Reservoir),
    Tank(Tank),
    Junction(Junction),
}

// helper methods for nodes to check if they are fixed head
// and to get the head pattern
impl Node {
    pub fn is_fixed(&self) -> bool {
        matches!(self.node_type, NodeType::Reservoir(_) | NodeType::Tank(_))
    }
    pub fn head_pattern(&self) -> Option<&str> {
        if let NodeType::Reservoir(reservoir) = &self.node_type {
            if reservoir.head_pattern.is_some() {
                return Some(reservoir.head_pattern.as_ref().unwrap());
            }
        }
        None
    }
    pub fn initial_head(&self) -> f64 {
        if self.is_fixed() {
            if let NodeType::Tank(tank) = &self.node_type {
                return self.elevation + tank.initial_level;
            }
            return self.elevation;
        }
        return 0.0;
    }
}

impl UnitConversion for Node {
    fn convert_to_standard(&mut self, options: &SimulationOptions) {
        // convert elevation to Feet
        self.elevation = self.elevation / options.unit_system.per_feet();

        match &mut self.node_type {
            NodeType::Reservoir(_reservoir) => (),
            NodeType::Tank(tank) => tank.convert_to_standard(options),
            NodeType::Junction(junction) => junction.convert_to_standard(options),
        }
    }

    fn convert_from_standard(&mut self, options: &SimulationOptions) {
        // convert elevation from Feet to the given unit system
        self.elevation = self.elevation * options.unit_system.per_feet();

        match &mut self.node_type {
            NodeType::Reservoir(_reservoir) => (),
            NodeType::Tank(tank) => tank.convert_from_standard(options),
            NodeType::Junction(junction) => junction.convert_from_standard(options),
        }
    }
}
