//! Simple `Control` rules that open/close links based on time, tank level or nodal pressure.

use crate::constants::{H_TOL, PSIperFT};
use crate::model::link::LinkStatus;
use crate::model::network::Network;
use crate::model::node::NodeType;
use crate::solver::state::SolverState;
use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize, Clone)]
pub enum ControlCondition {
    HighPressure { node_index: usize, target: f64 },
    LowPressure { node_index: usize, target: f64 },
    HighLevel { tank_index: usize, target: f64 },
    LowLevel { tank_index: usize, target: f64 },
    Time { seconds: usize },
    ClockTime { seconds: usize },
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Control {
    pub condition: ControlCondition,
    pub link_id: Box<str>,
    pub setting: Option<f64>,
    pub status: Option<LinkStatus>,
}

impl Control {
    pub fn is_active(
        &self,
        state: &SolverState,
        network: &Network,
        time: usize,
        clocktime: usize,
    ) -> bool {
        match &self.condition {
            ControlCondition::Time { seconds } => *seconds == time,
            ControlCondition::ClockTime { seconds } => *seconds == clocktime,
            ControlCondition::HighPressure { node_index, target } => {
                let node = &network.nodes[*node_index];
                let value = (state.heads[*node_index] + node.elevation) * PSIperFT; // convert head to pressure
                value - *target >= -H_TOL
            }
            ControlCondition::LowPressure { node_index, target } => {
                let node = &network.nodes[*node_index];
                let value = (state.heads[*node_index] + node.elevation) * PSIperFT; // convert head to pressure
                value - *target <= H_TOL
            }
            ControlCondition::HighLevel { tank_index, target } => {
                let node = &network.nodes[*tank_index];
                let NodeType::Tank(tank) = &node.node_type else {
                    return false;
                };
                // convert head to level
                let value = state.heads[*tank_index] - tank.elevation;
                value - *target >= -H_TOL
            }
            ControlCondition::LowLevel { tank_index, target } => {
                let node = &network.nodes[*tank_index];
                let NodeType::Tank(tank) = &node.node_type else {
                    return false;
                };
                // convert head to level
                let value = state.heads[*tank_index] - tank.elevation;
                value - *target <= H_TOL
            }
        }
    }

    pub fn activate(&self, state: &mut SolverState, network: &Network) -> bool {
        let link_index = network.link_map.get(&self.link_id).unwrap();

        if let Some(status) = self.status {
            let changed = state.statuses[*link_index] != status;
            state.statuses[*link_index] = status;
            return changed;
        }
        if let Some(setting) = self.setting {
            let changed = state.settings[*link_index] != setting;
            state.settings[*link_index] = setting;
            return changed;
        }
        false
    }
}
