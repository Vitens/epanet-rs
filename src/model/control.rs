//! Simple `Control` rules that open/close links based on time, tank level or nodal pressure.

use crate::constants::{H_TOL, L_TOL, PSIperFT};
use crate::model::link::{Link, LinkStatus, LinkType, LinkTrait};
use crate::model::network::Network;
use crate::model::node::NodeType;
use crate::model::options::SimulationOptions;
use crate::model::units::{Ft, UnitConversion, UnitSystem};
use crate::model::valve::ValveType;
use crate::solver::state::SolverState;
use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize, Clone)]
pub enum ControlCondition {
    /// Pressure target stored as feet of head above node elevation.
    HighPressure { node_index: usize, target: Ft },
    LowPressure { node_index: usize, target: Ft },
    HighLevel { tank_index: usize, target: Ft },
    LowLevel { tank_index: usize, target: Ft },
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

impl UnitConversion for ControlCondition {
    fn convert_to_standard(&mut self, options: &SimulationOptions) {
        match self {
            Self::HighPressure { target, .. } | Self::LowPressure { target, .. } => {
                *target /= options.pressure_units.per_feet();
            }
            Self::HighLevel { target, .. } | Self::LowLevel { target, .. } => {
                *target /= options.unit_system.per_feet();
            }
            Self::Time { .. } | Self::ClockTime { .. } => {}
        }
    }

    fn convert_from_standard(&mut self, options: &SimulationOptions) {
        match self {
            Self::HighPressure { target, .. } | Self::LowPressure { target, .. } => {
                *target *= options.pressure_units.per_feet();
            }
            Self::HighLevel { target, .. } | Self::LowLevel { target, .. } => {
                *target *= options.unit_system.per_feet();
            }
            Self::Time { .. } | Self::ClockTime { .. } => {}
        }
    }
}

impl UnitConversion for Control {
    fn convert_to_standard(&mut self, options: &SimulationOptions) {
        self.condition.convert_to_standard(options);
    }

    fn convert_from_standard(&mut self, options: &SimulationOptions) {
        self.condition.convert_from_standard(options);
    }
}

impl Control {
    /// Convert a link setting from user units to internal units.
    pub fn convert_setting_to_standard(&mut self, link: &Link, options: &SimulationOptions) {
        let Some(setting) = self.setting else {
            return;
        };

        if let LinkType::Valve(valve) = &link.link_type {
            self.setting = Some(match valve.valve_type {
                ValveType::PRV | ValveType::PSV | ValveType::PBV => {
                    let mut s = setting;
                    if options.unit_system == UnitSystem::US {
                        s /= PSIperFT;
                    }
                    s / options.unit_system.per_feet()
                }
                ValveType::FCV => setting / options.flow_units.per_cfs(),
                _ => setting,
            });
        }
    }

    /// Convert a link setting from internal units back to user units.
    pub fn convert_setting_from_standard(&mut self, link: &Link, options: &SimulationOptions) {
        let Some(setting) = self.setting else {
            return;
        };

        if let LinkType::Valve(valve) = &link.link_type {
            self.setting = Some(match valve.valve_type {
                ValveType::PRV | ValveType::PSV | ValveType::PBV => {
                    let mut s = setting * options.unit_system.per_feet();
                    if options.unit_system == UnitSystem::US {
                        s *= PSIperFT;
                    }
                    s
                }
                ValveType::FCV => setting * options.flow_units.per_cfs(),
                _ => setting,
            });
        }
    }

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
                let value = state.heads[*node_index] - node.elevation;
                value - *target >= -H_TOL
            }
            ControlCondition::LowPressure { node_index, target } => {
                let node = &network.nodes[*node_index];
                let value = state.heads[*node_index] - node.elevation;
                value - *target <= H_TOL
            }
            ControlCondition::HighLevel { tank_index, target } => {
                let node = &network.nodes[*tank_index];
                let NodeType::Tank(tank) = &node.node_type else {
                    return false;
                };
                let value = state.heads[*tank_index] - tank.elevation;
                value - *target >= -L_TOL
            }
            ControlCondition::LowLevel { tank_index, target } => {
                let node = &network.nodes[*tank_index];
                let NodeType::Tank(tank) = &node.node_type else {
                    return false;
                };
                let value = state.heads[*tank_index] - tank.elevation;
                value - *target <= L_TOL
            }
        }
    }

    pub fn activate(&self, state: &mut SolverState, network: &Network) -> bool {
        let link_index = network.link_map.get(&self.link_id).unwrap();
        let link = &network.links[*link_index];

        if let Some(status) = self.status {
            let changed = state.statuses[*link_index] != status;
            state.statuses[*link_index] = status;
            // reset the flow of the link to the initial flow
            state.flows[*link_index] = link.initial_flow();
            return changed;
        }
        if let Some(setting) = self.setting {
            let changed = state.settings[*link_index] != setting;
            state.settings[*link_index] = setting;
            // reset the flow of the link to the initial flow
            state.flows[*link_index] = link.initial_flow();
            return changed;
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::PSIperFT;
    use crate::model::options::{HeadlossFormula, SimulationOptions};
    use crate::model::units::FlowUnits;

    fn us_options() -> SimulationOptions {
        SimulationOptions::new(FlowUnits::CFS, HeadlossFormula::HazenWilliams)
    }

    fn si_options() -> SimulationOptions {
        SimulationOptions::new(FlowUnits::LPS, HeadlossFormula::HazenWilliams)
    }

    #[test]
    fn pressure_target_converts_to_feet_of_head() {
        let mut condition = ControlCondition::LowPressure {
            node_index: 0,
            target: 20.0,
        };
        let options = us_options();
        condition.convert_to_standard(&options);

        assert!((condition.target() - 20.0 / PSIperFT).abs() < 1e-10);
    }

    #[test]
    fn pressure_target_converts_from_feet_of_head() {
        let mut condition = ControlCondition::HighPressure {
            node_index: 0,
            target: 20.0 / PSIperFT,
        };
        let options = us_options();
        condition.convert_from_standard(&options);

        assert!((condition.target() - 20.0).abs() < 1e-10);
    }

    #[test]
    fn si_pressure_target_converts_to_feet_of_head() {
        let mut condition = ControlCondition::LowPressure {
            node_index: 0,
            target: 10.0,
        };
        let options = si_options();
        condition.convert_to_standard(&options);

        assert!((condition.target() - 10.0 / options.unit_system.per_feet()).abs() < 1e-10);
    }

    #[test]
    fn level_target_converts_to_feet() {
        let mut condition = ControlCondition::LowLevel {
            tank_index: 0,
            target: 5.0,
        };
        let options = si_options();
        condition.convert_to_standard(&options);

        assert!((condition.target() - 5.0 / options.unit_system.per_feet()).abs() < 1e-10);
    }

    impl ControlCondition {
        fn target(&self) -> f64 {
            match self {
                Self::HighPressure { target, .. }
                | Self::LowPressure { target, .. }
                | Self::HighLevel { target, .. }
                | Self::LowLevel { target, .. } => *target,
                _ => panic!("expected target condition"),
            }
        }
    }
}
