use crate::constants::{H_TOL, PSIperFT};
use crate::model::link::LinkStatus;
use crate::model::network::Network;
use crate::model::node::NodeType;
use crate::solver::state::SolverState;
use serde::{Deserialize, Serialize};

use strum;

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Rule {
    /// Rule ID
    pub id: Box<str>,
    // List of conditions
    pub conditions: Vec<RuleCondition>,
    // Actions
    pub actions: Vec<RuleAction>,
    // Rule priority
    pub priority: Option<usize>,
}

/// Node condition attributes
#[derive(Debug, Deserialize, Serialize, Clone, PartialEq, Eq)]
pub enum NodeAttribute {
    Demand,
    Head,
    Pressure,
}
/// Tank condition attributes
#[derive(Debug, Deserialize, Serialize, Clone, PartialEq, Eq)]
pub enum TankAttribute {
    Level,
    FillTime,
    DrainTime,
}
/// Link condition attributes
#[derive(Debug, Deserialize, Serialize, Clone, PartialEq, Eq)]
pub enum LinkAttribute {
    Flow,
    Status,
    Setting,
}
/// System condition attributes
#[derive(Debug, Deserialize, Serialize, Clone, PartialEq, Eq)]
pub enum SystemAttribute {
    Demand,
    Time,
    ClockTime,
}

/// Rule condition target
#[derive(Debug, Deserialize, Serialize, Clone, PartialEq, Eq)]
pub enum RuleConditionTarget {
    Node {
        id: Box<str>,
        attribute: NodeAttribute,
    },
    Tank {
        id: Box<str>,
        attribute: TankAttribute,
    },
    Link {
        id: Box<str>,
        attribute: LinkAttribute,
    },
    System {
        attribute: SystemAttribute,
    },
}
/// Rule condition value
#[derive(Debug, Deserialize, Serialize, Clone, PartialEq)]
pub enum ConditionValue {
    Number(f64),
    Status(LinkStatus),
    Time(usize),
}

#[derive(Debug, Deserialize, Serialize, Clone, PartialEq, Eq)]
pub enum RuleConditionOperator {
    Or,
    And,
}

#[derive(strum::EnumString, strum::Display, Debug, Deserialize, Serialize, Clone, PartialEq, Eq)]
pub enum ComparisonOperator {
    #[strum(serialize = "IS", serialize = "=")]
    Eq,
    #[strum(serialize = "NOT", serialize = "<>")]
    Ne,
    #[strum(serialize = "ABOVE", serialize = ">")]
    Gt,
    #[strum(serialize = "BELOW", serialize = "<")]
    Lt,
    #[strum(serialize = ">=")]
    Ge,
    #[strum(serialize = "<=")]
    Le,
}

impl ComparisonOperator {
    /// Compare two floating point values with a small tolerance to absorb numerical noise.
    fn compare_f64(&self, lhs: f64, rhs: f64) -> bool {
        let diff = lhs - rhs;
        match self {
            ComparisonOperator::Eq => diff.abs() <= H_TOL,
            ComparisonOperator::Ne => diff.abs() > H_TOL,
            ComparisonOperator::Gt => diff > H_TOL,
            ComparisonOperator::Lt => diff < -H_TOL,
            ComparisonOperator::Ge => diff >= -H_TOL,
            ComparisonOperator::Le => diff <= H_TOL,
        }
    }

    /// Compare two link statuses; only equality / inequality are meaningful.
    fn compare_status(&self, lhs: LinkStatus, rhs: LinkStatus) -> bool {
        match self {
            ComparisonOperator::Eq => lhs == rhs,
            ComparisonOperator::Ne => lhs != rhs,
            _ => false,
        }
    }
}

/// Condition
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct RuleCondition {
    pub operator: RuleConditionOperator,
    pub target: RuleConditionTarget,
    pub comparison: ComparisonOperator,
    pub value: ConditionValue,
}

impl RuleCondition {
    /// Evaluate whether the condition currently holds against the given solver state.
    /// `time` is the elapsed simulation time in seconds and `clocktime` is the time of day in seconds.
    pub fn is_active(
        &self,
        state: &SolverState,
        network: &Network,
        time: usize,
        clocktime: usize,
    ) -> bool {
        match &self.target {
            RuleConditionTarget::System { attribute } => {
                let value = match attribute {
                    SystemAttribute::Demand => state.demands.iter().sum::<f64>(),
                    SystemAttribute::Time => time as f64,
                    SystemAttribute::ClockTime => clocktime as f64,
                };
                let target = match self.value {
                    ConditionValue::Number(v) => v,
                    ConditionValue::Time(t) => t as f64,
                    ConditionValue::Status(_) => return false,
                };
                self.comparison.compare_f64(value, target)
            }
            RuleConditionTarget::Node { id, attribute } => {
                let Some(&idx) = network.node_map.get(id) else {
                    return false;
                };
                let node = &network.nodes[idx];
                let value = match attribute {
                    NodeAttribute::Demand => state.demands[idx],
                    NodeAttribute::Head => state.heads[idx],
                    NodeAttribute::Pressure => (state.heads[idx] - node.elevation) * PSIperFT,
                };
                let target = match self.value {
                    ConditionValue::Number(v) => v,
                    _ => return false,
                };
                self.comparison.compare_f64(value, target)
            }
            RuleConditionTarget::Tank { id, attribute } => {
                let Some(&idx) = network.node_map.get(id) else {
                    return false;
                };
                let node = &network.nodes[idx];
                let NodeType::Tank(tank) = &node.node_type else {
                    return false;
                };
                let level = state.heads[idx] - tank.elevation;
                let demand = state.demands[idx];
                // FILLTIME / DRAINTIME are expressed in hours by EPANET convention.
                let value = match attribute {
                    TankAttribute::Level => level,
                    TankAttribute::FillTime => {
                        tank.time_to_reach_level(level, tank.max_level, demand) as f64 / 3600.0
                    }
                    TankAttribute::DrainTime => {
                        tank.time_to_reach_level(level, tank.min_level, demand) as f64 / 3600.0
                    }
                };
                let target = match self.value {
                    ConditionValue::Number(v) => v,
                    _ => return false,
                };
                self.comparison.compare_f64(value, target)
            }
            RuleConditionTarget::Link { id, attribute } => {
                let Some(&idx) = network.link_map.get(id) else {
                    return false;
                };
                match attribute {
                    LinkAttribute::Flow => {
                        let target = match self.value {
                            ConditionValue::Number(v) => v,
                            _ => return false,
                        };
                        self.comparison.compare_f64(state.flows[idx].abs(), target)
                    }
                    LinkAttribute::Setting => {
                        let target = match self.value {
                            ConditionValue::Number(v) => v,
                            _ => return false,
                        };
                        self.comparison.compare_f64(state.settings[idx], target)
                    }
                    LinkAttribute::Status => {
                        let target = match self.value {
                            ConditionValue::Status(s) => s,
                            _ => return false,
                        };
                        self.comparison.compare_status(state.statuses[idx], target)
                    }
                }
            }
        }
    }
}

/// Action struct
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct RuleAction {
    /// Id of the link targeted by action
    pub link_id: Box<str>,
    /// Setting to apply (if set)
    pub setting: Option<f64>,
    /// Status to apply (if set)
    pub status: Option<LinkStatus>,
    /// Rule active by default (part of ELSE clause)
    pub default_active: bool,
}

