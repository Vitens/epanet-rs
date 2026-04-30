use crate::model::link::LinkStatus;
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
#[derive(Debug, Deserialize, Serialize, Clone)]
pub enum NodeAttribute {
    Demand,
    Head,
    Pressure,
}
/// Tank condition attributes
#[derive(Debug, Deserialize, Serialize, Clone)]
pub enum TankAttribute {
    Level,
    FillTime,
    DrainTime,
}
/// Link condition attributes
#[derive(Debug, Deserialize, Serialize, Clone)]
pub enum LinkAttribute {
    Flow,
    Status,
    Setting,
}
/// System condition attributes
#[derive(Debug, Deserialize, Serialize, Clone)]
pub enum SystemAttribute {
    Demand,
    Time,
    ClockTime,
}

/// Rule condition target
#[derive(Debug, Deserialize, Serialize, Clone)]
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
#[derive(Debug, Deserialize, Serialize, Clone)]
pub enum ConditionValue {
    Number(f64),
    Status(LinkStatus),
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub enum RuleConditionOperator {
    Or,
    And,
}

#[derive(strum::EnumString, strum::Display, Debug, Deserialize, Serialize, Clone)]
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

/// Condition
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct RuleCondition {
    pub operator: RuleConditionOperator,
    pub target: RuleConditionTarget,
    pub comparison: ComparisonOperator,
    pub value: ConditionValue,
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
