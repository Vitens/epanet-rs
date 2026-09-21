//! Python enums mirroring the corresponding native Rust enums 1:1.
//!
//! These are thin wrapper types (the core `epanet-rs` crate does not depend on
//! `pyo3`), but every variant name and ordering matches the Rust source
//! exactly so callers see the same vocabulary as the Rust API. Several
//! variant names are acronyms (`CFS`, `PSI`, `PRV`, ...), matching the
//! upstream Rust identifiers verbatim.
#![allow(clippy::upper_case_acronyms)]

use pyo3::prelude::*;

use epanet::model::link::LinkStatus as RsLinkStatus;
use epanet::model::options::{DemandModel as RsDemandModel, HeadlossFormula as RsHeadlossFormula};
use epanet::model::units::{
    FlowUnits as RsFlowUnits, PressureUnits as RsPressureUnits, UnitSystem as RsUnitSystem,
};
use epanet::model::valve::ValveType as RsValveType;

/// Discriminant for [`crate::native::types::Node::node_type`]. The associated
/// per-type data is available via `Node.junction` / `Node.tank` / `Node.reservoir`.
#[pyclass(eq, eq_int, module = "epanet_rs.native")]
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum NodeType {
    Junction,
    Reservoir,
    Tank,
}

/// Discriminant for [`crate::native::types::Link::link_type`]. The associated
/// per-type data is available via `Link.pipe` / `Link.pump` / `Link.valve`.
#[pyclass(eq, eq_int, module = "epanet_rs.native")]
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum LinkType {
    Pipe,
    Pump,
    Valve,
}

/// Mirrors `epanet_rs::model::valve::ValveType`.
#[pyclass(eq, eq_int, module = "epanet_rs.native")]
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum ValveType {
    PRV,
    PSV,
    PBV,
    FCV,
    TCV,
    PCV,
    GPV,
}

impl From<RsValveType> for ValveType {
    fn from(value: RsValveType) -> Self {
        match value {
            RsValveType::PRV => ValveType::PRV,
            RsValveType::PSV => ValveType::PSV,
            RsValveType::PBV => ValveType::PBV,
            RsValveType::FCV => ValveType::FCV,
            RsValveType::TCV => ValveType::TCV,
            RsValveType::PCV => ValveType::PCV,
            RsValveType::GPV => ValveType::GPV,
        }
    }
}

impl From<ValveType> for RsValveType {
    fn from(value: ValveType) -> Self {
        match value {
            ValveType::PRV => RsValveType::PRV,
            ValveType::PSV => RsValveType::PSV,
            ValveType::PBV => RsValveType::PBV,
            ValveType::FCV => RsValveType::FCV,
            ValveType::TCV => RsValveType::TCV,
            ValveType::PCV => RsValveType::PCV,
            ValveType::GPV => RsValveType::GPV,
        }
    }
}

/// Mirrors `epanet_rs::model::link::LinkStatus`.
#[pyclass(eq, eq_int, module = "epanet_rs.native")]
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum LinkStatus {
    Xhead,
    TempClosed,
    Closed,
    Open,
    Active,
    Xflow,
    XFCV,
    XPressure,
    FixedOpen,
    FixedClosed,
}

impl From<RsLinkStatus> for LinkStatus {
    fn from(value: RsLinkStatus) -> Self {
        match value {
            RsLinkStatus::Xhead => LinkStatus::Xhead,
            RsLinkStatus::TempClosed => LinkStatus::TempClosed,
            RsLinkStatus::Closed => LinkStatus::Closed,
            RsLinkStatus::Open => LinkStatus::Open,
            RsLinkStatus::Active => LinkStatus::Active,
            RsLinkStatus::Xflow => LinkStatus::Xflow,
            RsLinkStatus::XFCV => LinkStatus::XFCV,
            RsLinkStatus::XPressure => LinkStatus::XPressure,
            RsLinkStatus::FixedOpen => LinkStatus::FixedOpen,
            RsLinkStatus::FixedClosed => LinkStatus::FixedClosed,
        }
    }
}

impl From<LinkStatus> for RsLinkStatus {
    fn from(value: LinkStatus) -> Self {
        match value {
            LinkStatus::Xhead => RsLinkStatus::Xhead,
            LinkStatus::TempClosed => RsLinkStatus::TempClosed,
            LinkStatus::Closed => RsLinkStatus::Closed,
            LinkStatus::Open => RsLinkStatus::Open,
            LinkStatus::Active => RsLinkStatus::Active,
            LinkStatus::Xflow => RsLinkStatus::Xflow,
            LinkStatus::XFCV => RsLinkStatus::XFCV,
            LinkStatus::XPressure => RsLinkStatus::XPressure,
            LinkStatus::FixedOpen => RsLinkStatus::FixedOpen,
            LinkStatus::FixedClosed => RsLinkStatus::FixedClosed,
        }
    }
}

/// Mirrors `epanet_rs::model::options::HeadlossFormula`.
#[pyclass(eq, eq_int, module = "epanet_rs.native")]
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum HeadlossFormula {
    HazenWilliams,
    DarcyWeisbach,
    ChezyManning,
}

impl From<RsHeadlossFormula> for HeadlossFormula {
    fn from(value: RsHeadlossFormula) -> Self {
        match value {
            RsHeadlossFormula::HazenWilliams => HeadlossFormula::HazenWilliams,
            RsHeadlossFormula::DarcyWeisbach => HeadlossFormula::DarcyWeisbach,
            RsHeadlossFormula::ChezyManning => HeadlossFormula::ChezyManning,
        }
    }
}

impl From<HeadlossFormula> for RsHeadlossFormula {
    fn from(value: HeadlossFormula) -> Self {
        match value {
            HeadlossFormula::HazenWilliams => RsHeadlossFormula::HazenWilliams,
            HeadlossFormula::DarcyWeisbach => RsHeadlossFormula::DarcyWeisbach,
            HeadlossFormula::ChezyManning => RsHeadlossFormula::ChezyManning,
        }
    }
}

/// Mirrors `epanet_rs::model::options::DemandModel`.
#[pyclass(eq, eq_int, module = "epanet_rs.native")]
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum DemandModel {
    PDA,
    DDA,
}

impl From<RsDemandModel> for DemandModel {
    fn from(value: RsDemandModel) -> Self {
        match value {
            RsDemandModel::PDA => DemandModel::PDA,
            RsDemandModel::DDA => DemandModel::DDA,
        }
    }
}

/// Mirrors `epanet_rs::model::units::FlowUnits` (the full set, including `CMS`
/// which the EPANET-2.3-compatible `Project` mirror does not expose).
#[pyclass(eq, eq_int, module = "epanet_rs.native")]
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum FlowUnits {
    CFS,
    GPM,
    MGD,
    IMGD,
    AFD,
    LPS,
    LPM,
    MLD,
    CMS,
    CMH,
    CMD,
}

impl From<RsFlowUnits> for FlowUnits {
    fn from(value: RsFlowUnits) -> Self {
        match value {
            RsFlowUnits::CFS => FlowUnits::CFS,
            RsFlowUnits::GPM => FlowUnits::GPM,
            RsFlowUnits::MGD => FlowUnits::MGD,
            RsFlowUnits::IMGD => FlowUnits::IMGD,
            RsFlowUnits::AFD => FlowUnits::AFD,
            RsFlowUnits::LPS => FlowUnits::LPS,
            RsFlowUnits::LPM => FlowUnits::LPM,
            RsFlowUnits::MLD => FlowUnits::MLD,
            RsFlowUnits::CMS => FlowUnits::CMS,
            RsFlowUnits::CMH => FlowUnits::CMH,
            RsFlowUnits::CMD => FlowUnits::CMD,
        }
    }
}

impl From<FlowUnits> for RsFlowUnits {
    fn from(value: FlowUnits) -> Self {
        match value {
            FlowUnits::CFS => RsFlowUnits::CFS,
            FlowUnits::GPM => RsFlowUnits::GPM,
            FlowUnits::MGD => RsFlowUnits::MGD,
            FlowUnits::IMGD => RsFlowUnits::IMGD,
            FlowUnits::AFD => RsFlowUnits::AFD,
            FlowUnits::LPS => RsFlowUnits::LPS,
            FlowUnits::LPM => RsFlowUnits::LPM,
            FlowUnits::MLD => RsFlowUnits::MLD,
            FlowUnits::CMS => RsFlowUnits::CMS,
            FlowUnits::CMH => RsFlowUnits::CMH,
            FlowUnits::CMD => RsFlowUnits::CMD,
        }
    }
}

/// Mirrors `epanet_rs::model::units::UnitSystem`.
#[pyclass(eq, eq_int, module = "epanet_rs.native")]
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum UnitSystem {
    US,
    SI,
}

impl From<RsUnitSystem> for UnitSystem {
    fn from(value: RsUnitSystem) -> Self {
        match value {
            RsUnitSystem::US => UnitSystem::US,
            RsUnitSystem::SI => UnitSystem::SI,
        }
    }
}

/// Mirrors `epanet_rs::model::units::PressureUnits`.
#[pyclass(eq, eq_int, module = "epanet_rs.native")]
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum PressureUnits {
    PSI,
    KPA,
    METERS,
    FEET,
    BAR,
}

impl From<RsPressureUnits> for PressureUnits {
    fn from(value: RsPressureUnits) -> Self {
        match value {
            RsPressureUnits::PSI => PressureUnits::PSI,
            RsPressureUnits::KPA => PressureUnits::KPA,
            RsPressureUnits::METERS => PressureUnits::METERS,
            RsPressureUnits::FEET => PressureUnits::FEET,
            RsPressureUnits::BAR => PressureUnits::BAR,
        }
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<NodeType>()?;
    m.add_class::<LinkType>()?;
    m.add_class::<ValveType>()?;
    m.add_class::<LinkStatus>()?;
    m.add_class::<HeadlossFormula>()?;
    m.add_class::<DemandModel>()?;
    m.add_class::<FlowUnits>()?;
    m.add_class::<UnitSystem>()?;
    m.add_class::<PressureUnits>()?;
    Ok(())
}
