//! EPANET 2.3 enumerated types.
//!
//! Values match the official C header `epanet2_enums.h`.
//! See: <http://wateranalytics.org/EPANET/group___enumerations.html>

use strum::FromRepr;

pub const MISSING_VALUE: f64 = -1.0e10;

/// Types of network objects.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum ObjectType {
    Node = 0,
    Link = 1,
    TimePat = 2,
    Curve = 3,
    Control = 4,
    Rule = 5,
}

/// Types of objects to count.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum CountType {
    NodeCount = 0,
    TankCount = 1,
    LinkCount = 2,
    PatCount = 3,
    CurveCount = 4,
    ControlCount = 5,
    RuleCount = 6,
}

/// Node types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum NodeType {
    Junction = 0,
    Reservoir = 1,
    Tank = 2,
}

/// Link types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum LinkType {
    CVPipe = 0,
    Pipe = 1,
    Pump = 2,
    PRV = 3,
    PSV = 4,
    PBV = 5,
    FCV = 6,
    TCV = 7,
    GPV = 8,
    PCV = 9,
}

/// Types of pump curves.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum PumpType {
    ConstHp = 0,
    PowerFunc = 1,
    Custom = 2,
    NoCurve = 3,
}

/// Pump states.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum PumpStateType {
    XHead = 0,
    Closed = 2,
    Open = 3,
    XFlow = 5,
}

/// Types of data curves.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum CurveType {
    Volume = 0,
    Pump = 1,
    Effic = 2,
    HLoss = 3,
    Generic = 4,
    Valve = 5,
}

/// Types of water quality analyses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum QualityType {
    None = 0,
    Chem = 1,
    Age = 2,
    Trace = 3,
}

/// Water quality source types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum SourceType {
    Concen = 0,
    Mass = 1,
    SetPoint = 2,
    FlowPaced = 3,
}

/// Simple control types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum ControlType {
    LowLevel = 0,
    HiLevel = 1,
    Timer = 2,
    TimeOfDay = 3,
}

/// Head loss formulas.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum HeadLossType {
    HW = 0,
    DW = 1,
    CM = 2,
}

/// Node properties.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum NodeProperty {
    Elevation = 0,
    BaseDemand = 1,
    Pattern = 2,
    Emitter = 3,
    InitQual = 4,
    SourceQual = 5,
    SourcePat = 6,
    SourceType = 7,
    TankLevel = 8,
    Demand = 9,
    Head = 10,
    Pressure = 11,
    Quality = 12,
    SourceMass = 13,
    InitVolume = 14,
    MixModel = 15,
    MixZoneVol = 16,
    TankDiam = 17,
    MinVolume = 18,
    VolCurve = 19,
    MinLevel = 20,
    MaxLevel = 21,
    MixFraction = 22,
    TankKBulk = 23,
    TankVolume = 24,
    MaxVolume = 25,
    CanOverflow = 26,
    DemandDeficit = 27,
    NodeInControl = 28,
    EmitterFlow = 29,
    LeakageFlow = 30,
    DemandFlow = 31,
    FullDemand = 32,
}

/// Link properties.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum LinkProperty {
    Diameter = 0,
    Length = 1,
    Roughness = 2,
    MinorLoss = 3,
    InitStatus = 4,
    InitSetting = 5,
    KBulk = 6,
    KWall = 7,
    Flow = 8,
    Velocity = 9,
    HeadLoss = 10,
    Status = 11,
    Setting = 12,
    Energy = 13,
    LinkQual = 14,
    LinkPattern = 15,
    PumpState = 16,
    PumpEffic = 17,
    PumpPower = 18,
    PumpHCurve = 19,
    PumpECurve = 20,
    PumpECost = 21,
    PumpEPat = 22,
    LinkInControl = 23,
    GpvCurve = 24,
    PcvCurve = 25,
    LeakArea = 26,
    LeakExpan = 27,
    LinkLeakage = 28,
    ValveType = 29,
}

/// Link status.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum LinkStatusType {
    Closed = 0,
    Open = 1,
}

/// Time parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum TimeParameter {
    Duration = 0,
    HydStep = 1,
    QualStep = 2,
    PatternStep = 3,
    PatternStart = 4,
    ReportStep = 5,
    ReportStart = 6,
    RuleStep = 7,
    Statistic = 8,
    Periods = 9,
    StartTime = 10,
    HTime = 11,
    QTime = 12,
    HaltFlag = 13,
    NextEvent = 14,
    NextEventTank = 15,
}

/// Time step events.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum TimestepEvent {
    Report = 0,
    Hyd = 1,
    Wq = 2,
    TankEvent = 3,
    ControlEvent = 4,
}

/// Simulation options.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum SimOption {
    Trials = 0,
    Accuracy = 1,
    Tolerance = 2,
    EmitExpon = 3,
    DemandMult = 4,
    HeadError = 5,
    FlowChange = 6,
    HeadLossForm = 7,
    GlobalEffic = 8,
    GlobalPrice = 9,
    GlobalPattern = 10,
    DemandCharge = 11,
    SpGravity = 12,
    SpViscos = 13,
    Unbalanced = 14,
    CheckFreq = 15,
    MaxCheck = 16,
    DampLimit = 17,
    SpDiffus = 18,
    BulkOrder = 19,
    WallOrder = 20,
    TankOrder = 21,
    ConcenLimit = 22,
    DemandPattern = 23,
    EmitBackflow = 24,
    PressUnits = 25,
    StatusReport = 26,
}

/// Flow units.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum FlowUnits {
    Cfs = 0,
    Gpm = 1,
    Mgd = 2,
    Imgd = 3,
    Afd = 4,
    Lps = 5,
    Lpm = 6,
    Mld = 7,
    Cmh = 8,
    Cmd = 9,
    Cms = 10,
}

/// Pressure units.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum PressUnits {
    Psi = 0,
    Kpa = 1,
    Meters = 2,
    Bar = 3,
    Feet = 4,
}

/// Demand models.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum DemandModel {
    Dda = 0,
    Pda = 1,
}

/// Tank mixing models.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum MixingModel {
    Mix1 = 0,
    Mix2 = 1,
    Fifo = 2,
    Lifo = 3,
}

/// Reporting statistic choices.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum StatisticType {
    Series = 0,
    Average = 1,
    Minimum = 2,
    Maximum = 3,
    Range = 4,
}

/// Hydraulic initialization options.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum InitHydOption {
    NoSave = 0,
    Save = 1,
    InitFlow = 10,
    SaveAndInit = 11,
}

/// Deletion action codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum ActionCodeType {
    Unconditional = 0,
    Conditional = 1,
}

/// Analysis convergence statistics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum AnalysisStatistic {
    Iterations = 0,
    RelativeError = 1,
    MaxHeadError = 2,
    MaxFlowChange = 3,
    MassBalance = 4,
    DeficientNodes = 5,
    DemandReduction = 6,
    LeakageLoss = 7,
}

/// Status reporting levels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum StatusReport {
    NoReport = 0,
    NormalReport = 1,
    FullReport = 2,
}

/// Network objects used in rule-based controls.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum RuleObject {
    Node = 6,
    Link = 7,
    System = 8,
}

/// Object variables used in rule-based controls.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum RuleVariable {
    Demand = 0,
    Head = 1,
    Grade = 2,
    Level = 3,
    Pressure = 4,
    Flow = 5,
    Status = 6,
    Setting = 7,
    Power = 8,
    Time = 9,
    ClockTime = 10,
    FillTime = 11,
    DrainTime = 12,
}

/// Comparison operators used in rule-based controls.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum RuleOperator {
    Eq = 0,
    Ne = 1,
    Le = 2,
    Ge = 3,
    Lt = 4,
    Gt = 5,
    Is = 6,
    Not = 7,
    Below = 8,
    Above = 9,
}

/// Link status codes used in rule-based controls.
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromRepr)]
#[repr(i32)]
pub enum RuleStatus {
    IsOpen = 1,
    IsClosed = 2,
    IsActive = 3,
}
