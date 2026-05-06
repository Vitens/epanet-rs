//! EPANET toolkit numeric error and warning codes mapped to a strongly-typed [`ErrorCode`] enum.

use strum::FromRepr;
use thiserror::Error;

/// EPANET toolkit error/warning codes.
///
/// Values match the official EPANET 2.3 error numbering.
/// See: <http://wateranalytics.org/EPANET/group___error_codes.html>
///
/// Uses `#[repr(i32)]` so the enum is ABI-compatible with `c_int`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error, FromRepr)]
#[repr(i32)]
pub enum ErrorCode {
    #[error("Invalid Project Handle")]
    InvalidHandle = -1,

    #[error("No error")]
    Ok = 0,

    // --- System errors (1xx) ---
    #[error("Insufficient memory available")]
    InsufficientMemory = 101,
    #[error("No network data available")]
    NoNetworkData = 102,
    #[error("Hydraulic solver not opened")]
    HydraulicSolverNotOpened = 103,
    #[error("No hydraulics for water quality analysis")]
    NoHydraulicsForQuality = 104,
    #[error("Water quality solver not opened")]
    QualitySolverNotOpened = 105,
    #[error("No results saved to report on")]
    NoResultsSaved = 106,
    #[error("Hydraulics supplied from external file")]
    HydraulicsFromExternalFile = 107,
    #[error("Cannot use external file while hydraulics solver is open")]
    ExternalFileConflict = 108,
    #[error("Cannot solve network hydraulic equations")]
    CannotSolveHydraulics = 110,
    #[error("Cannot solve water quality transport equations")]
    CannotSolveQuality = 120,

    // --- Input / data errors (2xx) ---
    #[error("One or more errors in an input file")]
    InputError = 200,
    #[error("Syntax error")]
    SyntaxError = 201,
    #[error("Function call contains an illegal numeric value")]
    IllegalNumericValue = 202,
    #[error("Function call refers to an undefined node")]
    UndefinedNode = 203,
    #[error("Function call refers to an undefined link")]
    UndefinedLink = 204,
    #[error("Function call refers to an undefined time pattern")]
    UndefinedPattern = 205,
    #[error("Function call refers to an undefined curve")]
    UndefinedCurve = 206,
    #[error("Function call attempts to control a check valve pipe or a GPV valve")]
    IllegalValveControl = 207,
    #[error("Function call contains illegal PDA pressure limits")]
    IllegalPdaPressureLimits = 208,
    #[error("Function call contains an illegal node property value")]
    IllegalNodeProperty = 209,
    #[error("Function call contains an illegal link property value")]
    IllegalLinkProperty = 211,
    #[error("Function call refers to an undefined Trace Node")]
    UndefinedTraceNode = 212,
    #[error("Function call contains an invalid option value")]
    InvalidOptionValue = 213,
    #[error("Too many characters in a line of an input file")]
    LineTooLong = 214,
    #[error("Function call contains a duplicate ID label")]
    DuplicateId = 215,
    #[error("Function call refers to an undefined pump")]
    UndefinedPump = 216,
    #[error("Invalid pump energy data")]
    InvalidPumpEnergyData = 217,
    #[error("Illegal valve connection to tank node")]
    IllegalValveToTank = 219,
    #[error("Illegal valve connection to another valve")]
    IllegalValveToValve = 220,
    #[error("Mis-placed clause in rule-based control")]
    MisplacedRuleClause = 221,
    #[error("Link assigned same start and end nodes")]
    LinkSameStartEnd = 222,
    #[error("Not enough nodes in network")]
    NotEnoughNodes = 223,
    #[error("No tanks or reservoirs in network")]
    NoTanksOrReservoirs = 224,
    #[error("Invalid lower/upper levels for tank")]
    InvalidTankLevels = 225,
    #[error("No head curve or power rating for pump")]
    NoPumpCurveOrPower = 226,
    #[error("Invalid head curve for pump")]
    InvalidPumpCurve = 227,
    #[error("Nonincreasing x-values for curve")]
    NonincreasingCurveX = 230,
    #[error("No data provided for a curve")]
    NoCurveData = 231,
    #[error("No data provided for a pattern")]
    NoPatternData = 232,
    #[error("Network has unconnected nodes")]
    UnconnectedNodes = 233,
    #[error("Function call refers to nonexistent water quality source")]
    NonexistentQualitySource = 240,
    #[error("Function call refers to nonexistent control")]
    NonexistentControl = 241,
    #[error("Function call contains invalid format")]
    InvalidFormat = 250,
    #[error("Function call contains invalid parameter code")]
    InvalidParameterCode = 251,
    #[error("Function call refers to an invalid ID name")]
    InvalidIdName = 252,
    #[error("Function call refers to nonexistent demand category")]
    NonexistentDemandCategory = 253,
    #[error("Function call refers to node with no coordinates")]
    NodeNoCoordinates = 254,
    #[error("Function call refers to link with no vertices")]
    LinkNoVertices = 255,
    #[error("Function call refers to nonexistent rule")]
    NonexistentRule = 257,
    #[error("Function call refers to nonexistent rule clause")]
    NonexistentRuleClause = 258,
    #[error("Function call attempts to delete a node that still has links connected to it")]
    DeleteNodeWithLinks = 259,
    #[error("Function call attempts to delete node assigned as a Trace Node")]
    DeleteTraceNode = 260,
    #[error("Function call attempts to delete a node or link contained in a control")]
    DeleteNodeOrLinkInControl = 261,
    #[error("Function call attempts to modify network structure while a solver is open")]
    ModifyWhileSolverOpen = 262,
    #[error("Function call refers to node that is not a tank")]
    NotATank = 263,
    #[error("Function call refers to a link that is not a valve")]
    NotAValve = 264,
    #[error("An invalid section keyword was detected in an input file")]
    InvalidSectionKeyword = 299,

    // --- File errors (3xx) ---
    #[error("Identical file names used for different types of files")]
    IdenticalFileNames = 301,
    #[error("Cannot open input file")]
    CannotOpenInputFile = 302,
    #[error("Cannot open report file")]
    CannotOpenReportFile = 303,
    #[error("Cannot open output file")]
    CannotOpenOutputFile = 304,
    #[error("Cannot open hydraulics file")]
    CannotOpenHydraulicsFile = 305,
    #[error("Hydraulics file does not match network data")]
    HydraulicsFileMismatch = 306,
    #[error("Cannot read hydraulics file")]
    CannotReadHydraulicsFile = 307,
    #[error("Cannot save results to binary file")]
    CannotSaveBinaryResults = 308,
    #[error("Cannot save results to report file")]
    CannotSaveReportResults = 309,
}
