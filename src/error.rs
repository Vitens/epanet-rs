//! Error types returned by the input parser and the hydraulic solver.

use std::str::FromStr;
use thiserror::Error;

#[derive(Debug, Error)]
/// Errors related to the solver.
pub enum SolverError {
    /// The solver has not been initialized
    #[error("Solver not initialized")]
    NotInitialized,
    /// The solver has encountered a singular matrix, check connectivity at node '{node_id}'
    #[error("Singular matrix: check connectivity at node '{node_id}'")]
    SingularMatrix { node_id: Box<str> },
    /// The solver has reached the maximum number of iterations
    #[error("Maximum number of iterations reached: {max_trials}")]
    MaxIterations { max_trials: usize },
    /// The network topology has changed since the solver was created
    #[error(
        "Network topology has changed since the solver was created, please recreate the solver"
    )]
    StaleTopology,
    /// The network properties have changed since the solver state was created
    #[error(
        "Network properties have changed since the solver state was created, please update the solver state"
    )]
    StaleProperties,
    /// Failed to compute Cholesky factorization
    #[error("Failed to compute Cholesky factorization: {0}")]
    Factorization(String),
}

#[derive(Debug, Error)]
/// Errors related to input file handling.
pub enum InputError {
    /// Cannot open file
    #[error("Cannot open file '{path}': {source}")]
    FileOpen {
        path: String,
        source: std::io::Error,
    },

    /// General IO error while handling file
    #[error("IO error: {0}")]
    FileRead(#[from] std::io::Error),

    /// Node wit ID already exists
    #[error("Node {node_id} already exists")]
    NodeExists { node_id: Box<str> },

    /// Link with ID already exists
    #[error("Link {link_id} already exists")]
    LinkExists { link_id: Box<str> },

    /// Node with ID not found
    #[error("Node {node_id} not found")]
    NodeNotFound { node_id: Box<str> },

    /// Node is not a junction, raised when attempting to update a node that is not a junction
    #[error("Node {node_id} is not a junction")]
    NodeNotAJunction { node_id: Box<str> },

    /// Node is not a tank, raised when attempting to update a node that is not a tank
    #[error("Node {node_id} is not a tank")]
    NodeNotATank { node_id: Box<str> },

    /// Node is not a reservoir, raised when attempting to update a node that is not a reservoir
    #[error("Node {node_id} is not a reservoir")]
    NodeNotAReservoir { node_id: Box<str> },

    /// Link with ID not found
    #[error("Link {link_id} not found")]
    LinkNotFound { link_id: Box<str> },

    /// Link is not a pipe, raised when attempting to update a link that is not a pipe
    #[error("Link {link_id} is not a pipe")]
    LinkNotAPipe { link_id: Box<str> },

    /// Link is not a pump, raised when attempting to update a link that is not a pump
    #[error("Link {link_id} is not a pump")]
    LinkNotAPump { link_id: Box<str> },

    /// Link is not a valve, raised when attempting to update a link that is not a valve
    #[error("Link {link_id} is not a valve")]
    LinkNotAValve { link_id: Box<str> },

    /// Pattern with ID not found
    #[error("Pattern {pattern_id} not found")]
    PatternNotFound { pattern_id: Box<str> },

    /// Curve with ID not found
    #[error("Curve {curve_id} not found")]
    CurveNotFound { curve_id: Box<str> },

    /// Tank levels are invalid
    #[error(
        "Tank {tank_id} levels are invalid: initial_level={initial_level}, min_level={min_level}, max_level={max_level}"
    )]
    TankLevelsInvalid {
        tank_id: Box<str>,
        initial_level: f64,
        min_level: f64,
        max_level: f64,
    },

    /// General parse error
    #[error("{message}{}", format_suffix(.line, .context))]
    Parse {
        message: String,
        line: Option<usize>,
        context: Option<String>,
    },
}

fn format_suffix(line: &Option<usize>, context: &Option<String>) -> String {
    let mut s = String::new();
    if let Some(l) = line {
        s.push_str(&format!(" (line {})", l));
    }
    if let Some(ctx) = context {
        s.push_str(&format!(" [{}]", ctx));
    }
    s
}

impl InputError {
    /// Creates a parse error. This is the most common constructor.
    pub fn new(message: impl Into<String>) -> Self {
        Self::Parse {
            message: message.into(),
            line: None,
            context: None,
        }
    }

    pub fn file_open(path: impl Into<String>, source: std::io::Error) -> Self {
        Self::FileOpen {
            path: path.into(),
            source,
        }
    }

    pub fn with_line(mut self, line_num: usize) -> Self {
        if let Self::Parse { ref mut line, .. } = self {
            *line = Some(line_num);
        }
        self
    }

    pub fn with_context(mut self, context: impl Into<String>) -> Self {
        if let Self::Parse {
            context: ref mut ctx,
            ..
        } = self
        {
            *ctx = Some(context.into());
        }
        self
    }

    pub fn is_file_error(&self) -> bool {
        matches!(self, Self::FileOpen { .. } | Self::FileRead(_))
    }
}

/// Helper trait for converting Option to Result with error context
pub trait OptionExt<T> {
    fn ok_or_missing(self, field: &str) -> Result<T, InputError>;
}

impl<T> OptionExt<T> for Option<T> {
    fn ok_or_missing(self, field: &str) -> Result<T, InputError> {
        self.ok_or_else(|| InputError::new(format!("Missing required field: {}", field)))
    }
}

/// Helper trait for parsing with better error messages
pub trait ParseExt {
    fn parse_field<T: FromStr>(&self, field: &str) -> Result<T, InputError>;
}

impl ParseExt for &str {
    fn parse_field<T: FromStr>(&self, field: &str) -> Result<T, InputError> {
        self.parse::<T>()
            .map_err(|_| InputError::new(format!("Invalid {}: '{}'", field, self)))
    }
}
