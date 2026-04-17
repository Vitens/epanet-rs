use std::str::FromStr;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SolverError {
  #[error("Solver not initialized")]
  NotInitialized,
  #[error("Singular matrix: check connectivity at node '{node_id}'")]
  SingularMatrix { node_id: Box<str> },
  #[error("Maximum number of iterations reached: {max_trials}")]
  MaxIterations { max_trials: usize },
  #[error("Network topology has changed since the solver was created, please recreate the solver")]
  StaleTopology,
  #[error("Failed to compute Cholesky factorization: {0}")]
  Factorization(String),
}

#[derive(Debug, Error)]
pub enum InputError {
  #[error("Cannot open file '{path}': {source}")]
  FileOpen {
    path: String,
    source: std::io::Error,
  },

  #[error("IO error: {0}")]
  FileRead(#[from] std::io::Error),

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
    Self::Parse { message: message.into(), line: None, context: None }
  }

  pub fn file_open(path: impl Into<String>, source: std::io::Error) -> Self {
    Self::FileOpen { path: path.into(), source }
  }

  pub fn with_line(mut self, line_num: usize) -> Self {
    if let Self::Parse { ref mut line, .. } = self {
      *line = Some(line_num);
    }
    self
  }

  pub fn with_context(mut self, context: impl Into<String>) -> Self {
    if let Self::Parse { context: ref mut ctx, .. } = self {
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
    self.parse::<T>().map_err(|_| InputError::new(format!("Invalid {}: '{}'", field, self)))
  }
}