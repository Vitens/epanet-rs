use std::str::FromStr;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SolverError {
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
#[error("{message}{}", format_suffix(.line, .context))]
pub struct InputError {
  pub message: String,
  pub line: Option<usize>,
  pub context: Option<String>,
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
  pub fn new(message: impl Into<String>) -> Self {
    Self { message: message.into(), line: None, context: None }
  }

  pub fn with_line(mut self, line: usize) -> Self {
    self.line = Some(line);
    self
  }

  pub fn with_context(mut self, context: impl Into<String>) -> Self {
    self.context = Some(context.into());
    self
  }
}

impl From<std::io::Error> for InputError {
  fn from(err: std::io::Error) -> Self {
    InputError::new(format!("IO error: {}", err))
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