use std::str::FromStr;

#[derive(Debug)]
pub struct InputError {
  pub message: String,
  pub line: Option<usize>,
  pub context: Option<String>,
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

impl std::fmt::Display for InputError {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.message)?;
    if let Some(line) = self.line {
      write!(f, " (line {})", line)?;
    }
    if let Some(ctx) = &self.context {
      write!(f, " [{}]", ctx)?;
    }
    Ok(())
  }
}

impl std::error::Error for InputError {}

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