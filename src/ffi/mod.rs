//! C-compatible FFI layer implementing the EPANET 2.3 toolkit API.
//!
//! Every function of the official toolkit is exported. Functions whose
//! underlying feature is not supported by `epanet-rs` (water quality,
//! report/hydraulics files, rule evaluation, ...) return
//! [`error_codes::ErrorCode::NotImplemented`] and are marked with a `TODO`.

pub mod analysis_options;
pub mod controls;
pub mod curves;
pub mod demands;
pub mod enums;
pub mod error_codes;
pub mod hydraulic_solver;
pub mod links;
pub mod nodes;
pub mod patterns;
pub mod project;
pub mod quality_solver;
pub mod report;
pub mod rules;
pub mod util;

#[cfg(test)]
mod tests;
