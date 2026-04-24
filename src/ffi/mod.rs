//! C-compatible FFI layer implementing a subset of the EPANET 2.3 toolkit API.

pub mod project;
pub mod report;
pub mod nodes;
pub mod links;
pub mod curves;
pub mod patterns;
pub mod error_codes;
pub mod analysis_options;
pub mod hydraulic_solver;
pub mod quality_solver;
pub mod enums;