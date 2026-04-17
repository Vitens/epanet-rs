//! C-compatible FFI layer implementing the EPANET 2.3 toolkit API.
//!
//! Equivalent C header:
//! ```c
//! typedef struct Project *EN_Project;
//!
//! int EN_createproject(EN_Project *ph);
//! int EN_deleteproject(EN_Project ph);
//! int EN_open(EN_Project ph, const char *inpFile, const char *rptFile, const char *outFile);
//! ```

mod project;
mod report;
mod nodes;
mod links;
mod curves;
mod error_codes;
mod analysis_options;
mod hydraulic_solver;
pub mod enums;