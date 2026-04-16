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
mod error_codes;

pub use project::Project;
pub use crate::solver::simulation::Simulation;
pub use error_codes::ErrorCode;
