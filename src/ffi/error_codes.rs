/// EPANET toolkit error/warning codes.
///
/// Values match the official EPANET 2.3 error numbering.
/// Uses `#[repr(i32)]` so the enum is ABI-compatible with `c_int`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum ErrorCode {
    Ok                = 0,
    /// Input data error (e.g. malformed INP file)
    InputError        = 200,
    /// Specified file not found
    FileNotFound      = 302,
    /// Invalid project handle (null pointer)
    InvalidHandle     = 303,
}
