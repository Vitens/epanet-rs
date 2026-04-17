use crate::ffi::enums::CountType;
use crate::ffi::error_codes::ErrorCode;
use crate::ffi::project::{Project, get_simulation};

use crate::model::node::NodeType;

use std::os::raw::{c_char, c_int};

/// Retrieves the number of objects of a given type in a project.
#[unsafe(no_mangle)]
pub extern "C" fn EN_getcount(ph: *mut Project, object: c_int, out_count: *mut c_int) -> ErrorCode {
  let simulation = get_simulation!(ph);

  let count_type = match CountType::from_repr(object) {
    Some(ct) => ct,
    None => return ErrorCode::InvalidParameterCode,
  };

  let net = &simulation.network;
  let count = match count_type {
    CountType::NodeCount    => net.nodes.len(),
    CountType::TankCount    => net.nodes.iter().filter(|n| matches!(n.node_type, NodeType::Tank(_) | NodeType::Reservoir(_))).count(),
    CountType::LinkCount    => net.links.len(),
    CountType::PatCount     => net.patterns.len(),
    CountType::CurveCount   => net.curves.len(),
    CountType::ControlCount => net.controls.len(),
    CountType::RuleCount    => 0, // TODO: implement rule counting
  };

  unsafe { *out_count = count as c_int };
  ErrorCode::Ok
}

/// Retrieves the text of an error message given its error code.
///
/// Writes up to `maxLen` bytes (including the null terminator) into `errmsg`.
#[unsafe(no_mangle)]
pub extern "C" fn EN_geterror(errcode: c_int, errmsg: *mut c_char, max_len: c_int) -> ErrorCode {
    if errmsg.is_null() || max_len <= 0 {
        return ErrorCode::InvalidFormat;
    }

    let msg = match ErrorCode::from_repr(errcode) {
        Some(code) => code.to_string(),
        None => format!("Unknown error code: {}", errcode),
    };

    let buf = unsafe { std::slice::from_raw_parts_mut(errmsg as *mut u8, max_len as usize) };
    let bytes = msg.as_bytes();
    let copy_len = bytes.len().min(buf.len() - 1);
    buf[..copy_len].copy_from_slice(&bytes[..copy_len]);
    buf[copy_len] = 0;

    ErrorCode::Ok
}
