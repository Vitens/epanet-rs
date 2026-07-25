//! Helpers shared by the FFI layer for reading and writing C strings and out-parameters.

use std::ffi::CStr;
use std::os::raw::c_char;

/// Maximum number of characters in an ID name (EPANET `EN_MAXID`).
pub const MAX_ID: usize = 31;
/// Maximum number of characters in a message or title line (EPANET `EN_MAXMSG`).
pub const MAX_MSG: usize = 255;

/// Borrows a NUL-terminated C string as a `&str`.
///
/// Returns `None` for null pointers or non-UTF-8 content.
///
/// # Safety
///
/// `ptr` must either be null or point to a NUL-terminated C string that stays
/// valid for the duration of the borrow.
pub(crate) unsafe fn read_str<'a>(ptr: *const c_char) -> Option<&'a str> {
    if ptr.is_null() {
        return None;
    }
    unsafe { CStr::from_ptr(ptr) }.to_str().ok()
}

/// Writes `value` followed by a NUL terminator into `dest`, truncating the
/// value to `max_len` bytes.
///
/// # Safety
///
/// `dest` must either be null or point to a buffer of at least `max_len + 1`
/// writable bytes.
pub(crate) unsafe fn write_str(dest: *mut c_char, value: &str, max_len: usize) {
    if dest.is_null() {
        return;
    }
    let bytes = value.as_bytes();
    let mut len = bytes.len().min(max_len);
    // never split a multi-byte character in half
    while len > 0 && !value.is_char_boundary(len) {
        len -= 1;
    }
    unsafe {
        std::ptr::copy_nonoverlapping(bytes.as_ptr(), dest as *mut u8, len);
        *dest.add(len) = 0;
    }
}

/// Writes `value` to `ptr` unless it is null.
///
/// # Safety
///
/// `ptr` must either be null or point to writable memory for one `T`.
pub(crate) unsafe fn write_out<T>(ptr: *mut T, value: T) {
    if !ptr.is_null() {
        unsafe { *ptr = value };
    }
}
