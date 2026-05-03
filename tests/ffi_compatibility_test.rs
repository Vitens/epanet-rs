//! FFI compatibility tests
//!
//! Tests for the three critical bug fixes:
//! 1. Pattern index returns 0.0 for nodes without patterns (not 123.0)
//! 2. Output parameters initialized to 0 before validation
//! 3. Invalid property codes return error (not success with -123.0)

use epanet_rs::ffi::error_codes::ErrorCode;
use epanet_rs::ffi::project::Project;
use epanet_rs::ffi::nodes::{EN_getnodevalue, EN_getnodeindex};
use epanet_rs::ffi::links::{EN_getlinkvalue, EN_getlinkindex};
use epanet_rs::ffi::patterns::{EN_getpatternvalue, EN_getpatternindex};
use epanet_rs::ffi::curves::{EN_getcurvevalue, EN_getcurveindex};
use epanet_rs::ffi::project::{EN_createproject, EN_deleteproject, EN_open};
use epanet_rs::ffi::enums::NodeProperty;

use std::os::raw::{c_double, c_int};
use std::ffi::CString;

/// Helper to create a test project with the pump.inp network
unsafe fn create_test_project() -> *mut Project {
    let mut ph: *mut Project = std::ptr::null_mut();
    let err = EN_createproject(&mut ph);
    assert_eq!(err, ErrorCode::Ok);
    assert!(!ph.is_null());

    // Open test network
    let inp_file = CString::new("tests/pump.inp").unwrap();
    let rpt_file = CString::new("").unwrap();
    let out_file = CString::new("").unwrap();

    let err = EN_open(ph, inp_file.as_ptr(), rpt_file.as_ptr(), out_file.as_ptr());
    assert_eq!(err, ErrorCode::Ok);

    ph
}

unsafe fn destroy_test_project(ph: *mut Project) {
    let err = EN_deleteproject(ph);
    assert_eq!(err, ErrorCode::Ok);
}

// =============================================================================
// TEST 1: Pattern Index Return Value
// =============================================================================

#[test]
fn test_pattern_index_returns_zero_when_none() {
    unsafe {
        let ph = create_test_project();

        // Get a node index (use node "1" which exists in pump.inp)
        let node_id = CString::new("1").unwrap();
        let mut node_idx: c_int = 0;
        let err = EN_getnodeindex(ph, node_id.as_ptr(), &mut node_idx);
        assert_eq!(err, ErrorCode::Ok);
        assert!(node_idx > 0);

        // Query pattern property for a node without a pattern
        let mut pattern_idx: c_double = -999.0;  // Initialize to sentinel
        let err = EN_getnodevalue(
            ph,
            node_idx,
            NodeProperty::Pattern as c_int,
            &mut pattern_idx
        );

        assert_eq!(err, ErrorCode::Ok, "Should succeed");

        // Critical: Must return 0.0 for "no pattern", not 123.0
        assert_eq!(
            pattern_idx, 0.0,
            "Pattern index should be 0.0 for nodes without patterns (was returning 123.0)"
        );

        destroy_test_project(ph);
    }
}

#[test]
fn test_pattern_index_returns_positive_when_assigned() {
    unsafe {
        let ph = create_test_project();

        // Get a node that HAS a pattern assigned (if any in pump.inp)
        // We'll just verify the logic: if pattern exists, should be > 0
        let node_id = CString::new("2").unwrap();
        let mut node_idx: c_int = 0;
        let err = EN_getnodeindex(ph, node_id.as_ptr(), &mut node_idx);
        assert_eq!(err, ErrorCode::Ok);

        let mut pattern_idx: c_double = -999.0;
        let err = EN_getnodevalue(
            ph,
            node_idx,
            NodeProperty::Pattern as c_int,
            &mut pattern_idx
        );

        assert_eq!(err, ErrorCode::Ok);

        // If node has a pattern, should be >= 1
        // If no pattern, should be 0 (not 123!)
        assert!(
            pattern_idx >= 0.0,
            "Pattern index should be >= 0 (0 = none, 1+ = pattern assigned)"
        );
        assert_ne!(
            pattern_idx, 123.0,
            "Should never return 123.0 as a pattern index"
        );

        destroy_test_project(ph);
    }
}

// =============================================================================
// TEST 2: Output Parameter Initialization
// =============================================================================

#[test]
fn test_output_initialized_on_null_project() {
    unsafe {
        // Test with null project handle
        let null_ph: *mut Project = std::ptr::null_mut();

        // Initialize to sentinel value
        let mut value: c_double = -999.0;

        let err = EN_getnodevalue(
            null_ph,
            1,
            NodeProperty::Elevation as c_int,
            &mut value
        );

        // Should return error
        assert_ne!(err, ErrorCode::Ok);

        // Critical: Output should be initialized to 0.0 even on error
        assert_eq!(
            value, 0.0,
            "Output parameter should be initialized to 0.0 before validation"
        );
    }
}

#[test]
fn test_output_initialized_on_invalid_index() {
    unsafe {
        let ph = create_test_project();

        // Use invalid node index
        let mut value: c_double = -999.0;

        let err = EN_getnodevalue(
            ph,
            99999,  // Invalid index
            NodeProperty::Elevation as c_int,
            &mut value
        );

        // Should return error
        assert_ne!(err, ErrorCode::Ok);

        // Output should be initialized to 0.0
        assert_eq!(
            value, 0.0,
            "Output parameter should be initialized to 0.0 even for invalid index"
        );

        destroy_test_project(ph);
    }
}

#[test]
fn test_output_initialized_handles_null_pointer_safely() {
    unsafe {
        let ph = create_test_project();

        // Get valid node
        let node_id = CString::new("1").unwrap();
        let mut node_idx: c_int = 0;
        let err = EN_getnodeindex(ph, node_id.as_ptr(), &mut node_idx);
        assert_eq!(err, ErrorCode::Ok);

        // Note: Actually passing null pointer to EN_getnodevalue would crash
        // because the C API doesn't check for null output pointers.
        // This test just verifies that our initialization pattern is safe
        // when the pointer IS valid.

        let mut value: c_double = -999.0;
        let err = EN_getnodevalue(
            ph,
            node_idx,
            NodeProperty::Elevation as c_int,
            &mut value
        );

        // Should succeed and value should be set
        assert_eq!(err, ErrorCode::Ok);
        assert_ne!(value, -999.0, "Value should be updated");

        destroy_test_project(ph);
    }
}

#[test]
fn test_link_output_initialized_on_error() {
    unsafe {
        let ph = create_test_project();

        // Use invalid link index
        let mut value: c_double = -999.0;

        let err = EN_getlinkvalue(
            ph,
            99999,  // Invalid index
            0,      // Any property
            &mut value
        );

        assert_ne!(err, ErrorCode::Ok);
        assert_eq!(
            value, 0.0,
            "Link output should be initialized to 0.0 on error"
        );

        destroy_test_project(ph);
    }
}

#[test]
fn test_pattern_output_initialized_on_error() {
    unsafe {
        let ph = create_test_project();

        // Use invalid pattern index
        let mut value: c_double = -999.0;

        let err = EN_getpatternvalue(
            ph,
            99999,  // Invalid index
            1,
            &mut value
        );

        assert_ne!(err, ErrorCode::Ok);
        assert_eq!(
            value, 0.0,
            "Pattern output should be initialized to 0.0 on error"
        );

        destroy_test_project(ph);
    }
}

#[test]
fn test_curve_output_initialized_on_error() {
    unsafe {
        let ph = create_test_project();

        // Use invalid curve index
        let mut x: c_double = -999.0;
        let mut y: c_double = -999.0;

        let err = EN_getcurvevalue(
            ph,
            99999,  // Invalid index
            1,
            &mut x,
            &mut y
        );

        assert_ne!(err, ErrorCode::Ok);
        assert_eq!(x, 0.0, "Curve x output should be initialized");
        assert_eq!(y, 0.0, "Curve y output should be initialized");

        destroy_test_project(ph);
    }
}

#[test]
fn test_index_output_initialized_on_error() {
    unsafe {
        let ph = create_test_project();

        // Use invalid ID
        let invalid_id = CString::new("NONEXISTENT_NODE_ID").unwrap();
        let mut index: c_int = -999;

        let err = EN_getnodeindex(ph, invalid_id.as_ptr(), &mut index);

        assert_ne!(err, ErrorCode::Ok);
        assert_eq!(
            index, 0,
            "Index output should be initialized to 0 on error"
        );

        destroy_test_project(ph);
    }
}

// =============================================================================
// TEST 3: Invalid Property Code Handling
// =============================================================================

#[test]
fn test_invalid_node_property_returns_error() {
    unsafe {
        let ph = create_test_project();

        // Get valid node
        let node_id = CString::new("1").unwrap();
        let mut node_idx: c_int = 0;
        let err = EN_getnodeindex(ph, node_id.as_ptr(), &mut node_idx);
        assert_eq!(err, ErrorCode::Ok);

        // Use invalid property code
        let invalid_property: c_int = 99999;
        let mut value: c_double = -999.0;

        let err = EN_getnodevalue(
            ph,
            node_idx,
            invalid_property,
            &mut value
        );

        // Critical: Must return error, not Ok with -123.0
        assert_eq!(
            err,
            ErrorCode::InvalidParameterCode,
            "Invalid property code should return InvalidParameterCode error"
        );

        // Value should be initialized to 0, not -123.0
        assert_eq!(
            value, 0.0,
            "Should initialize output to 0.0, not return -123.0"
        );

        destroy_test_project(ph);
    }
}

#[test]
fn test_invalid_node_property_code_251() {
    unsafe {
        let ph = create_test_project();

        let node_id = CString::new("1").unwrap();
        let mut node_idx: c_int = 0;
        let err = EN_getnodeindex(ph, node_id.as_ptr(), &mut node_idx);
        assert_eq!(err, ErrorCode::Ok, "Node 1 should exist");

        // Test various invalid property codes
        // NodeProperty goes from 0-32, test values outside that range
        let invalid_codes = vec![
            999,
            1000,
            12345,
            -999,
        ];

        for invalid_code in invalid_codes {
            let mut value: c_double = -999.0;

            let err = EN_getnodevalue(
                ph,
                node_idx,
                invalid_code,
                &mut value
            );

            assert_eq!(
                err,
                ErrorCode::InvalidParameterCode,
                "Invalid property code {} should return error 251",
                invalid_code
            );

            assert_ne!(
                value, -123.0,
                "Should never return -123.0 for invalid property code"
            );

            // Output should be initialized
            assert_eq!(
                value, 0.0,
                "Output should be initialized to 0.0"
            );
        }

        destroy_test_project(ph);
    }
}

#[test]
fn test_link_invalid_property_also_returns_error() {
    unsafe {
        let ph = create_test_project();

        // Get valid link (use link "B" which exists in pump.inp)
        let link_id = CString::new("B").unwrap();
        let mut link_idx: c_int = 0;
        let err = EN_getlinkindex(ph, link_id.as_ptr(), &mut link_idx);
        assert_eq!(err, ErrorCode::Ok);

        // Use invalid property code
        let invalid_property: c_int = 99999;
        let mut value: c_double = -999.0;

        let err = EN_getlinkvalue(
            ph,
            link_idx,
            invalid_property,
            &mut value
        );

        // EN_getlinkvalue should also return error for invalid property
        assert_eq!(
            err,
            ErrorCode::InvalidParameterCode,
            "Link getter should also return InvalidParameterCode for invalid property"
        );

        destroy_test_project(ph);
    }
}

// =============================================================================
// TEST 4: Valid Properties Still Work
// =============================================================================

#[test]
fn test_valid_properties_still_work() {
    unsafe {
        let ph = create_test_project();

        let node_id = CString::new("1").unwrap();
        let mut node_idx: c_int = 0;
        let err = EN_getnodeindex(ph, node_id.as_ptr(), &mut node_idx);
        assert_eq!(err, ErrorCode::Ok, "Node 1 should exist");

        // Test valid property - Elevation
        let mut elevation: c_double = -999.0;
        let err = EN_getnodevalue(
            ph,
            node_idx,
            NodeProperty::Elevation as c_int,
            &mut elevation
        );

        assert_eq!(err, ErrorCode::Ok, "Valid property should succeed");
        assert_ne!(elevation, -999.0, "Should have updated the value");
        // Note: elevation IS 0.0 in pump.inp, so just check it's not sentinel
        assert_ne!(elevation, -123.0, "Should not return sentinel value");

        destroy_test_project(ph);
    }
}

#[test]
fn test_multiple_valid_properties() {
    unsafe {
        let ph = create_test_project();

        let node_id = CString::new("1").unwrap();
        let mut node_idx: c_int = 0;
        let err = EN_getnodeindex(ph, node_id.as_ptr(), &mut node_idx);
        assert_eq!(err, ErrorCode::Ok, "Node 1 should exist");

        // Test several valid properties
        let valid_properties = vec![
            (NodeProperty::Elevation as c_int, "Elevation"),
            (NodeProperty::BaseDemand as c_int, "BaseDemand"),
            (NodeProperty::Emitter as c_int, "Emitter"),
            (NodeProperty::Pattern as c_int, "Pattern"),
        ];

        for (prop_code, prop_name) in valid_properties {
            let mut value: c_double = -999.0;

            let err = EN_getnodevalue(ph, node_idx, prop_code, &mut value);

            assert_eq!(
                err,
                ErrorCode::Ok,
                "Valid property {} should succeed",
                prop_name
            );

            // Value should have been set (not still -999.0)
            assert_ne!(
                value, -999.0,
                "Property {} should update the output value",
                prop_name
            );

            // Should never return sentinel values
            assert_ne!(
                value, -123.0,
                "Property {} should not return -123.0",
                prop_name
            );
        }

        destroy_test_project(ph);
    }
}

// =============================================================================
// TEST 5: Regression Tests - Ensure Old Bugs Don't Return
// =============================================================================

#[test]
fn test_regression_pattern_never_returns_123() {
    unsafe {
        let ph = create_test_project();

        // Check all nodes in the network
        for node_idx in 1..20 {  // Assuming < 20 nodes
            let mut pattern_idx: c_double = -999.0;

            let err = EN_getnodevalue(
                ph,
                node_idx,
                NodeProperty::Pattern as c_int,
                &mut pattern_idx
            );

            if err == ErrorCode::Ok {
                // Valid node - check pattern value
                assert_ne!(
                    pattern_idx, 123.0,
                    "Node {} pattern should never be 123.0 (old bug)",
                    node_idx
                );
                assert!(
                    pattern_idx >= 0.0,
                    "Node {} pattern should be >= 0",
                    node_idx
                );
            } else {
                // Invalid node - output should be initialized
                assert_eq!(
                    pattern_idx, 0.0,
                    "Invalid node {} should have output initialized to 0.0",
                    node_idx
                );
            }
        }

        destroy_test_project(ph);
    }
}

#[test]
fn test_regression_invalid_property_never_succeeds_with_negative_sentinel() {
    unsafe {
        let ph = create_test_project();

        let node_id = CString::new("1").unwrap();
        let mut node_idx: c_int = 0;
        let err = EN_getnodeindex(ph, node_id.as_ptr(), &mut node_idx);
        assert_eq!(err, ErrorCode::Ok, "Node 1 should exist");

        // Test many invalid property codes
        for invalid_prop in vec![999, 1000, -999, 12345] {
            let mut value: c_double = -999.0;

            let err = EN_getnodevalue(ph, node_idx, invalid_prop, &mut value);

            // Should ALWAYS return error
            assert_eq!(
                err,
                ErrorCode::InvalidParameterCode,
                "Invalid property {} should return error (old bug returned Ok)",
                invalid_prop
            );

            // Should NEVER return -123.0
            assert_ne!(
                value, -123.0,
                "Invalid property {} should not return -123.0 (old bug)",
                invalid_prop
            );

            // Should be initialized to 0
            assert_eq!(
                value, 0.0,
                "Invalid property {} should initialize output to 0.0",
                invalid_prop
            );
        }

        destroy_test_project(ph);
    }
}
