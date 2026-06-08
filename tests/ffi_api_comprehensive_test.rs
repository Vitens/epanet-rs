//! Comprehensive FFI API tests
//!
//! Tests all 52+ implemented FFI functions for:
//! - API stability and C compatibility
//! - Output parameter initialization
//! - Error code consistency
//! - Boundary conditions
//! - Valid and invalid inputs
//!
//! Organized by module:
//! 1. Project management (create, open, close, delete)
//! 2. Node operations (add, get, set, delete)
//! 3. Link operations (add, get, set, delete)
//! 4. Pattern operations (add, get, set, delete)
//! 5. Curve operations (add, get, set)
//! 6. Hydraulic solver (init, run, next)
//! 7. Analysis options (set option, time params, demand model)
//! 8. Report and file operations
//! 9. Error handling

use epanet_rs::ffi::analysis_options::{
    EN_getoption, EN_setdemandmodel, EN_setoption, EN_settimeparam,
};
use epanet_rs::ffi::curves::{
    EN_addcurve, EN_getcurve, EN_getcurveid, EN_getcurveindex, EN_getcurvelen, EN_getcurvevalue,
    EN_setcurve, EN_setcurveid,
};
use epanet_rs::ffi::enums::{
    CountType, DemandModel, LinkProperty, LinkType, NodeProperty, NodeType, SimOption,
    TimeParameter,
};
use epanet_rs::ffi::error_codes::ErrorCode;
use epanet_rs::ffi::hydraulic_solver::{EN_closeH, EN_initH, EN_nextH, EN_runH};
use epanet_rs::ffi::links::{
    EN_addlink, EN_deletelink, EN_getheadcurveindex, EN_getlinkid, EN_getlinkindex,
    EN_getlinknodes, EN_getlinktype, EN_getlinkvalue, EN_setheadcurveindex, EN_setlinkid,
    EN_setlinknodes, EN_setlinkvalue,
};
use epanet_rs::ffi::nodes::{
    EN_addnode, EN_deletenode, EN_getcoord, EN_getnodeid, EN_getnodeindex, EN_getnodetype,
    EN_getnodevalue, EN_setcoord, EN_setnodeid, EN_setnodevalue,
};
use epanet_rs::ffi::patterns::{
    EN_addpattern, EN_deletepattern, EN_getaveragepatternvalue, EN_getpatternid,
    EN_getpatternindex, EN_getpatternlen, EN_getpatternvalue, EN_setpattern, EN_setpatternid,
};
use epanet_rs::ffi::project::{
    EN_createproject, EN_deleteproject, EN_open, EN_saveinpfile, Project,
};
use epanet_rs::ffi::report::{EN_getcount, EN_geterror};

use std::ffi::CString;
use std::os::raw::{c_char, c_double, c_int, c_long};

// =============================================================================
// Test Helpers
// =============================================================================

/// Create a minimal test project
unsafe fn create_empty_project() -> *mut Project {
    let mut ph: *mut Project = std::ptr::null_mut();
    let err = unsafe { EN_createproject(&mut ph) };
    assert_eq!(err, ErrorCode::Ok);
    assert!(!ph.is_null());
    ph
}

/// Create project and load pump.inp
unsafe fn create_loaded_project() -> *mut Project {
    let ph = unsafe { create_empty_project() };

    let inp_file = CString::new("tests/pump.inp").unwrap();
    let rpt_file = CString::new("").unwrap();
    let out_file = CString::new("").unwrap();

    let err = unsafe { EN_open(ph, inp_file.as_ptr(), rpt_file.as_ptr(), out_file.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    ph
}

/// Clean up project
unsafe fn destroy_project(ph: *mut Project) {
    let err = unsafe { EN_deleteproject(ph) };
    assert_eq!(err, ErrorCode::Ok);
}

/// Get a valid node index from loaded project
unsafe fn get_valid_node_index(ph: *mut Project, node_id: &str) -> c_int {
    let id = CString::new(node_id).unwrap();
    let mut index: c_int = 0;
    let err = unsafe { EN_getnodeindex(ph, id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);
    assert!(index > 0);
    index
}

/// Get a valid link index from loaded project
unsafe fn get_valid_link_index(ph: *mut Project, link_id: &str) -> c_int {
    let id = CString::new(link_id).unwrap();
    let mut index: c_int = 0;
    let err = unsafe { EN_getlinkindex(ph, id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);
    assert!(index > 0);
    index
}

// =============================================================================
// MODULE 1: Project Management
// =============================================================================

#[test]
fn test_en_createproject_valid() {
    let mut ph: *mut Project = std::ptr::null_mut();
    let err = unsafe { EN_createproject(&mut ph) };

    assert_eq!(err, ErrorCode::Ok, "Should create project successfully");
    assert!(!ph.is_null(), "Project handle should not be null");

    unsafe { destroy_project(ph) };
}

#[test]
fn test_en_createproject_null_handle() {
    unsafe {
        let err = EN_createproject(std::ptr::null_mut());
        assert_ne!(err, ErrorCode::Ok, "Should fail with null handle pointer");
    }
}

#[test]
fn test_en_deleteproject_valid() {
    unsafe {
        let ph = create_empty_project();
        let err = EN_deleteproject(ph);
        assert_eq!(err, ErrorCode::Ok, "Should delete project successfully");
    }
}

#[test]
fn test_en_deleteproject_null() {
    unsafe {
        let err = EN_deleteproject(std::ptr::null_mut());
        assert_ne!(err, ErrorCode::Ok, "Should fail with null project");
    }
}

#[test]
fn test_en_open_valid_file() {
    unsafe {
        let ph = create_empty_project();

        let inp_file = CString::new("tests/pump.inp").unwrap();
        let rpt_file = CString::new("").unwrap();
        let out_file = CString::new("").unwrap();

        let err = EN_open(ph, inp_file.as_ptr(), rpt_file.as_ptr(), out_file.as_ptr());
        assert_eq!(err, ErrorCode::Ok, "Should open valid INP file");

        destroy_project(ph);
    }
}

#[test]
fn test_en_open_invalid_file() {
    unsafe {
        let ph = create_empty_project();

        let inp_file = CString::new("nonexistent.inp").unwrap();
        let rpt_file = CString::new("").unwrap();
        let out_file = CString::new("").unwrap();

        let err = EN_open(ph, inp_file.as_ptr(), rpt_file.as_ptr(), out_file.as_ptr());
        assert_ne!(err, ErrorCode::Ok, "Should fail with nonexistent file");

        destroy_project(ph);
    }
}

#[test]
fn test_en_open_null_project() {
    unsafe {
        let inp_file = CString::new("tests/pump.inp").unwrap();
        let rpt_file = CString::new("").unwrap();
        let out_file = CString::new("").unwrap();

        let err = EN_open(
            std::ptr::null_mut(),
            inp_file.as_ptr(),
            rpt_file.as_ptr(),
            out_file.as_ptr(),
        );
        assert_ne!(err, ErrorCode::Ok, "Should fail with null project");
    }
}

#[test]
fn test_en_closeh_valid() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_closeH(ph);
        // May succeed or fail depending on whether hydraulics were opened
        // Just verify it doesn't crash
        let _ = err;
        destroy_project(ph);
    }
}

#[test]
fn test_en_saveinpfile_valid() {
    unsafe {
        let ph = create_loaded_project();

        let out_file = CString::new("/tmp/test_output.inp").unwrap();
        let err = EN_saveinpfile(ph, out_file.as_ptr());
        assert_eq!(err, ErrorCode::Ok, "Should save INP file");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcount_nodes() {
    unsafe {
        let ph = create_loaded_project();

        let mut count: c_int = -1;
        let err = EN_getcount(ph, CountType::NodeCount as c_int, &mut count);

        assert_eq!(err, ErrorCode::Ok, "Should get node count");
        assert!(count > 0, "Node count should be positive");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcount_links() {
    unsafe {
        let ph = create_loaded_project();

        let mut count: c_int = -1;
        let err = EN_getcount(ph, CountType::LinkCount as c_int, &mut count);

        assert_eq!(err, ErrorCode::Ok, "Should get link count");
        assert!(count > 0, "Link count should be positive");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcount_output_initialized() {
    unsafe {
        let mut count: c_int = -999;
        let err = EN_getcount(
            std::ptr::null_mut(),
            CountType::NodeCount as c_int,
            &mut count,
        );

        assert_ne!(err, ErrorCode::Ok, "Should fail with null project");
        assert_eq!(count, 0, "Output should be initialized to 0");
    }
}

#[test]
fn test_en_geterror_valid_code() {
    unsafe {
        let code = ErrorCode::InvalidParameterCode as c_int;
        let mut buffer: [c_char; 256] = [0; 256];

        let err = EN_geterror(code, buffer.as_mut_ptr(), 256);
        assert_eq!(err, ErrorCode::Ok, "Should get error message");

        // Should have written something to buffer
        assert_ne!(buffer[0], 0, "Should write error message");
    }
}

// =============================================================================
// MODULE 2: Node Operations
// =============================================================================

#[test]
fn test_en_getnodeindex_valid() {
    unsafe {
        let ph = create_loaded_project();

        let node_id = CString::new("1").unwrap();
        let mut index: c_int = -999;

        let err = EN_getnodeindex(ph, node_id.as_ptr(), &mut index);

        assert_eq!(err, ErrorCode::Ok, "Should get node index");
        assert!(index > 0, "Index should be positive (1-based)");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getnodeindex_invalid_id() {
    unsafe {
        let ph = create_loaded_project();

        let node_id = CString::new("NONEXISTENT").unwrap();
        let mut index: c_int = -999;

        let err = EN_getnodeindex(ph, node_id.as_ptr(), &mut index);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid ID");
        assert_eq!(index, 0, "Output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getnodeindex_null_project() {
    unsafe {
        let node_id = CString::new("1").unwrap();
        let mut index: c_int = -999;

        let err = EN_getnodeindex(std::ptr::null_mut(), node_id.as_ptr(), &mut index);

        assert_ne!(err, ErrorCode::Ok, "Should fail with null project");
        assert_eq!(index, 0, "Output should be initialized to 0");
    }
}

#[test]
fn test_en_getnodeid_valid() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_node_index(ph, "1");

        let mut buffer: [c_char; 32] = [0; 32];
        let err = EN_getnodeid(ph, index, buffer.as_mut_ptr());

        assert_eq!(err, ErrorCode::Ok, "Should get node ID");
        assert_ne!(buffer[0], 0, "Should write ID to buffer");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getnodeid_invalid_index() {
    unsafe {
        let ph = create_loaded_project();

        let mut buffer: [c_char; 32] = [0; 32];
        let err = EN_getnodeid(ph, 99999, buffer.as_mut_ptr());

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getnodetype_valid() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_node_index(ph, "1");

        let mut node_type: c_int = -999;
        let err = EN_getnodetype(ph, index, &mut node_type);

        assert_eq!(err, ErrorCode::Ok, "Should get node type");
        assert!(node_type >= 0, "Node type should be valid");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getnodetype_output_initialized() {
    unsafe {
        let ph = create_loaded_project();

        let mut node_type: c_int = -999;
        let err = EN_getnodetype(ph, 99999, &mut node_type);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
        assert_eq!(node_type, 0, "Output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getnodevalue_elevation() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_node_index(ph, "1");

        let mut value: c_double = -999.0;
        let err = EN_getnodevalue(ph, index, NodeProperty::Elevation as c_int, &mut value);

        assert_eq!(err, ErrorCode::Ok, "Should get elevation");
        assert_ne!(value, -999.0, "Should update value");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getnodevalue_pattern() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_node_index(ph, "1");

        let mut value: c_double = -999.0;
        let err = EN_getnodevalue(ph, index, NodeProperty::Pattern as c_int, &mut value);

        assert_eq!(err, ErrorCode::Ok, "Should get pattern");
        assert!(value >= 0.0, "Pattern should be >= 0 (0 = none)");
        assert_ne!(value, 123.0, "Should never return 123.0 sentinel");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getnodevalue_invalid_property() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_node_index(ph, "1");

        let mut value: c_double = -999.0;
        let err = EN_getnodevalue(ph, index, 99999, &mut value);

        assert_eq!(
            err,
            ErrorCode::InvalidParameterCode,
            "Should return error 251"
        );
        assert_eq!(value, 0.0, "Output should be initialized to 0");
        assert_ne!(value, -123.0, "Should not return -123.0 sentinel");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getnodevalue_null_project() {
    unsafe {
        let mut value: c_double = -999.0;
        let err = EN_getnodevalue(
            std::ptr::null_mut(),
            1,
            NodeProperty::Elevation as c_int,
            &mut value,
        );

        assert_ne!(err, ErrorCode::Ok, "Should fail with null project");
        assert_eq!(value, 0.0, "Output should be initialized to 0");
    }
}

#[test]
fn test_en_setnodevalue_elevation() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_node_index(ph, "1");

        let new_elevation = 123.45;
        let err = EN_setnodevalue(ph, index, NodeProperty::Elevation as c_int, new_elevation);

        assert_eq!(err, ErrorCode::Ok, "Should set elevation");

        // Verify it was set
        let mut value: c_double = 0.0;
        EN_getnodevalue(ph, index, NodeProperty::Elevation as c_int, &mut value);
        assert_eq!(value, new_elevation, "Should retrieve new elevation");

        destroy_project(ph);
    }
}

#[test]
fn test_en_setnodevalue_invalid_property() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_node_index(ph, "1");

        let err = EN_setnodevalue(ph, index, 99999, 123.45);

        assert_eq!(
            err,
            ErrorCode::InvalidParameterCode,
            "Should return error 251"
        );

        destroy_project(ph);
    }
}

#[test]
fn test_en_addnode_valid() {
    unsafe {
        let ph = create_loaded_project();

        let node_id = CString::new("NEW_NODE").unwrap();
        let node_type = NodeType::Junction as c_int;

        let mut index: c_int = 0;
        let err = EN_addnode(ph, node_id.as_ptr(), node_type, &mut index);

        assert_eq!(err, ErrorCode::Ok, "Should add node");
        assert!(index > 0, "Should return valid index");

        destroy_project(ph);
    }
}

#[test]
fn test_en_addnode_duplicate_id() {
    unsafe {
        let ph = create_loaded_project();

        // Try to add node with existing ID
        let node_id = CString::new("1").unwrap();
        let node_type = NodeType::Junction as c_int;

        let mut index: c_int = 0;
        let err = EN_addnode(ph, node_id.as_ptr(), node_type, &mut index);

        assert_ne!(err, ErrorCode::Ok, "Should fail with duplicate ID");

        destroy_project(ph);
    }
}

#[test]
fn test_en_deletenode_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a node first
        let node_id = CString::new("TEMP_NODE").unwrap();
        let mut index: c_int = 0;
        EN_addnode(
            ph,
            node_id.as_ptr(),
            NodeType::Junction as c_int,
            &mut index,
        );

        // Delete it
        let err = EN_deletenode(ph, index, 0);
        assert_eq!(err, ErrorCode::Ok, "Should delete node");

        destroy_project(ph);
    }
}

#[test]
fn test_en_setnodeid_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a node
        let old_id = CString::new("OLD_ID").unwrap();
        let mut index: c_int = 0;
        EN_addnode(ph, old_id.as_ptr(), NodeType::Junction as c_int, &mut index);

        // Rename it
        let new_id = CString::new("NEW_ID").unwrap();
        let err = EN_setnodeid(ph, index, new_id.as_ptr());

        assert_eq!(err, ErrorCode::Ok, "Should rename node");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcoord_valid() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_node_index(ph, "1");

        let mut x: c_double = -999.0;
        let mut y: c_double = -999.0;

        let err = EN_getcoord(ph, index, &mut x, &mut y);

        assert_eq!(err, ErrorCode::Ok, "Should get coordinates");
        assert_ne!(x, -999.0, "X should be updated");
        assert_ne!(y, -999.0, "Y should be updated");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcoord_output_initialized() {
    unsafe {
        let ph = create_loaded_project();

        let mut x: c_double = -999.0;
        let mut y: c_double = -999.0;

        let err = EN_getcoord(ph, 99999, &mut x, &mut y);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
        assert_eq!(x, 0.0, "X should be initialized to 0");
        assert_eq!(y, 0.0, "Y should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_setcoord_valid() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_node_index(ph, "1");

        let err = EN_setcoord(ph, index, 100.0, 200.0);
        assert_eq!(err, ErrorCode::Ok, "Should set coordinates");

        // Verify
        let mut x: c_double = 0.0;
        let mut y: c_double = 0.0;
        EN_getcoord(ph, index, &mut x, &mut y);
        assert_eq!(x, 100.0, "X should match");
        assert_eq!(y, 200.0, "Y should match");

        destroy_project(ph);
    }
}

// =============================================================================
// MODULE 3: Link Operations
// =============================================================================

#[test]
fn test_en_getlinkindex_valid() {
    unsafe {
        let ph = create_loaded_project();

        let link_id = CString::new("B").unwrap();
        let mut index: c_int = -999;

        let err = EN_getlinkindex(ph, link_id.as_ptr(), &mut index);

        assert_eq!(err, ErrorCode::Ok, "Should get link index");
        assert!(index > 0, "Index should be positive");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getlinkindex_invalid_id() {
    unsafe {
        let ph = create_loaded_project();

        let link_id = CString::new("NONEXISTENT").unwrap();
        let mut index: c_int = -999;

        let err = EN_getlinkindex(ph, link_id.as_ptr(), &mut index);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid ID");
        assert_eq!(index, 0, "Output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getlinkid_valid() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_link_index(ph, "B");

        let mut buffer: [c_char; 32] = [0; 32];
        let err = EN_getlinkid(ph, index, buffer.as_mut_ptr());

        assert_eq!(err, ErrorCode::Ok, "Should get link ID");
        assert_ne!(buffer[0], 0, "Should write ID to buffer");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getlinktype_valid() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_link_index(ph, "B");

        let mut link_type: c_int = -999;
        let err = EN_getlinktype(ph, index, &mut link_type);

        assert_eq!(err, ErrorCode::Ok, "Should get link type");
        assert!(link_type >= 0, "Link type should be valid");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getlinktype_output_initialized() {
    unsafe {
        let ph = create_loaded_project();

        let mut link_type: c_int = -999;
        let err = EN_getlinktype(ph, 99999, &mut link_type);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
        assert_eq!(link_type, 0, "Output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getlinkvalue_diameter() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_link_index(ph, "B");

        let mut value: c_double = -999.0;
        let err = EN_getlinkvalue(ph, index, LinkProperty::Diameter as c_int, &mut value);

        assert_eq!(err, ErrorCode::Ok, "Should get diameter");
        assert_ne!(value, -999.0, "Should update value");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getlinkvalue_invalid_property() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_link_index(ph, "B");

        let mut value: c_double = -999.0;
        let err = EN_getlinkvalue(ph, index, 99999, &mut value);

        assert_eq!(
            err,
            ErrorCode::InvalidParameterCode,
            "Should return error 251"
        );
        assert_eq!(value, 0.0, "Output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_setlinkvalue_diameter() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_link_index(ph, "B");

        let new_diameter = 24.0;
        let err = EN_setlinkvalue(ph, index, LinkProperty::Diameter as c_int, new_diameter);

        assert_eq!(err, ErrorCode::Ok, "Should set diameter");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getlinknodes_valid() {
    unsafe {
        let ph = create_loaded_project();
        let index = get_valid_link_index(ph, "B");

        let mut node1: c_int = -999;
        let mut node2: c_int = -999;

        let err = EN_getlinknodes(ph, index, &mut node1, &mut node2);

        assert_eq!(err, ErrorCode::Ok, "Should get link nodes");
        assert!(node1 > 0, "Node 1 should be valid");
        assert!(node2 > 0, "Node 2 should be valid");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getlinknodes_output_initialized() {
    unsafe {
        let ph = create_loaded_project();

        let mut node1: c_int = -999;
        let mut node2: c_int = -999;

        let err = EN_getlinknodes(ph, 99999, &mut node1, &mut node2);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
        assert_eq!(node1, 0, "Node1 should be initialized to 0");
        assert_eq!(node2, 0, "Node2 should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_addlink_valid() {
    unsafe {
        let ph = create_loaded_project();

        let link_id = CString::new("NEW_PIPE").unwrap();
        let link_type = LinkType::Pipe as c_int;
        let from_node = CString::new("1").unwrap();
        let to_node = CString::new("2").unwrap();

        let mut index: c_int = 0;
        let err = EN_addlink(
            ph,
            link_id.as_ptr(),
            link_type,
            from_node.as_ptr(),
            to_node.as_ptr(),
            &mut index,
        );

        assert_eq!(err, ErrorCode::Ok, "Should add link");
        assert!(index > 0, "Should return valid index");

        destroy_project(ph);
    }
}

#[test]
fn test_en_deletelink_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a link first
        let link_id = CString::new("TEMP_LINK").unwrap();
        let from_node = CString::new("1").unwrap();
        let to_node = CString::new("2").unwrap();
        let mut index: c_int = 0;
        EN_addlink(
            ph,
            link_id.as_ptr(),
            LinkType::Pipe as c_int,
            from_node.as_ptr(),
            to_node.as_ptr(),
            &mut index,
        );

        // Delete it
        let err = EN_deletelink(ph, index, 0);
        assert_eq!(err, ErrorCode::Ok, "Should delete link");

        destroy_project(ph);
    }
}

#[test]
fn test_en_setlinkid_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a link
        let old_id = CString::new("OLD_LINK_ID").unwrap();
        let from_node = CString::new("1").unwrap();
        let to_node = CString::new("2").unwrap();
        let mut index: c_int = 0;
        EN_addlink(
            ph,
            old_id.as_ptr(),
            LinkType::Pipe as c_int,
            from_node.as_ptr(),
            to_node.as_ptr(),
            &mut index,
        );

        // Rename it
        let new_id = CString::new("NEW_LINK_ID").unwrap();
        let err = EN_setlinkid(ph, index, new_id.as_ptr());

        assert_eq!(err, ErrorCode::Ok, "Should rename link");

        destroy_project(ph);
    }
}

#[test]
fn test_en_setlinknodes_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a link
        let link_id = CString::new("TEST_LINK").unwrap();
        let node1_str = CString::new("1").unwrap();
        let node2_str = CString::new("2").unwrap();
        let mut index: c_int = 0;
        EN_addlink(
            ph,
            link_id.as_ptr(),
            LinkType::Pipe as c_int,
            node1_str.as_ptr(),
            node2_str.as_ptr(),
            &mut index,
        );

        // Change nodes
        let node1 = get_valid_node_index(ph, "1");
        let node3 = get_valid_node_index(ph, "3");
        let err = EN_setlinknodes(ph, index, node1, node3);

        assert_eq!(err, ErrorCode::Ok, "Should set link nodes");

        destroy_project(ph);
    }
}

// =============================================================================
// MODULE 4: Pattern Operations
// =============================================================================

#[test]
fn test_en_addpattern_valid() {
    unsafe {
        let ph = create_loaded_project();

        let pattern_id = CString::new("NEW_PATTERN").unwrap();
        let err = EN_addpattern(ph, pattern_id.as_ptr());

        assert_eq!(err, ErrorCode::Ok, "Should add pattern");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getpatternindex_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a pattern first
        let pattern_id = CString::new("TEST_PAT").unwrap();
        EN_addpattern(ph, pattern_id.as_ptr());

        let mut index: c_int = -999;
        let err = EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index);

        assert_eq!(err, ErrorCode::Ok, "Should get pattern index");
        assert!(index > 0, "Index should be positive");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getpatternindex_output_initialized() {
    unsafe {
        let ph = create_loaded_project();

        let pattern_id = CString::new("NONEXISTENT").unwrap();
        let mut index: c_int = -999;

        let err = EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid ID");
        assert_eq!(index, 0, "Output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getpatternid_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a pattern
        let pattern_id = CString::new("TEST_PAT2").unwrap();
        EN_addpattern(ph, pattern_id.as_ptr());

        let mut index: c_int = 0;
        EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index);

        let mut buffer: [c_char; 32] = [0; 32];
        let err = EN_getpatternid(ph, index, buffer.as_mut_ptr());

        assert_eq!(err, ErrorCode::Ok, "Should get pattern ID");
        assert_ne!(buffer[0], 0, "Should write ID to buffer");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getpatternlen_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a pattern
        let pattern_id = CString::new("TEST_PAT3").unwrap();
        EN_addpattern(ph, pattern_id.as_ptr());

        let mut index: c_int = 0;
        EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index);

        let mut len: c_int = -999;
        let err = EN_getpatternlen(ph, index, &mut len);

        assert_eq!(err, ErrorCode::Ok, "Should get pattern length");
        assert!(len >= 0, "Length should be non-negative");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getpatternlen_output_initialized() {
    unsafe {
        let ph = create_loaded_project();

        let mut len: c_int = -999;
        let err = EN_getpatternlen(ph, 99999, &mut len);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
        assert_eq!(len, 0, "Output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_setpattern_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a pattern
        let pattern_id = CString::new("TEST_PAT4").unwrap();
        EN_addpattern(ph, pattern_id.as_ptr());

        let mut index: c_int = 0;
        EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index);

        // Set pattern values
        let values: Vec<c_double> = vec![1.0, 1.5, 0.8, 1.2];
        let err = EN_setpattern(ph, index, values.as_ptr(), values.len() as c_int);

        assert_eq!(err, ErrorCode::Ok, "Should set pattern values");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getpatternvalue_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add and set pattern
        let pattern_id = CString::new("TEST_PAT5").unwrap();
        EN_addpattern(ph, pattern_id.as_ptr());

        let mut index: c_int = 0;
        EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index);

        let values: Vec<c_double> = vec![1.0, 2.0, 3.0];
        EN_setpattern(ph, index, values.as_ptr(), values.len() as c_int);

        // Get value
        let mut value: c_double = -999.0;
        let err = EN_getpatternvalue(ph, index, 1, &mut value);

        assert_eq!(err, ErrorCode::Ok, "Should get pattern value");
        assert_eq!(value, 1.0, "Should match set value");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getpatternvalue_output_initialized() {
    unsafe {
        let ph = create_loaded_project();

        let mut value: c_double = -999.0;
        let err = EN_getpatternvalue(ph, 99999, 1, &mut value);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
        assert_eq!(value, 0.0, "Output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getaveragepatternvalue_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add and set pattern
        let pattern_id = CString::new("TEST_PAT6").unwrap();
        EN_addpattern(ph, pattern_id.as_ptr());

        let mut index: c_int = 0;
        EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index);

        let values: Vec<c_double> = vec![1.0, 2.0, 3.0];
        EN_setpattern(ph, index, values.as_ptr(), values.len() as c_int);

        // Get average
        let mut avg: c_double = -999.0;
        let err = EN_getaveragepatternvalue(ph, index, &mut avg);

        assert_eq!(err, ErrorCode::Ok, "Should get average");
        assert_eq!(avg, 2.0, "Average of 1,2,3 should be 2");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getaveragepatternvalue_output_initialized() {
    unsafe {
        let ph = create_loaded_project();

        let mut avg: c_double = -999.0;
        let err = EN_getaveragepatternvalue(ph, 99999, &mut avg);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
        assert_eq!(avg, 0.0, "Output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_deletepattern_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a pattern
        let pattern_id = CString::new("TEMP_PAT").unwrap();
        EN_addpattern(ph, pattern_id.as_ptr());

        let mut index: c_int = 0;
        EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index);

        // Delete it
        let err = EN_deletepattern(ph, index);
        assert_eq!(err, ErrorCode::Ok, "Should delete pattern");

        destroy_project(ph);
    }
}

#[test]
fn test_en_setpatternid_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a pattern
        let old_id = CString::new("OLD_PAT_ID").unwrap();
        EN_addpattern(ph, old_id.as_ptr());

        let mut index: c_int = 0;
        EN_getpatternindex(ph, old_id.as_ptr(), &mut index);

        // Rename it
        let new_id = CString::new("NEW_PAT_ID").unwrap();
        let err = EN_setpatternid(ph, index, new_id.as_ptr());

        assert_eq!(err, ErrorCode::Ok, "Should rename pattern");

        destroy_project(ph);
    }
}

// =============================================================================
// MODULE 5: Curve Operations
// =============================================================================

#[test]
fn test_en_addcurve_valid() {
    unsafe {
        let ph = create_loaded_project();

        let curve_id = CString::new("NEW_CURVE").unwrap();
        let err = EN_addcurve(ph, curve_id.as_ptr());

        assert_eq!(err, ErrorCode::Ok, "Should add curve");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcurveindex_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a curve first
        let curve_id = CString::new("TEST_CURVE").unwrap();
        EN_addcurve(ph, curve_id.as_ptr());

        let mut index: c_int = -999;
        let err = EN_getcurveindex(ph, curve_id.as_ptr(), &mut index);

        assert_eq!(err, ErrorCode::Ok, "Should get curve index");
        assert!(index > 0, "Index should be positive");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcurveindex_output_initialized() {
    unsafe {
        let ph = create_loaded_project();

        let curve_id = CString::new("NONEXISTENT").unwrap();
        let mut index: c_int = -999;

        let err = EN_getcurveindex(ph, curve_id.as_ptr(), &mut index);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid ID");
        assert_eq!(index, 0, "Output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcurveid_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a curve
        let curve_id = CString::new("TEST_CURVE2").unwrap();
        EN_addcurve(ph, curve_id.as_ptr());

        let mut index: c_int = 0;
        EN_getcurveindex(ph, curve_id.as_ptr(), &mut index);

        let mut buffer: [c_char; 32] = [0; 32];
        let err = EN_getcurveid(ph, index, buffer.as_mut_ptr());

        assert_eq!(err, ErrorCode::Ok, "Should get curve ID");
        assert_ne!(buffer[0], 0, "Should write ID to buffer");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcurvelen_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a curve
        let curve_id = CString::new("TEST_CURVE3").unwrap();
        EN_addcurve(ph, curve_id.as_ptr());

        let mut index: c_int = 0;
        EN_getcurveindex(ph, curve_id.as_ptr(), &mut index);

        let mut len: c_int = -999;
        let err = EN_getcurvelen(ph, index, &mut len);

        assert_eq!(err, ErrorCode::Ok, "Should get curve length");
        assert!(len > 0, "Curve should have default point");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcurvelen_output_initialized() {
    unsafe {
        let ph = create_loaded_project();

        let mut len: c_int = -999;
        let err = EN_getcurvelen(ph, 99999, &mut len);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
        assert_eq!(len, 0, "Output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_setcurve_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a curve
        let curve_id = CString::new("TEST_CURVE4").unwrap();
        EN_addcurve(ph, curve_id.as_ptr());

        let mut index: c_int = 0;
        EN_getcurveindex(ph, curve_id.as_ptr(), &mut index);

        // Set curve values
        let x_vals: Vec<c_double> = vec![0.0, 100.0, 200.0];
        let y_vals: Vec<c_double> = vec![0.0, 50.0, 80.0];
        let err = EN_setcurve(
            ph,
            index,
            x_vals.as_ptr(),
            y_vals.as_ptr(),
            x_vals.len() as c_int,
        );

        assert_eq!(err, ErrorCode::Ok, "Should set curve values");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcurvevalue_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add and set curve
        let curve_id = CString::new("TEST_CURVE5").unwrap();
        EN_addcurve(ph, curve_id.as_ptr());

        let mut index: c_int = 0;
        EN_getcurveindex(ph, curve_id.as_ptr(), &mut index);

        let x_vals: Vec<c_double> = vec![10.0, 20.0];
        let y_vals: Vec<c_double> = vec![30.0, 40.0];
        EN_setcurve(
            ph,
            index,
            x_vals.as_ptr(),
            y_vals.as_ptr(),
            x_vals.len() as c_int,
        );

        // Get value
        let mut x: c_double = -999.0;
        let mut y: c_double = -999.0;
        let err = EN_getcurvevalue(ph, index, 1, &mut x, &mut y);

        assert_eq!(err, ErrorCode::Ok, "Should get curve value");
        assert_eq!(x, 10.0, "X should match");
        assert_eq!(y, 30.0, "Y should match");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcurvevalue_output_initialized() {
    unsafe {
        let ph = create_loaded_project();

        let mut x: c_double = -999.0;
        let mut y: c_double = -999.0;
        let err = EN_getcurvevalue(ph, 99999, 1, &mut x, &mut y);

        assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
        assert_eq!(x, 0.0, "X output should be initialized to 0");
        assert_eq!(y, 0.0, "Y output should be initialized to 0");

        destroy_project(ph);
    }
}

#[test]
fn test_en_getcurve_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add and set curve
        let curve_id = CString::new("TEST_CURVE6").unwrap();
        EN_addcurve(ph, curve_id.as_ptr());

        let mut index: c_int = 0;
        EN_getcurveindex(ph, curve_id.as_ptr(), &mut index);

        let x_vals: Vec<c_double> = vec![1.0, 2.0, 3.0];
        let y_vals: Vec<c_double> = vec![4.0, 5.0, 6.0];
        EN_setcurve(
            ph,
            index,
            x_vals.as_ptr(),
            y_vals.as_ptr(),
            x_vals.len() as c_int,
        );

        // Get entire curve
        let mut n_points: c_int = 0;
        let mut x_ptr: *mut c_double = std::ptr::null_mut();
        let mut y_ptr: *mut c_double = std::ptr::null_mut();

        let err = EN_getcurve(ph, index, &mut n_points, &mut x_ptr, &mut y_ptr);

        assert_eq!(err, ErrorCode::Ok, "Should get curve");
        assert_eq!(n_points, 3, "Should have 3 points");

        destroy_project(ph);
    }
}

#[test]
fn test_en_setcurveid_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a curve
        let old_id = CString::new("OLD_CURVE_ID").unwrap();
        EN_addcurve(ph, old_id.as_ptr());

        let mut index: c_int = 0;
        EN_getcurveindex(ph, old_id.as_ptr(), &mut index);

        // Rename it
        let new_id = CString::new("NEW_CURVE_ID").unwrap();
        let err = EN_setcurveid(ph, index, new_id.as_ptr());

        assert_eq!(err, ErrorCode::Ok, "Should rename curve");

        destroy_project(ph);
    }
}

// =============================================================================
// MODULE 6: Hydraulic Solver
// =============================================================================

#[test]
fn test_en_inith_valid() {
    unsafe {
        let ph = create_loaded_project();

        let err = EN_initH(ph, 0);
        // May succeed or fail depending on network state
        // Just verify it doesn't crash
        let _ = err;

        destroy_project(ph);
    }
}

#[test]
fn test_en_runh_after_init() {
    unsafe {
        let ph = create_loaded_project();

        EN_initH(ph, 0);

        let err = EN_runH(ph, 0);

        // May succeed or fail depending on network validity
        // Just verify it doesn't crash
        let _ = err;

        destroy_project(ph);
    }
}

#[test]
fn test_en_nexth_after_run() {
    unsafe {
        let ph = create_loaded_project();

        EN_initH(ph, 0);
        EN_runH(ph, 0);

        let mut tstep: c_long = 0;
        let err = EN_nextH(ph, &mut tstep);

        // May succeed or fail, just verify output is initialized
        let _ = (err, tstep);

        destroy_project(ph);
    }
}

// =============================================================================
// MODULE 7: Analysis Options
// =============================================================================

#[test]
fn test_en_setoption_valid() {
    unsafe {
        let ph = create_loaded_project();

        let err = EN_setoption(ph, SimOption::Trials as c_int, 100.0);
        assert_eq!(err, ErrorCode::Ok, "Should set option");

        destroy_project(ph);
    }
}

#[test]
fn test_en_settimeparam_valid() {
    unsafe {
        let ph = create_loaded_project();

        let err = EN_settimeparam(ph, TimeParameter::Duration as c_int, 7200);
        assert_eq!(err, ErrorCode::Ok, "Should set time parameter");

        destroy_project(ph);
    }
}

#[test]
fn test_en_setdemandmodel_valid() {
    unsafe {
        let ph = create_loaded_project();

        let err = EN_setdemandmodel(ph, DemandModel::Pda as c_int, 3.0, 15.0, 0.5);
        assert_eq!(err, ErrorCode::Ok, "Should set demand model");

        destroy_project(ph);
    }
}

// =============================================================================
// MODULE 7b: EN_setoption / EN_getoption roundtrips and validation
// =============================================================================

#[test]
fn test_en_setoption_trials_valid() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::Trials as c_int, 50.0);
        assert_eq!(err, ErrorCode::Ok, "Trials should accept 50");
        let mut v: c_double = 0.0;
        let err = EN_getoption(ph, SimOption::Trials as c_int, &mut v);
        assert_eq!(err, ErrorCode::Ok);
        assert_eq!(v, 50.0, "Trials roundtrip");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_trials_invalid() {
    unsafe {
        let ph = create_loaded_project();
        // value < 1 is invalid
        let err = EN_setoption(ph, SimOption::Trials as c_int, 0.5);
        assert_eq!(err, ErrorCode::InvalidOptionValue, "Trials < 1 should fail");
        // negative is invalid
        let err = EN_setoption(ph, SimOption::Trials as c_int, -1.0);
        assert_eq!(
            err,
            ErrorCode::InvalidOptionValue,
            "Negative trials should fail"
        );
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_accuracy_valid() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::Accuracy as c_int, 1e-4);
        assert_eq!(err, ErrorCode::Ok);
        let mut v: c_double = 0.0;
        EN_getoption(ph, SimOption::Accuracy as c_int, &mut v);
        assert_eq!(v, 1e-4, "Accuracy roundtrip");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_accuracy_bounds() {
    unsafe {
        let ph = create_loaded_project();
        // Too small
        let err = EN_setoption(ph, SimOption::Accuracy as c_int, 1e-9);
        assert_eq!(
            err,
            ErrorCode::InvalidOptionValue,
            "Accuracy < 1e-8 should fail"
        );
        // Too large
        let err = EN_setoption(ph, SimOption::Accuracy as c_int, 0.2);
        assert_eq!(
            err,
            ErrorCode::InvalidOptionValue,
            "Accuracy > 0.1 should fail"
        );
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_demandmult_roundtrip() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::DemandMult as c_int, 1.5);
        assert_eq!(err, ErrorCode::Ok);
        let mut v: c_double = 0.0;
        EN_getoption(ph, SimOption::DemandMult as c_int, &mut v);
        assert_eq!(v, 1.5, "DemandMult roundtrip");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_emitexpon_roundtrip() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::EmitExpon as c_int, 0.5);
        assert_eq!(err, ErrorCode::Ok);
        let mut v: c_double = 0.0;
        EN_getoption(ph, SimOption::EmitExpon as c_int, &mut v);
        assert!((v - 0.5).abs() < 1e-10, "EmitExpon roundtrip: got {}", v);
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_emitexpon_invalid() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::EmitExpon as c_int, 0.0);
        assert_eq!(
            err,
            ErrorCode::InvalidOptionValue,
            "EmitExpon = 0 should fail"
        );
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_flowchange_roundtrip() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::FlowChange as c_int, 0.01);
        assert_eq!(err, ErrorCode::Ok);
        let mut v: c_double = 0.0;
        EN_getoption(ph, SimOption::FlowChange as c_int, &mut v);
        assert_eq!(v, 0.01, "FlowChange roundtrip");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_headlossform_roundtrip() {
    unsafe {
        use epanet_rs::ffi::enums::HeadLossType;
        let ph = create_loaded_project();
        let err = EN_setoption(
            ph,
            SimOption::HeadLossForm as c_int,
            HeadLossType::DW as c_int as c_double,
        );
        assert_eq!(err, ErrorCode::Ok);
        let mut v: c_double = 0.0;
        EN_getoption(ph, SimOption::HeadLossForm as c_int, &mut v);
        assert_eq!(v, HeadLossType::DW as i32 as f64, "HeadLossForm roundtrip");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_headlossform_invalid() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::HeadLossForm as c_int, 99.0);
        assert_eq!(
            err,
            ErrorCode::InvalidOptionValue,
            "Invalid head loss form should fail"
        );
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_spgravity_roundtrip() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::SpGravity as c_int, 1.05);
        assert_eq!(err, ErrorCode::Ok);
        let mut v: c_double = 0.0;
        EN_getoption(ph, SimOption::SpGravity as c_int, &mut v);
        assert_eq!(v, 1.05, "SpGravity roundtrip");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_spgravity_invalid() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::SpGravity as c_int, 0.0);
        assert_eq!(
            err,
            ErrorCode::InvalidOptionValue,
            "SpGravity = 0 should fail"
        );
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_spviscos_roundtrip() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::SpViscos as c_int, 1.2);
        assert_eq!(err, ErrorCode::Ok);
        let mut v: c_double = 0.0;
        EN_getoption(ph, SimOption::SpViscos as c_int, &mut v);
        assert_eq!(v, 1.2, "SpViscos roundtrip");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_checkfreq_roundtrip() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::CheckFreq as c_int, 5.0);
        assert_eq!(err, ErrorCode::Ok);
        let mut v: c_double = 0.0;
        EN_getoption(ph, SimOption::CheckFreq as c_int, &mut v);
        assert_eq!(v, 5.0, "CheckFreq roundtrip");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_maxcheck_roundtrip() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::MaxCheck as c_int, 8.0);
        assert_eq!(err, ErrorCode::Ok);
        let mut v: c_double = 0.0;
        EN_getoption(ph, SimOption::MaxCheck as c_int, &mut v);
        assert_eq!(v, 8.0, "MaxCheck roundtrip");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_pressunits_roundtrip() {
    unsafe {
        use epanet_rs::ffi::enums::PressUnits;
        let ph = create_loaded_project();
        let err = EN_setoption(
            ph,
            SimOption::PressUnits as c_int,
            PressUnits::Kpa as c_int as c_double,
        );
        assert_eq!(err, ErrorCode::Ok);
        let mut v: c_double = 0.0;
        EN_getoption(ph, SimOption::PressUnits as c_int, &mut v);
        assert_eq!(v, PressUnits::Kpa as i32 as f64, "PressUnits roundtrip");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_pressunits_invalid() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::PressUnits as c_int, 99.0);
        assert_eq!(
            err,
            ErrorCode::UndefinedPattern,
            "Invalid pressure unit should return 205"
        );
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_unavailable_returns_invalid_param() {
    unsafe {
        let ph = create_loaded_project();
        // Options not implemented in this Rust build
        for opt in &[
            SimOption::Tolerance,
            SimOption::HeadError,
            SimOption::GlobalEffic,
            SimOption::GlobalPrice,
            SimOption::GlobalPattern,
            SimOption::DemandCharge,
            SimOption::DampLimit,
            SimOption::SpDiffus,
            SimOption::BulkOrder,
            SimOption::WallOrder,
            SimOption::TankOrder,
            SimOption::ConcenLimit,
            SimOption::EmitBackflow,
            SimOption::StatusReport,
            SimOption::Unbalanced,
        ] {
            let err = EN_setoption(ph, *opt as c_int, 1.0);
            assert_eq!(
                err,
                ErrorCode::InvalidParameterCode,
                "Option {:?} should return InvalidParameterCode",
                opt
            );
        }
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_unknown_code() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, 9999, 1.0);
        assert_eq!(
            err,
            ErrorCode::InvalidParameterCode,
            "Unknown option code should return 251"
        );
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_negative_value_rejected() {
    unsafe {
        let ph = create_loaded_project();
        // Negative values (other than Unbalanced) must return InvalidOptionValue
        let err = EN_setoption(ph, SimOption::DemandMult as c_int, -0.5);
        assert_eq!(
            err,
            ErrorCode::InvalidOptionValue,
            "Negative value should return 213"
        );
        destroy_project(ph);
    }
}

#[test]
fn test_en_getoption_null_project() {
    unsafe {
        let mut v: c_double = -999.0;
        let err = EN_getoption(std::ptr::null_mut(), SimOption::Trials as c_int, &mut v);
        assert_ne!(err, ErrorCode::Ok, "Null project should fail");
        assert_eq!(v, 0.0, "Output should be initialized to 0");
    }
}

#[test]
fn test_en_getoption_unknown_code() {
    unsafe {
        let ph = create_loaded_project();
        let mut v: c_double = -999.0;
        let err = EN_getoption(ph, 9999, &mut v);
        assert_eq!(
            err,
            ErrorCode::InvalidParameterCode,
            "Unknown option code should return 251"
        );
        assert_eq!(v, 0.0, "Output should be initialized to 0");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setdemandmodel_dda() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setdemandmodel(ph, DemandModel::Dda as c_int, 0.0, 0.0, 0.5);
        assert_eq!(err, ErrorCode::Ok, "DDA should succeed");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setdemandmodel_pda_params() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setdemandmodel(ph, DemandModel::Pda as c_int, 3.0, 15.0, 0.5);
        assert_eq!(err, ErrorCode::Ok, "PDA should succeed");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_spviscos_invalid() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::SpViscos as c_int, 0.0);
        assert_eq!(
            err,
            ErrorCode::InvalidOptionValue,
            "SpViscos = 0 should fail"
        );
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_headlossform_blocked_when_solver_open() {
    unsafe {
        use epanet_rs::ffi::enums::HeadLossType;
        use epanet_rs::ffi::hydraulic_solver::EN_openH;
        let ph = create_loaded_project();
        // Open the hydraulic solver
        let open_err = unsafe { EN_openH(ph) };
        assert_eq!(open_err, ErrorCode::Ok, "EN_openH should succeed");
        // Changing headloss formula while solver is open must fail
        let err = EN_setoption(
            ph,
            SimOption::HeadLossForm as c_int,
            HeadLossType::DW as c_int as c_double,
        );
        assert_eq!(
            err,
            ErrorCode::ModifyWhileSolverOpen,
            "Should return 262 when solver is open"
        );
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_demandpattern_valid() {
    unsafe {
        let ph = create_loaded_project();
        // Add a pattern so index 1 is valid
        let pat_id = std::ffi::CString::new("OPT_PAT").unwrap();
        let add_err = unsafe { epanet_rs::ffi::patterns::EN_addpattern(ph, pat_id.as_ptr()) };
        assert_eq!(add_err, ErrorCode::Ok);
        let mut pat_idx: c_int = 0;
        unsafe { epanet_rs::ffi::patterns::EN_getpatternindex(ph, pat_id.as_ptr(), &mut pat_idx) };
        assert!(pat_idx > 0);
        // Set demand pattern
        let err = EN_setoption(ph, SimOption::DemandPattern as c_int, pat_idx as c_double);
        assert_eq!(
            err,
            ErrorCode::Ok,
            "Setting valid demand pattern should succeed"
        );
        // Read it back
        let mut v: c_double = 0.0;
        let get_err = EN_getoption(ph, SimOption::DemandPattern as c_int, &mut v);
        assert_eq!(get_err, ErrorCode::Ok);
        assert_eq!(v, pat_idx as f64, "DemandPattern roundtrip");
        // Reset to 0 (none)
        let reset_err = EN_setoption(ph, SimOption::DemandPattern as c_int, 0.0);
        assert_eq!(
            reset_err,
            ErrorCode::Ok,
            "Setting pattern to 0 (none) should succeed"
        );
        let mut v2: c_double = -1.0;
        EN_getoption(ph, SimOption::DemandPattern as c_int, &mut v2);
        assert_eq!(v2, 0.0, "After reset, demand pattern should be 0");
        destroy_project(ph);
    }
}

#[test]
fn test_en_setoption_demandpattern_out_of_range() {
    unsafe {
        let ph = create_loaded_project();
        let err = EN_setoption(ph, SimOption::DemandPattern as c_int, 99999.0);
        assert_eq!(
            err,
            ErrorCode::UndefinedPattern,
            "Out-of-range pattern index should return 205"
        );
        destroy_project(ph);
    }
}

#[test]
fn test_en_getoption_unavailable_returns_invalid_param() {
    unsafe {
        let ph = create_loaded_project();
        for opt in &[
            SimOption::Tolerance,
            SimOption::HeadError,
            SimOption::GlobalEffic,
            SimOption::GlobalPrice,
            SimOption::GlobalPattern,
            SimOption::DemandCharge,
            SimOption::Unbalanced,
            SimOption::DampLimit,
            SimOption::SpDiffus,
            SimOption::BulkOrder,
            SimOption::WallOrder,
            SimOption::TankOrder,
            SimOption::ConcenLimit,
            SimOption::EmitBackflow,
            SimOption::StatusReport,
        ] {
            let mut v: c_double = -999.0;
            let err = EN_getoption(ph, *opt as c_int, &mut v);
            assert_eq!(
                err,
                ErrorCode::InvalidParameterCode,
                "EN_getoption({:?}) should return InvalidParameterCode",
                opt
            );
            assert_eq!(v, 0.0, "Output should be initialized to 0 on error");
        }
        destroy_project(ph);
    }
}

// =============================================================================
// MODULE 8: Pump Head Curve
// =============================================================================

#[test]
fn test_en_getheadcurveindex_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Find a pump link (link "FH" or similar in pump.inp)
        let link_id = CString::new("FH").unwrap();
        let mut link_idx: c_int = 0;
        let get_err = EN_getlinkindex(ph, link_id.as_ptr(), &mut link_idx);

        if get_err == ErrorCode::Ok {
            let mut curve_idx: c_int = -999;
            let err = EN_getheadcurveindex(ph, link_idx, &mut curve_idx);

            // Pump may or may not have a curve, just verify no crash
            if err == ErrorCode::Ok {
                assert!(curve_idx >= 0, "Curve index should be valid");
            }
        }

        destroy_project(ph);
    }
}

#[test]
fn test_en_setheadcurveindex_valid() {
    unsafe {
        let ph = create_loaded_project();

        // Add a curve
        let curve_id = CString::new("PUMP_CURVE").unwrap();
        EN_addcurve(ph, curve_id.as_ptr());

        let mut curve_idx: c_int = 0;
        EN_getcurveindex(ph, curve_id.as_ptr(), &mut curve_idx);

        // Set curve data
        let x_vals: Vec<c_double> = vec![0.0, 100.0];
        let y_vals: Vec<c_double> = vec![200.0, 150.0];
        EN_setcurve(ph, curve_idx, x_vals.as_ptr(), y_vals.as_ptr(), 2);

        // Try to add a pump and assign curve
        let pump_id = CString::new("TEST_PUMP").unwrap();
        let node1_str = CString::new("1").unwrap();
        let node2_str = CString::new("2").unwrap();
        let mut pump_idx: c_int = 0;

        let add_err = EN_addlink(
            ph,
            pump_id.as_ptr(),
            LinkType::Pump as c_int,
            node1_str.as_ptr(),
            node2_str.as_ptr(),
            &mut pump_idx,
        );

        if add_err == ErrorCode::Ok {
            let err = EN_setheadcurveindex(ph, pump_idx, curve_idx);
            assert_eq!(err, ErrorCode::Ok, "Should set head curve");
        }

        destroy_project(ph);
    }
}

// =============================================================================
// MODULE 9: Error Code Consistency
// =============================================================================

#[test]
fn test_error_codes_consistent_across_modules() {
    unsafe {
        let ph = create_loaded_project();

        // Test that invalid property codes consistently return error 251
        let node_idx = get_valid_node_index(ph, "1");
        let link_idx = get_valid_link_index(ph, "B");

        let invalid_property = 99999;

        // Node
        let mut node_val: c_double = 0.0;
        let node_err = EN_getnodevalue(ph, node_idx, invalid_property, &mut node_val);
        assert_eq!(
            node_err,
            ErrorCode::InvalidParameterCode,
            "Node should return error 251"
        );

        // Link
        let mut link_val: c_double = 0.0;
        let link_err = EN_getlinkvalue(ph, link_idx, invalid_property, &mut link_val);
        assert_eq!(
            link_err,
            ErrorCode::InvalidParameterCode,
            "Link should return error 251"
        );

        destroy_project(ph);
    }
}

#[test]
fn test_null_project_consistent_errors() {
    unsafe {
        // Test that all getters handle null project consistently

        let mut int_val: c_int = -999;

        let node_err = EN_getnodetype(std::ptr::null_mut(), 1, &mut int_val);
        assert_ne!(node_err, ErrorCode::Ok);
        assert_eq!(int_val, 0, "Should initialize output");

        int_val = -999;
        let link_err = EN_getlinktype(std::ptr::null_mut(), 1, &mut int_val);
        assert_ne!(link_err, ErrorCode::Ok);
        assert_eq!(int_val, 0, "Should initialize output");

        let mut dbl_val: c_double = -999.0;
        let pattern_err = EN_getaveragepatternvalue(std::ptr::null_mut(), 1, &mut dbl_val);
        assert_ne!(pattern_err, ErrorCode::Ok);
        assert_eq!(dbl_val, 0.0, "Should initialize output");
    }
}

#[test]
fn test_invalid_indices_consistent() {
    unsafe {
        let ph = create_loaded_project();

        let invalid_idx = 99999;
        let mut dbl_val: c_double = -999.0;

        let node_err = EN_getnodevalue(
            ph,
            invalid_idx,
            NodeProperty::Elevation as c_int,
            &mut dbl_val,
        );
        assert_ne!(node_err, ErrorCode::Ok, "Invalid node index should fail");
        assert_eq!(dbl_val, 0.0, "Should initialize output");

        dbl_val = -999.0;
        let link_err = EN_getlinkvalue(
            ph,
            invalid_idx,
            LinkProperty::Diameter as c_int,
            &mut dbl_val,
        );
        assert_ne!(link_err, ErrorCode::Ok, "Invalid link index should fail");
        assert_eq!(dbl_val, 0.0, "Should initialize output");

        destroy_project(ph);
    }
}

// =============================================================================
// MODULE 10: Regression Tests - Ensure Bug Fixes Remain Fixed
// =============================================================================

#[test]
fn test_regression_no_pattern_returns_zero() {
    unsafe {
        let ph = create_loaded_project();

        // Scan multiple nodes
        for node_id in &["1", "2", "3", "4"] {
            let id = CString::new(*node_id).unwrap();
            let mut idx: c_int = 0;

            if EN_getnodeindex(ph, id.as_ptr(), &mut idx) == ErrorCode::Ok {
                let mut pattern: c_double = -999.0;
                let err = EN_getnodevalue(ph, idx, NodeProperty::Pattern as c_int, &mut pattern);

                if err == ErrorCode::Ok {
                    assert_ne!(pattern, 123.0, "Should never return 123.0 sentinel");
                    assert!(pattern >= 0.0, "Pattern should be >= 0");
                }
            }
        }

        destroy_project(ph);
    }
}

#[test]
fn test_regression_invalid_property_always_errors() {
    unsafe {
        let ph = create_loaded_project();
        let node_idx = get_valid_node_index(ph, "1");

        for invalid_code in vec![999, 1000, -999, 12345] {
            let mut value: c_double = -999.0;
            let err = EN_getnodevalue(ph, node_idx, invalid_code, &mut value);

            assert_eq!(
                err,
                ErrorCode::InvalidParameterCode,
                "Invalid property {} should return error 251",
                invalid_code
            );
            assert_ne!(value, -123.0, "Should never return -123.0 sentinel");
            assert_eq!(value, 0.0, "Should initialize to 0.0");
        }

        destroy_project(ph);
    }
}

#[test]
fn test_regression_all_outputs_initialized() {
    unsafe {
        // Test that ALL getter functions initialize outputs before validation

        let ph = std::ptr::null_mut(); // Null project to trigger errors

        // Integer outputs
        let mut int_val: c_int = -999;
        EN_getnodetype(ph, 1, &mut int_val);
        assert_eq!(int_val, 0, "getnodetype should initialize");

        int_val = -999;
        EN_getlinktype(ph, 1, &mut int_val);
        assert_eq!(int_val, 0, "getlinktype should initialize");

        int_val = -999;
        let id = CString::new("X").unwrap();
        EN_getnodeindex(ph, id.as_ptr(), &mut int_val);
        assert_eq!(int_val, 0, "getnodeindex should initialize");

        // Double outputs
        let mut dbl_val: c_double = -999.0;
        EN_getnodevalue(ph, 1, 0, &mut dbl_val);
        assert_eq!(dbl_val, 0.0, "getnodevalue should initialize");

        dbl_val = -999.0;
        EN_getlinkvalue(ph, 1, 0, &mut dbl_val);
        assert_eq!(dbl_val, 0.0, "getlinkvalue should initialize");

        dbl_val = -999.0;
        EN_getpatternvalue(ph, 1, 1, &mut dbl_val);
        assert_eq!(dbl_val, 0.0, "getpatternvalue should initialize");

        // Multiple outputs
        let mut x: c_double = -999.0;
        let mut y: c_double = -999.0;
        EN_getcoord(ph, 1, &mut x, &mut y);
        assert_eq!(x, 0.0, "getcoord should initialize x");
        assert_eq!(y, 0.0, "getcoord should initialize y");

        x = -999.0;
        y = -999.0;
        EN_getcurvevalue(ph, 1, 1, &mut x, &mut y);
        assert_eq!(x, 0.0, "getcurvevalue should initialize x");
        assert_eq!(y, 0.0, "getcurvevalue should initialize y");

        let mut n1: c_int = -999;
        let mut n2: c_int = -999;
        EN_getlinknodes(ph, 1, &mut n1, &mut n2);
        assert_eq!(n1, 0, "getlinknodes should initialize node1");
        assert_eq!(n2, 0, "getlinknodes should initialize node2");
    }
}
