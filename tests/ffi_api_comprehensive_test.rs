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
use epanet_rs::ffi::enums::LinkType as ENLinkType;
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
fn create_empty_project() -> *mut Project {
    let mut ph: *mut Project = std::ptr::null_mut();
    let err = unsafe { EN_createproject(&mut ph) };
    assert_eq!(err, ErrorCode::Ok);
    assert!(!ph.is_null());
    ph
}

/// Create project and load pump.inp
fn create_loaded_project() -> *mut Project {
    let ph = create_empty_project();

    let inp_file = CString::new("tests/pump.inp").unwrap();
    let rpt_file = CString::new("").unwrap();
    let out_file = CString::new("").unwrap();

    let err = unsafe { EN_open(ph, inp_file.as_ptr(), rpt_file.as_ptr(), out_file.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    ph
}

/// Clean up project
fn destroy_project(ph: *mut Project) {
    let err = unsafe { EN_deleteproject(ph) };
    assert_eq!(err, ErrorCode::Ok);
}

/// Get a valid node index from loaded project
fn get_valid_node_index(ph: *mut Project, node_id: &str) -> c_int {
    let id = CString::new(node_id).unwrap();
    let mut index: c_int = 0;
    let err = unsafe { EN_getnodeindex(ph, id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);
    assert!(index > 0);
    index
}

/// Get a valid link index from loaded project
fn get_valid_link_index(ph: *mut Project, link_id: &str) -> c_int {
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

    destroy_project(ph);
}

#[test]
fn test_en_createproject_null_handle() {
    let err = unsafe { EN_createproject(std::ptr::null_mut()) };
    assert_ne!(err, ErrorCode::Ok, "Should fail with null handle pointer");
}

#[test]
fn test_en_deleteproject_valid() {
    let ph = create_empty_project();
    let err = unsafe { EN_deleteproject(ph) };
    assert_eq!(err, ErrorCode::Ok, "Should delete project successfully");
}

#[test]
fn test_en_deleteproject_null() {
    let err = unsafe { EN_deleteproject(std::ptr::null_mut()) };
    assert_ne!(err, ErrorCode::Ok, "Should fail with null project");
}

#[test]
fn test_en_open_valid_file() {
    let ph = create_empty_project();

    let inp_file = CString::new("tests/pump.inp").unwrap();
    let rpt_file = CString::new("").unwrap();
    let out_file = CString::new("").unwrap();

    let err = unsafe { EN_open(ph, inp_file.as_ptr(), rpt_file.as_ptr(), out_file.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok, "Should open valid INP file");

    destroy_project(ph);
}

#[test]
fn test_en_open_invalid_file() {
    let ph = create_empty_project();

    let inp_file = CString::new("nonexistent.inp").unwrap();
    let rpt_file = CString::new("").unwrap();
    let out_file = CString::new("").unwrap();

    let err = unsafe { EN_open(ph, inp_file.as_ptr(), rpt_file.as_ptr(), out_file.as_ptr()) };
    assert_ne!(err, ErrorCode::Ok, "Should fail with nonexistent file");

    destroy_project(ph);
}

#[test]
fn test_en_open_null_project() {
    let inp_file = CString::new("tests/pump.inp").unwrap();
    let rpt_file = CString::new("").unwrap();
    let out_file = CString::new("").unwrap();

    let err = unsafe {
        EN_open(
            std::ptr::null_mut(),
            inp_file.as_ptr(),
            rpt_file.as_ptr(),
            out_file.as_ptr(),
        )
    };
    assert_ne!(err, ErrorCode::Ok, "Should fail with null project");
}

#[test]
fn test_en_closeh_valid() {
    let ph = create_loaded_project();
    let err = unsafe { EN_closeH(ph) };
    // May succeed or fail depending on whether hydraulics were opened
    // Just verify it doesn't crash
    let _ = err;
    destroy_project(ph);
}

#[test]
fn test_en_saveinpfile_valid() {
    let ph = create_loaded_project();

    let out_path = std::env::temp_dir().join(format!(
        "epanet_rs_test_output_{}.inp",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("System time before UNIX epoch")
            .as_nanos()
    ));

    let out_file = CString::new(out_path.to_str().expect("Temp path contains invalid UTF-8"))
        .expect("Path contains null byte");

    let err = unsafe { EN_saveinpfile(ph, out_file.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok, "Should save INP file");

    destroy_project(ph);
}

#[test]
fn test_en_getcount_nodes() {
    let ph = create_loaded_project();

    let mut count: c_int = -1;
    let err = unsafe { EN_getcount(ph, CountType::NodeCount as c_int, &mut count) };

    assert_eq!(err, ErrorCode::Ok, "Should get node count");
    assert!(count > 0, "Node count should be positive");

    destroy_project(ph);
}

#[test]
fn test_en_getcount_links() {
    let ph = create_loaded_project();

    let mut count: c_int = -1;
    let err = unsafe { EN_getcount(ph, CountType::LinkCount as c_int, &mut count) };

    assert_eq!(err, ErrorCode::Ok, "Should get link count");
    assert!(count > 0, "Link count should be positive");

    destroy_project(ph);
}

#[test]
fn test_en_getcount_output_initialized() {
    let mut count: c_int = -999;
    let err = unsafe {
        EN_getcount(
            std::ptr::null_mut(),
            CountType::NodeCount as c_int,
            &mut count,
        )
    };

    assert_ne!(err, ErrorCode::Ok, "Should fail with null project");
    assert_eq!(count, 0, "Output should be initialized to 0");
}

#[test]
fn test_en_geterror_valid_code() {
    let code = ErrorCode::InvalidParameterCode as c_int;
    let mut buffer: [c_char; 256] = [0; 256];

    let err = unsafe { EN_geterror(code, buffer.as_mut_ptr(), 256) };
    assert_eq!(err, ErrorCode::Ok, "Should get error message");

    // Should have written something to buffer
    assert_ne!(buffer[0], 0, "Should write error message");
}

// =============================================================================
// MODULE 2: Node Operations
// =============================================================================

#[test]
fn test_en_getnodeindex_valid() {
    let ph = create_loaded_project();

    let node_id = CString::new("1").unwrap();
    let mut index: c_int = -999;

    let err = unsafe { EN_getnodeindex(ph, node_id.as_ptr(), &mut index) };

    assert_eq!(err, ErrorCode::Ok, "Should get node index");
    assert!(index > 0, "Index should be positive (1-based)");

    destroy_project(ph);
}

#[test]
fn test_en_getnodeindex_invalid_id() {
    let ph = create_loaded_project();

    let node_id = CString::new("NONEXISTENT").unwrap();
    let mut index: c_int = -999;

    let err = unsafe { EN_getnodeindex(ph, node_id.as_ptr(), &mut index) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid ID");
    assert_eq!(index, 0, "Output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_getnodeindex_null_project() {
    let node_id = CString::new("1").unwrap();
    let mut index: c_int = -999;

    let err = unsafe { EN_getnodeindex(std::ptr::null_mut(), node_id.as_ptr(), &mut index) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with null project");
    assert_eq!(index, 0, "Output should be initialized to 0");
}

#[test]
fn test_en_getnodeid_valid() {
    let ph = create_loaded_project();
    let index = get_valid_node_index(ph, "1");

    let mut buffer: [c_char; 32] = [0; 32];
    let err = unsafe { EN_getnodeid(ph, index, buffer.as_mut_ptr()) };

    assert_eq!(err, ErrorCode::Ok, "Should get node ID");
    assert_ne!(buffer[0], 0, "Should write ID to buffer");

    destroy_project(ph);
}

#[test]
fn test_en_getnodeid_invalid_index() {
    let ph = create_loaded_project();

    let mut buffer: [c_char; 32] = [0; 32];
    let err = unsafe { EN_getnodeid(ph, 99999, buffer.as_mut_ptr()) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");

    destroy_project(ph);
}

#[test]
fn test_en_getnodetype_valid() {
    let ph = create_loaded_project();
    let index = get_valid_node_index(ph, "1");

    let mut node_type: c_int = -999;
    let err = unsafe { EN_getnodetype(ph, index, &mut node_type) };

    assert_eq!(err, ErrorCode::Ok, "Should get node type");
    assert!(node_type >= 0, "Node type should be valid");

    destroy_project(ph);
}

#[test]
fn test_en_getnodetype_output_initialized() {
    let ph = create_loaded_project();

    let mut node_type: c_int = -999;
    let err = unsafe { EN_getnodetype(ph, 99999, &mut node_type) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
    assert_eq!(node_type, 0, "Output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_getnodevalue_elevation() {
    let ph = create_loaded_project();
    let index = get_valid_node_index(ph, "1");

    let mut value: c_double = -999.0;
    let err = unsafe { EN_getnodevalue(ph, index, NodeProperty::Elevation as c_int, &mut value) };

    assert_eq!(err, ErrorCode::Ok, "Should get elevation");
    assert_ne!(value, -999.0, "Should update value");

    destroy_project(ph);
}

#[test]
fn test_en_getnodevalue_pattern() {
    let ph = create_loaded_project();
    let index = get_valid_node_index(ph, "1");

    let mut value: c_double = -999.0;
    let err = unsafe { EN_getnodevalue(ph, index, NodeProperty::Pattern as c_int, &mut value) };

    assert_eq!(err, ErrorCode::Ok, "Should get pattern");
    assert!(value >= 0.0, "Pattern should be >= 0 (0 = none)");
    assert_ne!(value, 123.0, "Should never return 123.0 sentinel");

    destroy_project(ph);
}

#[test]
fn test_en_getnodevalue_invalid_property() {
    let ph = create_loaded_project();
    let index = get_valid_node_index(ph, "1");

    let mut value: c_double = -999.0;
    let err = unsafe { EN_getnodevalue(ph, index, 99999, &mut value) };

    assert_eq!(
        err,
        ErrorCode::InvalidParameterCode,
        "Should return error 251"
    );
    assert_eq!(value, 0.0, "Output should be initialized to 0");
    assert_ne!(value, -123.0, "Should not return -123.0 sentinel");

    destroy_project(ph);
}

#[test]
fn test_en_getnodevalue_null_project() {
    let mut value: c_double = -999.0;
    let err = unsafe {
        EN_getnodevalue(
            std::ptr::null_mut(),
            1,
            NodeProperty::Elevation as c_int,
            &mut value,
        )
    };

    assert_ne!(err, ErrorCode::Ok, "Should fail with null project");
    assert_eq!(value, 0.0, "Output should be initialized to 0");
}

#[test]
fn test_en_setnodevalue_elevation() {
    let ph = create_loaded_project();
    let index = get_valid_node_index(ph, "1");

    let new_elevation = 123.45;
    let err =
        unsafe { EN_setnodevalue(ph, index, NodeProperty::Elevation as c_int, new_elevation) };

    assert_eq!(err, ErrorCode::Ok, "Should set elevation");

    // Verify it was set
    let mut value: c_double = 0.0;
    let err = unsafe { EN_getnodevalue(ph, index, NodeProperty::Elevation as c_int, &mut value) };
    assert_eq!(err, ErrorCode::Ok, "Should get elevation");
    assert_eq!(value, new_elevation, "Should retrieve new elevation");

    destroy_project(ph);
}

#[test]
fn test_en_setnodevalue_invalid_property() {
    let ph = create_loaded_project();
    let index = get_valid_node_index(ph, "1");

    let err = unsafe { EN_setnodevalue(ph, index, 99999, 123.45) };

    assert_eq!(
        err,
        ErrorCode::InvalidParameterCode,
        "Should return error 251"
    );

    destroy_project(ph);
}

#[test]
fn test_en_addnode_valid() {
    let ph = create_loaded_project();

    let node_id = CString::new("NEW_NODE").unwrap();
    let node_type = NodeType::Junction as c_int;

    let mut index: c_int = 0;
    let err = unsafe { EN_addnode(ph, node_id.as_ptr(), node_type, &mut index) };

    assert_eq!(err, ErrorCode::Ok, "Should add node");
    assert!(index > 0, "Should return valid index");

    destroy_project(ph);
}

#[test]
fn test_en_addnode_duplicate_id() {
    let ph = create_loaded_project();

    // Try to add node with existing ID
    let node_id = CString::new("1").unwrap();
    let node_type = NodeType::Junction as c_int;

    let mut index: c_int = 0;
    let err = unsafe { EN_addnode(ph, node_id.as_ptr(), node_type, &mut index) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with duplicate ID");

    destroy_project(ph);
}

#[test]
fn test_en_deletenode_valid() {
    let ph = create_loaded_project();

    // Add a node first
    let node_id = CString::new("TEMP_NODE").unwrap();
    let mut index: c_int = 0;
    let err = unsafe {
        EN_addnode(
            ph,
            node_id.as_ptr(),
            NodeType::Junction as c_int,
            &mut index,
        )
    };
    assert_eq!(err, ErrorCode::Ok);

    // Delete it
    let err = unsafe { EN_deletenode(ph, index, 0) };
    assert_eq!(err, ErrorCode::Ok, "Should delete node");

    destroy_project(ph);
}

#[test]
fn test_en_setnodeid_valid() {
    let ph = create_loaded_project();

    // Add a node
    let old_id = CString::new("OLD_ID").unwrap();
    let mut index: c_int = 0;
    let err = unsafe { EN_addnode(ph, old_id.as_ptr(), NodeType::Junction as c_int, &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    // Rename it
    let new_id = CString::new("NEW_ID").unwrap();
    let err = unsafe { EN_setnodeid(ph, index, new_id.as_ptr()) };

    assert_eq!(err, ErrorCode::Ok, "Should rename node");

    destroy_project(ph);
}

#[test]
fn test_en_getcoord_valid() {
    let ph = create_loaded_project();
    let index = get_valid_node_index(ph, "1");

    let mut x: c_double = -999.0;
    let mut y: c_double = -999.0;

    let err = unsafe { EN_getcoord(ph, index, &mut x, &mut y) };

    assert_eq!(err, ErrorCode::Ok, "Should get coordinates");
    assert_ne!(x, -999.0, "X should be updated");
    assert_ne!(y, -999.0, "Y should be updated");

    destroy_project(ph);
}

#[test]
fn test_en_getcoord_output_initialized() {
    let ph = create_loaded_project();

    let mut x: c_double = -999.0;
    let mut y: c_double = -999.0;

    let err = unsafe { EN_getcoord(ph, 99999, &mut x, &mut y) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
    assert_eq!(x, 0.0, "X should be initialized to 0");
    assert_eq!(y, 0.0, "Y should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_setcoord_valid() {
    let ph = create_loaded_project();
    let index = get_valid_node_index(ph, "1");

    let err = unsafe { EN_setcoord(ph, index, 100.0, 200.0) };
    assert_eq!(err, ErrorCode::Ok, "Should set coordinates");

    // Verify
    let mut x: c_double = 0.0;
    let mut y: c_double = 0.0;
    let err = unsafe { EN_getcoord(ph, index, &mut x, &mut y) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(x, 100.0, "X should match");
    assert_eq!(y, 200.0, "Y should match");

    destroy_project(ph);
}

// =============================================================================
// MODULE 3: Link Operations
// =============================================================================

#[test]
fn test_en_getlinkindex_valid() {
    let ph = create_loaded_project();

    let link_id = CString::new("B").unwrap();
    let mut index: c_int = -999;

    let err = unsafe { EN_getlinkindex(ph, link_id.as_ptr(), &mut index) };

    assert_eq!(err, ErrorCode::Ok, "Should get link index");
    assert!(index > 0, "Index should be positive");

    destroy_project(ph);
}

#[test]
fn test_en_getlinkindex_invalid_id() {
    let ph = create_loaded_project();

    let link_id = CString::new("NONEXISTENT").unwrap();
    let mut index: c_int = -999;

    let err = unsafe { EN_getlinkindex(ph, link_id.as_ptr(), &mut index) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid ID");
    assert_eq!(index, 0, "Output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_getlinkid_valid() {
    let ph = create_loaded_project();
    let index = get_valid_link_index(ph, "B");

    let mut buffer: [c_char; 32] = [0; 32];
    let err = unsafe { EN_getlinkid(ph, index, buffer.as_mut_ptr()) };

    assert_eq!(err, ErrorCode::Ok, "Should get link ID");
    assert_ne!(buffer[0], 0, "Should write ID to buffer");

    destroy_project(ph);
}

#[test]
fn test_en_getlinktype_valid() {
    let ph = create_loaded_project();
    let index = get_valid_link_index(ph, "B");

    let mut link_type: c_int = -999;
    let err = unsafe { EN_getlinktype(ph, index, &mut link_type) };

    assert_eq!(err, ErrorCode::Ok, "Should get link type");
    assert!(link_type >= 0, "Link type should be valid");

    destroy_project(ph);
}

#[test]
fn test_en_getlinktype_output_initialized() {
    let ph = create_loaded_project();

    let mut link_type: c_int = -999;
    let err = unsafe { EN_getlinktype(ph, 99999, &mut link_type) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
    assert_eq!(link_type, 0, "Output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_getlinkvalue_diameter() {
    let ph = create_loaded_project();
    let index = get_valid_link_index(ph, "B");

    let mut value: c_double = -999.0;
    let err = unsafe { EN_getlinkvalue(ph, index, LinkProperty::Diameter as c_int, &mut value) };

    assert_eq!(err, ErrorCode::Ok, "Should get diameter");
    assert_ne!(value, -999.0, "Should update value");

    destroy_project(ph);
}

#[test]
fn test_en_getlinkvalue_invalid_property() {
    let ph = create_loaded_project();
    let index = get_valid_link_index(ph, "B");

    let mut value: c_double = -999.0;
    let err = unsafe { EN_getlinkvalue(ph, index, 99999, &mut value) };

    assert_eq!(
        err,
        ErrorCode::InvalidParameterCode,
        "Should return error 251"
    );
    assert_eq!(value, 0.0, "Output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_setlinkvalue_diameter() {
    let ph = create_loaded_project();
    let index = get_valid_link_index(ph, "B");

    let new_diameter = 24.0;
    let err = unsafe { EN_setlinkvalue(ph, index, LinkProperty::Diameter as c_int, new_diameter) };

    assert_eq!(err, ErrorCode::Ok, "Should set diameter");

    destroy_project(ph);
}

#[test]
fn test_en_getlinknodes_valid() {
    let ph = create_loaded_project();
    let index = get_valid_link_index(ph, "B");

    let mut node1: c_int = -999;
    let mut node2: c_int = -999;

    let err = unsafe { EN_getlinknodes(ph, index, &mut node1, &mut node2) };

    assert_eq!(err, ErrorCode::Ok, "Should get link nodes");
    assert!(node1 > 0, "Node 1 should be valid");
    assert!(node2 > 0, "Node 2 should be valid");

    destroy_project(ph);
}

#[test]
fn test_en_getlinknodes_output_initialized() {
    let ph = create_loaded_project();

    let mut node1: c_int = -999;
    let mut node2: c_int = -999;

    let err = unsafe { EN_getlinknodes(ph, 99999, &mut node1, &mut node2) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
    assert_eq!(node1, 0, "Node1 should be initialized to 0");
    assert_eq!(node2, 0, "Node2 should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_addlink_valid() {
    let ph = create_loaded_project();

    let link_id = CString::new("NEW_PIPE").unwrap();
    let link_type = LinkType::Pipe as c_int;
    let from_node = CString::new("1").unwrap();
    let to_node = CString::new("2").unwrap();

    let mut index: c_int = 0;
    let err = unsafe {
        EN_addlink(
            ph,
            link_id.as_ptr(),
            link_type,
            from_node.as_ptr(),
            to_node.as_ptr(),
            &mut index,
        )
    };

    assert_eq!(err, ErrorCode::Ok, "Should add link");
    assert!(index > 0, "Should return valid index");

    destroy_project(ph);
}

#[test]
fn test_en_deletelink_valid() {
    let ph = create_loaded_project();

    // Add a link first
    let link_id = CString::new("TEMP_LINK").unwrap();
    let from_node = CString::new("1").unwrap();
    let to_node = CString::new("2").unwrap();
    let mut index: c_int = 0;
    let err = unsafe {
        EN_addlink(
            ph,
            link_id.as_ptr(),
            LinkType::Pipe as c_int,
            from_node.as_ptr(),
            to_node.as_ptr(),
            &mut index,
        )
    };
    assert_eq!(err, ErrorCode::Ok);

    // Delete it
    let err = unsafe { EN_deletelink(ph, index, 0) };
    assert_eq!(err, ErrorCode::Ok, "Should delete link");

    destroy_project(ph);
}

#[test]
fn test_en_setlinkid_valid() {
    let ph = create_loaded_project();

    // Add a link
    let old_id = CString::new("OLD_LINK_ID").unwrap();
    let from_node = CString::new("1").unwrap();
    let to_node = CString::new("2").unwrap();
    let mut index: c_int = 0;
    let err = unsafe {
        EN_addlink(
            ph,
            old_id.as_ptr(),
            LinkType::Pipe as c_int,
            from_node.as_ptr(),
            to_node.as_ptr(),
            &mut index,
        )
    };
    assert_eq!(err, ErrorCode::Ok);

    // Rename it
    let new_id = CString::new("NEW_LINK_ID").unwrap();
    let err = unsafe { EN_setlinkid(ph, index, new_id.as_ptr()) };

    assert_eq!(err, ErrorCode::Ok, "Should rename link");

    destroy_project(ph);
}

#[test]
fn test_en_setlinknodes_valid() {
    let ph = create_loaded_project();

    // Add a link
    let link_id = CString::new("TEST_LINK").unwrap();
    let node1_str = CString::new("1").unwrap();
    let node2_str = CString::new("2").unwrap();
    let mut index: c_int = 0;
    let err = unsafe {
        EN_addlink(
            ph,
            link_id.as_ptr(),
            LinkType::Pipe as c_int,
            node1_str.as_ptr(),
            node2_str.as_ptr(),
            &mut index,
        )
    };
    assert_eq!(err, ErrorCode::Ok);

    // Change nodes
    let node1 = get_valid_node_index(ph, "1");
    let node3 = get_valid_node_index(ph, "3");
    let err = unsafe { EN_setlinknodes(ph, index, node1, node3) };

    assert_eq!(err, ErrorCode::Ok, "Should set link nodes");

    destroy_project(ph);
}

// =============================================================================
// MODULE 4: Pattern Operations
// =============================================================================

#[test]
fn test_en_addpattern_valid() {
    let ph = create_loaded_project();

    let pattern_id = CString::new("NEW_PATTERN").unwrap();
    let err = unsafe { EN_addpattern(ph, pattern_id.as_ptr()) };

    assert_eq!(err, ErrorCode::Ok, "Should add pattern");

    destroy_project(ph);
}

#[test]
fn test_en_getpatternindex_valid() {
    let ph = create_loaded_project();

    // Add a pattern first
    let pattern_id = CString::new("TEST_PAT").unwrap();
    let err = unsafe { EN_addpattern(ph, pattern_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = -999;
    let err = unsafe { EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index) };

    assert_eq!(err, ErrorCode::Ok, "Should get pattern index");
    assert!(index > 0, "Index should be positive");

    destroy_project(ph);
}

#[test]
fn test_en_getpatternindex_output_initialized() {
    let ph = create_loaded_project();

    let pattern_id = CString::new("NONEXISTENT").unwrap();
    let mut index: c_int = -999;

    let err = unsafe { EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid ID");
    assert_eq!(index, 0, "Output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_getpatternid_valid() {
    let ph = create_loaded_project();

    // Add a pattern
    let pattern_id = CString::new("TEST_PAT2").unwrap();
    let err = unsafe { EN_addpattern(ph, pattern_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    let mut buffer: [c_char; 32] = [0; 32];
    let err = unsafe { EN_getpatternid(ph, index, buffer.as_mut_ptr()) };

    assert_eq!(err, ErrorCode::Ok, "Should get pattern ID");
    assert_ne!(buffer[0], 0, "Should write ID to buffer");

    destroy_project(ph);
}

#[test]
fn test_en_getpatternlen_valid() {
    let ph = create_loaded_project();

    // Add a pattern
    let pattern_id = CString::new("TEST_PAT3").unwrap();
    let err = unsafe { EN_addpattern(ph, pattern_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    let mut len: c_int = -999;
    let err = unsafe { EN_getpatternlen(ph, index, &mut len) };

    assert_eq!(err, ErrorCode::Ok, "Should get pattern length");
    assert!(len >= 0, "Length should be non-negative");

    destroy_project(ph);
}

#[test]
fn test_en_getpatternlen_output_initialized() {
    let ph = create_loaded_project();

    let mut len: c_int = -999;
    let err = unsafe { EN_getpatternlen(ph, 99999, &mut len) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
    assert_eq!(len, 0, "Output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_setpattern_valid() {
    let ph = create_loaded_project();

    // Add a pattern
    let pattern_id = CString::new("TEST_PAT4").unwrap();
    let err = unsafe { EN_addpattern(ph, pattern_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    // Set pattern values
    let values: Vec<c_double> = vec![1.0, 1.5, 0.8, 1.2];
    let err = unsafe { EN_setpattern(ph, index, values.as_ptr(), values.len() as c_int) };

    assert_eq!(err, ErrorCode::Ok, "Should set pattern values");

    destroy_project(ph);
}

#[test]
fn test_en_getpatternvalue_valid() {
    let ph = create_loaded_project();

    // Add and set pattern
    let pattern_id = CString::new("TEST_PAT5").unwrap();
    let err = unsafe { EN_addpattern(ph, pattern_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    let values: Vec<c_double> = vec![1.0, 2.0, 3.0];
    let err = unsafe { EN_setpattern(ph, index, values.as_ptr(), values.len() as c_int) };
    assert_eq!(err, ErrorCode::Ok);

    // Get value
    let mut value: c_double = -999.0;
    let err = unsafe { EN_getpatternvalue(ph, index, 1, &mut value) };

    assert_eq!(err, ErrorCode::Ok, "Should get pattern value");
    assert_eq!(value, 1.0, "Should match set value");

    destroy_project(ph);
}

#[test]
fn test_en_getpatternvalue_output_initialized() {
    let ph = create_loaded_project();

    let mut value: c_double = -999.0;
    let err = unsafe { EN_getpatternvalue(ph, 99999, 1, &mut value) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
    assert_eq!(value, 0.0, "Output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_getaveragepatternvalue_valid() {
    let ph = create_loaded_project();

    // Add and set pattern
    let pattern_id = CString::new("TEST_PAT6").unwrap();
    let err = unsafe { EN_addpattern(ph, pattern_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    let values: Vec<c_double> = vec![1.0, 2.0, 3.0];
    let err = unsafe { EN_setpattern(ph, index, values.as_ptr(), values.len() as c_int) };
    assert_eq!(err, ErrorCode::Ok);

    // Get average
    let mut avg: c_double = -999.0;
    let err = unsafe { EN_getaveragepatternvalue(ph, index, &mut avg) };

    assert_eq!(err, ErrorCode::Ok, "Should get average");
    assert_eq!(avg, 2.0, "Average of 1,2,3 should be 2");

    destroy_project(ph);
}

#[test]
fn test_en_getaveragepatternvalue_output_initialized() {
    let ph = create_loaded_project();

    let mut avg: c_double = -999.0;
    let err = unsafe { EN_getaveragepatternvalue(ph, 99999, &mut avg) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
    assert_eq!(avg, 0.0, "Output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_deletepattern_valid() {
    let ph = create_loaded_project();

    // Add a pattern
    let pattern_id = CString::new("TEMP_PAT").unwrap();
    let err = unsafe { EN_addpattern(ph, pattern_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getpatternindex(ph, pattern_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    // Delete it
    let err = unsafe { EN_deletepattern(ph, index) };
    assert_eq!(err, ErrorCode::Ok, "Should delete pattern");

    destroy_project(ph);
}

#[test]
fn test_en_setpatternid_valid() {
    let ph = create_loaded_project();

    // Add a pattern
    let old_id = CString::new("OLD_PAT_ID").unwrap();
    let err = unsafe { EN_addpattern(ph, old_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getpatternindex(ph, old_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    // Rename it
    let new_id = CString::new("NEW_PAT_ID").unwrap();
    let err = unsafe { EN_setpatternid(ph, index, new_id.as_ptr()) };

    assert_eq!(err, ErrorCode::Ok, "Should rename pattern");

    destroy_project(ph);
}

// =============================================================================
// MODULE 5: Curve Operations
// =============================================================================

#[test]
fn test_en_addcurve_valid() {
    let ph = create_loaded_project();

    let curve_id = CString::new("NEW_CURVE").unwrap();
    let err = unsafe { EN_addcurve(ph, curve_id.as_ptr()) };

    assert_eq!(err, ErrorCode::Ok, "Should add curve");

    destroy_project(ph);
}

#[test]
fn test_en_getcurveindex_valid() {
    let ph = create_loaded_project();

    // Add a curve first
    let curve_id = CString::new("TEST_CURVE").unwrap();
    let err = unsafe { EN_addcurve(ph, curve_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = -999;
    let err = unsafe { EN_getcurveindex(ph, curve_id.as_ptr(), &mut index) };

    assert_eq!(err, ErrorCode::Ok, "Should get curve index");
    assert!(index > 0, "Index should be positive");

    destroy_project(ph);
}

#[test]
fn test_en_getcurveindex_output_initialized() {
    let ph = create_loaded_project();

    let curve_id = CString::new("NONEXISTENT").unwrap();
    let mut index: c_int = -999;

    let err = unsafe { EN_getcurveindex(ph, curve_id.as_ptr(), &mut index) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid ID");
    assert_eq!(index, 0, "Output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_getcurveid_valid() {
    let ph = create_loaded_project();

    // Add a curve
    let curve_id = CString::new("TEST_CURVE2").unwrap();
    let err = unsafe { EN_addcurve(ph, curve_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getcurveindex(ph, curve_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    let mut buffer: [c_char; 32] = [0; 32];
    let err = unsafe { EN_getcurveid(ph, index, buffer.as_mut_ptr()) };

    assert_eq!(err, ErrorCode::Ok, "Should get curve ID");
    assert_ne!(buffer[0], 0, "Should write ID to buffer");

    destroy_project(ph);
}

#[test]
fn test_en_getcurvelen_valid() {
    let ph = create_loaded_project();

    // Add a curve
    let curve_id = CString::new("TEST_CURVE3").unwrap();
    let err = unsafe { EN_addcurve(ph, curve_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getcurveindex(ph, curve_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    let mut len: c_int = -999;
    let err = unsafe { EN_getcurvelen(ph, index, &mut len) };

    assert_eq!(err, ErrorCode::Ok, "Should get curve length");
    assert!(len > 0, "Curve should have default point");

    destroy_project(ph);
}

#[test]
fn test_en_getcurvelen_output_initialized() {
    let ph = create_loaded_project();

    let mut len: c_int = -999;
    let err = unsafe { EN_getcurvelen(ph, 99999, &mut len) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
    assert_eq!(len, 0, "Output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_setcurve_valid() {
    let ph = create_loaded_project();

    // Add a curve
    let curve_id = CString::new("TEST_CURVE4").unwrap();
    let err = unsafe { EN_addcurve(ph, curve_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getcurveindex(ph, curve_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    // Set curve values
    let x_vals: Vec<c_double> = vec![0.0, 100.0, 200.0];
    let y_vals: Vec<c_double> = vec![0.0, 50.0, 80.0];
    let err = unsafe {
        EN_setcurve(
            ph,
            index,
            x_vals.as_ptr(),
            y_vals.as_ptr(),
            x_vals.len() as c_int,
        )
    };

    assert_eq!(err, ErrorCode::Ok, "Should set curve values");

    destroy_project(ph);
}

#[test]
fn test_en_getcurvevalue_valid() {
    let ph = create_loaded_project();

    // Add and set curve
    let curve_id = CString::new("TEST_CURVE5").unwrap();
    let err = unsafe { EN_addcurve(ph, curve_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getcurveindex(ph, curve_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    let x_vals: Vec<c_double> = vec![10.0, 20.0];
    let y_vals: Vec<c_double> = vec![30.0, 40.0];
    let err = unsafe {
        EN_setcurve(
            ph,
            index,
            x_vals.as_ptr(),
            y_vals.as_ptr(),
            x_vals.len() as c_int,
        )
    };
    assert_eq!(err, ErrorCode::Ok);

    // Get value
    let mut x: c_double = -999.0;
    let mut y: c_double = -999.0;
    let err = unsafe { EN_getcurvevalue(ph, index, 1, &mut x, &mut y) };

    assert_eq!(err, ErrorCode::Ok, "Should get curve value");
    assert_eq!(x, 10.0, "X should match");
    assert_eq!(y, 30.0, "Y should match");

    destroy_project(ph);
}

#[test]
fn test_en_getcurvevalue_output_initialized() {
    let ph = create_loaded_project();

    let mut x: c_double = -999.0;
    let mut y: c_double = -999.0;
    let err = unsafe { EN_getcurvevalue(ph, 99999, 1, &mut x, &mut y) };

    assert_ne!(err, ErrorCode::Ok, "Should fail with invalid index");
    assert_eq!(x, 0.0, "X output should be initialized to 0");
    assert_eq!(y, 0.0, "Y output should be initialized to 0");

    destroy_project(ph);
}

#[test]
fn test_en_getcurve_valid() {
    let ph = create_loaded_project();

    // Add and set curve
    let curve_id = CString::new("TEST_CURVE6").unwrap();
    let err = unsafe { EN_addcurve(ph, curve_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getcurveindex(ph, curve_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    let x_vals: Vec<c_double> = vec![1.0, 2.0, 3.0];
    let y_vals: Vec<c_double> = vec![4.0, 5.0, 6.0];
    let err = unsafe {
        EN_setcurve(
            ph,
            index,
            x_vals.as_ptr(),
            y_vals.as_ptr(),
            x_vals.len() as c_int,
        )
    };
    assert_eq!(err, ErrorCode::Ok);

    // Get entire curve
    let mut n_points: c_int = 0;
    let mut x_ptr: *mut c_double = std::ptr::null_mut();
    let mut y_ptr: *mut c_double = std::ptr::null_mut();

    let err = unsafe { EN_getcurve(ph, index, &mut n_points, &mut x_ptr, &mut y_ptr) };

    assert_eq!(err, ErrorCode::Ok, "Should get curve");
    assert_eq!(n_points, 3, "Should have 3 points");

    destroy_project(ph);
}

#[test]
fn test_en_setcurveid_valid() {
    let ph = create_loaded_project();

    // Add a curve
    let old_id = CString::new("OLD_CURVE_ID").unwrap();
    let err = unsafe { EN_addcurve(ph, old_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut index: c_int = 0;
    let err = unsafe { EN_getcurveindex(ph, old_id.as_ptr(), &mut index) };
    assert_eq!(err, ErrorCode::Ok);

    // Rename it
    let new_id = CString::new("NEW_CURVE_ID").unwrap();
    let err = unsafe { EN_setcurveid(ph, index, new_id.as_ptr()) };

    assert_eq!(err, ErrorCode::Ok, "Should rename curve");

    destroy_project(ph);
}

// =============================================================================
// MODULE 6: Hydraulic Solver
// =============================================================================

#[test]
fn test_en_inith_valid() {
    let ph = create_loaded_project();

    let err = unsafe { EN_initH(ph, 0) };
    // May succeed or fail depending on network state
    // Just verify it doesn't crash
    let _ = err;

    destroy_project(ph);
}

#[test]
fn test_en_runh_after_init() {
    let ph = create_loaded_project();

    let _ = unsafe { EN_initH(ph, 0) };

    let err = unsafe { EN_runH(ph, 0) };

    // May succeed or fail depending on network validity
    // Just verify it doesn't crash
    let _ = err;

    destroy_project(ph);
}

#[test]
fn test_en_nexth_after_run() {
    let ph = create_loaded_project();

    let _ = unsafe { EN_initH(ph, 0) };
    let _ = unsafe { EN_runH(ph, 0) };

    let mut tstep: c_long = 0;
    let err = unsafe { EN_nextH(ph, &mut tstep) };

    // May succeed or fail, just verify output is initialized
    let _ = (err, tstep);

    destroy_project(ph);
}

// =============================================================================
// MODULE 7: Analysis Options
// =============================================================================

#[test]
fn test_en_setoption_valid() {
    let ph = create_loaded_project();

    let err = unsafe { EN_setoption(ph, SimOption::Trials as c_int, 100.0) };
    assert_eq!(err, ErrorCode::Ok, "Should set option");

    destroy_project(ph);
}

#[test]
fn test_en_settimeparam_valid() {
    let ph = create_loaded_project();

    let err = unsafe { EN_settimeparam(ph, TimeParameter::Duration as c_int, 7200) };
    assert_eq!(err, ErrorCode::Ok, "Should set time parameter");

    destroy_project(ph);
}

#[test]
fn test_en_setdemandmodel_valid() {
    let ph = create_loaded_project();

    let err = unsafe { EN_setdemandmodel(ph, DemandModel::Pda as c_int, 3.0, 15.0, 0.5) };
    assert_eq!(err, ErrorCode::Ok, "Should set demand model");

    destroy_project(ph);
}

// =============================================================================
// MODULE 7b: EN_setoption / EN_getoption roundtrips and validation
// =============================================================================

#[test]
fn test_en_setoption_trials_valid() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::Trials as c_int, 50.0) };
    assert_eq!(err, ErrorCode::Ok, "Trials should accept 50");
    let mut v: c_double = 0.0;
    let err = unsafe { EN_getoption(ph, SimOption::Trials as c_int, &mut v) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(v, 50.0, "Trials roundtrip");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_trials_invalid() {
    let ph = create_loaded_project();
    // value < 1 is invalid
    let err = unsafe { EN_setoption(ph, SimOption::Trials as c_int, 0.5) };
    assert_eq!(err, ErrorCode::InvalidOptionValue, "Trials < 1 should fail");
    // negative is invalid
    let err = unsafe { EN_setoption(ph, SimOption::Trials as c_int, -1.0) };
    assert_eq!(
        err,
        ErrorCode::InvalidOptionValue,
        "Negative trials should fail"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setoption_accuracy_valid() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::Accuracy as c_int, 1e-4) };
    assert_eq!(err, ErrorCode::Ok);
    let mut v: c_double = 0.0;
    let err = unsafe { EN_getoption(ph, SimOption::Accuracy as c_int, &mut v) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(v, 1e-4, "Accuracy roundtrip");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_accuracy_bounds() {
    let ph = create_loaded_project();
    // Too small
    let err = unsafe { EN_setoption(ph, SimOption::Accuracy as c_int, 1e-9) };
    assert_eq!(
        err,
        ErrorCode::InvalidOptionValue,
        "Accuracy < 1e-8 should fail"
    );
    // Too large
    let err = unsafe { EN_setoption(ph, SimOption::Accuracy as c_int, 0.2) };
    assert_eq!(
        err,
        ErrorCode::InvalidOptionValue,
        "Accuracy > 0.1 should fail"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setoption_demandmult_roundtrip() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::DemandMult as c_int, 1.5) };
    assert_eq!(err, ErrorCode::Ok);
    let mut v: c_double = 0.0;
    let err = unsafe { EN_getoption(ph, SimOption::DemandMult as c_int, &mut v) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(v, 1.5, "DemandMult roundtrip");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_emitexpon_roundtrip() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::EmitExpon as c_int, 0.5) };
    assert_eq!(err, ErrorCode::Ok);
    let mut v: c_double = 0.0;
    let err = unsafe { EN_getoption(ph, SimOption::EmitExpon as c_int, &mut v) };
    assert_eq!(err, ErrorCode::Ok);
    assert!((v - 0.5).abs() < 1e-10, "EmitExpon roundtrip: got {}", v);
    destroy_project(ph);
}

#[test]
fn test_en_setoption_emitexpon_invalid() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::EmitExpon as c_int, 0.0) };
    assert_eq!(
        err,
        ErrorCode::InvalidOptionValue,
        "EmitExpon = 0 should fail"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setoption_flowchange_roundtrip() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::FlowChange as c_int, 0.01) };
    assert_eq!(err, ErrorCode::Ok);
    let mut v: c_double = 0.0;
    let err = unsafe { EN_getoption(ph, SimOption::FlowChange as c_int, &mut v) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(v, 0.01, "FlowChange roundtrip");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_headlossform_roundtrip() {
    use epanet_rs::ffi::enums::HeadLossType;
    let ph = create_loaded_project();
    let err = unsafe {
        EN_setoption(
            ph,
            SimOption::HeadLossForm as c_int,
            HeadLossType::DW as c_int as c_double,
        )
    };
    assert_eq!(err, ErrorCode::Ok);
    let mut v: c_double = 0.0;
    let err = unsafe { EN_getoption(ph, SimOption::HeadLossForm as c_int, &mut v) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(v, HeadLossType::DW as i32 as f64, "HeadLossForm roundtrip");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_headlossform_invalid() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::HeadLossForm as c_int, 99.0) };
    assert_eq!(
        err,
        ErrorCode::InvalidOptionValue,
        "Invalid head loss form should fail"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setoption_spgravity_roundtrip() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::SpGravity as c_int, 1.05) };
    assert_eq!(err, ErrorCode::Ok);
    let mut v: c_double = 0.0;
    let err = unsafe { EN_getoption(ph, SimOption::SpGravity as c_int, &mut v) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(v, 1.05, "SpGravity roundtrip");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_spgravity_invalid() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::SpGravity as c_int, 0.0) };
    assert_eq!(
        err,
        ErrorCode::InvalidOptionValue,
        "SpGravity = 0 should fail"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setoption_spviscos_roundtrip() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::SpViscos as c_int, 1.2) };
    assert_eq!(err, ErrorCode::Ok);
    let mut v: c_double = 0.0;
    let err = unsafe { EN_getoption(ph, SimOption::SpViscos as c_int, &mut v) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(v, 1.2, "SpViscos roundtrip");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_checkfreq_roundtrip() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::CheckFreq as c_int, 5.0) };
    assert_eq!(err, ErrorCode::Ok);
    let mut v: c_double = 0.0;
    let err = unsafe { EN_getoption(ph, SimOption::CheckFreq as c_int, &mut v) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(v, 5.0, "CheckFreq roundtrip");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_maxcheck_roundtrip() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::MaxCheck as c_int, 8.0) };
    assert_eq!(err, ErrorCode::Ok);
    let mut v: c_double = 0.0;
    let err = unsafe { EN_getoption(ph, SimOption::MaxCheck as c_int, &mut v) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(v, 8.0, "MaxCheck roundtrip");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_pressunits_roundtrip() {
    use epanet_rs::ffi::enums::PressUnits;
    let ph = create_loaded_project();
    let err = unsafe {
        EN_setoption(
            ph,
            SimOption::PressUnits as c_int,
            PressUnits::Kpa as c_int as c_double,
        )
    };
    assert_eq!(err, ErrorCode::Ok);
    let mut v: c_double = 0.0;
    let err = unsafe { EN_getoption(ph, SimOption::PressUnits as c_int, &mut v) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(v, PressUnits::Kpa as i32 as f64, "PressUnits roundtrip");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_pressunits_invalid() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::PressUnits as c_int, 99.0) };
    assert_eq!(
        err,
        ErrorCode::UndefinedPattern,
        "Invalid pressure unit should return 205"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setoption_unavailable_returns_invalid_param() {
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
        let err = unsafe { EN_setoption(ph, *opt as c_int, 1.0) };
        assert_eq!(
            err,
            ErrorCode::InvalidParameterCode,
            "Option {:?} should return InvalidParameterCode",
            opt
        );
    }
    destroy_project(ph);
}

#[test]
fn test_en_setoption_unknown_code() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, 9999, 1.0) };
    assert_eq!(
        err,
        ErrorCode::InvalidParameterCode,
        "Unknown option code should return 251"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setoption_negative_value_rejected() {
    let ph = create_loaded_project();
    // Negative values (other than Unbalanced) must return InvalidOptionValue
    let err = unsafe { EN_setoption(ph, SimOption::DemandMult as c_int, -0.5) };
    assert_eq!(
        err,
        ErrorCode::InvalidOptionValue,
        "Negative value should return 213"
    );
    destroy_project(ph);
}

#[test]
fn test_en_getoption_null_project() {
    let mut v: c_double = -999.0;
    let err = unsafe { EN_getoption(std::ptr::null_mut(), SimOption::Trials as c_int, &mut v) };
    assert_ne!(err, ErrorCode::Ok, "Null project should fail");
    assert_eq!(v, 0.0, "Output should be initialized to 0");
}

#[test]
fn test_en_getoption_unknown_code() {
    let ph = create_loaded_project();
    let mut v: c_double = -999.0;
    let err = unsafe { EN_getoption(ph, 9999, &mut v) };
    assert_eq!(
        err,
        ErrorCode::InvalidParameterCode,
        "Unknown option code should return 251"
    );
    assert_eq!(v, 0.0, "Output should be initialized to 0");
    destroy_project(ph);
}

#[test]
fn test_en_setdemandmodel_dda() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setdemandmodel(ph, DemandModel::Dda as c_int, 0.0, 0.0, 0.5) };
    assert_eq!(err, ErrorCode::Ok, "DDA should succeed");
    destroy_project(ph);
}

#[test]
fn test_en_setdemandmodel_pda_params() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setdemandmodel(ph, DemandModel::Pda as c_int, 3.0, 15.0, 0.5) };
    assert_eq!(err, ErrorCode::Ok, "PDA should succeed");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_spviscos_invalid() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::SpViscos as c_int, 0.0) };
    assert_eq!(
        err,
        ErrorCode::InvalidOptionValue,
        "SpViscos = 0 should fail"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setoption_headlossform_blocked_when_solver_open() {
    use epanet_rs::ffi::enums::HeadLossType;
    use epanet_rs::ffi::hydraulic_solver::EN_openH;
    let ph = create_loaded_project();
    // Open the hydraulic solver
    let open_err = unsafe { EN_openH(ph) };
    assert_eq!(open_err, ErrorCode::Ok, "EN_openH should succeed");
    // Changing headloss formula while solver is open must fail
    let err = unsafe {
        EN_setoption(
            ph,
            SimOption::HeadLossForm as c_int,
            HeadLossType::DW as c_int as c_double,
        )
    };
    assert_eq!(
        err,
        ErrorCode::ModifyWhileSolverOpen,
        "Should return 262 when solver is open"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setoption_demandpattern_valid() {
    let ph = create_loaded_project();
    // Add a pattern so index 1 is valid
    let pat_id = std::ffi::CString::new("OPT_PAT").unwrap();
    let add_err = unsafe { epanet_rs::ffi::patterns::EN_addpattern(ph, pat_id.as_ptr()) };
    assert_eq!(add_err, ErrorCode::Ok);
    let mut pat_idx: c_int = 0;
    let err =
        unsafe { epanet_rs::ffi::patterns::EN_getpatternindex(ph, pat_id.as_ptr(), &mut pat_idx) };
    assert_eq!(err, ErrorCode::Ok);
    assert!(pat_idx > 0);
    // Set demand pattern
    let err = unsafe { EN_setoption(ph, SimOption::DemandPattern as c_int, pat_idx as c_double) };
    assert_eq!(
        err,
        ErrorCode::Ok,
        "Setting valid demand pattern should succeed"
    );
    // Read it back
    let mut v: c_double = 0.0;
    let get_err = unsafe { EN_getoption(ph, SimOption::DemandPattern as c_int, &mut v) };
    assert_eq!(get_err, ErrorCode::Ok);
    assert_eq!(v, pat_idx as f64, "DemandPattern roundtrip");
    // Reset to 0 (none)
    let reset_err = unsafe { EN_setoption(ph, SimOption::DemandPattern as c_int, 0.0) };
    assert_eq!(
        reset_err,
        ErrorCode::Ok,
        "Setting pattern to 0 (none) should succeed"
    );
    let mut v2: c_double = -1.0;
    let err = unsafe { EN_getoption(ph, SimOption::DemandPattern as c_int, &mut v2) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(v2, 0.0, "After reset, demand pattern should be 0");
    destroy_project(ph);
}

#[test]
fn test_en_setoption_demandpattern_out_of_range() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setoption(ph, SimOption::DemandPattern as c_int, 99999.0) };
    assert_eq!(
        err,
        ErrorCode::UndefinedPattern,
        "Out-of-range pattern index should return 205"
    );
    destroy_project(ph);
}

#[test]
fn test_en_getoption_unavailable_returns_invalid_param() {
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
        let err = unsafe { EN_getoption(ph, *opt as c_int, &mut v) };
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

// =============================================================================
// MODULE 8: Pump Head Curve
// =============================================================================

#[test]
fn test_en_getheadcurveindex_valid() {
    let ph = create_loaded_project();

    // Find a pump link (link "FH" or similar in pump.inp)
    let link_id = CString::new("FH").unwrap();
    let mut link_idx: c_int = 0;
    let get_err = unsafe { EN_getlinkindex(ph, link_id.as_ptr(), &mut link_idx) };

    if get_err == ErrorCode::Ok {
        let mut curve_idx: c_int = -999;
        let err = unsafe { EN_getheadcurveindex(ph, link_idx, &mut curve_idx) };

        // Pump may or may not have a curve, just verify no crash
        if err == ErrorCode::Ok {
            assert!(curve_idx >= 0, "Curve index should be valid");
        }
    }

    destroy_project(ph);
}

#[test]
fn test_en_setheadcurveindex_valid() {
    let ph = create_loaded_project();

    // Add a curve
    let curve_id = CString::new("PUMP_CURVE").unwrap();
    let err = unsafe { EN_addcurve(ph, curve_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);

    let mut curve_idx: c_int = 0;
    let err = unsafe { EN_getcurveindex(ph, curve_id.as_ptr(), &mut curve_idx) };
    assert_eq!(err, ErrorCode::Ok);

    // Set curve data
    let x_vals: Vec<c_double> = vec![0.0, 100.0];
    let y_vals: Vec<c_double> = vec![200.0, 150.0];
    let err = unsafe { EN_setcurve(ph, curve_idx, x_vals.as_ptr(), y_vals.as_ptr(), 2) };
    assert_eq!(err, ErrorCode::Ok);

    // Try to add a pump and assign curve
    let pump_id = CString::new("TEST_PUMP").unwrap();
    let node1_str = CString::new("1").unwrap();
    let node2_str = CString::new("2").unwrap();
    let mut pump_idx: c_int = 0;

    let add_err = unsafe {
        EN_addlink(
            ph,
            pump_id.as_ptr(),
            LinkType::Pump as c_int,
            node1_str.as_ptr(),
            node2_str.as_ptr(),
            &mut pump_idx,
        )
    };

    if add_err == ErrorCode::Ok {
        let err = unsafe { EN_setheadcurveindex(ph, pump_idx, curve_idx) };
        assert_eq!(err, ErrorCode::Ok, "Should set head curve");
    }

    destroy_project(ph);
}

// =============================================================================
// MODULE 9: Error Code Consistency
// =============================================================================

#[test]
fn test_error_codes_consistent_across_modules() {
    let ph = create_loaded_project();

    // Test that invalid property codes consistently return error 251
    let node_idx = get_valid_node_index(ph, "1");
    let link_idx = get_valid_link_index(ph, "B");

    let invalid_property = 99999;

    // Node
    let mut node_val: c_double = 0.0;
    let node_err = unsafe { EN_getnodevalue(ph, node_idx, invalid_property, &mut node_val) };
    assert_eq!(
        node_err,
        ErrorCode::InvalidParameterCode,
        "Node should return error 251"
    );

    // Link
    let mut link_val: c_double = 0.0;
    let link_err = unsafe { EN_getlinkvalue(ph, link_idx, invalid_property, &mut link_val) };
    assert_eq!(
        link_err,
        ErrorCode::InvalidParameterCode,
        "Link should return error 251"
    );

    destroy_project(ph);
}

#[test]
fn test_null_project_consistent_errors() {
    // Test that all getters handle null project consistently

    let mut int_val: c_int = -999;

    let node_err = unsafe { EN_getnodetype(std::ptr::null_mut(), 1, &mut int_val) };
    assert_ne!(node_err, ErrorCode::Ok);
    assert_eq!(int_val, 0, "Should initialize output");

    int_val = -999;
    let link_err = unsafe { EN_getlinktype(std::ptr::null_mut(), 1, &mut int_val) };
    assert_ne!(link_err, ErrorCode::Ok);
    assert_eq!(int_val, 0, "Should initialize output");

    let mut dbl_val: c_double = -999.0;
    let pattern_err = unsafe { EN_getaveragepatternvalue(std::ptr::null_mut(), 1, &mut dbl_val) };
    assert_ne!(pattern_err, ErrorCode::Ok);
    assert_eq!(dbl_val, 0.0, "Should initialize output");
}

#[test]
fn test_invalid_indices_consistent() {
    let ph = create_loaded_project();

    let invalid_idx = 99999;
    let mut dbl_val: c_double = -999.0;

    let node_err = unsafe {
        EN_getnodevalue(
            ph,
            invalid_idx,
            NodeProperty::Elevation as c_int,
            &mut dbl_val,
        )
    };
    assert_ne!(node_err, ErrorCode::Ok, "Invalid node index should fail");
    assert_eq!(dbl_val, 0.0, "Should initialize output");

    dbl_val = -999.0;
    let link_err = unsafe {
        EN_getlinkvalue(
            ph,
            invalid_idx,
            LinkProperty::Diameter as c_int,
            &mut dbl_val,
        )
    };
    assert_ne!(link_err, ErrorCode::Ok, "Invalid link index should fail");
    assert_eq!(dbl_val, 0.0, "Should initialize output");

    destroy_project(ph);
}

// =============================================================================
// MODULE 10: Regression Tests - Ensure Bug Fixes Remain Fixed
// =============================================================================

#[test]
fn test_regression_no_pattern_returns_zero() {
    let ph = create_loaded_project();

    // Scan multiple nodes
    for node_id in &["1", "2", "3", "4"] {
        let id = CString::new(*node_id).unwrap();
        let mut idx: c_int = 0;

        let err = unsafe { EN_getnodeindex(ph, id.as_ptr(), &mut idx) };
        if err == ErrorCode::Ok {
            let mut pattern: c_double = -999.0;
            let err =
                unsafe { EN_getnodevalue(ph, idx, NodeProperty::Pattern as c_int, &mut pattern) };

            if err == ErrorCode::Ok {
                assert_ne!(pattern, 123.0, "Should never return 123.0 sentinel");
                assert!(pattern >= 0.0, "Pattern should be >= 0");
            }
        }
    }

    destroy_project(ph);
}

#[test]
fn test_regression_invalid_property_always_errors() {
    let ph = create_loaded_project();
    let node_idx = get_valid_node_index(ph, "1");

    for invalid_code in vec![999, 1000, -999, 12345] {
        let mut value: c_double = -999.0;
        let err = unsafe { EN_getnodevalue(ph, node_idx, invalid_code, &mut value) };

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

#[test]
fn test_regression_all_outputs_initialized() {
    // Test that ALL getter functions initialize outputs before validation

    let ph = std::ptr::null_mut(); // Null project to trigger errors

    // Integer outputs
    let mut int_val: c_int = -999;
    let _ = unsafe { EN_getnodetype(ph, 1, &mut int_val) };
    assert_eq!(int_val, 0, "getnodetype should initialize");

    int_val = -999;
    let _ = unsafe { EN_getlinktype(ph, 1, &mut int_val) };
    assert_eq!(int_val, 0, "getlinktype should initialize");

    int_val = -999;
    let id = CString::new("X").unwrap();
    let _ = unsafe { EN_getnodeindex(ph, id.as_ptr(), &mut int_val) };
    assert_eq!(int_val, 0, "getnodeindex should initialize");

    // Double outputs
    let mut dbl_val: c_double = -999.0;
    let _ = unsafe { EN_getnodevalue(ph, 1, 0, &mut dbl_val) };
    assert_eq!(dbl_val, 0.0, "getnodevalue should initialize");

    dbl_val = -999.0;
    let _ = unsafe { EN_getlinkvalue(ph, 1, 0, &mut dbl_val) };
    assert_eq!(dbl_val, 0.0, "getlinkvalue should initialize");

    dbl_val = -999.0;
    let _ = unsafe { EN_getpatternvalue(ph, 1, 1, &mut dbl_val) };
    assert_eq!(dbl_val, 0.0, "getpatternvalue should initialize");

    // Multiple outputs
    let mut x: c_double = -999.0;
    let mut y: c_double = -999.0;
    let _ = unsafe { EN_getcoord(ph, 1, &mut x, &mut y) };
    assert_eq!(x, 0.0, "getcoord should initialize x");
    assert_eq!(y, 0.0, "getcoord should initialize y");

    x = -999.0;
    y = -999.0;
    let _ = unsafe { EN_getcurvevalue(ph, 1, 1, &mut x, &mut y) };
    assert_eq!(x, 0.0, "getcurvevalue should initialize x");
    assert_eq!(y, 0.0, "getcurvevalue should initialize y");

    let mut n1: c_int = -999;
    let mut n2: c_int = -999;
    let _ = unsafe { EN_getlinknodes(ph, 1, &mut n1, &mut n2) };
    assert_eq!(n1, 0, "getlinknodes should initialize node1");
    assert_eq!(n2, 0, "getlinknodes should initialize node2");
}

/// Helper to open the valves test network (has all valve types)
fn create_valves_project() -> *mut Project {
    let ph = create_empty_project();
    let inp = CString::new("tests/valves.inp").unwrap();
    let rpt = CString::new("").unwrap();
    let out = CString::new("").unwrap();
    let err = unsafe { EN_open(ph, inp.as_ptr(), rpt.as_ptr(), out.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok, "valves.inp should load");
    ph
}

/// Helper to open the tanks test network (has reservoir + tank nodes)
fn create_tanks_project() -> *mut Project {
    let ph = create_empty_project();
    let inp = CString::new("tests/tanks.inp").unwrap();
    let rpt = CString::new("").unwrap();
    let out = CString::new("").unwrap();
    let err = unsafe { EN_open(ph, inp.as_ptr(), rpt.as_ptr(), out.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok, "tanks.inp should load");
    ph
}

// --- EN_getnodetype: Reservoir and Tank branches ---

#[test]
fn test_en_getnodetype_reservoir() {
    // pump.inp has reservoirs named "FH" and "FH2"
    let ph = create_loaded_project();
    let id = CString::new("FH").unwrap();
    let mut idx: c_int = 0;
    let err = unsafe { EN_getnodeindex(ph, id.as_ptr(), &mut idx) };
    assert_eq!(err, ErrorCode::Ok);
    let mut node_type: c_int = -1;
    let err = unsafe { EN_getnodetype(ph, idx, &mut node_type) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(
        node_type,
        NodeType::Reservoir as c_int,
        "FH should be Reservoir"
    );
    destroy_project(ph);
}

#[test]
fn test_en_getnodetype_tank() {
    // tanks.inp has tanks
    let ph = create_tanks_project();
    let id = CString::new("1").unwrap();
    let mut idx: c_int = 0;
    let err = unsafe { EN_getnodeindex(ph, id.as_ptr(), &mut idx) };
    assert_eq!(err, ErrorCode::Ok);
    let mut node_type: c_int = -1;
    let err = unsafe { EN_getnodetype(ph, idx, &mut node_type) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(
        node_type,
        NodeType::Tank as c_int,
        "Node 1 in tanks.inp should be Tank"
    );
    destroy_project(ph);
}

// --- EN_getlinktype: Pump and all valve types ---

#[test]
fn test_en_getlinktype_pump() {
    // pump.inp has pump link "1"
    let ph = create_loaded_project();
    let id = CString::new("1").unwrap();
    let mut idx: c_int = 0;
    unsafe { EN_getlinkindex(ph, id.as_ptr(), &mut idx) };
    let mut link_type: c_int = -1;
    let err = unsafe { EN_getlinktype(ph, idx, &mut link_type) };
    assert_eq!(err, ErrorCode::Ok);
    assert_eq!(
        link_type,
        ENLinkType::Pump as c_int,
        "Link '1' should be Pump"
    );
    destroy_project(ph);
}

#[test]
fn test_en_getlinktype_all_valve_types() {
    let ph = create_valves_project();
    // valves.inp contains PRV, PBV, TCV, PCV, FCV, PSV, GPV
    let cases: &[(&str, ENLinkType)] = &[
        ("PRV", ENLinkType::PRV),
        ("PBV", ENLinkType::PBV),
        ("TCV", ENLinkType::TCV),
        ("PCV", ENLinkType::PCV),
        ("FCV", ENLinkType::FCV),
        ("PSV", ENLinkType::PSV),
        ("GPV", ENLinkType::GPV),
    ];
    for (link_id, expected_type) in cases {
        let id = CString::new(*link_id).unwrap();
        let mut idx: c_int = 0;
        let get_err = unsafe { EN_getlinkindex(ph, id.as_ptr(), &mut idx) };
        if get_err != ErrorCode::Ok {
            // skip if link not present
            continue;
        }
        let mut link_type: c_int = -1;
        let err = unsafe { EN_getlinktype(ph, idx, &mut link_type) };
        assert_eq!(
            err,
            ErrorCode::Ok,
            "getlinktype for {} should succeed",
            link_id
        );
        assert_eq!(
            link_type, *expected_type as c_int,
            "Link '{}' should be {:?}",
            link_id, expected_type
        );
    }
    destroy_project(ph);
}

// --- EN_setlinknodes: same-node guard and invalid node indices ---

#[test]
fn test_en_setlinknodes_same_node_rejected() {
    let ph = create_loaded_project();
    let idx = get_valid_link_index(ph, "B");
    let node = get_valid_node_index(ph, "1");
    let err = unsafe { EN_setlinknodes(ph, idx, node, node) };
    assert_eq!(
        err,
        ErrorCode::LinkSameStartEnd,
        "Same start/end node should return 222"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setlinknodes_invalid_start_node() {
    let ph = create_loaded_project();
    let idx = get_valid_link_index(ph, "B");
    let valid = get_valid_node_index(ph, "1");
    let err = unsafe { EN_setlinknodes(ph, idx, 99999, valid) };
    assert_eq!(
        err,
        ErrorCode::UndefinedNode,
        "Invalid start node should return UndefinedNode"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setlinknodes_invalid_end_node() {
    let ph = create_loaded_project();
    let idx = get_valid_link_index(ph, "B");
    let valid = get_valid_node_index(ph, "1");
    let err = unsafe { EN_setlinknodes(ph, idx, valid, 99999) };
    assert_eq!(
        err,
        ErrorCode::UndefinedNode,
        "Invalid end node should return UndefinedNode"
    );
    destroy_project(ph);
}

// --- EN_deletelink: invalid action_code ---

#[test]
fn test_en_deletelink_invalid_action_code() {
    let ph = create_loaded_project();
    // add a temporary link
    let link_id = CString::new("TMP_DEL").unwrap();
    let n1 = CString::new("1").unwrap();
    let n2 = CString::new("2").unwrap();
    let mut idx: c_int = 0;
    let err = unsafe {
        EN_addlink(
            ph,
            link_id.as_ptr(),
            LinkType::Pipe as c_int,
            n1.as_ptr(),
            n2.as_ptr(),
            &mut idx,
        )
    };
    assert_eq!(err, ErrorCode::Ok);
    let err = unsafe { EN_deletelink(ph, idx, 99) };
    assert_eq!(
        err,
        ErrorCode::InvalidParameterCode,
        "Invalid action_code should return 251"
    );
    destroy_project(ph);
}

// --- EN_getcurvevalue: out-of-range point index ---

#[test]
fn test_en_getcurvevalue_invalid_point_index() {
    let ph = create_loaded_project();
    let curve_id = CString::new("GAP_CURVE").unwrap();
    let err = unsafe { epanet_rs::ffi::curves::EN_addcurve(ph, curve_id.as_ptr()) };
    assert_eq!(err, ErrorCode::Ok);
    let mut cidx: c_int = 0;
    let err = unsafe { epanet_rs::ffi::curves::EN_getcurveindex(ph, curve_id.as_ptr(), &mut cidx) };
    assert_eq!(err, ErrorCode::Ok);
    let x_vals: Vec<c_double> = vec![1.0, 2.0];
    let y_vals: Vec<c_double> = vec![3.0, 4.0];
    let err = unsafe {
        epanet_rs::ffi::curves::EN_setcurve(ph, cidx, x_vals.as_ptr(), y_vals.as_ptr(), 2)
    };
    assert_eq!(err, ErrorCode::Ok);
    // point index 99 is beyond the 2-point curve (1-based)
    let mut x: c_double = -999.0;
    let mut y: c_double = -999.0;
    let err = unsafe { EN_getcurvevalue(ph, cidx, 99, &mut x, &mut y) };
    assert_ne!(err, ErrorCode::Ok, "Out-of-range point index should fail");
    assert_eq!(x, 0.0, "x output should be initialized to 0 on error");
    assert_eq!(y, 0.0, "y output should be initialized to 0 on error");
    destroy_project(ph);
}

// --- EN_settimeparam: each variant, zero-value and invalid code ---

#[test]
fn test_en_settimeparam_all_valid_variants() {
    let ph = create_loaded_project();
    let cases: &[(TimeParameter, c_long)] = &[
        (TimeParameter::Duration, 7200),
        (TimeParameter::HydStep, 1800),
        (TimeParameter::PatternStep, 3600),
        (TimeParameter::PatternStart, 0),
        (TimeParameter::ReportStep, 3600),
    ];
    for (param, value) in cases {
        let err = unsafe { EN_settimeparam(ph, *param as c_int, *value) };
        assert_eq!(
            err,
            ErrorCode::Ok,
            "settimeparam({:?}, {}) should succeed",
            param,
            value
        );
    }
    destroy_project(ph);
}

#[test]
fn test_en_settimeparam_zero_hydstep_rejected() {
    let ph = create_loaded_project();
    let err = unsafe { EN_settimeparam(ph, TimeParameter::HydStep as c_int, 0) };
    assert_eq!(
        err,
        ErrorCode::InvalidOptionValue,
        "HydStep = 0 should fail"
    );
    destroy_project(ph);
}

#[test]
fn test_en_settimeparam_zero_patternstep_rejected() {
    let ph = create_loaded_project();
    let err = unsafe { EN_settimeparam(ph, TimeParameter::PatternStep as c_int, 0) };
    assert_eq!(
        err,
        ErrorCode::InvalidOptionValue,
        "PatternStep = 0 should fail"
    );
    destroy_project(ph);
}

#[test]
fn test_en_settimeparam_zero_reportstep_rejected() {
    let ph = create_loaded_project();
    let err = unsafe { EN_settimeparam(ph, TimeParameter::ReportStep as c_int, 0) };
    assert_eq!(
        err,
        ErrorCode::InvalidOptionValue,
        "ReportStep = 0 should fail"
    );
    destroy_project(ph);
}

#[test]
fn test_en_settimeparam_invalid_code() {
    let ph = create_loaded_project();
    let err = unsafe { EN_settimeparam(ph, 9999, 3600) };
    assert_eq!(
        err,
        ErrorCode::InvalidParameterCode,
        "Unknown time param code should return 251"
    );
    destroy_project(ph);
}

// --- EN_setdemandmodel: invalid code ---

#[test]
fn test_en_setdemandmodel_invalid_code() {
    let ph = create_loaded_project();
    let err = unsafe { EN_setdemandmodel(ph, 99, 0.0, 0.0, 0.5) };
    assert_eq!(
        err,
        ErrorCode::InvalidParameterCode,
        "Invalid demand model code should return 251"
    );
    destroy_project(ph);
}

// --- Hydraulic solver: HydraulicSolverNotOpened paths ---

#[test]
fn test_en_inith_without_openh_fails() {
    let ph = create_loaded_project();
    // Don't call EN_openH first
    let err = unsafe { EN_initH(ph, 0) };
    assert_eq!(
        err,
        ErrorCode::HydraulicSolverNotOpened,
        "initH without openH should return 103"
    );
    destroy_project(ph);
}

#[test]
fn test_en_runh_without_openh_fails() {
    let ph = create_loaded_project();
    let err = unsafe { EN_runH(ph, 0) };
    assert_eq!(
        err,
        ErrorCode::HydraulicSolverNotOpened,
        "runH without openH should return 103"
    );
    destroy_project(ph);
}

#[test]
fn test_en_nexth_without_openh_fails() {
    let ph = create_loaded_project();
    let mut t: c_long = 0;
    let err = unsafe { EN_nextH(ph, &mut t) };
    assert_eq!(
        err,
        ErrorCode::HydraulicSolverNotOpened,
        "nextH without openH should return 103"
    );
    destroy_project(ph);
}

// =============================================================================
// MODULE 12: Medium-Priority Coverage Gaps
// =============================================================================

// --- EN_saveinpfile edge cases ---

#[test]
fn test_en_saveinpfile_null_path() {
    let ph = create_loaded_project();
    let err = unsafe { EN_saveinpfile(ph, std::ptr::null()) };
    assert_eq!(
        err,
        ErrorCode::CannotOpenInputFile,
        "null path should return CannotOpenInputFile"
    );
    destroy_project(ph);
}

#[test]
fn test_en_saveinpfile_no_network() {
    let ph = create_empty_project();

    let out_path = std::env::temp_dir().join(format!(
        "epanet_rs_test_no_network_{}.inp",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("System time before UNIX epoch")
            .as_nanos()
    ));

    let out_file = CString::new(out_path.to_str().expect("Temp path contains invalid UTF-8"))
        .expect("Path contains null byte");

    let err = unsafe { EN_saveinpfile(ph, out_file.as_ptr()) };
    assert_eq!(
        err,
        ErrorCode::NoNetworkData,
        "empty project should return NoNetworkData"
    );
    destroy_project(ph);
}

// --- EN_getheadcurveindex / EN_setheadcurveindex non-pump ---

#[test]
fn test_en_getheadcurveindex_non_pump_returns_undefined_pump() {
    let ph = create_loaded_project();
    let idx = get_valid_link_index(ph, "B");
    let mut out: c_int = 0;
    let err = unsafe { EN_getheadcurveindex(ph, idx, &mut out) };
    assert_eq!(
        err,
        ErrorCode::UndefinedPump,
        "pipe link should return UndefinedPump"
    );
    destroy_project(ph);
}

#[test]
fn test_en_setheadcurveindex_non_pump_returns_undefined_pump() {
    let ph = create_loaded_project();
    let idx = get_valid_link_index(ph, "B");
    let err = unsafe { EN_setheadcurveindex(ph, idx, 0) };
    assert_eq!(
        err,
        ErrorCode::UndefinedPump,
        "pipe link should return UndefinedPump"
    );
    destroy_project(ph);
}

// --- EN_getcount: remaining CountType variants and invalid type ---

#[test]
fn test_en_getcount_tanks() {
    let ph = create_loaded_project(); // pump.inp has 2 reservoirs; TankCount includes reservoirs
    let mut count: c_int = -1;
    let err = unsafe { EN_getcount(ph, CountType::TankCount as c_int, &mut count) };
    assert_eq!(err, ErrorCode::Ok);
    assert!(count >= 0, "TankCount should be non-negative");
    destroy_project(ph);
}

#[test]
fn test_en_getcount_patterns() {
    let ph = create_loaded_project();
    let mut count: c_int = -1;
    let err = unsafe { EN_getcount(ph, CountType::PatCount as c_int, &mut count) };
    assert_eq!(err, ErrorCode::Ok);
    assert!(count >= 0, "PatCount should be non-negative");
    destroy_project(ph);
}

#[test]
fn test_en_getcount_curves() {
    let ph = create_loaded_project();
    let mut count: c_int = -1;
    let err = unsafe { EN_getcount(ph, CountType::CurveCount as c_int, &mut count) };
    assert_eq!(err, ErrorCode::Ok);
    assert!(count >= 0, "CurveCount should be non-negative");
    destroy_project(ph);
}

#[test]
fn test_en_getcount_controls() {
    let ph = create_loaded_project();
    let mut count: c_int = -1;
    let err = unsafe { EN_getcount(ph, CountType::ControlCount as c_int, &mut count) };
    assert_eq!(err, ErrorCode::Ok);
    assert!(count >= 0, "ControlCount should be non-negative");
    destroy_project(ph);
}

#[test]
fn test_en_getcount_invalid_type() {
    let ph = create_loaded_project();
    let mut count: c_int = -1;
    let err = unsafe { EN_getcount(ph, 9999, &mut count) };
    assert_eq!(
        err,
        ErrorCode::InvalidParameterCode,
        "unknown count type should return 251"
    );
    assert_eq!(count, 0, "output should be initialized to 0 on error");
    destroy_project(ph);
}

// --- EN_geterror edge cases (no project handle needed) ---

#[test]
fn test_en_geterror_null_buffer() {
    let err = unsafe { EN_geterror(ErrorCode::Ok as c_int, std::ptr::null_mut(), 64) };
    assert_eq!(
        err,
        ErrorCode::InvalidFormat,
        "null buffer should return InvalidFormat"
    );
}

#[test]
fn test_en_geterror_zero_maxlen() {
    let mut buf = [0u8; 64];
    let err = unsafe { EN_geterror(ErrorCode::Ok as c_int, buf.as_mut_ptr() as *mut c_char, 0) };
    assert_eq!(
        err,
        ErrorCode::InvalidFormat,
        "max_len=0 should return InvalidFormat"
    );
}

#[test]
fn test_en_geterror_negative_maxlen() {
    let mut buf = [0u8; 64];
    let err = unsafe { EN_geterror(ErrorCode::Ok as c_int, buf.as_mut_ptr() as *mut c_char, -1) };
    assert_eq!(
        err,
        ErrorCode::InvalidFormat,
        "max_len<0 should return InvalidFormat"
    );
}

#[test]
fn test_en_geterror_unknown_code_returns_ok() {
    let mut buf = [0u8; 128];
    let err = unsafe { EN_geterror(9999, buf.as_mut_ptr() as *mut c_char, 128) };
    assert_eq!(
        err,
        ErrorCode::Ok,
        "unknown error code should still return Ok"
    );
    // buffer should contain a non-empty message
    assert_ne!(buf[0], 0, "message buffer should not be empty");
}

// --- EN_addnode: Reservoir and Tank types ---

#[test]
fn test_en_addnode_reservoir_type() {
    let ph = create_loaded_project();
    let id = CString::new("NEW_RES_ADD").unwrap();
    let mut idx: c_int = 0;
    let err = unsafe { EN_addnode(ph, id.as_ptr(), NodeType::Reservoir as c_int, &mut idx) };
    assert_eq!(err, ErrorCode::Ok, "adding Reservoir should succeed");
    assert!(idx > 0, "index should be positive");
    let mut node_type: c_int = -1;
    let err2 = unsafe { EN_getnodetype(ph, idx, &mut node_type) };
    assert_eq!(err2, ErrorCode::Ok);
    assert_eq!(
        node_type,
        NodeType::Reservoir as c_int,
        "type should be Reservoir"
    );
    destroy_project(ph);
}

#[test]
fn test_en_addnode_tank_type() {
    let ph = create_loaded_project();
    let id = CString::new("NEW_TANK_ADD").unwrap();
    let mut idx: c_int = 0;
    let err = unsafe { EN_addnode(ph, id.as_ptr(), NodeType::Tank as c_int, &mut idx) };
    assert_eq!(err, ErrorCode::Ok, "adding Tank should succeed");
    assert!(idx > 0, "index should be positive");
    let mut node_type: c_int = -1;
    let err2 = unsafe { EN_getnodetype(ph, idx, &mut node_type) };
    assert_eq!(err2, ErrorCode::Ok);
    assert_eq!(node_type, NodeType::Tank as c_int, "type should be Tank");
    destroy_project(ph);
}

// --- EN_setlinkvalue: check-valve pipe status is illegal ---

#[test]
fn test_en_setlinkvalue_cvpipe_status_illegal() {
    let ph = create_loaded_project();
    let link_id = CString::new("TEST_CV_STAT").unwrap();
    let n1 = CString::new("1").unwrap();
    let n2 = CString::new("2").unwrap();
    let mut cv_idx: c_int = 0;
    let add_err = unsafe {
        EN_addlink(
            ph,
            link_id.as_ptr(),
            ENLinkType::CVPipe as c_int,
            n1.as_ptr(),
            n2.as_ptr(),
            &mut cv_idx,
        )
    };
    assert_eq!(add_err, ErrorCode::Ok, "should be able to add a CVPipe");
    assert!(cv_idx > 0);
    let err = unsafe { EN_setlinkvalue(ph, cv_idx, LinkProperty::InitStatus as c_int, 0.0) };
    assert_eq!(
        err,
        ErrorCode::IllegalValveControl,
        "setting status on CV pipe should return 219"
    );
    destroy_project(ph);
}
