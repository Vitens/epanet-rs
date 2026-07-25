//! Smoke tests exercising the toolkit API the way a C caller would.

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_double, c_int, c_long};
use std::ptr;

use crate::ffi::analysis_options::*;
use crate::ffi::controls::*;
use crate::ffi::curves::*;
use crate::ffi::demands::*;
use crate::ffi::enums::{
    ControlType, CountType, CurveType, LinkProperty, LinkType, NodeProperty, TimeParameter,
};
use crate::ffi::error_codes::ErrorCode;
use crate::ffi::hydraulic_solver::*;
use crate::ffi::links::*;
use crate::ffi::nodes::*;
use crate::ffi::patterns::*;
use crate::ffi::project::*;
use crate::ffi::quality_solver::*;
use crate::ffi::report::*;
use crate::ffi::rules::*;
use crate::ffi::util::MAX_ID;

/// A project handle that closes itself at the end of a test.
struct TestProject(*mut Project);

impl TestProject {
    fn open(path: &str) -> Self {
        let mut ph: *mut Project = ptr::null_mut();
        assert_eq!(unsafe { EN_createproject(&mut ph) }, ErrorCode::Ok);
        let inp = CString::new(path).unwrap();
        assert_eq!(
            unsafe { EN_open(ph, inp.as_ptr(), ptr::null(), ptr::null()) },
            ErrorCode::Ok
        );
        Self(ph)
    }

    fn handle(&self) -> *mut Project {
        self.0
    }
}

impl Drop for TestProject {
    fn drop(&mut self) {
        unsafe {
            EN_close(self.0);
            EN_deleteproject(self.0);
        }
    }
}

fn id_buffer() -> Vec<c_char> {
    vec![0; MAX_ID + 1]
}

fn buffer_to_string(buffer: &[c_char]) -> String {
    unsafe { CStr::from_ptr(buffer.as_ptr()) }
        .to_string_lossy()
        .into_owned()
}

fn count(ph: *mut Project, what: CountType) -> c_int {
    let mut value = 0;
    assert_eq!(
        unsafe { EN_getcount(ph, what as c_int, &mut value) },
        ErrorCode::Ok
    );
    value
}

fn node_index(ph: *mut Project, id: &str) -> c_int {
    let id = CString::new(id).unwrap();
    let mut index = 0;
    assert_eq!(
        unsafe { EN_getnodeindex(ph, id.as_ptr(), &mut index) },
        ErrorCode::Ok
    );
    index
}

fn link_index(ph: *mut Project, id: &str) -> c_int {
    let id = CString::new(id).unwrap();
    let mut index = 0;
    assert_eq!(
        unsafe { EN_getlinkindex(ph, id.as_ptr(), &mut index) },
        ErrorCode::Ok
    );
    index
}

fn node_value(ph: *mut Project, index: c_int, property: NodeProperty) -> c_double {
    let mut value = 0.0;
    assert_eq!(
        unsafe { EN_getnodevalue(ph, index, property as c_int, &mut value) },
        ErrorCode::Ok
    );
    value
}

fn link_value(ph: *mut Project, index: c_int, property: LinkProperty) -> c_double {
    let mut value = 0.0;
    assert_eq!(
        unsafe { EN_getlinkvalue(ph, index, property as c_int, &mut value) },
        ErrorCode::Ok
    );
    value
}

#[test]
fn version_matches_epanet_23() {
    let mut version = 0;
    assert_eq!(unsafe { EN_getversion(&mut version) }, ErrorCode::Ok);
    assert_eq!(version, 20300);
}

#[test]
fn open_reports_counts_and_ids() {
    let project = TestProject::open("tests/pump.inp");
    let ph = project.handle();

    assert_eq!(count(ph, CountType::NodeCount), 9);
    assert_eq!(count(ph, CountType::LinkCount), 9);
    assert_eq!(count(ph, CountType::CurveCount), 1);
    assert_eq!(count(ph, CountType::ControlCount), 0);
    assert_eq!(count(ph, CountType::RuleCount), 0);

    let index = node_index(ph, "FH");
    let mut buffer = id_buffer();
    assert_eq!(
        unsafe { EN_getnodeid(ph, index, buffer.as_mut_ptr()) },
        ErrorCode::Ok
    );
    assert_eq!(buffer_to_string(&buffer), "FH");

    let mut node_type = 0;
    assert_eq!(
        unsafe { EN_getnodetype(ph, index, &mut node_type) },
        ErrorCode::Ok
    );
    assert_eq!(node_type, crate::ffi::enums::NodeType::Reservoir as c_int);
}

#[test]
fn solve_hydraulics_and_read_results() {
    let project = TestProject::open("tests/pump.inp");
    let ph = project.handle();

    assert_eq!(unsafe { EN_solveH(ph) }, ErrorCode::Ok);

    let junction = node_index(ph, "3");
    let pressure = node_value(ph, junction, NodeProperty::Pressure);
    let head = node_value(ph, junction, NodeProperty::Head);
    assert!(
        pressure > 0.0,
        "expected a positive pressure, got {pressure}"
    );
    assert!(head > 0.0);

    let pipe = link_index(ph, "B");
    let flow = link_value(ph, pipe, LinkProperty::Flow);
    assert!(flow.abs() > 0.0);
}

#[test]
fn step_through_a_simulation() {
    let project = TestProject::open("tests/2tanks.inp");
    let ph = project.handle();

    assert_eq!(unsafe { EN_openH(ph) }, ErrorCode::Ok);
    assert_eq!(unsafe { EN_initH(ph, 0) }, ErrorCode::Ok);

    let mut steps = 0;
    loop {
        let mut now: c_long = 0;
        assert_eq!(unsafe { EN_runH(ph, &mut now) }, ErrorCode::Ok);
        steps += 1;

        let mut step: c_long = 0;
        assert_eq!(unsafe { EN_nextH(ph, &mut step) }, ErrorCode::Ok);
        if step == 0 {
            break;
        }
        assert!(steps < 1000, "simulation did not terminate");
    }

    assert!(steps > 1);
    assert_eq!(unsafe { EN_closeH(ph) }, ErrorCode::Ok);
}

#[test]
fn node_and_link_properties_round_trip() {
    let project = TestProject::open("tests/pump.inp");
    let ph = project.handle();

    let junction = node_index(ph, "3");
    assert_eq!(
        unsafe { EN_setnodevalue(ph, junction, NodeProperty::Elevation as c_int, 12.5) },
        ErrorCode::Ok
    );
    assert!((node_value(ph, junction, NodeProperty::Elevation) - 12.5).abs() < 1e-9);

    let pipe = link_index(ph, "B");
    for (property, value) in [
        (LinkProperty::Diameter, 8.0),
        (LinkProperty::Length, 250.0),
        (LinkProperty::Roughness, 120.0),
        (LinkProperty::MinorLoss, 0.5),
    ] {
        assert_eq!(
            unsafe { EN_setlinkvalue(ph, pipe, property as c_int, value) },
            ErrorCode::Ok
        );
        assert!(
            (link_value(ph, pipe, property) - value).abs() < 1e-6,
            "{property:?} did not round-trip"
        );
    }

    // one call for all pipe properties at once
    assert_eq!(
        unsafe { EN_setpipedata(ph, pipe, 300.0, 10.0, 110.0, 0.0) },
        ErrorCode::Ok
    );
    assert!((link_value(ph, pipe, LinkProperty::Length) - 300.0).abs() < 1e-6);
    assert!((link_value(ph, pipe, LinkProperty::Diameter) - 10.0).abs() < 1e-6);
}

#[test]
fn link_type_conversions() {
    let project = TestProject::open("tests/pump.inp");
    let ph = project.handle();

    let mut index = link_index(ph, "B");

    // a check valve pipe keeps its index
    let original_index = index;
    assert_eq!(
        unsafe { EN_setlinktype(ph, &mut index, LinkType::CVPipe as c_int, 0) },
        ErrorCode::Ok
    );
    assert_eq!(index, original_index);
    let mut link_type = 0;
    assert_eq!(
        unsafe { EN_getlinktype(ph, index, &mut link_type) },
        ErrorCode::Ok
    );
    assert_eq!(link_type, LinkType::CVPipe as c_int);

    // turning it into a valve recreates the link, so the index is updated
    assert_eq!(
        unsafe { EN_setlinktype(ph, &mut index, LinkType::TCV as c_int, 0) },
        ErrorCode::Ok
    );
    assert_eq!(
        unsafe { EN_getlinktype(ph, index, &mut link_type) },
        ErrorCode::Ok
    );
    assert_eq!(link_type, LinkType::TCV as c_int);
    assert_eq!(index, link_index(ph, "B"));

    // changing the kind of valve is an in-place update
    assert_eq!(
        unsafe { EN_setlinktype(ph, &mut index, LinkType::PRV as c_int, 0) },
        ErrorCode::Ok
    );
    assert_eq!(
        unsafe { EN_getlinkvalue(ph, index, LinkProperty::ValveType as c_int, &mut 0.0) },
        ErrorCode::Ok
    );
    assert_eq!(
        unsafe { EN_getlinktype(ph, index, &mut link_type) },
        ErrorCode::Ok
    );
    assert_eq!(link_type, LinkType::PRV as c_int);
}

#[test]
fn build_a_network_from_scratch() {
    let mut ph: *mut Project = ptr::null_mut();
    assert_eq!(unsafe { EN_createproject(&mut ph) }, ErrorCode::Ok);
    assert_eq!(
        unsafe {
            EN_init(
                ph,
                ptr::null(),
                ptr::null(),
                crate::ffi::enums::FlowUnits::Lps as c_int,
                crate::ffi::enums::HeadLossType::DW as c_int,
            )
        },
        ErrorCode::Ok
    );

    let reservoir = CString::new("R1").unwrap();
    let junction = CString::new("J1").unwrap();
    let pipe = CString::new("P1").unwrap();

    let mut index = 0;
    assert_eq!(
        unsafe {
            EN_addnode(
                ph,
                reservoir.as_ptr(),
                crate::ffi::enums::NodeType::Reservoir as c_int,
                &mut index,
            )
        },
        ErrorCode::Ok
    );
    let reservoir_index = index;
    assert_eq!(
        unsafe {
            EN_addnode(
                ph,
                junction.as_ptr(),
                crate::ffi::enums::NodeType::Junction as c_int,
                &mut index,
            )
        },
        ErrorCode::Ok
    );
    let junction_index = index;

    assert_eq!(
        unsafe {
            EN_addlink(
                ph,
                pipe.as_ptr(),
                LinkType::Pipe as c_int,
                reservoir.as_ptr(),
                junction.as_ptr(),
                &mut index,
            )
        },
        ErrorCode::Ok
    );
    let pipe_index = index;

    assert_eq!(
        unsafe { EN_setnodevalue(ph, reservoir_index, NodeProperty::Elevation as c_int, 100.0) },
        ErrorCode::Ok
    );
    assert_eq!(
        unsafe { EN_setjuncdata(ph, junction_index, 0.0, 10.0, ptr::null()) },
        ErrorCode::Ok
    );
    assert_eq!(
        unsafe { EN_setpipedata(ph, pipe_index, 1000.0, 200.0, 0.1, 0.0) },
        ErrorCode::Ok
    );

    assert_eq!(unsafe { EN_solveH(ph) }, ErrorCode::Ok);

    let mut flow = 0.0;
    assert_eq!(
        unsafe { EN_getlinkvalue(ph, pipe_index, LinkProperty::Flow as c_int, &mut flow) },
        ErrorCode::Ok
    );
    assert!(
        (flow - 10.0).abs() < 1e-3,
        "expected the demand to be met, got {flow}"
    );

    assert_eq!(unsafe { EN_deleteproject(ph) }, ErrorCode::Ok);
}

#[test]
fn patterns_round_trip() {
    let project = TestProject::open("tests/pump.inp");
    let ph = project.handle();

    let id = CString::new("PAT1").unwrap();
    assert_eq!(unsafe { EN_addpattern(ph, id.as_ptr()) }, ErrorCode::Ok);

    let mut index = 0;
    assert_eq!(
        unsafe { EN_getpatternindex(ph, id.as_ptr(), &mut index) },
        ErrorCode::Ok
    );

    let multipliers = [0.5, 1.0, 1.5, 2.0];
    assert_eq!(
        unsafe { EN_setpattern(ph, index, multipliers.as_ptr(), multipliers.len() as c_int) },
        ErrorCode::Ok
    );

    let mut length = 0;
    assert_eq!(
        unsafe { EN_getpatternlen(ph, index, &mut length) },
        ErrorCode::Ok
    );
    assert_eq!(length, 4);

    assert_eq!(
        unsafe { EN_setpatternvalue(ph, index, 2, 1.25) },
        ErrorCode::Ok
    );
    let mut value = 0.0;
    assert_eq!(
        unsafe { EN_getpatternvalue(ph, index, 2, &mut value) },
        ErrorCode::Ok
    );
    assert!((value - 1.25).abs() < 1e-9);

    // periods are 1-based, so 0 and len + 1 are out of range
    assert_eq!(
        unsafe { EN_getpatternvalue(ph, index, 0, &mut value) },
        ErrorCode::InvalidParameterCode
    );
    assert_eq!(
        unsafe { EN_getpatternvalue(ph, index, 5, &mut value) },
        ErrorCode::InvalidParameterCode
    );

    let mut average = 0.0;
    assert_eq!(
        unsafe { EN_getaveragepatternvalue(ph, index, &mut average) },
        ErrorCode::Ok
    );
    assert!((average - (0.5 + 1.25 + 1.5 + 2.0) / 4.0).abs() < 1e-9);

    assert_eq!(unsafe { EN_deletepattern(ph, index) }, ErrorCode::Ok);
}

#[test]
fn curves_round_trip() {
    let project = TestProject::open("tests/pump.inp");
    let ph = project.handle();

    let id = CString::new("CURVE1").unwrap();
    assert_eq!(unsafe { EN_addcurve(ph, id.as_ptr()) }, ErrorCode::Ok);

    let mut index = 0;
    assert_eq!(
        unsafe { EN_getcurveindex(ph, id.as_ptr(), &mut index) },
        ErrorCode::Ok
    );

    let x = [0.0, 100.0, 200.0];
    let y = [50.0, 40.0, 10.0];
    assert_eq!(
        unsafe { EN_setcurve(ph, index, x.as_ptr(), y.as_ptr(), x.len() as c_int) },
        ErrorCode::Ok
    );

    let mut buffer = id_buffer();
    let mut n_points = 0;
    let mut out_x = [0.0; 3];
    let mut out_y = [0.0; 3];
    assert_eq!(
        unsafe {
            EN_getcurve(
                ph,
                index,
                buffer.as_mut_ptr(),
                &mut n_points,
                out_x.as_mut_ptr(),
                out_y.as_mut_ptr(),
            )
        },
        ErrorCode::Ok
    );
    assert_eq!(buffer_to_string(&buffer), "CURVE1");
    assert_eq!(n_points, 3);
    assert_eq!(out_x, x);
    assert_eq!(out_y, y);

    assert_eq!(
        unsafe { EN_setcurvevalue(ph, index, 2, 120.0, 35.0) },
        ErrorCode::Ok
    );
    let (mut px, mut py) = (0.0, 0.0);
    assert_eq!(
        unsafe { EN_getcurvevalue(ph, index, 2, &mut px, &mut py) },
        ErrorCode::Ok
    );
    assert_eq!((px, py), (120.0, 35.0));

    // an unused curve has no type of its own, the pump curve is a pump curve
    let mut curve_type = 0;
    assert_eq!(
        unsafe { EN_getcurvetype(ph, index, &mut curve_type) },
        ErrorCode::Ok
    );
    assert_eq!(curve_type, CurveType::Generic as c_int);

    let pump_curve = CString::new("P1").unwrap();
    let mut pump_curve_index = 0;
    assert_eq!(
        unsafe { EN_getcurveindex(ph, pump_curve.as_ptr(), &mut pump_curve_index) },
        ErrorCode::Ok
    );
    assert_eq!(
        unsafe { EN_getcurvetype(ph, pump_curve_index, &mut curve_type) },
        ErrorCode::Ok
    );
    assert_eq!(curve_type, CurveType::Pump as c_int);

    assert_eq!(unsafe { EN_deletecurve(ph, index) }, ErrorCode::Ok);
}

#[test]
fn demands_round_trip() {
    let project = TestProject::open("tests/pump.inp");
    let ph = project.handle();

    let junction = node_index(ph, "3");

    let mut num_demands = 0;
    assert_eq!(
        unsafe { EN_getnumdemands(ph, junction, &mut num_demands) },
        ErrorCode::Ok
    );
    assert_eq!(num_demands, 1);

    let mut base = 0.0;
    assert_eq!(
        unsafe { EN_getbasedemand(ph, junction, 1, &mut base) },
        ErrorCode::Ok
    );
    assert!((base - 1.0).abs() < 1e-9);

    let name = CString::new("fire flow").unwrap();
    assert_eq!(
        unsafe { EN_adddemand(ph, junction, 3.0, ptr::null(), name.as_ptr()) },
        ErrorCode::Ok
    );
    assert_eq!(
        unsafe { EN_getnumdemands(ph, junction, &mut num_demands) },
        ErrorCode::Ok
    );
    assert_eq!(num_demands, 2);

    let mut index = 0;
    assert_eq!(
        unsafe { EN_getdemandindex(ph, junction, name.as_ptr(), &mut index) },
        ErrorCode::Ok
    );
    assert_eq!(index, 2);

    assert_eq!(
        unsafe { EN_setbasedemand(ph, junction, index, 4.0) },
        ErrorCode::Ok
    );
    assert_eq!(
        unsafe { EN_getbasedemand(ph, junction, index, &mut base) },
        ErrorCode::Ok
    );
    assert!((base - 4.0).abs() < 1e-9);

    // the junction's total demand is the sum of both categories
    assert_eq!(unsafe { EN_solveH(ph) }, ErrorCode::Ok);
    let demand = node_value(ph, junction, NodeProperty::Demand);
    assert!((demand - 5.0).abs() < 1e-3, "expected 1 + 4, got {demand}");

    assert_eq!(
        unsafe { EN_deletedemand(ph, junction, index) },
        ErrorCode::Ok
    );
    assert_eq!(
        unsafe { EN_getnumdemands(ph, junction, &mut num_demands) },
        ErrorCode::Ok
    );
    assert_eq!(num_demands, 1);
}

#[test]
fn controls_round_trip() {
    let project = TestProject::open("tests/2tanks.inp");
    let ph = project.handle();

    let tank = node_index(ph, "1");
    let pipe = link_index(ph, "1");

    let controls_before = count(ph, CountType::ControlCount);

    let mut index = 0;
    assert_eq!(
        unsafe {
            EN_addcontrol(
                ph,
                ControlType::LowLevel as c_int,
                pipe,
                0.0,
                tank,
                3.5,
                &mut index,
            )
        },
        ErrorCode::Ok
    );
    assert_eq!(count(ph, CountType::ControlCount), controls_before + 1);

    let mut control_type = 0;
    let mut control_link = 0;
    let mut setting = 0.0;
    let mut control_node = 0;
    let mut level = 0.0;
    assert_eq!(
        unsafe {
            EN_getcontrol(
                ph,
                index,
                &mut control_type,
                &mut control_link,
                &mut setting,
                &mut control_node,
                &mut level,
            )
        },
        ErrorCode::Ok
    );
    assert_eq!(control_type, ControlType::LowLevel as c_int);
    assert_eq!(control_link, pipe);
    assert_eq!(control_node, tank);
    assert_eq!(setting, 0.0);
    assert!((level - 3.5).abs() < 1e-6);

    // a timer control watches no node at all
    assert_eq!(
        unsafe { EN_setcontrol(ph, index, ControlType::Timer as c_int, pipe, 1.0, 0, 7200.0,) },
        ErrorCode::Ok
    );
    assert_eq!(
        unsafe {
            EN_getcontrol(
                ph,
                index,
                &mut control_type,
                &mut control_link,
                &mut setting,
                &mut control_node,
                &mut level,
            )
        },
        ErrorCode::Ok
    );
    assert_eq!(control_type, ControlType::Timer as c_int);
    assert_eq!(control_node, 0);
    assert_eq!(setting, 1.0);
    assert!((level - 7200.0).abs() < 1e-6);

    assert_eq!(unsafe { EN_deletecontrol(ph, index) }, ErrorCode::Ok);
    assert_eq!(count(ph, CountType::ControlCount), controls_before);
}

#[test]
fn options_and_time_parameters_round_trip() {
    let project = TestProject::open("tests/pump.inp");
    let ph = project.handle();

    let mut units = 0;
    assert_eq!(unsafe { EN_getflowunits(ph, &mut units) }, ErrorCode::Ok);
    assert_eq!(units, crate::ffi::enums::FlowUnits::Cfs as c_int);

    assert_eq!(
        unsafe { EN_settimeparam(ph, TimeParameter::Duration as c_int, 7200) },
        ErrorCode::Ok
    );
    let mut duration: c_long = 0;
    assert_eq!(
        unsafe { EN_gettimeparam(ph, TimeParameter::Duration as c_int, &mut duration) },
        ErrorCode::Ok
    );
    assert_eq!(duration, 7200);

    let mut trials = 0.0;
    assert_eq!(
        unsafe {
            EN_getoption(
                ph,
                crate::ffi::enums::SimOption::Trials as c_int,
                &mut trials,
            )
        },
        ErrorCode::Ok
    );
    assert_eq!(trials, 40.0);

    assert_eq!(
        unsafe { EN_setoption(ph, crate::ffi::enums::SimOption::Trials as c_int, 25.0) },
        ErrorCode::Ok
    );
    assert_eq!(
        unsafe {
            EN_getoption(
                ph,
                crate::ffi::enums::SimOption::Trials as c_int,
                &mut trials,
            )
        },
        ErrorCode::Ok
    );
    assert_eq!(trials, 25.0);
}

#[test]
fn invalid_arguments_are_rejected() {
    let project = TestProject::open("tests/pump.inp");
    let ph = project.handle();

    let mut value = 0.0;
    assert_eq!(
        unsafe { EN_getnodevalue(ph, 999, NodeProperty::Elevation as c_int, &mut value) },
        ErrorCode::UndefinedNode
    );
    assert_eq!(
        unsafe { EN_getlinkvalue(ph, 999, LinkProperty::Diameter as c_int, &mut value) },
        ErrorCode::UndefinedLink
    );
    assert_eq!(
        unsafe { EN_getnodevalue(ph, 1, 9999, &mut value) },
        ErrorCode::InvalidParameterCode
    );

    let mut index = 0;
    let unknown = CString::new("does-not-exist").unwrap();
    assert_eq!(
        unsafe { EN_getnodeindex(ph, unknown.as_ptr(), &mut index) },
        ErrorCode::UndefinedNode
    );

    // a null handle is caught before anything is dereferenced
    assert_eq!(
        unsafe { EN_getnodeindex(ptr::null_mut(), unknown.as_ptr(), &mut index) },
        ErrorCode::InvalidHandle
    );
}

#[test]
fn unsupported_features_report_not_implemented() {
    let project = TestProject::open("tests/pump.inp");
    let ph = project.handle();

    assert_eq!(unsafe { EN_solveQ(ph) }, ErrorCode::NotImplemented);
    assert_eq!(unsafe { EN_openQ(ph) }, ErrorCode::NotImplemented);

    let rule = CString::new("RULE 1\nIF SYSTEM TIME = 0\nTHEN LINK B STATUS IS CLOSED").unwrap();
    assert_eq!(
        unsafe { EN_addrule(ph, rule.as_ptr()) },
        ErrorCode::NotImplemented
    );

    let mut enabled = 0;
    assert_eq!(
        unsafe { EN_getruleenabled(ph, 1, &mut enabled) },
        ErrorCode::NotImplemented
    );

    assert_eq!(
        unsafe { EN_setcurvetype(ph, 1, CurveType::Pump as c_int) },
        ErrorCode::NotImplemented
    );
}

#[test]
fn error_messages_are_written_back() {
    let mut buffer = vec![0 as c_char; 256];
    assert_eq!(
        unsafe { EN_geterror(ErrorCode::UndefinedNode as c_int, buffer.as_mut_ptr(), 255) },
        ErrorCode::Ok
    );
    assert!(buffer_to_string(&buffer).contains("undefined node"));
}
