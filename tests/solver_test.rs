//! Integration test for the hydraulic solver using pump.inp

use epanet_rs::model::network::Network;
use epanet_rs::solver::simulation::Simulation;
use epanet_rs::solver::result::SolverResult;

fn verify_heads_and_flows(network: &Network, result: &SolverResult, expected_heads: &Vec<(&str, f64)>, expected_flows: &Vec<(&str, f64)>) {

  let result_length = result.heads.len();
  // Verify heads
  for (node_id, expected_head) in expected_heads {
    let idx = *network.node_map.get(*node_id).expect(&format!("Node {} not found", node_id));
    let actual_head = result.heads[result_length-1][idx];
    assert!(
      (actual_head - expected_head).abs() < 0.01,
      "Head mismatch for node {}: expected {:.2}, got {:.2}",
      node_id, expected_head, actual_head
    );
  }

  // Verify flows
  for (link_id, expected_flow) in expected_flows {
    let idx = *network.link_map.get(*link_id).expect(&format!("Link {} not found", link_id));
    let actual_flow = result.flows[result_length-1][idx];
    assert!(
      (actual_flow - expected_flow).abs() < 0.01,
      "Flow mismatch for link {}: expected {:.2}, got {:.2}",
      link_id, expected_flow, actual_flow
    );
  }
}

/// Test solving pump.inp and verify exact head and flow values
#[test]
fn test_solve_pump_network() {
    let mut network = Network::default();
    network.read_inp("tests/pump.inp").expect("Failed to load pump.inp");

    let mut simulation = Simulation::new(&network);
    let result = simulation.solve_hydraulics(false);

    // Expected heads (in feet)
    let expected_heads: Vec<(&str, f64)> = vec![
        ("1", 166.00),
        ("2", 164.35),
        ("3", 164.61),
        ("4", 163.76),
        ("5", 163.67),
        ("6", 163.05),
        ("7", 162.95),
        ("FH", 100.00),
        ("FH2", 85.00),
    ];

    // Expected flows (in CFS)
    let expected_flows: Vec<(&str, f64)> = vec![
        ("B", 4.29),
        ("C", 4.71),
        ("D", 2.71),
        ("E", 3.29),
        ("F", 1.00),
        ("G", 3.00),
        ("H", 1.00),
        ("I", 0.00),
        ("1", 10.00),  // pump
    ];

    verify_heads_and_flows(&network, &result, &expected_heads, &expected_flows);
}

/// Test solving valve.inp and verify exact head and flow values
#[test]
fn test_solve_valve_network() {
  let mut network = Network::default();
  network.read_inp("tests/valves.inp").expect("Failed to load valves.inp");

  let mut simulation = Simulation::new(&network);
  let result = simulation.solve_hydraulics(false);

  // Expected heads (in CFS)
  let expected_flows: Vec<(&str, f64)> = vec![
    ("2", 1.000000),
    ("3", 1.000000),
    ("4", 1.000000),
    ("5", 5.000000),
    ("6", 50.000001),
    ("7", 50.000001),
    ("8", 76.813932),
    ("9", 76.813789),
    ("PRV", 1.000000),
    ("PBV", 1.000000),
    ("TCV", 1.000000),
    ("PCV", 5.000000),
    ("FCV", 50.000001),
    ("PSV", 76.813931),
    ("FCV-WARN", 50.0000),
  ];

  let expected_heads: Vec<(&str, f64)> = vec![
    ("2", 100.000000),
    ("3", 100.000000),
    ("4", 100.000000),
    ("5", 99.999998),
    ("PRV-o", 23.078698),
    ("PBV-o", 76.921302),
    ("TCV-o", 97.482999),
    ("PCV-o", 49.659997),
    ("9", 99.999985),
    ("10", 0.000003),
    ("7", 56.157397),
    ("8", 0.296530),
    ("1", 100.000000),
    ("6", 0.000000),
    ("13", 100.000),
    ("14", 100.000),
    ("15", 90.000),
    ("16", 49.659997)
  ];

  verify_heads_and_flows(&network, &result, &expected_heads, &expected_flows);
}

#[test]
fn test_solve_tanks_network() {
  let mut network = Network::default();
  network.read_inp("tests/tanks.inp").expect("Failed to load tanks.inp");

  let mut simulation = Simulation::new(&network);
  let result = simulation.solve_hydraulics(false);

  let expected_heads: Vec<(&str, f64)> = vec![
    ("1", 15.00),
    ("3", 25.00),
    ("4", 10.00),
    ("5", 5.00),
    ("6", 5.00),
  ];

  let expected_flows: Vec<(&str, f64)> = vec![
    ("1", 0.00),
    ("2", 0.00),
    ("3", 8.58),
    ("4", 15.52),
    ("5", 29.73),
  ];

  verify_heads_and_flows(&network, &result, &expected_heads, &expected_flows);
}

/// Test solving 2tanks.inp and verify exact head and flow values
#[test]
fn test_solve_2tanks_controls_network() {
  let mut network = Network::default();
  network.read_inp("tests/2tanks-controls.inp").expect("Failed to load 2tanks-controls.inp");

  let mut simulation = Simulation::new(&network);
  let result = simulation.solve_hydraulics(false);

  let expected_heads: Vec<(&str, f64)> = vec![
    ("1", 5.00),
    ("3", 4.00),
    ("2", 3.86)
  ];
  let expected_flows: Vec<(&str, f64)> = vec![
    ("1", 0.00),
    ("2", 1.00)
  ];
  verify_heads_and_flows(&network, &result, &expected_heads, &expected_flows);
}

#[test]
fn test_solve_emitters_network() {
  let mut network = Network::default();
  network.read_inp("tests/emitter.inp").expect("Failed to load emitters.inp");

  let mut simulation = Simulation::new(&network);
  let result = simulation.solve_hydraulics(false);

  let expected_heads: Vec<(&str, f64)> = vec![
    ("1", 10.00),
    ("2", 6.41),
    ("3", 4.85),
    ("4", 3.13),
    ("5", 2.84),
  ];
  let expected_flows: Vec<(&str, f64)> = vec![
    ("1", 8.61),
    ("2", 2.56),
    ("3", 2.70),
    ("4", 1.03),
    ("5", -4.00),
  ];

  verify_heads_and_flows(&network, &result, &expected_heads, &expected_flows);
}

#[test]
fn test_solve_pda_network() {
  let mut network = Network::default();
  network.read_inp("tests/pda.inp").expect("Failed to load pda.inp");

  let mut simulation = Simulation::new(&network);
  let result = simulation.solve_hydraulics(false);

  let expected_flows: Vec<(&str, f64)> = vec![
    ("A", 35.994202),
    ("B", 15.244768),
    ("C", 16.814336),
    ("D", 9.042474),
    ("E", 11.350634),
    ("F", 3.866947),
    ("G", 8.785734),
    ("H", 1.075225),
  ];

  let expected_heads: Vec<(&str, f64)> = vec![
    ("1", 1.548499),
    ("2", 1.510046),
    ("3", 1.516428),
    ("4", 1.497855),
    ("5", 1.495327),
    ("6", 1.486298),
    ("7", 0.028894),
    ("FH", 100.000000)
  ];

  verify_heads_and_flows(&network, &result, &expected_heads, &expected_flows);
}
