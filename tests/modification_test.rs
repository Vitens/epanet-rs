use epanet_rs::model::network::{Network};
use epanet_rs::model::link::LinkStatus;
use epanet_rs::model::options::HeadlossFormula;
use epanet_rs::model::units::FlowUnits;
use epanet_rs::simulation::Simulation;
use epanet_rs::model::network::modify::*;


#[test]
fn test_network_creation() {

  let mut network = Network::new(FlowUnits::CFS, HeadlossFormula::DarcyWeisbach);

  network.add_reservoir("FH", &ReservoirData { elevation: 100.0, head_pattern: None, coordinates: None, }).unwrap();

  let junctions = vec![
    ("1", 0.0, 1.0),
    ("2", 0.0, 2.0),
    ("3", 0.0, 1.0),
    ("4", 0.0, 2.0),
    ("5", 0.0, 1.0),
    ("6", 0.0, 2.0),
    ("7", 0.0, 1.0),
  ];

  // add the junctions to the network
  for (id, elevation, basedemand) in junctions {
    network.add_junction(id, &JunctionData { elevation, basedemand, emitter_coefficient: 0.0, pattern: None, coordinates: None, }).unwrap();
  }

  // add the pipes to the network
  let pipes = vec![
    ("A", "FH", "1", 100.0, 12.0, 0.6),
    ("B", "1", "3", 100.0, 12.0, 0.6),
    ("C", "1", "2", 100.0, 12.0, 0.6),
    ("D", "2", "4", 100.0, 12.0, 0.6),
    ("E", "3", "4", 100.0, 12.0, 0.6),
    ("F", "4", "5", 100.0, 12.0, 0.6),
    ("G", "4", "6", 100.0, 12.0, 0.6),
    ("H", "6", "7", 100.0, 12.0, 0.6),
  ];

  for (id, start_node, end_node, length, diameter, roughness) in pipes {
    network.add_pipe(id, &PipeData { start_node: start_node.into(), end_node: end_node.into(), length, diameter, roughness, minor_loss: 0.0, check_valve: false, initial_status: LinkStatus::Open, vertices: None, }).unwrap();
  }

  assert_eq!(network.nodes.len(), 8);
  assert_eq!(network.links.len(), 8);

  let mut simulation = Simulation::new(network);
  let result = simulation.solve_hydraulics(false).unwrap();

  // verify initial solution
  assert_eq!(result.heads[0][0], 100.00);
  assert!((result.heads[0][1] - 95.50).abs() < 0.01);

  // update the reservoir elevation 
  simulation.network.update_reservoir("FH", &ReservoirUpdate { elevation: Some(110.0), ..Default::default()}).unwrap();

  // verify updated solution
  let result = simulation.solve_hydraulics(false).unwrap();
  assert_eq!(result.heads[0][0], 110.00);
  assert!((result.heads[0][1] - 105.50).abs() < 0.01);

  // remove node "7"
  simulation.network.remove_node("7", false).unwrap();
  // check that the network now has 7 nodes and 7 links
  assert_eq!(simulation.network.nodes.len(), 7);
  assert_eq!(simulation.network.links.len(), 7);

  let result = simulation.solve_hydraulics(false).unwrap();
  assert!((result.heads[0][1] - 106.35).abs() < 0.01);

  // switch link "G" to a valve 

}