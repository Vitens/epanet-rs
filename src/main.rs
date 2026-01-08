mod input;
mod network;
mod solver;

use std::time::Instant;


fn main() {
  let start_time = Instant::now();
  let mut network = network::Network::from_inp("luke.inp").unwrap();
  let end_time = Instant::now();
  println!("Load time: {:?}", end_time.duration_since(start_time));
  println!("Loaded network with {} nodes and {} links", network.nodes.len(), network.links.len());

  let start_time = Instant::now();
  network.solve().unwrap();

  // check the mass balance
  let total_demand: f64 = network.nodes.iter().map(|n| n.demand).sum();
  let reservoir_supply: f64 = network.links.iter().map(|l| {
    if !matches!(network.nodes[l.end_node].node_type, network::NodeType::Junction { .. }) {
      -l.result.flow
    } 
    else if !matches!(network.nodes[l.start_node].node_type, network::NodeType::Junction { .. }) {
      l.result.flow
    }
    else {
      0.0
    }
  }).sum();

  
  println!("Mass balance: demand = {:.4}, supply = {:.4}, error = {:.2e}", total_demand, reservoir_supply, (total_demand - reservoir_supply).abs());

  // print the heads of the nodes
  for node in network.nodes.iter() {
    println!("Node {} head = {:.4}", node.id, node.result.head);
  }

  // print the flows of the links
  for link in network.links.iter() {
    println!("Link {} flow = {:.4}", link.id, link.result.flow);
  }


  let end_time = Instant::now();
  println!("Solve time: {:?}", end_time.duration_since(start_time));
}