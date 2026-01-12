mod input;
mod network;
mod solver;

use std::{env, time::Instant};

const BANNER: [&str; 6] = [r"  _____ ____   _    _   _ _____ _____     ____  ____  ", r" | ____|  _ \ / \  | \ | | ____|_   _|   |  _ \/ ___| ", r" |  _| | |_) / _ \ |  \| |  _|   | |_____| |_) \___ \ ", r" | |___|  __/ ___ \| |\  | |___  | |_____|  _ < ___) |", r" |_____|_| /_/   \_\_| \_|_____| |_|     |_| \_\____/ ", r"                                                      "];
fn main() {
  let start_time = Instant::now();
  let file = env::args().nth(1).unwrap_or("grid100.inp".to_string());
  let parallel = env::args().nth(2).unwrap_or("false".to_string()) == "true";

  println!("{}", BANNER.join("\n"));
  println!("Loading network from file: {}", file);

  let mut network = network::Network::default();
  network.read_inp(&file.as_str()).expect("Failed to load network");
  let end_time = Instant::now();
  println!("Loaded network with {} nodes and {} links", network.nodes.len(), network.links.len());
  println!("Network loaded in {:?}", end_time.duration_since(start_time));

  let start_time = Instant::now();
  let solver = solver::HydraulicSolver::new(&network);
  solver.run(parallel);

  let end_time = Instant::now();

  // // print the heads of the nodes
  // for node in network.nodes.iter() {
  //   println!("Node {} head = {:.4}", node.id, node.result.head * 1.0/network::UCF_H);
  // }

  // // print the flows of the links
  // for link in network.links.iter() {
  //   println!("Link {} flow = {:.4}", link.id, link.result.flow * 1.0/network::UCF_Q);
  // }


  println!("Solve time: {:?}", end_time.duration_since(start_time));
}