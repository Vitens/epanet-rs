mod input;
mod model;
mod solver;

use std::{env, time::Instant};

use model::network::Network;
use solver::HydraulicSolver;

const BANNER: [&str; 6] = [r"  _____ ____   _    _   _ _____ _____     ____  ____  ", r" | ____|  _ \ / \  | \ | | ____|_   _|   |  _ \/ ___| ", r" |  _| | |_) / _ \ |  \| |  _|   | |_____| |_) \___ \ ", r" | |___|  __/ ___ \| |\  | |___  | |_____|  _ < ___) |", r" |_____|_| /_/   \_\_| \_|_____| |_|     |_| \_\____/ ", r"                                                      "];
fn main() {
  let start_time = Instant::now();
  let file = env::args().nth(1).unwrap_or("networks/luke.inp".to_string());
  let parallel = env::args().nth(2).unwrap_or("false".to_string()) == "true";

  println!("{}", BANNER.join("\n"));
  println!("Loading network from file: {}", file);

  let mut network = Network::default();
  network.read_inp(&file.as_str()).expect("Failed to load network");
  let end_time = Instant::now();
  println!("Loaded network with {} nodes and {} links", network.nodes.len(), network.links.len());
  println!("Network loaded in {:?}", end_time.duration_since(start_time));

  let solver = HydraulicSolver::new(&network);
  let result = solver.run(parallel);
  println!("Solver finished in {:?}", end_time.duration_since(start_time));
  println!("Result: {:?}", result.heads[0]);
}