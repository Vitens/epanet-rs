mod input;
mod model;
mod solver;
mod constants;
mod output;

use std::{env, time::Instant};
use clap::Parser;

use model::network::Network;
use solver::HydraulicSolver;

const BANNER: [&str; 6] = [r"  _____ ____   _    _   _ _____ _____     ____  ____  ", r" | ____|  _ \ / \  | \ | | ____|_   _|   |  _ \/ ___| ", r" |  _| | |_) / _ \ |  \| |  _|   | |_____| |_) \___ \ ", r" | |___|  __/ ___ \| |\  | |___  | |_____|  _ < ___) |", r" |_____|_| /_/   \_\_| \_|_____| |_|     |_| \_\____/ ", r"                                                      "];

#[derive(Parser, Debug)]
#[command(author="Abel Heinsbroek (Vitens N.V.)", version = "0.1.0", about, long_about = "A very fast, modern and safe re-implementation of the EPANET2 hydraulic solver, written in Rust")]
struct Args {
  input_file: String,
  output_file: Option<String>,
  #[arg(short, long)]
  parallel: bool,
  #[arg(short, long)]
  verbose: bool,
}

fn main() {

  let args = Args::parse();

  let input_file = args.input_file;
  let parallel = args.parallel;
  let output_file = args.output_file;
  let verbose = args.verbose;

  let start_time = Instant::now();

  println!("{}", BANNER.join("\n"));
  println!("Loading network from file: {}", input_file);

  let mut network = Network::default();
  network.read_inp(&input_file.as_str()).expect("Failed to load network");
  let end_time = Instant::now();
  println!("Loaded network with {} nodes and {} links", network.nodes.len(), network.links.len());
  println!("Network loaded in {:?}", end_time.duration_since(start_time));

  let start_time = Instant::now();
  let solver = HydraulicSolver::new(&network);
  let result = solver.run(parallel, verbose);
  let end_time = Instant::now();
  println!("Solver finished in {:?}", end_time.duration_since(start_time));

  if let Some(output_file) = output_file {
    let start_time = Instant::now();
    network.write_results(&result, &output_file).expect("Failed to write results");
    let end_time = Instant::now();
    println!("Results written in {:?}", end_time.duration_since(start_time));
  }
}