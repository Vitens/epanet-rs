use std::time::Instant;
use clap::Parser;

use epanet_rs::model::network::Network;
use epanet_rs::solver::HydraulicSolver;

const BANNER: [&str; 6] = [r"  _____ ____   _    _   _ _____ _____     ____  ____  ", r" | ____|  _ \ / \  | \ | | ____|_   _|   |  _ \/ ___| ", r" |  _| | |_) / _ \ |  \| |  _|   | |_____| |_) \___ \ ", r" | |___|  __/ ___ \| |\  | |___  | |_____|  _ < ___) |", r" |_____|_| /_/   \_\_| \_|_____| |_|     |_| \_\____/ ", r"                                                      "];
#[derive(Parser, Debug)]
#[command(author="Abel Heinsbroek (Vitens N.V.)", version = "0.1.0", about, long_about = "A very fast, modern and safe re-implementation of the EPANET2 hydraulic solver, written in Rust")]
struct Args {
  input_file: String,
  output_file: Option<String>,
  network_file: Option<String>,
  #[arg(short, long)]
  parallel: bool,
  #[arg(short, long)]
  verbose: bool,
  #[arg(long)]
  print_results: bool,
}

fn main() {

  let args = Args::parse();

  let input_file = args.input_file;
  let parallel = args.parallel;
  let output_file = args.output_file;
  let verbose = args.verbose;
  let print_results = args.print_results;
  let network_file = args.network_file;

  let start_time = Instant::now();

  println!("{}", BANNER.join("\n"));
  println!("Loading network from file: {}", input_file);

  let mut network = Network::default();
  network.read_file(&input_file.as_str()).expect("Failed to load network");
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

  if print_results {
    println!("Results:");
    println!("=== Heads:");
    for (i, node) in network.nodes.iter().enumerate() {
      println!("Node {}: {:.2}", node.id, result.heads[0][i]);
    }
    println!("=== Flows:");
    for (i, link) in network.links.iter().enumerate() {
      println!("Link {}: {:.2}", link.id, result.flows[0][i]);
    }
  }

  if let Some(network_file) = network_file {
    let start_time = Instant::now();
    network.save_network(network_file.as_str()).expect("Failed to save network");
    let end_time = Instant::now();
    println!("Network saved in {:?}", end_time.duration_since(start_time));
  }
}