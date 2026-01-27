use std::process::Command;
use std::process::Stdio;
use std::time::Instant;

use clap::{Parser, Subcommand};

use simplelog::{info, warn, error, debug, LevelFilter, TerminalMode, ColorChoice, Config, TermLogger};

use epanet_rs::model::network::Network;
use epanet_rs::solver::HydraulicSolver;
use epanet_rs::utils::binfile::read_outfile;

const BANNER: [&str; 6] = [
  r"  _____ ____   _    _   _ _____ _____     ____  ____  ",
  r" | ____|  _ \ / \  | \ | | ____|_   _|   |  _ \/ ___| ",
  r" |  _| | |_) / _ \ |  \| |  _|   | |_____| |_) \___ \ ",
  r" | |___|  __/ ___ \| |\  | |___  | |_____|  _ < ___) |",
  r" |_____|_| /_/   \_\_| \_|_____| |_|     |_| \_\____/ ",
  r"                                                      "
];

#[derive(Parser, Debug)]
#[command(
  author = "Abel Heinsbroek (Vitens N.V.)",
  version = "0.1.0",
  about = "A very fast, modern and safe re-implementation of the EPANET2 hydraulic solver, written in Rust"
)]
struct Cli {
  #[command(subcommand)]
  command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
  /// Run the hydraulic solver on a network
  Run {
    /// Input file (EPANET .inp format)
    input_file: String,
    /// Output file for results (.rpt format)
    output_file: Option<String>,
    /// Run timesteps in parallel (experimental)
    #[arg(short, long)]
    parallel: bool,
    /// Print verbose output during solving
    #[arg(short, long)]
    verbose: bool,
    /// Print results to stdout
    #[arg(long)]
    print_results: bool,
    /// Suppress all output except for errors
    #[arg(long)]
    quiet: bool,
  },
  /// Convert a network file to a different format
  Convert {
    /// Input file (EPANET .inp format)
    input_file: String,
    /// Output file (.json or .msgpack/.mpk)
    output_file: String,
  },
  /// Validate a network file against EPANET results
  Validate {
    /// Input file to validate
    input_file: String,
    /// Number of decimal places for comparison (default: 2)
    #[arg(short, long, default_value = "2")]
    precision: u32,
  },
}

fn main() -> Result<(), String> {
  let cli = Cli::parse();

  // Determine log level based on command
  let log_level = match &cli.command {
    Commands::Run { quiet, verbose, .. } => {
      if *quiet { LevelFilter::Error }
      else if *verbose { LevelFilter::Debug }
      else { LevelFilter::Info }
    }
    _ => LevelFilter::Info,
  };

  // Initialize the logger with colors
  TermLogger::init(
    log_level,
    Config::default(),
    TerminalMode::Mixed,
    ColorChoice::Auto,
  ).expect("Failed to initialize logger");

  // Run the command
  match cli.command {
    Commands::Run { input_file, output_file, parallel, verbose, print_results, quiet } => {
      run_solver(&input_file, output_file.as_deref(), parallel, verbose, print_results, quiet);
      Ok(())
    }
    Commands::Convert { input_file, output_file } => {
      convert_network(&input_file, &output_file);
      Ok(())
    }
    Commands::Validate { input_file, precision } => {
      if validate_network(&input_file, precision) {
        Ok(())
      } else {
        Err("Validation failed".to_string())
      }

    }
  }
}

/// Run the hydraulic solver on a network
fn run_solver(input_file: &str, output_file: Option<&str>, parallel: bool, verbose: bool, print_results: bool, quiet: bool) {
  let start_time = Instant::now();

  if !quiet {
    println!("{}", BANNER.join("\n"));
  }
  info!("Loading network from file: {}", input_file);

  let mut network = Network::default();
  network.read_file(input_file).expect("Failed to load network");
  let end_time = Instant::now();

  info!("Loaded network with {} nodes and {} links", network.nodes.len(), network.links.len());
  debug!("Network loaded in {:?}", end_time.duration_since(start_time));

  let start_time = Instant::now();
  let solver = HydraulicSolver::new(&network);
  let result = solver.run(parallel, verbose);
  let end_time = Instant::now();
  info!("Solver finished in {:?}", end_time.duration_since(start_time));

  if let Some(output_file) = output_file {
    let start_time = Instant::now();
    network.write_results(&result, output_file).expect("Failed to write results");
    let end_time = Instant::now();
    info!("Results written to {} in {:?}", output_file, end_time.duration_since(start_time));
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
}

/// Convert a network file to a different format
fn convert_network(input_file: &str, output_file: &str) {
  let start_time = Instant::now();

  info!("Loading network from file: {}", input_file);
  let mut network = Network::default();
  network.read_file(input_file).expect("Failed to load network");
  
  let load_time = Instant::now();
  info!("Loaded network with {} nodes and {} links in {:?}", 
    network.nodes.len(), network.links.len(), load_time.duration_since(start_time));

  info!("Converting to: {}", output_file);
  network.save_network(output_file).expect("Failed to save network");
  
  let end_time = Instant::now();
  info!("Network saved in {:?}", end_time.duration_since(load_time));
}

/// Check if two values are equal within a tolerance based on decimal precision
/// For precision=2, tolerance is 0.005 (half of 0.01)
fn values_equal(a: f64, b: f64, precision: u32) -> bool {
  let tolerance = 0.5 * 10_f64.powi(-(precision as i32));
  (a - b).abs() <= tolerance
}

/// Validate the results of a network with EPANET
fn validate_network(input_file: &str, precision: u32) -> bool {
  let tolerance = 0.5 * 10_f64.powi(-(precision as i32));
  info!("Loading network from file: {}", input_file);
  info!("Using precision: {} decimal places (tolerance: {})", precision, tolerance);

  // check if the input file is a .inp file
  if !input_file.ends_with(".inp") {
    error!("Input file must be a .inp file");
    return false;
  }

  let mut network = Network::default();
  network.read_file(input_file).expect("Failed to load network");
  let solver = HydraulicSolver::new(&network);
  let rs_result = solver.run(false, false);

  info!("Running EPANET to validate results");
  
  let mut epanet_process = match Command::new("runepanet")
    .arg(input_file)
    .arg("/dev/null")
    .arg("/tmp/validate.out")
    .stdout(Stdio::null())
    .spawn() {
      Ok(p) => p,
      Err(e) => {
        error!("Failed to run EPANET: {}. Make sure 'runepanet' is installed and in PATH.", e);
        return false;
      }
    };

  let epanet_result = epanet_process.wait().expect("Failed to wait for EPANET");
  if !epanet_result.success() {
    error!("EPANET failed to run");
    return false;
  }

  let epanet_results = read_outfile("/tmp/validate.out");

  let mut head_mismatches = 0;
  let mut flow_mismatches = 0;

  // compare the heads
  for i in 0..rs_result.heads.len() {
    for j in 0..rs_result.heads[i].len() {
      if !values_equal(rs_result.heads[i][j], epanet_results.heads[i][j], precision) {
        if head_mismatches < 5 {
          let diff = (rs_result.heads[i][j] - epanet_results.heads[i][j]).abs();
          warn!("Head mismatch at node '{}' in period {} (RS vs EPA): {} vs {} (diff: {:.6})", 
            epanet_results.node_ids[j], i, rs_result.heads[i][j], epanet_results.heads[i][j], diff);
        }
        head_mismatches += 1;
      }
    }
  }

  // compare the flows
  for i in 0..rs_result.flows.len() {
    for j in 0..rs_result.flows[i].len() {
      if !values_equal(rs_result.flows[i][j], epanet_results.flows[i][j], precision) {
        if flow_mismatches < 5 {
          let diff = (rs_result.flows[i][j] - epanet_results.flows[i][j]).abs();
          warn!("Flow mismatch at link '{}' in period {} (RS vs EPA): {} vs {} (diff: {:.6})", 
            epanet_results.link_ids[j], i, rs_result.flows[i][j], epanet_results.flows[i][j], diff);
        }
        flow_mismatches += 1;
      }
    }
  }

  if head_mismatches > 0 || flow_mismatches > 0 {
    error!("Validation <on-red><b> FAILED </> : {} head mismatches, {} flow mismatches", head_mismatches, flow_mismatches);
    return false;
  } else {
    info!("Validation <on-green><b> PASSED! </>");
    return true;
  }
}
