use std::time::Instant;

use clap::{Parser, Subcommand};

use simplelog::format_description;
use simplelog::{
    ColorChoice, ConfigBuilder, LevelFilter, TermLogger, TerminalMode, debug, error, info,
};

use epanet_rs::model::network::Network;
use epanet_rs::simulation::Simulation;
use epanet_rs::utils::validate_epanet::validate_with_epanet;

const BANNER: [&str; 6] = [
    r"  _____ ____   _    _   _ _____ _____     ____  ____  ",
    r" | ____|  _ \ / \  | \ | | ____|_   _|   |  _ \/ ___| ",
    r" |  _| | |_) / _ \ |  \| |  _|   | |_____| |_) \___ \ ",
    r" | |___|  __/ ___ \| |\  | |___  | |_____|  _ < ___) |",
    r" |_____|_| /_/   \_\_| \_|_____| |_|     |_| \_\____/ ",
    r"                                                      ",
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
        /// Input file (.json, .msgpack or .inp format)
        input_file: String,
        /// Output file (.json, .msgpack or .mpk format)
        output_file: String,
    },
    /// Validate a network file against EPANET results
    Validate {
        /// Input file to validate
        input_file: String,
        /// Relative tolerance (default: 0.001 = 0.1%)
        #[arg(short, long, default_value = "0.001")]
        rtol: f64,
        /// Absolute tolerance (default: 0.01)
        #[arg(short, long, default_value = "0.01")]
        atol: f64,
        /// Run in parallel
        #[arg(short = 'P', long, default_value = "false")]
        parallel: bool,
    },
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();

    // Determine log level based on command
    let log_level = match &cli.command {
        Commands::Run { quiet, verbose, .. } => {
            if *quiet {
                LevelFilter::Error
            } else if *verbose {
                LevelFilter::Debug
            } else {
                LevelFilter::Info
            }
        }
        _ => LevelFilter::Info,
    };

    let logconfig = ConfigBuilder::new()
        .set_time_format_custom(format_description!(
            "[hour]:[minute]:[second].[subsecond digits:3]"
        ))
        .set_location_level(LevelFilter::Trace)
        .set_target_level(LevelFilter::Trace)
        .build();

    // Initialize the logger with colors
    TermLogger::init(log_level, logconfig, TerminalMode::Mixed, ColorChoice::Auto)
        .expect("Failed to initialize logger");

    // Run the command
    match cli.command {
        Commands::Run {
            input_file,
            output_file,
            parallel,
            print_results,
            quiet,
            ..
        } => {
            run_solver(
                &input_file,
                output_file.as_deref(),
                parallel,
                print_results,
                quiet,
            );
            Ok(())
        }
        Commands::Convert {
            input_file,
            output_file,
        } => {
            convert_network(&input_file, &output_file);
            Ok(())
        }
        Commands::Validate {
            input_file,
            rtol,
            atol,
            parallel,
        } => {
            if validate_network(&input_file, rtol, atol, parallel) {
                Ok(())
            } else {
                Err("Validation failed".to_string())
            }
        }
    }
}

/// Run the hydraulic solver on a network
fn run_solver(
    input_file: &str,
    output_file: Option<&str>,
    parallel: bool,
    print_results: bool,
    quiet: bool,
) {
    let start_time = Instant::now();

    if !quiet {
        println!("{}", BANNER.join("\n"));
    }
    info!("Loading network from file: {}", input_file);

    let mut network = Network::default();
    network.read_file(input_file).unwrap_or_else(|e| {
        error!("Failed to load network: {}", e);
        std::process::exit(1);
    });
    let end_time = Instant::now();

    info!(
        "Loaded network with {} nodes and {} links",
        network.nodes.len(),
        network.links.len()
    );
    debug!(
        "Network loaded in {:?}",
        end_time.duration_since(start_time)
    );

    let start_time = Instant::now();
    let mut simulation = Simulation::new(network);

    simulation.initialize_hydraulics().unwrap_or_else(|e| {
        error!("Failed to initialize simulation: {}", e);
        std::process::exit(1);
    });

    let result = simulation.solve_hydraulics(parallel).unwrap_or_else(|e| {
        error!("Solver failed: {}", e);
        std::process::exit(1);
    });
    let end_time = Instant::now();
    info!(
        "Solver finished in {:?}",
        end_time.duration_since(start_time)
    );

    if let Some(output_file) = output_file {
        let start_time = Instant::now();
        simulation
            .network
            .write_results(&result, output_file)
            .expect("Failed to write results");
        let end_time = Instant::now();
        info!(
            "Results written to {} in {:?}",
            output_file,
            end_time.duration_since(start_time)
        );
    }

    if print_results {
        println!("Results:");
        println!("=== Heads:");
        for (i, node) in simulation.network.nodes.iter().enumerate() {
            println!(
                "Node {}: {:.2}",
                node.id,
                result.heads[result.heads.len() - 1][i]
            );
        }
        println!("=== Flows:");
        for (i, link) in simulation.network.links.iter().enumerate() {
            println!(
                "Link {}: {:.2}",
                link.id,
                result.flows[result.flows.len() - 1][i]
            );
        }
    }
}

/// Convert a network file to a different format
fn convert_network(input_file: &str, output_file: &str) {
    let start_time = Instant::now();

    info!("Loading network from file: {}", input_file);
    let mut network = Network::default();
    network.read_file(input_file).unwrap_or_else(|e| {
        error!("Failed to load network: {}", e);
        std::process::exit(1);
    });

    let load_time = Instant::now();
    info!(
        "Loaded network with {} nodes and {} links in {:?}",
        network.nodes.len(),
        network.links.len(),
        load_time.duration_since(start_time)
    );

    info!("Converting to: {}", output_file);
    network
        .save_network(output_file)
        .expect("Failed to save network");

    let end_time = Instant::now();
    info!("Network saved in {:?}", end_time.duration_since(load_time));
}

/// Validate a network file against EPANET results
fn validate_network(input_file: &str, rtol: f64, atol: f64, parallel: bool) -> bool {
    validate_with_epanet(input_file, rtol, atol, parallel)
}
