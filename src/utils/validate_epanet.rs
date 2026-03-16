use std::env;
use std::path::PathBuf;

use std::process::Command;
use std::process::Stdio;

use simplelog::{info, warn, error};

use crate::model::network::Network;
use crate::solver::solver::HydraulicSolver;
use crate::utils::binfile::read_outfile;

pub fn validate_with_epanet(input_file: &str, rtol: f64, atol: f64, parallel: bool) -> bool {
  info!("Loading network from file: {}", input_file);
  info!("Using tolerance: rtol={} ({}%), atol={}", rtol, rtol * 100.0, atol);

  // check if the input file is a .inp file
  if !input_file.ends_with(".inp") {
    error!("Input file must be a .inp file");
    return false;
  }

  let mut network = Network::default();
  network.read_file(input_file).unwrap_or_else(|e| {
    error!("Failed to load network: {}", e);
    std::process::exit(1);
  });
  let mut solver = HydraulicSolver::new(&network);
  // set the validate_with_epanet flag to true to match epanet timestep behaviour
  solver.skip_timesteps = false;
  let rs_result = solver.run(parallel);

  info!("Running EPANET to validate results");

  // null file for EPANET output (NUL on Windows, /dev/null on Unix)
  let null_file = if cfg!(windows) { "NUL" } else { "/dev/null" };
  // temporary file for EPANET output
  let temp_file = PathBuf::from(env::temp_dir()).join("validate.out");
  let temp_file_str = temp_file.to_str().unwrap();
  
  let mut epanet_process = match Command::new("runepanet")
    .arg(input_file)
    .arg(null_file)
    .arg(temp_file_str)
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

  let epanet_results = read_outfile(temp_file_str);

  let mut head_mismatches = 0;
  let mut flow_mismatches = 0;
  let mut demand_mismatches = 0;

  // compare the heads
  for i in 0..rs_result.heads.len() {
    for j in 0..rs_result.heads[i].len() {
      if !values_equal(rs_result.heads[i][j], epanet_results.heads[i][j], rtol, atol) {
        if head_mismatches < 5 {
          let diff = (rs_result.heads[i][j] - epanet_results.heads[i][j]).abs();
          let allowed = atol + rtol * epanet_results.heads[i][j].abs();
          warn!("Head mismatch at node '{}' in period {} (RS vs EPA): {} vs {} (diff: {:.4}, allowed: {:.4})", 
            epanet_results.node_ids[j], i, rs_result.heads[i][j], epanet_results.heads[i][j], diff, allowed);
        }
        head_mismatches += 1;
      }
    }
  }

  // compare the flows
  for i in 0..rs_result.flows.len() {
    for j in 0..rs_result.flows[i].len() {
      if !values_equal(rs_result.flows[i][j], epanet_results.flows[i][j], rtol, atol) {
        let diff = (rs_result.flows[i][j] - epanet_results.flows[i][j]).abs();
        let allowed = atol + rtol * epanet_results.flows[i][j].abs();
        if flow_mismatches < 5 {
          warn!("Flow mismatch at link '{}' in period {} (RS vs EPA): {} vs {} (diff: {:.4}, allowed: {:.4})", 
            epanet_results.link_ids[j], i, rs_result.flows[i][j], epanet_results.flows[i][j], diff, allowed);
        }
        flow_mismatches += 1;
      }
    }
  }

  // compare the node demands
  for i in 0..rs_result.demands.len() {
    for j in 0..rs_result.demands[i].len() {
      if !values_equal(rs_result.demands[i][j], epanet_results.demands[i][j], rtol, atol) {
        // skip results for fixed_head nodes
        if network.nodes[j].is_fixed() {
          continue;
        }
        let diff = (rs_result.demands[i][j] - epanet_results.demands[i][j]).abs();
        let allowed = atol + rtol * epanet_results.demands[i][j].abs();
        if demand_mismatches < 5 {
          warn!("Demand mismatch at node '{}' in period {} (RS vs EPA): {} vs {} (diff: {:.4}, allowed: {:.4})", 
            epanet_results.node_ids[j], i, rs_result.demands[i][j], epanet_results.demands[i][j], diff, allowed);
        }
        demand_mismatches += 1;
      }
    }
  }

  if head_mismatches > 0 || flow_mismatches > 0 || demand_mismatches > 0 {
    error!(
      "Validation <on-red><b> FAILED </> : {} head mismatches, {} flow mismatches, {} demand mismatches",
      head_mismatches,
      flow_mismatches,
      demand_mismatches
    );
    return false;
  } else {
    info!("Validation <on-green><b> PASSED! </>");
    return true;
  }
}

/// Check if two values are equal using hybrid tolerance (like numpy.isclose)
/// |a - b| <= atol + rtol * max(|a|, |b|)
/// - atol: absolute tolerance (dominates for small values near zero)
/// - rtol: relative tolerance (dominates for larger values)
fn values_equal(a: f64, b: f64, rtol: f64, atol: f64) -> bool {
  let diff = (a - b).abs();
  let max_val = a.abs().max(b.abs());
  diff <= atol + rtol * max_val
}