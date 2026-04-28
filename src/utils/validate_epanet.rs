//! Test helpers that run the reference EPANET2 binary and compare its output against `epanet-rs`.

use std::env;

use std::cmp::Ordering;
use std::process::Command;
use std::process::Stdio;

use simplelog::{error, info, warn};

use crate::simulation::Simulation;
use crate::utils::binfile::read_outfile;

pub fn validate_with_epanet(
    input_file: &str,
    rtol: f64,
    atol: f64,
    max_mismatches: usize,
    parallel: bool,
) -> bool {
    info!("Loading network from file: {}", input_file);
    info!(
        "Using tolerance: rtol={} ({}%), atol={}",
        rtol,
        rtol * 100.0,
        atol
    );

    // check if the input file is a .inp file
    if !input_file.ends_with(".inp") {
        error!("Input file must be a .inp file");
        return false;
    }

    let mut simulation = Simulation::from_file(input_file).unwrap_or_else(|e| {
        error!("Failed to load network: {}", e);
        std::process::exit(1);
    });
    simulation.skip_timesteps = false;
    let rs_result = match simulation.solve_hydraulics(parallel) {
        Ok(result) => result,
        Err(e) => {
            error!("Solver failed: {}", e);
            return false;
        }
    };

    info!("Running EPANET to validate results");

    // null file for EPANET output (NUL on Windows, /dev/null on Unix)
    let null_file = if cfg!(windows) { "NUL" } else { "/dev/null" };
    // temporary file for EPANET output
    let temp_file = env::temp_dir().join("validate.out");
    let temp_file_str = temp_file.to_str().unwrap();

    let mut epanet_process = match Command::new("runepanet")
        .arg(input_file)
        .arg(null_file)
        .arg(temp_file_str)
        .stdout(Stdio::null())
        .spawn()
    {
        Ok(p) => p,
        Err(e) => {
            error!(
                "Failed to run EPANET: {}. Make sure 'runepanet' is installed and in PATH.",
                e
            );
            return false;
        }
    };

    let epanet_result = epanet_process.wait().expect("Failed to wait for EPANET");
    if !epanet_result.success() {
        error!("EPANET failed to run");
        return false;
    }

    let epanet_results = read_outfile(temp_file_str);

    let mut head_mismatches = Vec::new();
    let mut flow_mismatches = Vec::new();
    let mut demand_mismatches = Vec::new();

    // compare the heads
    for i in 0..rs_result.heads.len() {
        for j in 0..rs_result.heads[i].len() {
            if !values_equal(
                rs_result.heads[i][j],
                epanet_results.heads[i][j],
                rtol,
                atol,
            ) {
                let diff = (rs_result.heads[i][j] - epanet_results.heads[i][j]).abs();
                let allowed = atol + rtol * epanet_results.heads[i][j].abs();
                head_mismatches.push(ValidationMismatch {
                    id: &epanet_results.node_ids[j],
                    period: i,
                    rs_value: rs_result.heads[i][j],
                    epanet_value: epanet_results.heads[i][j],
                    diff,
                    allowed,
                });
            }
        }
    }

    // compare the flows
    for i in 0..rs_result.flows.len() {
        for j in 0..rs_result.flows[i].len() {
            if !values_equal(
                rs_result.flows[i][j],
                epanet_results.flows[i][j],
                rtol,
                atol,
            ) {
                let diff = (rs_result.flows[i][j] - epanet_results.flows[i][j]).abs();
                let allowed = atol + rtol * epanet_results.flows[i][j].abs();
                flow_mismatches.push(ValidationMismatch {
                    id: &epanet_results.link_ids[j],
                    period: i,
                    rs_value: rs_result.flows[i][j],
                    epanet_value: epanet_results.flows[i][j],
                    diff,
                    allowed,
                });
            }
        }
    }

    // compare the node demands
    for i in 0..rs_result.demands.len() {
        for j in 0..rs_result.demands[i].len() {
            if !values_equal(
                rs_result.demands[i][j],
                epanet_results.demands[i][j],
                rtol,
                atol,
            ) {
                // skip results for fixed_head nodes
                if simulation.network.nodes[j].is_fixed() {
                    continue;
                }
                let diff = (rs_result.demands[i][j] - epanet_results.demands[i][j]).abs();
                let allowed = atol + rtol * epanet_results.demands[i][j].abs();
                demand_mismatches.push(ValidationMismatch {
                    id: &epanet_results.node_ids[j],
                    period: i,
                    rs_value: rs_result.demands[i][j],
                    epanet_value: epanet_results.demands[i][j],
                    diff,
                    allowed,
                });
            }
        }
    }

    log_top_mismatches("Head", "node", &mut head_mismatches, max_mismatches);
    log_top_mismatches("Flow", "link", &mut flow_mismatches, max_mismatches);
    log_top_mismatches("Demand", "node", &mut demand_mismatches, max_mismatches);

    if !head_mismatches.is_empty() || !flow_mismatches.is_empty() || !demand_mismatches.is_empty() {
        error!(
            "Validation <on-red><b> FAILED </> : {} head mismatches, {} flow mismatches, {} demand mismatches",
            head_mismatches.len(),
            flow_mismatches.len(),
            demand_mismatches.len()
        );
        false
    } else {
        info!("Validation <on-green><b> PASSED! </>");
        true
    }
}

struct ValidationMismatch<'a> {
    id: &'a str,
    period: usize,
    rs_value: f64,
    epanet_value: f64,
    diff: f64,
    allowed: f64,
}

fn log_top_mismatches(
    mismatch_type: &str,
    entity_type: &str,
    mismatches: &mut [ValidationMismatch],
    max_mismatches: usize,
) {
    mismatches.sort_by(|a, b| b.diff.partial_cmp(&a.diff).unwrap_or(Ordering::Equal));

    for mismatch_item in mismatches.iter().take(max_mismatches) {
        warn!(
            "{} mismatch at {} '{}' in period {} (RS vs EPA): {} vs {} (diff: {:.4}, allowed: {:.4})",
            mismatch_type,
            entity_type,
            mismatch_item.id,
            mismatch_item.period,
            mismatch_item.rs_value,
            mismatch_item.epanet_value,
            mismatch_item.diff,
            mismatch_item.allowed
        );
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
