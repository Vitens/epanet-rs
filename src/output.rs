use std::fs::File;
use std::io::BufWriter;
use crate::model::network::Network;
use crate::solver::SolverResult;
use std::collections::HashMap;
use serde::Serialize;
use rmp_serde::{Serializer};

#[derive(Serialize)]
struct JsonOutput {
  nodes: Vec<String>,
  links: Vec<String>,
  heads: Vec<Vec<f64>>,
  flows: Vec<Vec<f64>>,
}


const DIGITS: usize = 3;

// helper function to round to a given number of digits (prevent JSON file bloat due to floating point precision)
fn round_to_digits(value: f64, digits: usize) -> f64 {
  let factor = 10.0_f64.powi(digits as i32);
  (value * factor).round() / factor
}

impl Network {

  pub fn write_results(&self, results: &SolverResult, file: &str) -> Result<(), String> {
    // get file extension
    let file_extension = file.split('.').last().unwrap();

    let file = File::create(file).map_err(|e| format!("Failed to create output file: {}", e))?;
    let writer = BufWriter::new(file);

    let nodes = self.nodes.iter().map(|n| n.id.to_string()).collect();
    let links = self.links.iter().map(|l| l.id.to_string()).collect();

    let heads = results.heads.iter().map(|h| h.iter().map(|h| round_to_digits(*h, DIGITS)).collect()).collect();
    let flows = results.flows.iter().map(|f| f.iter().map(|f| round_to_digits(*f, DIGITS)).collect()).collect();

    let output = JsonOutput { nodes, links, heads, flows };
    if file_extension == "json" {
      serde_json::to_writer(writer, &output).map_err(|e| format!("Failed to write results to file: {}", e))?;
    } else if file_extension == "mpk" || file_extension == "msgpack" {
      let mut serializer = Serializer::new(writer);
      output.serialize(&mut serializer).map_err(|e| format!("Failed to write results to file: {}", e))?;
    } else {
      return Err(format!("Unsupported file extension: {}", file_extension));
    }

    Ok(())

  }
}

