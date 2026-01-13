use std::fs::File;
use std::io::{BufReader, BufRead};

use crate::model::network::Network;
use crate::model::node::{Node, NodeType};
use crate::model::link::{Link, LinkType};
use crate::model::curve::Curve;
use crate::model::pattern::Pattern;
use crate::model::junction::Junction;
use crate::model::reservoir::Reservoir;
use crate::model::pipe::{Pipe, PipeStatus};
use crate::model::valve::{Valve, ValveType};
// use crate::model::tank::Tank;
use crate::model::pump::Pump;

#[derive(Debug)]
enum ReadState {
  Junctions,
  Pipes,
  Reservoirs,
  Demands,
  Valves,
  Pumps,
  Curves,
  Patterns,
  None,
}

impl Network {
  /// Read a network from an INP file.
  pub fn read_inp(&mut self, inp: &str) -> Result<(), String> {

    // set the initial state to none
    let mut state = ReadState::None;

    // open the INP file
    let file = File::open(inp).or_else(|e| Err(format!("Failed to open file: {}: {}", inp, e)))?;

    // create a new reader
    let mut reader = BufReader::new(file);
    let mut line_buffer = String::with_capacity(512);

    // iterate over the lines in the file
    while reader.read_line(&mut line_buffer).unwrap() > 0 {
      let line = line_buffer.trim();

      if line.starts_with(";") || line.is_empty() {
        // skip comment and empty lines
      }
      // if the line starts with [, it is a new section
      else if line.starts_with("[") {
        state = match line {
          "[JUNCTIONS]" => ReadState::Junctions,
          "[RESERVOIRS]" => ReadState::Reservoirs,
          "[PIPES]" => ReadState::Pipes,
          "[DEMANDS]" => ReadState::Demands,
          "[VALVES]" => ReadState::Valves,
          "[PUMPS]" => ReadState::Pumps,
          "[CURVES]" => ReadState::Curves,
          "[PATTERNS]" => ReadState::Patterns,
          _ => ReadState::None,
        }
      }
      else {

        match state {
          ReadState::Junctions => { 
            let junction = self.read_junction(line);
            self.add_node(junction).unwrap(); 
          }
          ReadState::Valves => { 
            let valve = self.read_valve(line);
            self.add_link(valve).unwrap(); 
          }
          ReadState::Pipes => { 
            let pipe = self.read_pipe(line);
            self.add_link(pipe).unwrap(); 
          }
          ReadState::Reservoirs => { 
            let reservoir = self.read_reservoir(line);
            self.add_node(reservoir).unwrap(); 
          }
          ReadState::Pumps => { 
            let pump = self.read_pump(line);
            self.add_link(pump).unwrap(); 
          }
          ReadState::Curves => { 
            self.read_curve(line);
          }
          ReadState::Demands => { 
            self.read_demand(line); 
          }
          ReadState::Patterns => { 
            self.read_pattern(line); 
          }
          ReadState::None => {
            // skip unknown state
          }
        }
      }
      // clear the line buffer
      line_buffer.clear();
    }
    Ok(())
  }

  /// Read a junction from a parts iterator
  fn read_junction(&mut self, line: &str) -> Node {
    let mut parts = line.split_whitespace();
    // read the junction data
    let id = parts.next().unwrap().into();
    // read the elevation (optional, default 0.0)
    let elevation = parts.next().unwrap().parse::<f64>().unwrap_or(0.0);
    // read the demand (optional, default 0.0)
    let demand = parts.next().and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);

    Node {
      id,
      elevation,
      node_type: NodeType::Junction(Junction { basedemand: demand, pattern: None }),
    }
  }

  /// Read a reservoir from a parts iterator
  fn read_reservoir(&mut self, line: &str) -> Node {
    let mut parts = line.split_whitespace();

    let id = parts.next().unwrap().into();
    // read the elevation
    let elevation = parts.next().unwrap().parse::<f64>().unwrap();
    // add the node to the network
    Node {
      id,
      elevation,
      node_type: NodeType::Reservoir(Reservoir { head_pattern: None }),
    }
  }

  /// Read a valve from a parts iterator
  fn read_valve(&mut self, line: &str) -> Link {
    let mut parts = line.split_whitespace();

    let id = parts.next().unwrap().into();
    // read the start node
    let start_node: Box<str> = parts.next().unwrap().into();
    // read the end node
    let end_node: Box<str> = parts.next().unwrap().into();
    // read the diameter
    let diameter = parts.next().unwrap().parse::<f64>().unwrap();
    // create the link
    let valve_type = match parts.next().unwrap().to_uppercase().as_str() {
      "GPV" => ValveType::GPV,
      "TCV" => ValveType::TCV,
      "FCV" => ValveType::FCV,
      "PRV" => ValveType::PRV,
      "PSV" => ValveType::PSV,
      "PBV" => ValveType::PBV,
      "PCV" => ValveType::PCV,
      _ => panic!("Invalid valve type: {}", parts.next().unwrap()),
    };
    // read the setting
    let setting = parts.next().unwrap().parse::<f64>().unwrap_or(0.0);
    // read the minor loss
    let minor_loss = parts.next().unwrap().parse::<f64>().unwrap_or(0.0);

    // read the curve ID (optional, default none)
    let curve: Option<Box<str>> = parts.next().map(|s| s.to_string().into_boxed_str());

    let start_node_index = *self.node_map.get(&start_node).unwrap();
    let end_node_index = *self.node_map.get(&end_node).unwrap();

    Link {
      id,
      start_node: start_node_index,
      end_node: end_node_index,
      link_type: LinkType::Valve(Valve { diameter, setting, curve, valve_type }),
      minor_loss: minor_loss,
    }
  }

  /// Read a pipe from a parts iterator
  fn read_pipe(&self, line: &str) -> Link {
    let mut parts = line.split_whitespace();

    let id = parts.next().unwrap().into();
    // read the start node
    let start_node: Box<str> = parts.next().unwrap().into();
    // read the end node
    let end_node: Box<str> = parts.next().unwrap().into();
    // read the length
    let length = parts.next().unwrap().parse::<f64>().unwrap();
    // read the diameter
    let diameter = parts.next().unwrap().parse::<f64>().unwrap();
    // read the roughness
    let roughness = parts.next().unwrap().parse::<f64>().unwrap();

    // create the link
    let minor_loss = parts.next().unwrap().parse::<f64>().unwrap_or(0.0);
    // convert minor los
    let minor_loss = 0.02517 * minor_loss / diameter.powi(2) / diameter.powi(2);

    let start_node_index = *self.node_map.get(&start_node).unwrap();
    let end_node_index = *self.node_map.get(&end_node).unwrap();

    Link {
      id,
      start_node: start_node_index,
      end_node: end_node_index,
      minor_loss: minor_loss,
      link_type: LinkType::Pipe(Pipe { diameter, length, roughness, minor_loss, status: PipeStatus::Open }),
    }
  }
  /// Read a pump from a parts iterator
  fn read_pump(&self, line: &str) -> Link {
    let mut parts = line.split_whitespace();

    let id : Box<str> = parts.next().unwrap().into();
    // read the start node
    let start_node: Box<str> = parts.next().unwrap().into();
    // read the end node
    let end_node: Box<str> = parts.next().unwrap().into();

    // read the parameters
    let mut parameters = parts.skip(3);

    let mut speed = 1.0;
    let mut head_curve = "";
    let mut power = 0.0;

    // create the link
    let start_node_index = *self.node_map.get(&start_node).unwrap();
    let end_node_index = *self.node_map.get(&end_node).unwrap();

    while let Some(parameter) = parameters.next() {
      if parameter == ";" {
        continue;
      }
      if let Some(value) = parameters.next() {
        match parameter {
          "SPEED" => speed = value.parse::<f64>().unwrap(),
          "HEAD" => head_curve = value.trim(),
          "POWER" => power = value.parse::<f64>().unwrap(),
          _ => continue,
        }
      }
    }
    let head_curve = head_curve.to_string().into_boxed_str();

    Link {
      id,
      start_node: start_node_index,
      end_node: end_node_index,
      minor_loss: 0.0,
      link_type: LinkType::Pump(Pump { speed, head_curve, power }),
    }
  }
  /// Read a curve from a parts iterator
  /// Appends a point to the curve if it already exists, otherwise creates a new curve
  fn read_curve(&mut self, line: &str) {
    let mut parts = line.split_whitespace();

    let id: Box<str> = parts.next().unwrap().into();
    let x = parts.next().unwrap().parse::<f64>().unwrap();
    let y = parts.next().unwrap().parse::<f64>().unwrap();
    // create curve if it does not exist
    if !self.curves.contains_key(&id) {
      self.curves.insert(id.clone(), Curve { id, x: vec![x], y: vec![y] });
    }
    // otherwise, add the point to the curve
    else {
      let curve = self.curves.get_mut(&id).unwrap();
      // check if the x value is in ascending order
      if x < *curve.x.last().unwrap() {
        panic!("X values must be in ascending order for curve {}", id);
      }
      // add the point to the curve
      curve.x.push(x);
      curve.y.push(y);
    }
  }
  /// Read a pattern from a parts iterator
  fn read_pattern(&mut self, line: &str) {
    let mut parts = line.split_whitespace();
    let id: Box<str> = parts.next().unwrap().into();
    let multipliers = parts.map(|s| s.parse::<f64>().unwrap()).collect();
    if !self.patterns.contains_key(&id) {
      self.patterns.insert(id.clone(), Pattern { id, multipliers });
    }
    else {
      let pattern = self.patterns.get_mut(&id).unwrap();
      pattern.multipliers.extend(multipliers);
    }
  }
  
  /// Read a demand from a parts iterator and set the basedemand for the junction
  fn read_demand(&mut self, line: &str) {
    let mut parts = line.split_whitespace();

    let id: Box<str> = parts.next().unwrap().into();
    let demand = parts.next().unwrap().parse::<f64>().unwrap();
    let node_index = *self.node_map.get(&id).unwrap();
    let node = &mut self.nodes[node_index];
    match &mut node.node_type {
      NodeType::Junction(junction) => {
        junction.basedemand = demand;
      }
      _ => panic!("Demand can only be set for junctions"),
    }
  }
}

