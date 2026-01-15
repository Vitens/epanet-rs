use std::fs::File;
use std::io::{BufReader, BufRead};
use std::str::FromStr;

use crate::model::network::Network;
use crate::model::node::{Node, NodeType};
use crate::model::link::{Link, LinkType, LinkStatus};
use crate::model::curve::Curve;
use crate::model::pattern::Pattern;
use crate::model::junction::Junction;
use crate::model::reservoir::Reservoir;
use crate::model::pipe::Pipe;
use crate::model::pump::Pump;
use crate::model::valve::{Valve, ValveType};
use crate::model::units::{FlowUnits, UnitSystem, PressureUnits, UnitConversion};
use crate::model::options::*;
// use crate::model::tank::Tank;

#[derive(Debug)]
enum ReadState {
  Junctions,
  Pipes,
  Reservoirs,
  Demands,
  Valves,
  Pumps,
  Curves,
  Options,
  Patterns,
  Status,
  Times,
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
          "[OPTIONS]" => ReadState::Options,
          "[STATUS]" => ReadState::Status,
          "[TIMES]" => ReadState::Times,
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
          },
          ReadState::Options => {
            self.read_options(line);
          }
          ReadState::Times => {
            self.read_times(line);
          }
          ReadState::Status => {
            self.read_status(line);
          }
        }
      }
      // clear the line buffer
      line_buffer.clear();
    }
    // convert units
    self.convert_units();
    // update pump statistics

    for pump in self.links.iter_mut() {
      if let LinkType::Pump(pump) = &mut pump.link_type {
        let head_curve_statistics = self.curves.get(&pump.head_curve).unwrap().head_curve_statistics();
        pump.head_curve_statistics = Some(head_curve_statistics);
      }
    }

    Ok(())
  }

  /// convert units
  fn convert_units(&mut self) {
    // convert the nodes to standard units (US standard) and CFS
    for node in self.nodes.iter_mut() {
      node.convert_units(&self.options.flow_units, &self.options.unit_system, false);
    }
    // convert the links to standard units (US standard) and CFS
    for link in self.links.iter_mut() {
      link.convert_units(&self.options.flow_units, &self.options.unit_system, false);
    }

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
    // read the pattern (optional, default none)
    let mut pattern: Option<Box<str>> = parts.next().map(|s| s.into());
    // if pattern is not None and is ";", set it to None
    if pattern.is_some() && pattern.clone().unwrap().contains(";") {
      pattern = None;
    }

    Node {
      id,
      elevation,
      node_type: NodeType::Junction(Junction { basedemand: demand, pattern: pattern }),
    }
  }

  /// Read a reservoir from a parts iterator
  fn read_reservoir(&mut self, line: &str) -> Node {
    let mut parts = line.split_whitespace();

    let id = parts.next().unwrap().into();
    // read the elevation
    let elevation = parts.next().unwrap().parse::<f64>().unwrap();
    // add the node to the network
    let mut pattern: Option<Box<str>> = parts.next().map(|s| s.into());
    // if pattern is not None and is ";", set it to None
    if pattern.is_some() && pattern.clone().unwrap().contains(";") {
      pattern = None;
    }

    Node {
      id,
      elevation,
      node_type: NodeType::Reservoir(Reservoir { head_pattern: pattern }),
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
      initial_status: LinkStatus::Active,
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
    // check if the pipe has a status
    let mut status = LinkStatus::Open;
    let mut check_valve = false;

    if let Some(status_str) = parts.next() {
      if status_str.to_uppercase().as_str() == "CV" {
        check_valve = true;
      }
      if status_str.to_uppercase().as_str() == "CLOSED" {
        status = LinkStatus::Closed;
      }
    }

    let start_node_index = *self.node_map.get(&start_node).unwrap();
    let end_node_index = *self.node_map.get(&end_node).unwrap();
    let headloss_formula = HeadlossFormula::HazenWilliams;

    Link {
      id,
      start_node: start_node_index,
      end_node: end_node_index,
      minor_loss: minor_loss,
      link_type: LinkType::Pipe(Pipe { diameter, length, roughness, minor_loss, check_valve, headloss_formula }),
      initial_status: status,
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
    let mut parameters = parts;

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
      link_type: LinkType::Pump(Pump { speed, head_curve, power, head_curve_statistics: None}),
      initial_status: LinkStatus::Open,
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
      self.patterns.insert(id.clone(), Pattern { multipliers });
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

    // get pattern
    let mut pattern: Option<Box<str>> = parts.next().map(|s| s.into());
    if pattern.is_some() && pattern.clone().unwrap().contains(";") {
      pattern = None;
    }

    let node_index = *self.node_map.get(&id).unwrap();
    let node = &mut self.nodes[node_index];
    match &mut node.node_type {
      NodeType::Junction(junction) => {
        junction.basedemand = demand;
        junction.pattern = pattern;
      }
      _ => panic!("Demand can only be set for junctions"),
    }
  }
  /// Read the options from a parts iterator
  fn read_options(&mut self, line: &str) {
    let mut parts = line.split_whitespace(); // TODO: Implement option reading
    
    let option = parts.next().unwrap().trim().to_uppercase();
    let value = parts.next().unwrap();

    match option.as_str() {
      "UNITS" => {
        // set the flow units
        self.options.flow_units = FlowUnits::from_str(value).unwrap();
        
        // set the unit system based on the flow units
        self.options.unit_system = match self.options.flow_units {
          FlowUnits::CFS | FlowUnits::GMP | FlowUnits::MGD | FlowUnits::IMGD | FlowUnits::AFD => UnitSystem::US,
          FlowUnits::LPS | FlowUnits::LPM | FlowUnits::MLD | FlowUnits::CMS | FlowUnits::CMH | FlowUnits::CMD => UnitSystem::SI,
        };
        // set the default pressure units based on the unit system
        self.options.pressure_units = match self.options.unit_system {
          UnitSystem::US => PressureUnits::FEET,
          UnitSystem::SI => PressureUnits::METERS,
        }

      }
      "HEADLOSS" => {
        self.options.headloss_formula = match value.to_uppercase().as_str() {
          "H-W" => HeadlossFormula::HazenWilliams,
          "D-W" => HeadlossFormula::DarcyWeisbach,
          "C-M" => HeadlossFormula::ChezyManning,
          _ => panic!("Invalid headloss formula: {}", value),
        };
        // apply the headloss formula to the pipes
        for link in self.links.iter_mut() {
          if let LinkType::Pipe(pipe) = &mut link.link_type {
            pipe.headloss_formula = self.options.headloss_formula;
          }
        }
      },
      "PRESSURE" => {
        self.options.pressure_units = PressureUnits::from_str(value).unwrap();
      }
      "TRIALS" => {
        self.options.max_trials = value.parse::<usize>().unwrap();
      },
      "ACCURACY" => {
        self.options.accuracy = value.parse::<f64>().unwrap();
      },
      "CHECKFREQ" => {
        self.options.check_frequency = value.parse::<usize>().unwrap();
      },
      "MAXCHECK" => {
        self.options.max_check = value.parse::<usize>().unwrap();
      },
      _ => ()
    }
  }
  fn read_times(&mut self, line: &str) {
    let mut parts = line.split_whitespace();
    // read the time option name
    let mut time_option = parts.next().unwrap().to_uppercase();

    if time_option == "STATISTIC" {
      return;
    }

    // if the next part is not a number, append it to the time option
    let mut duration = parts.next().unwrap();

    // if the duration is not a number and does not contain ":", append it to the time option
    if !duration.parse::<usize>().is_ok() && !duration.contains(":") {
      time_option += " ";
      time_option += &duration.to_uppercase();
      duration = parts.next().unwrap();
    }

    let mut time_units = parts.next().unwrap_or("HOURS").to_uppercase();
    // remove last "S" if it exists
    if time_units.ends_with("S") {
      time_units.pop();
    }

    let mut seconds : usize;

    // if ":" in duration, split into hours and minutes
    if duration.contains(":") {
      let mut time_parts = duration.split(":");
      let hours = time_parts.next().unwrap().parse::<usize>().unwrap();
      let minutes = time_parts.next().unwrap().parse::<usize>().unwrap();
      seconds = hours * 3600 + minutes * 60;

      if time_units == "PM" {
        // add 12 hours to the seconds
        seconds += 12 * 3600;
      }

    } else {
      let duration_value = duration.parse::<usize>().unwrap();

      seconds = match time_units.as_str() {
        "HOUR" => duration_value * 3600,
        "MINUTE" => duration_value * 60,
        "MIN" => duration_value * 60,
        "SECOND" => duration_value,
        "SEC" => duration_value,
        "DAY" => duration_value * 86400,
        "AM" => duration_value * 3600,    // AM is the same as hours
        "PM" => duration_value * 3600 + 12 * 3600, // add 12 hours to the seconds 
        _ => panic!("Invalid time units: {}", time_units),
      };
    }
    // assign the duration to the time options
    match time_option.as_str() {
      "DURATION" => self.options.time_options.duration = seconds,
      "HYDRAULIC TIMESTEP" => self.options.time_options.hydraulic_timestep = seconds,
      "PATTERN TIMESTEP" => self.options.time_options.pattern_timestep = seconds,
      "PATTERN START" => self.options.time_options.pattern_start = seconds,
      "START CLOCKTIME" => self.options.time_options.start_clocktime = seconds,
      _ => ()
    }
  }
  fn read_status(&mut self, line: &str) {
    let mut parts = line.split_whitespace();
    let id : &str = parts.next().unwrap().into();
    let status : &str = parts.next().unwrap().into();

    // get the corresponding link type
    let link = &mut self.links[*self.link_map.get(id).unwrap()];

    if let Ok(setting) = status.parse::<f64>() {
      match &mut link.link_type {
        LinkType::Valve(valve) => valve.setting = setting,
        LinkType::Pump(pump) => pump.speed = setting,
        _ => panic!("Status can only be set for valves and pumps"),
      }
    } else {
      link.initial_status = LinkStatus::from_str(status);
    }
  }
}

