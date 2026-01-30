use std::fs::File;
use std::io::{BufReader, BufRead};
use std::str::FromStr;
use std::sync::Arc;

use crate::model::network::Network;
use crate::model::node::{Node, NodeType};
use crate::model::link::{Link, LinkType, LinkStatus};
use crate::model::curve::{Curve, HeadCurve};
use crate::model::pattern::Pattern;
use crate::model::junction::Junction;
use crate::model::reservoir::Reservoir;
use crate::model::tank::Tank;
use crate::model::pipe::Pipe;
use crate::model::pump::Pump;
use crate::model::valve::{Valve, ValveType};
use crate::model::control::{Control, ControlCondition};
use crate::model::units::{FlowUnits, UnitSystem, PressureUnits, UnitConversion};
use crate::model::options::*;
use crate::error::*;

/// Error type for input parsing operations
#[derive(Debug)]
enum ReadState {
  Junctions,
  Pipes,
  Reservoirs,
  Tanks,
  Demands,
  Valves,
  Pumps,
  Curves,
  Options,
  Patterns,
  Status,
  Times,
  Rules,
  Controls,
  None,
}

impl Network {

  pub fn read_file(&mut self, file: &str) -> Result<(), InputError> {
    let file_extension = file.split('.').last()
      .ok_or_else(|| InputError::new(format!("Cannot determine file extension for: {}", file)))?;
    
    match file_extension {
      "inp" => self.read_inp(file),
      "json" => self.read_json(file),
      "mpk" | "msgpack" => self.read_msgpack(file),
      _ => Err(InputError::new(format!("Unsupported file extension: {}", file_extension))),
    }
  }

  pub fn read_msgpack(&mut self, msgpack: &str) -> Result<(), InputError> {
    let file = File::open(msgpack)
      .map_err(|e| InputError::new(format!("Failed to open file '{}': {}", msgpack, e)))?;
    let reader = BufReader::new(file);
    let network: Network = rmp_serde::from_read(reader)
      .map_err(|e| InputError::new(format!("Failed to parse msgpack file '{}': {}", msgpack, e)))?;
    *self = network;
    Ok(())
  }

  /// Read a network from a JSON file using serde_json.
  pub fn read_json(&mut self, json: &str) -> Result<(), InputError> {
    let file = File::open(json)
      .map_err(|e| InputError::new(format!("Failed to open file '{}': {}", json, e)))?;
    let reader = BufReader::new(file);
    let network: Network = serde_json::from_reader(reader)
      .map_err(|e| InputError::new(format!("Failed to parse JSON file '{}': {}", json, e)))?;
    *self = network;
    Ok(())
  }
  /// Read a network from an INP file.
  pub fn read_inp(&mut self, inp: &str) -> Result<(), InputError> {
    // set the initial state to none
    let mut state = ReadState::None;
    let mut line_number: usize = 0;

    // open the INP file
    let file = File::open(inp)
      .map_err(|e| InputError::new(format!("Failed to open file '{}': {}", inp, e)))?;

    // create a new reader
    let mut reader = BufReader::new(file);
    let mut line_buffer = String::with_capacity(512);

    // iterate over the lines in the file
    while reader.read_line(&mut line_buffer)? > 0 {
      line_number += 1;
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
          "[TANKS]" => ReadState::Tanks,
          "[CURVES]" => ReadState::Curves,
          "[PATTERNS]" => ReadState::Patterns,
          "[OPTIONS]" => ReadState::Options,
          "[STATUS]" => ReadState::Status,
          "[TIMES]" => ReadState::Times,
          "[RULES]" => ReadState::Rules,
          "[CONTROLS]" => ReadState::Controls,
          _ => ReadState::None,
        }
      }
      else {
        let result: Result<(), InputError> = match state {
          ReadState::Junctions => { 
            let junction = self.read_junction(line)?;
            self.add_node(junction).map_err(|e| InputError::new(e))
          }
          ReadState::Valves => { 
            let valve = self.read_valve(line)?;
            self.add_link(valve).map_err(|e| InputError::new(e))
          }
          ReadState::Pipes => { 
            let pipe = self.read_pipe(line)?;
            self.add_link(pipe).map_err(|e| InputError::new(e))
          }
          ReadState::Tanks => { 
            let tank = self.read_tank(line)?;
            self.add_node(tank).map_err(|e| InputError::new(e))
          }
          ReadState::Reservoirs => { 
            let reservoir = self.read_reservoir(line)?;
            self.add_node(reservoir).map_err(|e| InputError::new(e))
          }
          ReadState::Pumps => { 
            let pump = self.read_pump(line)?;
            self.add_link(pump).map_err(|e| InputError::new(e))
          }
          ReadState::Curves => { 
            self.read_curve(line)
          }
          ReadState::Demands => { 
            self.read_demand(line)
          }
          ReadState::Patterns => { 
            self.read_pattern(line)
          }
          ReadState::None => {
            // skip unknown state
            Ok(())
          },
          ReadState::Options => {
            self.read_options(line)
          }
          ReadState::Times => {
            self.read_times(line)
          }
          ReadState::Status => {
            self.read_status(line)
          }
          ReadState::Rules => {
            Err(InputError::new("Rules are not supported yet"))
          }
          ReadState::Controls => {
            self.read_control(line)
          }
        };
        
        // Add line number context to any error
        result.map_err(|e| e.with_line(line_number).with_context(line.to_string()))?;
      }
      // clear the line buffer
      line_buffer.clear();
    }
    // convert units
    self.convert_units();
    self.update_links()?;
    Ok(())
  }

  fn update_links(&mut self) -> Result<(), InputError> {
    // update link specific statistics (TODO: move to separate method)
    for (link_index, link) in self.links.iter_mut().enumerate() {
      if let LinkType::Pump(pump) = &mut link.link_type {
        // get the head curve
        if let Some(head_curve_id) = &pump.head_curve_id {
          let curve = self.curves.get(head_curve_id)
            .ok_or_else(|| InputError::new(format!("Head curve '{}' not found for pump", head_curve_id)))?;
          // assign the head curve to the pump and convert units to standard units (US standard) and CFS
          pump.head_curve = Some(HeadCurve::new(curve, &self.options.flow_units, &self.options.unit_system));
        }
      }
      if let LinkType::Valve(valve) = &mut link.link_type {
        // if the valve is a PSV, add the elevation of the start node to the setting
        if valve.valve_type == ValveType::PSV {
          valve.setting += self.nodes[link.start_node].elevation;
          self.contains_pressure_control_valve = true;
        }
        // if the valve is a PRV, add the elevation of the end node from the setting
        if valve.valve_type == ValveType::PRV {
          valve.setting += self.nodes[link.end_node].elevation;
          self.contains_pressure_control_valve = true;
        }
        // assign the valve curve to the valve
        if let Some(curve_id) = &valve.curve_id {
          let curve = self.curves.get(curve_id)
            .ok_or_else(|| InputError::new(format!("Curve '{}' not found for valve", curve_id)))?;
          valve.valve_curve = Some(Arc::new(curve.clone()));
        }
      }
      if let LinkType::Pipe(pipe) = &mut link.link_type {
        // convert the minor loss to a coefficient
        pipe.minor_loss = 0.02517 * pipe.minor_loss / pipe.diameter.powi(4);
      }
      // check if the link is connected to a tank
      if let NodeType::Tank(tank) = &mut self.nodes[link.end_node].node_type {
        tank.links_to.push(link_index);
      }
      if let NodeType::Tank(tank) = &mut self.nodes[link.start_node].node_type {
        tank.links_from.push(link_index);
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
  fn read_junction(&mut self, line: &str) -> Result<Node, InputError> {
    let mut parts = parse_line(line);
    let id = parts.next().ok_or_missing("junction id")?.into();
    let elevation = parts.next().unwrap_or("0").parse_field::<f64>("elevation")?;
    let demand = parts.next().and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
    let pattern: Option<Box<str>> = parts.next().map(|s| s.into());

    Ok(Node {
      id,
      elevation,
      node_type: NodeType::Junction(Junction { basedemand: demand, pattern }),
    })
  }

  /// Read a reservoir from a parts iterator
  fn read_reservoir(&mut self, line: &str) -> Result<Node, InputError> {
    let mut parts = parse_line(line);
    let id = parts.next().ok_or_missing("reservoir id")?.into();
    let elevation = parts.next().ok_or_missing("elevation")?.parse_field::<f64>("elevation")?;
    let pattern: Option<Box<str>> = parts.next().map(|s| s.into());

    Ok(Node {
      id,
      elevation,
      node_type: NodeType::Reservoir(Reservoir { head_pattern: pattern }),
    })
  }

  fn read_tank(&mut self, line: &str) -> Result<Node, InputError> {
    let mut parts = parse_line(line);
    let id = parts.next().ok_or_missing("tank id")?.into();
    let elevation = parts.next().ok_or_missing("elevation")?.parse_field::<f64>("elevation")?;
    let initial_level = parts.next().ok_or_missing("initial level")?.parse_field::<f64>("initial level")?;
    let min_level = parts.next().ok_or_missing("min level")?.parse_field::<f64>("min level")?;
    let max_level = parts.next().ok_or_missing("max level")?.parse_field::<f64>("max level")?;
    let diameter = parts.next().ok_or_missing("diameter")?.parse_field::<f64>("diameter")?;
    let min_volume = parts.next()
      .map(|v| v.parse_field::<f64>("min volume"))
      .transpose()?
      .unwrap_or(0.0);

    // read the volume curve ID (optional, default none. if asterisk, no volume curve)
    let volume_curve_id: Option<Box<str>> = parts.next()
      .filter(|&v| v != "*")
      .map(|s| s.to_string().into_boxed_str());
    
    // read the overflow (optional, default false)
    let overflow = parts.next()
      .map(|s| s.to_uppercase() == "YES")
      .unwrap_or(false);

    Ok(Node {
      id,
      elevation,
      node_type: NodeType::Tank(Tank { elevation, initial_level, min_level, max_level, diameter, min_volume, volume_curve_id, overflow, volume_curve: None, links_to: Vec::new(), links_from: Vec::new() }),
    })
  }

  /// Read a valve from a parts iterator
  fn read_valve(&mut self, line: &str) -> Result<Link, InputError> {
    let mut parts = parse_line(line);

    let id = parts.next().ok_or_missing("valve id")?.into();
    // read the start node
    let start_node: Box<str> = parts.next().ok_or_missing("start node")?.into();
    // read the end node
    let end_node: Box<str> = parts.next().ok_or_missing("end node")?.into();
    // read the diameter
    let diameter = parts.next().ok_or_missing("diameter")?.parse_field::<f64>("diameter")?;
    // create the link
    let valve_type_str = parts.next().ok_or_missing("valve type")?;
    let valve_type = match valve_type_str.to_uppercase().as_str() {
      "GPV" => ValveType::GPV,
      "TCV" => ValveType::TCV,
      "FCV" => ValveType::FCV,
      "PRV" => ValveType::PRV,
      "PSV" => ValveType::PSV,
      "PBV" => ValveType::PBV,
      "PCV" => ValveType::PCV,
      _ => return Err(InputError::new(format!("Invalid valve type: {}", valve_type_str))),
    };

    // read the setting or curve ID
    let (curve_id, setting) = if valve_type == ValveType::GPV {
      let curve_id = parts.next().map(|s| s.to_string().into_boxed_str());
      (curve_id, 0.0)
    } else {
      let setting = parts.next()
        .map(|s| s.parse::<f64>().unwrap_or(0.0))
        .unwrap_or(0.0);
      (None, setting)
    };

    // read the minor loss
    let minor_loss = parts.next()
      .map(|s| s.parse::<f64>().unwrap_or(0.0))
      .unwrap_or(0.0);

    // read the PCV curve ID (optional, default none)
    let fcv_curve_id: Option<Box<str>> = parts.next().map(|s| s.to_string().into_boxed_str());
    // overwrite curve ID with FCV curve ID if present
    let curve_id = fcv_curve_id.or(curve_id);

    let start_node_index = *self.node_map.get(&start_node)
      .ok_or_else(|| InputError::new(format!("Start node '{}' not found for valve '{}'", start_node, id)))?;
    let end_node_index = *self.node_map.get(&end_node)
      .ok_or_else(|| InputError::new(format!("End node '{}' not found for valve '{}'", end_node, id)))?;

    Ok(Link {
      id,
      start_node: start_node_index,
      end_node: end_node_index,
      start_node_id: start_node,
      end_node_id: end_node,
      link_type: LinkType::Valve(Valve { diameter, setting, curve_id, valve_type, minor_loss, valve_curve: None }),
      initial_status: LinkStatus::Active,
    })
  }

  /// Read a pipe from a parts iterator
  fn read_pipe(&self, line: &str) -> Result<Link, InputError> {
    let mut parts = parse_line(line);

    let id = parts.next().ok_or_missing("pipe id")?.into();
    // read the start node
    let start_node: Box<str> = parts.next().ok_or_missing("start node")?.into();
    // read the end node
    let end_node: Box<str> = parts.next().ok_or_missing("end node")?.into();
    // read the length
    let length = parts.next().ok_or_missing("length")?.parse_field::<f64>("length")?;
    // read the diameter
    let diameter = parts.next().ok_or_missing("diameter")?.parse_field::<f64>("diameter")?;
    // read the roughness
    let roughness = parts.next().ok_or_missing("roughness")?.parse_field::<f64>("roughness")?;

    // create the link
    let minor_loss = parts.next()
      .map(|s| s.parse::<f64>().unwrap_or(0.0))
      .unwrap_or(0.0);
    
    // check if the pipe has a status
    let mut status = LinkStatus::Open;
    let mut check_valve = false;

    if let Some(status_str) = parts.next() {
      match status_str.to_uppercase().as_str() {
        "CV" => check_valve = true,
        "CLOSED" => status = LinkStatus::Closed,
        _ => {}
      }
    }

    let start_node_index = *self.node_map.get(&start_node)
      .ok_or_else(|| InputError::new(format!("Start node '{}' not found for pipe '{}'", start_node, id)))?;
    let end_node_index = *self.node_map.get(&end_node)
      .ok_or_else(|| InputError::new(format!("End node '{}' not found for pipe '{}'", end_node, id)))?;
    let headloss_formula = HeadlossFormula::HazenWilliams;

    Ok(Link {
      id,
      start_node: start_node_index,
      end_node: end_node_index,
      start_node_id: start_node,
      end_node_id: end_node,
      link_type: LinkType::Pipe(Pipe { diameter, length, roughness, minor_loss, check_valve, headloss_formula }),
      initial_status: status,
    })
  }

  /// Read a pump from a parts iterator
  fn read_pump(&self, line: &str) -> Result<Link, InputError> {
    let mut parts = parse_line(line);

    let id: Box<str> = parts.next().ok_or_missing("pump id")?.into();
    // read the start node
    let start_node: Box<str> = parts.next().ok_or_missing("start node")?.into();
    // read the end node
    let end_node: Box<str> = parts.next().ok_or_missing("end node")?.into();

    // read the parameters
    let mut parameters = parts;

    let mut speed = 1.0;
    let mut head_curve_id = None;
    let mut power = 0.0;

    // create the link
    let start_node_index = *self.node_map.get(&start_node)
      .ok_or_else(|| InputError::new(format!("Start node '{}' not found for pump '{}'", start_node, id)))?;
    let end_node_index = *self.node_map.get(&end_node)
      .ok_or_else(|| InputError::new(format!("End node '{}' not found for pump '{}'", end_node, id)))?;

    while let Some(parameter) = parameters.next() {
      if parameter == ";" {
        continue;
      }
      if let Some(value) = parameters.next() {
        match parameter {
          "SPEED" => speed = value.parse::<f64>().unwrap_or(1.0),
          "HEAD" => head_curve_id = Some(value.trim().into()),
          "POWER" => power = value.parse::<f64>().unwrap_or(0.0),
          _ => continue,
        }
      }
    }

    Ok(Link {
      id,
      start_node_id: start_node,
      end_node_id: end_node,
      start_node: start_node_index,
      end_node: end_node_index,
      link_type: LinkType::Pump(Pump { speed, head_curve_id, power, head_curve: None }),
      initial_status: LinkStatus::Open,
    })
  }
  /// Read a curve from a parts iterator
  /// Appends a point to the curve if it already exists, otherwise creates a new curve
  fn read_curve(&mut self, line: &str) -> Result<(), InputError> {
    let mut parts = parse_line(line);

    let id: Box<str> = parts.next().ok_or_missing("curve id")?.into();
    let x = parts.next().ok_or_missing("x value")?.parse_field::<f64>("x value")?;
    let y = parts.next().ok_or_missing("y value")?.parse_field::<f64>("y value")?;
    
    // create curve if it does not exist
    if !self.curves.contains_key(&id) {
      self.curves.insert(id.clone(), Curve { id, x: vec![x], y: vec![y]});
    }
    // otherwise, add the point to the curve
    else {
      let curve = self.curves.get_mut(&id)
        .ok_or_else(|| InputError::new(format!("Curve '{}' not found", id)))?;
      // check if the x value is in ascending order
      if let Some(last_x) = curve.x.last() {
        if x < *last_x {
          return Err(InputError::new(format!("X values must be in ascending order for curve '{}': {} < {}", id, x, last_x)));
        }
      }
      // add the point to the curve
      curve.x.push(x);
      curve.y.push(y);
    }
    Ok(())
  }

  /// Read a pattern from a parts iterator
  fn read_pattern(&mut self, line: &str) -> Result<(), InputError> {
    let mut parts = parse_line(line);
    let id: Box<str> = parts.next().ok_or_missing("pattern id")?.into();
    
    let multipliers: Result<Vec<f64>, InputError> = parts
      .map(|s| s.parse_field::<f64>("multiplier"))
      .collect();
    let multipliers = multipliers?;
    
    if !self.patterns.contains_key(&id) {
      self.patterns.insert(id.clone(), Pattern { multipliers });
    } else {
      let pattern = self.patterns.get_mut(&id)
        .ok_or_else(|| InputError::new(format!("Pattern '{}' not found", id)))?;
      pattern.multipliers.extend(multipliers);
    }
    Ok(())
  }
  
  /// Read a demand from a parts iterator and set the basedemand for the junction
  fn read_demand(&mut self, line: &str) -> Result<(), InputError> {
    let mut parts = parse_line(line);
    let id: Box<str> = parts.next().ok_or_missing("node id")?.into();
    let demand = parts.next().ok_or_missing("demand")?.parse_field::<f64>("demand")?;
    let pattern: Option<Box<str>> = parts.next().map(|s| s.into());

    let node_index = *self.node_map.get(&id)
      .ok_or_else(|| InputError::new(format!("Node '{}' not found for demand", id)))?;
    let node = &mut self.nodes[node_index];
    
    match &mut node.node_type {
      NodeType::Junction(junction) => {
        junction.basedemand = demand;
        junction.pattern = pattern;
        Ok(())
      }
      _ => Err(InputError::new(format!("Demand can only be set for junctions, but '{}' is not a junction", id))),
    }
  }
  /// Read the options from a parts iterator
  fn read_options(&mut self, line: &str) -> Result<(), InputError> {
    let mut parts = parse_line(line);
    
    let option = parts.next().ok_or_missing("option name")?.trim().to_uppercase();
    let value = parts.next().ok_or_missing("option value")?;

    match option.as_str() {
      "UNITS" => {
        // set the flow units
        self.options.flow_units = FlowUnits::from_str(value)
          .map_err(|_| InputError::new(format!("Invalid flow units: {}", value)))?;
        
        // set the unit system based on the flow units
        self.options.unit_system = match self.options.flow_units {
          FlowUnits::CFS | FlowUnits::GPM | FlowUnits::MGD | FlowUnits::IMGD | FlowUnits::AFD => UnitSystem::US,
          FlowUnits::LPS | FlowUnits::LPM | FlowUnits::MLD | FlowUnits::CMS | FlowUnits::CMH | FlowUnits::CMD => UnitSystem::SI,
        };
        // set the default pressure units based on the unit system
        self.options.pressure_units = match self.options.unit_system {
          UnitSystem::US => PressureUnits::FEET,
          UnitSystem::SI => PressureUnits::METERS,
        }
      }
      "DEMAND" => {
        let next_part = value.trim().to_uppercase();
        if next_part == "MULTIPLIER" {
          self.options.demand_multiplier = parts.next()
            .ok_or_missing("demand multiplier value")?
            .parse_field::<f64>("demand multiplier")?;
        } else if let Some(model) = parts.next() {
          if model.trim().to_uppercase() == "PDA" {
            return Err(InputError::new("PDA demand model is not supported yet"));
          }
        }
      },
      "HEADLOSS" => {
        self.options.headloss_formula = match value.to_uppercase().as_str() {
          "H-W" => HeadlossFormula::HazenWilliams,
          "D-W" => HeadlossFormula::DarcyWeisbach,
          "C-M" => HeadlossFormula::ChezyManning,
          _ => return Err(InputError::new(format!("Invalid headloss formula: {}", value))),
        };
        // apply the headloss formula to the pipes
        for link in self.links.iter_mut() {
          if let LinkType::Pipe(pipe) = &mut link.link_type {
            pipe.headloss_formula = self.options.headloss_formula;
          }
        }
      },
      "PRESSURE" => {
        // Handle "Pressure Exponent" as a separate option (skip if not a valid pressure unit)
        if let Ok(units) = PressureUnits::from_str(value) {
          self.options.pressure_units = units;
        }
        // Otherwise it's likely "Pressure Exponent" or similar - ignore
      }
      "TRIALS" => {
        self.options.max_trials = value.parse_field::<usize>("trials")?;
      },
      "ACCURACY" => {
        self.options.accuracy = value.parse_field::<f64>("accuracy")?;
      },
      "CHECKFREQ" => {
        self.options.check_frequency = value.parse_field::<usize>("check frequency")?;
      },
      "MAXCHECK" => {
        self.options.max_check = value.parse_field::<usize>("max check")?;
      },
      "FLOWCHANGE" => {
        self.options.max_flow_change = Some(value.parse_field::<f64>("flow change")?);
      },
      "PATTERN" => {
        let pattern: Box<str> = value.into();
        if self.patterns.contains_key(&pattern) {
          // set the default pattern to all junctions without a pattern
          for node in self.nodes.iter_mut() {
            if let NodeType::Junction(junction) = &mut node.node_type {
              if junction.pattern.is_none() {
                junction.pattern = Some(value.into());
              }
            }
          }
        }
      }
      _ => ()
    }
    Ok(())
  }

  fn read_times(&mut self, line: &str) -> Result<(), InputError> {
    let mut parts = parse_line(line);
    // read the time option name
    let mut time_option = parts.next().ok_or_missing("time option")?.to_uppercase();

    if time_option == "STATISTIC" {
      return Ok(());
    }

    // if the next part is not a number, append it to the time option
    let mut duration = parts.next().ok_or_missing("duration")?;

    // if the duration is not a number and does not contain ":", append it to the time option
    if duration.parse::<usize>().is_err() && !duration.contains(":") {
      time_option += " ";
      time_option += &duration.to_uppercase();
      duration = parts.next().ok_or_missing("duration value")?;
    }

    let time_units = parts.next();
    let seconds = parse_time(duration, time_units)?;

    // assign the duration to the time options
    match time_option.as_str() {
      "DURATION" => self.options.time_options.duration = seconds,
      "HYDRAULIC TIMESTEP" => self.options.time_options.hydraulic_timestep = seconds,
      "PATTERN TIMESTEP" => self.options.time_options.pattern_timestep = seconds,
      "REPORT TIMESTEP" => self.options.time_options.report_timestep = seconds,
      "PATTERN START" => self.options.time_options.pattern_start = seconds,
      "START CLOCKTIME" => self.options.time_options.start_clocktime = seconds,
      _ => ()
    }
    Ok(())
  }

  fn read_control(&mut self, line: &str) -> Result<(), InputError> {
    let mut parts = parse_line(line);

    // Skip "LINK" keyword
    let first = parts.next().ok_or_missing("control keyword")?.to_uppercase();
    if first != "LINK" {
      return Ok(());
    }

    // Read link id
    let link_id: Box<str> = parts.next().ok_or_missing("link id")?.into();

    // Read status or setting
    let status_or_setting = parts.next().ok_or_missing("status or setting")?;
    let (status, setting) = if let Ok(value) = status_or_setting.parse::<f64>() {
      (None, Some(value))
    } else {
      (Some(LinkStatus::from_str(status_or_setting, false)), None)
    };

    // Skip "IF" keyword
    parts.next();

    // Read condition type
    let condition_type = parts.next().ok_or_missing("condition type")?.to_uppercase();

    let condition = match condition_type.as_str() {
      "NODE" => {
        let node_id: Box<str> = parts.next().ok_or_missing("node id")?.into();
        let above_below = parts.next().ok_or_missing("ABOVE or BELOW")?;
        let above = above_below.to_uppercase() == "ABOVE";
        let value = parts.next().ok_or_missing("pressure value")?.parse_field::<f64>("pressure value")?;
        ControlCondition::Pressure { node_id, above, below: value }
      }
      "TIME" => {
        let time_str = parts.next().ok_or_missing("time value")?;
        let seconds = parse_time(time_str, parts.next())?;
        ControlCondition::Time { seconds }
      }
      "CLOCKTIME" => {
        let time_str = parts.next().ok_or_missing("clock time value")?;
        let seconds = parse_time(time_str, parts.next())?;
        ControlCondition::ClockTime { seconds }
      }
      _ => return Err(InputError::new(format!("Invalid control condition type: {}", condition_type))),
    };

    self.controls.push(Control {
      condition,
      link_id,
      setting,
      status,
    });
    Ok(())
  }

  fn read_status(&mut self, line: &str) -> Result<(), InputError> {
    let mut parts = parse_line(line);
    let id: &str = parts.next().ok_or_missing("link id")?;
    let status: &str = parts.next().ok_or_missing("status")?;

    // get the corresponding link type
    let link_index = *self.link_map.get(id)
      .ok_or_else(|| InputError::new(format!("Link '{}' not found for status", id)))?;
    let link = &mut self.links[link_index];

    if let Ok(setting) = status.parse::<f64>() {
      match &mut link.link_type {
        LinkType::Valve(valve) => valve.setting = setting,
        LinkType::Pump(pump) => pump.speed = setting,
        _ => return Err(InputError::new(format!("Status/setting can only be set for valves and pumps, not link '{}'", id))),
      }
    } else {
      // if the link is a valve, set the status as fixed open or fixed closed
      let is_valve = matches!(link.link_type, LinkType::Valve(_));
      link.initial_status = LinkStatus::from_str(status, is_valve);
    }
    Ok(())
  }
}

/// Strip comments (after ';') from a line and split into whitespace-separated parts
fn parse_line(line: &str) -> std::str::SplitWhitespace<'_> {
  line.split(';').next().unwrap_or("").trim().split_whitespace()
}

/// Parse a time string with optional time unit suffix.
/// Supports formats:
/// - "HH:MM" with optional AM/PM suffix
/// - Numeric value with unit (HOUR, MINUTE, MIN, SECOND, SEC, DAY, AM, PM)
/// - If no unit is provided, defaults to HOURS
fn parse_time(time_str: &str, unit_or_suffix: Option<&str>) -> Result<usize, InputError> {
  let seconds: usize;

  // If time contains ":", parse as HH:MM format
  if time_str.contains(":") {
    let mut time_parts = time_str.split(":");
    let hours = time_parts.next()
      .ok_or_else(|| InputError::new("Missing hours in time"))?
      .parse::<usize>()
      .map_err(|_| InputError::new(format!("Invalid hours value in time: {}", time_str)))?;
    let minutes = time_parts.next()
      .ok_or_else(|| InputError::new("Missing minutes in time"))?
      .parse::<usize>()
      .map_err(|_| InputError::new(format!("Invalid minutes value in time: {}", time_str)))?;
    seconds = hours * 3600 + minutes * 60 + if let Some(suffix) = unit_or_suffix {
      if suffix.to_uppercase() == "PM" { 12 * 3600 } else { 0 }
    } else { 0 };
  } else {
    // Parse as numeric value with optional time unit
    let value = time_str.parse::<usize>()
      .map_err(|_| InputError::new(format!("Invalid time value: {}", time_str)))?;
    let mut unit = unit_or_suffix.unwrap_or("HOURS").to_uppercase();
    
    // Remove trailing "S" for singular form
    if unit.ends_with("S") && unit != "HOURS" {
      unit.pop();
    }

    seconds = match unit.as_str() {
      "HOUR" | "HOURS" => value * 3600,
      "MINUTE" | "MIN" => value * 60,
      "SECOND" | "SEC" => value,
      "DAY" => value * 86400,
      "AM" => value * 3600,
      "PM" => value * 3600 + 12 * 3600,
      _ => return Err(InputError::new(format!("Invalid time unit: {}", unit))),
    };
  }

  Ok(seconds)
}

#[cfg(test)]
mod tests {
  use super::*;

  /// Helper to create a network with default options for testing.
  /// If `with_nodes` is true, adds two default nodes: "N1" (junction) and "N2" (junction).
  fn test_network(with_nodes: bool) -> Network {
    let mut network = Network::default();
    if with_nodes {
      network.add_node(Node {
        id: "N1".into(),
        elevation: 0.0,
        node_type: NodeType::Junction(Junction { basedemand: 0.0, pattern: None }),
      }).unwrap();
      network.add_node(Node {
        id: "N2".into(),
        elevation: 0.0,
        node_type: NodeType::Junction(Junction { basedemand: 0.0, pattern: None }),
      }).unwrap();
    }
    network
  }

  // ==================== Junction Tests ====================

  #[test]
  fn test_read_junction_basic() {
    let mut network = test_network(false);
    let node = network.read_junction("J1  100.5  25.0").unwrap();
    
    assert_eq!(&*node.id, "J1");
    assert_eq!(node.elevation, 100.5);
    
    let NodeType::Junction(junction) = &node.node_type else {
      panic!("Expected Junction node type");
    };

    assert_eq!(junction.basedemand, 25.0);
    assert!(junction.pattern.is_none());
  }

  #[test]
  fn test_read_junction_with_pattern() {
    let mut network = test_network(false);
    let node = network.read_junction("J2  50.0  100.0  PAT1").unwrap();
    
    assert_eq!(&*node.id, "J2");
    assert_eq!(node.elevation, 50.0);
    
    let NodeType::Junction(junction) = &node.node_type else {
      panic!("Expected Junction node type");
    };

    assert_eq!(junction.basedemand, 100.0);
    assert_eq!(junction.pattern.as_deref(), Some("PAT1"));
  }

  #[test]
  fn test_read_junction_with_comment() {
    let mut network = test_network(false);
    // Pattern field contains semicolon (comment marker) - should be ignored
    let node = network.read_junction("J3  75.0  50.0  ;comment").unwrap();
    
    let NodeType::Junction(junction) = &node.node_type else {
      panic!("Expected Junction node type");
    };

    assert!(junction.pattern.is_none());
  }

  #[test]
  fn test_read_junction_missing_demand() {
    let mut network = test_network(false);
    // Only ID and elevation provided - demand should default to 0.0
    let node = network.read_junction("J4  200.0").unwrap();

    let NodeType::Junction(junction) = &node.node_type else {
      panic!("Expected Junction node type");
    };

    assert_eq!(junction.basedemand, 0.0);
    
  }
  // ==================== Valve Tests ====================

  #[test]
  fn test_read_valve_basic() {
    let mut network = test_network(true);
    let valve = network.read_valve("V1  N1  N2 12.0 PRV 50.0  100.0").unwrap();
    
    assert_eq!(&*valve.id, "V1");
    let LinkType::Valve(valve) = &valve.link_type else {
      panic!("Expected Valve link type");
    };

    assert_eq!(valve.valve_type, ValveType::PRV);
    assert_eq!(valve.diameter, 12.0);
    assert_eq!(valve.setting, 50.0);
    assert_eq!(valve.minor_loss, 100.0);
  }

  #[test]
  fn test_read_valve_gpv() {
    let mut network = test_network(true);
    let valve = network.read_valve("V2  N1  N2 12.0 GPV GPV_CURVE  100.0").unwrap();
    
    assert_eq!(&*valve.id, "V2");
    let LinkType::Valve(valve) = &valve.link_type else {
      panic!("Expected Valve link type");
    };
    assert_eq!(valve.valve_type, ValveType::GPV);
    assert_eq!(valve.curve_id.as_deref(), Some("GPV_CURVE"));
  }

  #[test]
  fn test_read_valve_fcv() {
    let mut network = test_network(true);
    let valve = network.read_valve("V3  N1  N2 12.0 FCV 50.0  100.0 FCV_CURVE").unwrap();
    
    assert_eq!(&*valve.id, "V3");
    let LinkType::Valve(valve) = &valve.link_type else {
      panic!("Expected Valve link type");
    };
    assert_eq!(valve.valve_type, ValveType::FCV);
    assert_eq!(valve.diameter, 12.0);
    assert_eq!(valve.setting, 50.0);
    assert_eq!(valve.minor_loss, 100.0);
    assert_eq!(valve.curve_id.as_deref(), Some("FCV_CURVE"));
  }

  // ==================== Reservoir Tests ====================

  #[test]
  fn test_read_reservoir_basic() {
    let mut network = test_network(false);
    let node = network.read_reservoir("RES1  150.0").unwrap();
    
    assert_eq!(&*node.id, "RES1");
    assert_eq!(node.elevation, 150.0);
    
    let NodeType::Reservoir(reservoir) = &node.node_type else {
      panic!("Expected Reservoir node type");
    };

    assert!(reservoir.head_pattern.is_none());
  }

  #[test]
  fn test_read_reservoir_with_pattern() {
    let mut network = test_network(false);
    let node = network.read_reservoir("RES2  200.0  HEADPAT; comment").unwrap();
    
    let NodeType::Reservoir(reservoir) = &node.node_type else {
      panic!("Expected Reservoir node type");
    };

    assert_eq!(reservoir.head_pattern.as_deref(), Some("HEADPAT"));
  }

  // ==================== Pipe Tests ====================

  #[test]
  fn test_read_pipe_basic() {
    let network = test_network(true);

    let link = network.read_pipe("P1  N2  N1  1000.0  12.0  100.0  0.0").unwrap();
    
    assert_eq!(&*link.id, "P1");
    assert_eq!(link.start_node, 1);
    assert_eq!(link.end_node, 0);
    assert_eq!(link.initial_status, LinkStatus::Open);
    
    let LinkType::Pipe(pipe) = &link.link_type else {
      panic!("Expected Pipe link type");
    };

    assert_eq!(pipe.length, 1000.0);
    assert_eq!(pipe.diameter, 12.0);
    assert_eq!(pipe.roughness, 100.0);
    assert!(!pipe.check_valve);
  }

  #[test]
  fn test_read_pipe_with_check_valve() {
    let network = test_network(true);

    let link = network.read_pipe("P2  N1  N2  500.0  8.0  120.0  0.0  CV").unwrap();
    
    let LinkType::Pipe(pipe) = &link.link_type else {
      panic!("Expected Pipe link type");
    };

    assert!(pipe.check_valve);
  }

  #[test]
  fn test_read_pipe_closed() {
    let network = test_network(true);

    let link = network.read_pipe("P3  N1  N2  200.0  6.0  110.0  0.0  CLOSED").unwrap();
    
    assert_eq!(link.initial_status, LinkStatus::Closed);
  }
  // ==================== Tank Tests ====================

  #[test]
  fn test_read_tank_basic() {
    let mut network = test_network(false);
    let node = network.read_tank("T1  100  15  5  25  120  0  *  Yes").unwrap();

    let NodeType::Tank(tank) = &node.node_type else {
      panic!("Expected Tank node type");
    };

    assert_eq!(node.elevation, 100.0);
    assert_eq!(tank.initial_level, 15.0);
    assert_eq!(tank.min_level, 5.0);
    assert_eq!(tank.max_level, 25.0);
    assert_eq!(tank.diameter, 120.0);
    assert_eq!(tank.min_volume, 0.0);
    assert!(tank.volume_curve_id.is_none());
    assert!(tank.overflow);
  }
  #[test]
  fn test_read_tank_with_volume_curve() {
    let mut network = test_network(false);
    let node = network.read_tank("T2  100  15  5  25  120  0  VOLCURVE").unwrap();
    
    let NodeType::Tank(tank) = &node.node_type else {
      panic!("Expected Tank node type");
    };
    assert_eq!(tank.volume_curve_id.as_deref(), Some("VOLCURVE"));
    assert!(!tank.overflow);
  }

  // ==================== Pump Tests ====================

  #[test]
  fn test_read_pump_with_head_curve() {
    let network = test_network(true);

    let link = network.read_pump("PUMP1  N1  N2  HEAD CURVE1").unwrap();
    
    assert_eq!(&*link.id, "PUMP1");
    assert_eq!(link.initial_status, LinkStatus::Open);
    
    let LinkType::Pump(pump) = &link.link_type else {
      panic!("Expected Pump link type");
    };

    assert_eq!(pump.head_curve_id, Some("CURVE1".into()));
    assert_eq!(pump.speed, 1.0); // default speed
  }

  #[test]
  fn test_read_pump_with_speed() {
    let network = test_network(true);

    let link = network.read_pump("PUMP2  N1  N2  HEAD C1  SPEED 1.5").unwrap();
    
    let LinkType::Pump(pump) = &link.link_type else {
      panic!("Expected Pump link type");
    };

    assert_eq!(pump.speed, 1.5);
    assert_eq!(pump.head_curve_id, Some("C1".into()));
  }

  #[test]
  fn test_read_pump_with_power() {
    let network = test_network(true);
    let link = network.read_pump("PUMP3  N1  N2  HEAD C1  POWER 100.0").unwrap();
    
    let LinkType::Pump(pump) = &link.link_type else {
      panic!("Expected Pump link type");
    };
    assert_eq!(pump.power, 100.0);
  }

  // ==================== Curve Tests ====================

  #[test]
  fn test_read_curve_single_point() {
    let mut network = test_network(false);
    network.read_curve("CURVE1  100.0  50.0").unwrap();
    
    let curve = network.curves.get("CURVE1").unwrap();
    assert_eq!(curve.x, vec![100.0]);
    assert_eq!(curve.y, vec![50.0]);
  }

  #[test]
  fn test_read_curve_multiple_points() {
    let mut network = test_network(false);
    network.read_curve("CURVE2  0.0  100.0").unwrap();
    network.read_curve("CURVE2  50.0  75.0").unwrap();
    network.read_curve("CURVE2  100.0  25.0").unwrap();
    
    let curve = network.curves.get("CURVE2").unwrap();
    assert_eq!(curve.x, vec![0.0, 50.0, 100.0]);
    assert_eq!(curve.y, vec![100.0, 75.0, 25.0]);
  }

  // ==================== Pattern Tests ====================

  #[test]
  fn test_read_pattern_single_line() {
    let mut network = test_network(false);
    network.read_pattern("PAT1  1.0  1.2  0.8  1.1").unwrap();
    
    let pattern = network.patterns.get("PAT1").unwrap();
    assert_eq!(pattern.multipliers, vec![1.0, 1.2, 0.8, 1.1]);
  }

  #[test]
  fn test_read_pattern_multiple_lines() {
    let mut network = test_network(false);
    network.read_pattern("PAT2  1.0  1.5").unwrap();
    network.read_pattern("PAT2  2.0  0.5").unwrap();
    
    let pattern = network.patterns.get("PAT2").unwrap();
    assert_eq!(pattern.multipliers, vec![1.0, 1.5, 2.0, 0.5]);
  }

  // ==================== Options Tests ====================

  #[test]
  fn test_read_options_units_lps() {
    let mut network = test_network(false);
    network.read_options("UNITS  LPS").unwrap();
    
    assert_eq!(network.options.flow_units, FlowUnits::LPS);
    assert_eq!(network.options.unit_system, UnitSystem::SI);
    assert_eq!(network.options.pressure_units, PressureUnits::METERS);
  }

  #[test]
  fn test_read_options_units_cfs() {
    let mut network = test_network(false);
    network.read_options("UNITS  CFS").unwrap();
    
    assert_eq!(network.options.flow_units, FlowUnits::CFS);
    assert_eq!(network.options.unit_system, UnitSystem::US);
    assert_eq!(network.options.pressure_units, PressureUnits::FEET);
  }

  #[test]
  fn test_read_options_headloss() {
    let mut network = test_network(false);
    network.read_options("HEADLOSS  D-W").unwrap();
    
    assert_eq!(network.options.headloss_formula, HeadlossFormula::DarcyWeisbach);
  }

  #[test]
  fn test_read_options_trials() {
    let mut network = test_network(false);
    network.read_options("TRIALS  100").unwrap();
    
    assert_eq!(network.options.max_trials, 100);
  }

  #[test]
  fn test_read_options_accuracy() {
    let mut network = test_network(false);
    network.read_options("ACCURACY  0.0001").unwrap();
    
    assert!((network.options.accuracy - 0.0001).abs() < 1e-10);
  }

  // ==================== Times Tests ====================

  #[test]
  fn test_read_times_duration_hours() {
    let mut network = test_network(false);
    network.read_times("DURATION  24  HOURS").unwrap();
    
    assert_eq!(network.options.time_options.duration, 24 * 3600);
  }

  #[test]
  fn test_read_times_duration_colon_format() {
    let mut network = test_network(false);
    network.read_times("DURATION  12:30").unwrap();
    
    assert_eq!(network.options.time_options.duration, 12 * 3600 + 30 * 60);
  }

  #[test]
  fn test_read_times_hydraulic_timestep() {
    let mut network = test_network(false);
    network.read_times("HYDRAULIC TIMESTEP  1:00").unwrap();
    
    assert_eq!(network.options.time_options.hydraulic_timestep, 3600);
  }

  #[test]
  fn test_read_times_pattern_timestep() {
    let mut network = test_network(false);
    network.read_times("PATTERN TIMESTEP  2  HOURS").unwrap();
    
    assert_eq!(network.options.time_options.pattern_timestep, 2 * 3600);
  }

  #[test]
  fn test_read_times_start_clocktime_am() {
    let mut network = test_network(false);
    network.read_times("START CLOCKTIME  6  AM").unwrap();
    
    assert_eq!(network.options.time_options.start_clocktime, 6 * 3600);
  }

  #[test]
  fn test_read_times_start_clocktime_pm() {
    let mut network = test_network(false);
    network.read_times("START CLOCKTIME  6  PM").unwrap();
    
    assert_eq!(network.options.time_options.start_clocktime, 18 * 3600);
  }

  #[test]
  fn test_read_times_duration_minutes() {
    let mut network = test_network(false);
    network.read_times("DURATION  90  MINUTES").unwrap();
    
    assert_eq!(network.options.time_options.duration, 90 * 60);
  }

  #[test]
  fn test_read_times_duration_min() {
    let mut network = test_network(false);
    network.read_times("DURATION  30  MIN").unwrap();
    
    assert_eq!(network.options.time_options.duration, 30 * 60);
  }

  #[test]
  fn test_read_times_duration_seconds() {
    let mut network = test_network(false);
    network.read_times("DURATION  3600  SECONDS").unwrap();
    
    assert_eq!(network.options.time_options.duration, 3600);
  }

  #[test]
  fn test_read_times_duration_sec() {
    let mut network = test_network(false);
    network.read_times("DURATION  1800  SEC").unwrap();
    
    assert_eq!(network.options.time_options.duration, 1800);
  }

  #[test]
  fn test_read_times_duration_days() {
    let mut network = test_network(false);
    network.read_times("DURATION  2  DAYS").unwrap();
    
    assert_eq!(network.options.time_options.duration, 2 * 86400);
  }

  #[test]
  fn test_read_times_duration_day() {
    let mut network = test_network(false);
    network.read_times("DURATION  1  DAY").unwrap();
    
    assert_eq!(network.options.time_options.duration, 86400);
  }

  // ==================== Control Tests ====================
  #[test]
  fn test_pressure_control_above() {
    let mut network = test_network(false);
    network.read_control("LINK 12 CLOSED IF NODE 23 BELOW 20").unwrap();

    assert_eq!(network.controls.len(), 1);
    let control = network.controls.get(0).unwrap();
    assert_eq!(control.link_id, "12".into());
    assert_eq!(control.setting, None);
    assert_eq!(control.status, Some(LinkStatus::Closed));
    let ControlCondition::Pressure { node_id, above, below } = &control.condition else {
      panic!("Expected Pressure control condition");
    };
    assert_eq!(*node_id, "23".into());
    assert!(!above);
    assert_eq!(*below, 20.0);
  }
  #[test]
  fn test_pressure_control_setting() {
    let mut network = test_network(false);
    network.read_control("LINK 12 1.5 IF NODE 23 ABOVE 20").unwrap();

    assert_eq!(network.controls.len(), 1);
    let control = network.controls.get(0).unwrap();
    assert_eq!(control.link_id, "12".into());
    assert_eq!(control.setting, Some(1.5));
  }

  #[test]
  fn test_time_control() {
    let mut network = test_network(false);
    network.read_control("LINK 12 CLOSED IF TIME 1:15").unwrap();

    assert_eq!(network.controls.len(), 1);
    let control = network.controls.get(0).unwrap();
    assert_eq!(control.link_id, "12".into());
    let ControlCondition::Time { seconds } = &control.condition else {
      panic!("Expected Time control condition");
    };
    assert_eq!(*seconds, 1 * 3600 + 15 * 60);
  }
  #[test]
  fn test_clock_time_control() {
    let mut network = test_network(false);
    network.read_control("LINK 12 CLOSED IF CLOCKTIME 1:15 PM").unwrap();

    assert_eq!(network.controls.len(), 1);
    let control = network.controls.get(0).unwrap();
    assert_eq!(control.link_id, "12".into());
    let ControlCondition::ClockTime { seconds } = &control.condition else {
      panic!("Expected ClockTime control condition");
    };
    assert_eq!(*seconds, 13 * 3600 + 15 * 60);
  }

}

