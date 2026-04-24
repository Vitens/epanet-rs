//! Serializer that writes a [`Network`] back out as an EPANET `.inp` text file.

use std::io::{BufWriter, Write};
use std::fs::File;

use crate::model::network::Network;
use crate::model::node::NodeType;
use crate::model::link::{LinkType, LinkStatus};
use crate::model::control::ControlCondition;
use crate::model::valve::ValveType;
use crate::model::options::HeadlossFormula;
use crate::utils::time::seconds_to_hhmm;


fn write_line(buffer: &mut String, line: &str) {
  buffer.push_str(line);
  buffer.push_str("\n");
}

fn write_section(buffer: &mut String, section: &str, columns: &[(&str, usize)]) {

  write_line(buffer, ""); // empty line
  write_line(buffer, &format!("[{}]", section));


  let mut header = String::from(";");
  let mut separator = String::from(";");
  for (i, (name, width)) in columns.iter().enumerate() {
    let length = if i > 0 { width + 1 } else { *width };
    header.push_str(&format!("{:<length$}", name));
    separator.push_str(&"-".repeat(length));
  }
  write_line(buffer, &header);
  write_line(buffer, &separator);
}

pub fn write_inp(network: &Network, mut writer: BufWriter<File>) -> Result<(), String> {

  // write the network to the inp file
  let mut buffer = String::new();

  // write the title
  if let Some(title) = &network.title {
    write_line(&mut buffer, &format!("[TITLE]\n{}", title));
  }

  // write the junctions
  write_section(&mut buffer, "JUNCTIONS", &[("ID", 10), ("Elev", 12), ("Demand", 12), ("Pattern", 12)]);

  for node in network.nodes.iter() {
    if let NodeType::Junction(junction) = &node.node_type {
      write_line(&mut buffer, &format!("{:<10} {:<12} {:<12} {}", node.id, node.elevation, junction.basedemand, junction.pattern.as_deref().unwrap_or("")));
    }
  }

  // write the reservoirs
  write_section(&mut buffer, "RESERVOIRS", &[("ID", 10), ("Head", 15), ("Pattern", 12)]);
  for node in network.nodes.iter() {
    if let NodeType::Reservoir(reservoir) = &node.node_type {
      write_line(&mut buffer, &format!("{:<10} {:<12} {}", node.id, node.elevation, reservoir.head_pattern.as_deref().unwrap_or("")));
    }
  }
  
  // write the tanks
  write_section(&mut buffer, "TANKS", &[("ID", 10), ("Elevation", 15), ("InitLevel", 15), ("MinLevel", 15), ("MaxLevel", 15), ("Diameter", 15), ("MinVol", 15), ("VolCurve", 15), ("Overflow", 12)]);
  for node in network.nodes.iter() {
    if let NodeType::Tank(tank) = &node.node_type {
      write_line(&mut buffer, &format!("{:<10} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}", node.id, node.elevation, tank.initial_level, tank.min_level, tank.max_level, tank.diameter, tank.min_volume, tank.volume_curve_id.as_deref().unwrap_or(""), tank.overflow));
    }
  }
  
  // write the pipes
  write_section(&mut buffer, "PIPES", &[("ID", 10), ("Node1", 12), ("Node2", 12), ("Length", 12), ("Diameter", 12), ("Roughness", 12), ("MinorLoss", 12), ("Status", 12)]);
  for link in network.links.iter() {
    if let LinkType::Pipe(pipe) = &link.link_type {

      let mut status = match link.initial_status {
        LinkStatus::Open => "Open",
        LinkStatus::Closed => "Closed",
        _ => ""
      };
      if pipe.check_valve {
        status = "CV";
      }


      write_line(&mut buffer, &format!("{:<10} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}", link.id, link.start_node_id, link.end_node_id, pipe.length, pipe.diameter, pipe.roughness, pipe.minor_loss, status));
    }
  }

  // write the pumps
  write_section(&mut buffer, "PUMPS", &[("ID", 10), ("Node1", 12), ("Node2", 12), ("Parameters", 12)]);
  for link in network.links.iter() {
    if let LinkType::Pump(pump) = &link.link_type {
      let mut parameters = String::new();

      if pump.head_curve_id.is_some() {
        parameters.push_str(&format!("HEAD {}", pump.head_curve_id.as_deref().unwrap()));
      }
      if pump.speed != 1.0 {
        parameters.push_str(&format!("SPEED {}", pump.speed));
      }
      if pump.power != 0.0 {
        parameters.push_str(&format!("POWER {}", pump.power));
      }

      write_line(&mut buffer, &format!("{:<10} {:<12} {:<12} {}", link.id, link.start_node_id, link.end_node_id, parameters));
    }
  }

  // write the valves
  write_section(&mut buffer, "VALVES", &[("ID", 10), ("Node1", 12), ("Node2", 12), ("Diameter", 12), ("Type", 12), ("Setting", 12), ("MinorLoss", 12), ("PCV Curve", 12)]);
  for link in network.links.iter() {


    if let LinkType::Valve(valve) = &link.link_type {

      let curve = if valve.valve_type == ValveType::PCV {
        valve.curve_id.as_deref().unwrap_or("")
      } else {
        ""
      };

      let setting = if valve.valve_type == ValveType::GPV {
        valve.curve_id.as_deref().unwrap_or("")
      } else {
        &format!("{}", valve.setting)
      };

      let valve_type_str = valve.valve_type.to_string();
      write_line(&mut buffer, &format!("{:<10} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}", link.id, link.start_node_id, link.end_node_id, valve.diameter, valve_type_str, setting, valve.minor_loss, curve));
    }
  }

  // write demands
  write_section(&mut buffer, "DEMANDS", &[("Junction", 10), ("Demand", 12), ("Pattern", 12)]);

  for node in network.nodes.iter() {
    if let NodeType::Junction(junction) = &node.node_type {
      write_line(&mut buffer, &format!("{:<10} {:<12} {:<12}", node.id, junction.basedemand, junction.pattern.as_deref().unwrap_or("")));
    }
  }
  // write emitters
  write_section(&mut buffer, "EMITTERS", &[("Junction", 10), ("Coefficient", 12)]);
  for node in network.nodes.iter() {
    if let NodeType::Junction(junction) = &node.node_type {
      if junction.emitter_coefficient > 0.0 {
        write_line(&mut buffer, &format!("{:<10} {:<12}", node.id, junction.emitter_coefficient));
      }
    }
  }

  // write curves
  write_section(&mut buffer, "CURVES", &[("ID", 10), ("X-Value", 12), ("Y-Value", 12)]);
  for curve in network.curves.iter() {
    for (x, y) in curve.x.iter().zip(curve.y.iter()) {
      write_line(&mut buffer, &format!("{:<10} {:<12} {:<12}", curve.id, x, y));
    }
    write_line(&mut buffer, ";");
  }

  // write patterns
  write_section(&mut buffer, "PATTERNS", &[("ID", 10), ("Multipliers", 12)]);
  for pattern in network.patterns.iter() {
    for multiplier in pattern.multipliers.iter() {
      write_line(&mut buffer, &format!("{:<10} {:<12}", pattern.id, multiplier));
    }
    write_line(&mut buffer, ";");
  }

  // write controls
  write_line(&mut buffer, "\n[CONTROLS]");

  for control in network.controls.iter() {

    let mut control_str = String::from(format!("LINK {}", control.link_id));

    if let Some(status) = control.status {
      control_str.push_str(&format!(" {}", status));
    }
    else {
      control_str.push_str(&format!(" {}", control.setting.unwrap_or(0.0)));
    }

    match &control.condition {
      ControlCondition::Time { seconds } => {
        control_str.push_str(&format!(" AT TIME {}", seconds));
      }
      ControlCondition::ClockTime { seconds } => {
        control_str.push_str(&format!(" AT CLOCKTIME {}", seconds));
      }
      ControlCondition::HighPressure { node_index, target } => {
        control_str.push_str(&format!(" IF NODE {} ABOVE {}", node_index, target));
      }
      ControlCondition::LowPressure { node_index, target } => {
        control_str.push_str(&format!(" IF NODE {} BELOW {}", node_index, target));
      }
      ControlCondition::HighLevel { tank_index, target } => {
        control_str.push_str(&format!(" IF NODE {} ABOVE {}", tank_index, target));
      }
      ControlCondition::LowLevel { tank_index, target } => {
        control_str.push_str(&format!(" IF NODE {} BELOW {}", tank_index, target));
      }
    }
    write_line(&mut buffer, &control_str);
  }

  // write times
  write_line(&mut buffer, "\n[TIMES]");

  let opts = &network.options.time_options;

  write_line(&mut buffer, &format!("DURATION {}", seconds_to_hhmm(opts.duration)));
  write_line(&mut buffer, &format!("HYDRAULIC TIMESTEP {}", seconds_to_hhmm(opts.hydraulic_timestep)));
  write_line(&mut buffer, &format!("PATTERN TIMESTEP {}", seconds_to_hhmm(opts.pattern_timestep)));
  write_line(&mut buffer, &format!("REPORT TIMESTEP {}", seconds_to_hhmm(opts.report_timestep)));
  write_line(&mut buffer, &format!("PATTERN START {}", seconds_to_hhmm(opts.pattern_start)));
  write_line(&mut buffer, &format!("START CLOCKTIME {}", seconds_to_hhmm(opts.start_clocktime)));

  // write options
  write_line(&mut buffer, "\n[OPTIONS]");

  write_line(&mut buffer, &format!("Units {}", network.options.flow_units));

  let hl_model = match network.options.headloss_formula {
    HeadlossFormula::HazenWilliams => "H-W",
    HeadlossFormula::DarcyWeisbach => "D-W",
    HeadlossFormula::ChezyManning => "C-M",
  };
  write_line(&mut buffer, &format!("Headloss {}", hl_model));
  write_line(&mut buffer, &format!("DEMAND MODEL {}", network.options.demand_model));
  write_line(&mut buffer, &format!("Specific Gravity {}", network.options.specific_gravity));
  write_line(&mut buffer, &format!("Viscosity {}", network.options.viscosity));
  write_line(&mut buffer, &format!("Trials {}", network.options.max_trials));
  write_line(&mut buffer, &format!("Accuracy {}", network.options.accuracy));
  write_line(&mut buffer, &format!("Pattern {}", network.options.pattern.as_deref().unwrap_or("")));
  write_line(&mut buffer, &format!("Demand Multiplier {}", network.options.demand_multiplier));
  write_line(&mut buffer, &format!("Emitter Exponent {}", network.options.emitter_exponent));

  // write coordinates
  write_section(&mut buffer, "COORDINATES", &[("ID", 10), ("X", 12), ("Y", 12)]);
  for node in network.nodes.iter() {
    if let Some((x, y)) = node.coordinates {
      write_line(&mut buffer, &format!("{:<10} {:<12} {:<12}", node.id, x, y));
    }
  }

  write_section(&mut buffer, "VERTICES", &[("Link", 10), ("X", 12), ("Y", 12)]);
  for link in network.links.iter() {
    if let Some(vertices) = &link.vertices {
      for (x, y) in vertices.iter() {
        write_line(&mut buffer, &format!("{:<10} {:<12} {:<12}", link.id, x, y));
      }
    }
  }

  writer.write(buffer.as_bytes()).map_err(|e| format!("Failed to write network to file: {}", e))?;

  Ok(())
}
