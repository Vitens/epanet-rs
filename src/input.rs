use std::fs::File;
use std::io::{BufReader, BufRead};

use crate::network::*;

#[derive(Debug)]
enum ReadState {
  Junctions,
  Pipes,
  Reservoirs,
  Demands,
  Valves,
  Pumps,
  Curves,
  None,
}


impl Network {
  /// Read a network from an INP file.
  pub fn from_inp(inp: &str) -> Result<Network, String> {

    // create a new network
    let mut network = Network::default();

    network.nodes.reserve(25000);
    network.links.reserve(25000);

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
          _ => ReadState::None,
        }
      }
      else {
        let mut parts = line.split_whitespace();

        match state {
          ReadState::Junctions => {
            // read the junction data
            let id = parts.next().unwrap().into();
            // read the elevation (optional, default 0.0)
            let elevation = parts.next().unwrap().parse::<f64>().unwrap_or(0.0) * UCF_H;
            // read the demand (optional, default 0.0)
            let demand = parts.next().and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0) * UCF_Q;

            // add the node to the network
            network.add_node(Node {
              id,
              elevation,
              demand: 0.0,
              node_type: NodeType::Junction { basedemand: demand },
              result: NodeResult::default(),
            }).unwrap();
          }
          ReadState::Valves => {
            let id = parts.next().unwrap().into();
            // read the start node
            let start_node: Box<str> = parts.next().unwrap().into();
            // read the end node
            let end_node: Box<str> = parts.next().unwrap().into();
            // read the diameter
            let diameter = parts.next().unwrap().parse::<f64>().unwrap() * UCF_D;
            // create the link
            let valve_type = match parts.next().unwrap().to_uppercase().as_str() {
              "GPV" => ValveType::GPV,
              "TCV" => ValveType::TCV,
              "FCV" => ValveType::FCV,
              "PRV" => ValveType::PRV,
              "PSV" => ValveType::PSV,
              "PBV" => ValveType::PBV,
              "PCV" => ValveType::PCV,
              _ => return Err(format!("Invalid valve type: {}", parts.next().unwrap())),
            };
            // read the setting
            let setting = parts.next().unwrap().parse::<f64>().unwrap_or(0.0);
            // read the minor loss
            let minor_loss = parts.next().unwrap().parse::<f64>().unwrap_or(0.0);

            // read the curve ID (optional, default none)
            let curve_id = if let Some(curve_id) = parts.next() {
              Some(curve_id.to_string().into())
            } else {
              None
            };

            let start_node_index = *network.node_map.get(&start_node).unwrap();
            let end_node_index = *network.node_map.get(&end_node).unwrap();

            let _ = network.add_link(Link {
              id,
              start_node: start_node_index,
              end_node: end_node_index,
              link_type: LinkType::Valve { diameter, setting, curve: curve_id, valve_type },
              minor_loss: minor_loss,
              resistance: 0.0,
              result: LinkResult::default(),
              csc_index: CSCIndex::default(),
              g_inv: 0.0,
              y: 0.0,
            });
          }
          ReadState::Pipes => {
            let id = parts.next().unwrap().into();
            // read the start node
            let start_node: Box<str> = parts.next().unwrap().into();
            // read the end node
            let end_node: Box<str> = parts.next().unwrap().into();
            // read the length
            let length = parts.next().unwrap().parse::<f64>().unwrap() * UCF_H;
            // read the diameter
            let diameter = parts.next().unwrap().parse::<f64>().unwrap() * UCF_D;
            // read the roughness
            let roughness = parts.next().unwrap().parse::<f64>().unwrap();
            let roughness = 100.0;

            // create the link
            let minor_loss = parts.next().unwrap().parse::<f64>().unwrap_or(0.0);
            // convert minor los
            let minor_loss = 0.02517 * minor_loss / diameter.powi(2) / diameter.powi(2);

            let start_node_index = *network.node_map.get(&start_node).unwrap();
            let end_node_index = *network.node_map.get(&end_node).unwrap();

            let _ = network.add_link(Link {
              id,
              start_node: start_node_index,
              end_node: end_node_index,
              resistance: 0.0,
              minor_loss: minor_loss,
              link_type: LinkType::Pipe { diameter, length, roughness },
              result: LinkResult::default(),
              csc_index: CSCIndex::default(),
              g_inv: 0.0,
              y: 0.0
            });
          }
          ReadState::Reservoirs => {
            let id = parts.next().unwrap().into();
            // read the elevation
            let elevation = parts.next().unwrap().parse::<f64>().unwrap() * UCF_H;
            // add the node to the network
            let _ = network.add_node(Node {
              id,
              elevation,
              demand: 0.0,
              node_type: NodeType::Reservoir,
              result: NodeResult::default(),
            });
          }
          ReadState::Pumps => {

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
            let start_node_index = *network.node_map.get(&start_node).unwrap();
            let end_node_index = *network.node_map.get(&end_node).unwrap();

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

            let _ = network.add_link(Link {
              id,
              start_node: start_node_index,
              end_node: end_node_index,
              minor_loss: 0.0,
              resistance: 0.0,
              link_type: LinkType::Pump { speed, head_curve, power },
              result: LinkResult::default(),
              csc_index: CSCIndex::default(),
              g_inv: 0.0,
              y: 0.0,
            });

          }
          ReadState::Curves => {
            let id: Box<str> = parts.next().unwrap().into();
            let x = parts.next().unwrap().parse::<f64>().unwrap();
            let y = parts.next().unwrap().parse::<f64>().unwrap();
            // create curve if it does not exist
            if !network.curves.contains_key(&id) {
              network.curves.insert(id.clone(), Curve { id, x: vec![x], y: vec![y] });
            }
            // otherwise, add the point to the curve
            else {
              let curve = network.curves.get_mut(&id).unwrap();
              // check if the x value is in ascending order
              if x < *curve.x.last().unwrap() {
                return Err(format!("X values must be in ascending order for curve {}", id));
              }
              // add the point to the curve
              curve.x.push(x);
              curve.y.push(y);
            }
          }
          ReadState::Demands => {

            let id : Box<str> = parts.next().unwrap().into();
            // read the demand
            let demand = parts.next().unwrap().parse::<f64>().unwrap() * UCF_Q;
            // add the demand to the network
            let node_index = *network.node_map.get(&id).unwrap();
            // update basedemand for the node
            if let NodeType::Junction { basedemand } = &mut network.nodes[node_index].node_type {
              *basedemand = demand;
            }
          }
          ReadState::None => {
            // skip unknown state
          }
        }
      }
      // clear the line buffer
      line_buffer.clear();
    }

    Ok(network)
  }
}
