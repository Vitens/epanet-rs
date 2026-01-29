use hashbrown::HashMap;

use crate::model::link::{Link, LinkType};
use crate::model::node::{Node, NodeType};
use crate::model::curve::Curve;
use crate::model::valve::ValveType;
use crate::model::pattern::Pattern;
use crate::model::options::SimulationOptions;
use crate::model::control::Control;

use serde::{Deserialize, Deserializer, Serialize};

#[derive(Default, Serialize)]
pub struct Network {
    /// Simulation options
    pub options: SimulationOptions,
    /// Nodes in the network
    pub nodes: Vec<Node>,
    /// Links in the network
    pub links: Vec<Link>,

    /// Curves in the network
    pub curves: HashMap<Box<str>, Curve>,

    /// Patterns in the network
    pub patterns: HashMap<Box<str>, Pattern>,

    /// Controls in the network
    pub controls: Vec<Control>,

    // Skip serialization of the node_map and link_map to avoid bloating the file size with unnecessary data.
    #[serde(skip)]
    pub node_map: HashMap<Box<str>, usize>,
    #[serde(skip)]
    pub link_map: HashMap<Box<str>, usize>,
    #[serde(skip)]
    pub contains_pressure_control_valve: bool,
}
impl Network {
  /// Check if the network has any tanks
  pub fn has_tanks(&self) -> bool {
    self.nodes.iter().any(|n| matches!(n.node_type, NodeType::Tank(_)))
  }
}

/// Deserialize a network from a JSON or MessagePack file and build the node_map and link_map.
impl<'de> Deserialize<'de> for Network {
  fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
  where
    D: Deserializer<'de>,
  {
    // Deserialize the network data
    #[derive(Deserialize)]
    struct NetworkData {
      options: SimulationOptions,
      nodes: Vec<Node>,
      links: Vec<Link>,
      curves: HashMap<Box<str>, Curve>,
      patterns: HashMap<Box<str>, Pattern>,
      controls: Vec<Control>,
    }

    let mut data = NetworkData::deserialize(deserializer)?;

    // Build node_map from nodes
    let node_map: HashMap<Box<str>, usize> = data.nodes
      .iter()
      .enumerate()
      .map(|(i, n)| (n.id.clone(), i))
      .collect();

    // Build link_map from links
    let link_map: HashMap<Box<str>, usize> = data.links
      .iter()
      .enumerate()
      .map(|(i, l)| (l.id.clone(), i))
      .collect();

    let mut contains_pressure_control_valve = false;
    // Update link start and end node indices, and check if the network contains a pressure control valve
    for link in data.links.iter_mut() {
      link.start_node = *node_map.get(&link.start_node_id).unwrap();
      link.end_node = *node_map.get(&link.end_node_id).unwrap();
      if let LinkType::Valve(valve) = &mut link.link_type {
        if valve.valve_type == ValveType::PSV || valve.valve_type == ValveType::PRV {
          contains_pressure_control_valve = true;
        }
      }
    }

    Ok(Network {
      options: data.options,
      nodes: data.nodes,
      links: data.links,
      curves: data.curves,
      patterns: data.patterns,
      controls: data.controls,
      node_map,
      link_map,
      contains_pressure_control_valve,
    })
  }
}

/// Network methods to add nodes and links
impl Network {
  pub fn add_node(&mut self, node: Node) -> Result<(), String> {
    if self.node_map.contains_key(&node.id) {
      return Err(format!("Node {} already exists", node.id));
    }
    self.node_map.insert(node.id.clone(), self.nodes.len());
    self.nodes.push(node);
    Ok(())
  }
  pub fn add_link(&mut self, link: Link) -> Result<(), String> {
    if self.link_map.contains_key(&link.id) {
      return Err(format!("Link {} already exists", link.id));
    }
    self.link_map.insert(link.id.clone(), self.links.len());
    self.links.push(link);
    Ok(())
  }
}