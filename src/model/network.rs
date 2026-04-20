use hashbrown::HashMap;

use crate::model::link::{Link, LinkType};
use crate::model::node::{Node, NodeType};
use crate::model::curve::Curve;
use crate::model::valve::ValveType;
use crate::model::pattern::Pattern;
use crate::model::options::SimulationOptions;
use crate::model::control::{Control, ControlCondition};
use crate::model::units::UnitConversion;

use crate::error::InputError;

use serde::{Deserialize, Deserializer, Serialize};

#[derive(Default, Serialize, Clone)]
pub struct Network {
    /// Title of the network
    pub title: Option<Box<str>>,
    /// Simulation options
    pub options: SimulationOptions,
    /// Nodes in the network
    pub nodes: Vec<Node>,
    /// Links in the network
    pub links: Vec<Link>,

    /// Curves in the network
    pub curves: Vec<Curve>,

    /// Patterns in the network
    pub patterns: Vec<Pattern>,

    /// Controls in the network
    pub controls: Vec<Control>,

    // Skip serialization of the node_map and link_map to avoid bloating the file size with unnecessary data.
    #[serde(skip)]
    pub node_map: HashMap<Box<str>, usize>,
    #[serde(skip)]
    pub link_map: HashMap<Box<str>, usize>,
    #[serde(skip)]
    pub curve_map: HashMap<Box<str>, usize>,
    #[serde(skip)]
    pub pattern_map: HashMap<Box<str>, usize>,
    #[serde(skip)]
    pub contains_pressure_control_valve: bool,
    #[serde(skip)]
    // Network version, incremented when the network topology is modified
    pub network_version: i32,
}
impl Network {
  /// Check if the network has any tanks
  pub fn has_tanks(&self) -> bool {
    self.nodes.iter().any(|n| matches!(n.node_type, NodeType::Tank(_)))
  }
  pub fn has_pressure_controls(&self) -> bool {
    self.controls.iter().any(|c| matches!(c.condition, ControlCondition::LowPressure { .. } | ControlCondition::HighPressure { .. }))
  }
  pub fn has_quality(&self) -> bool {
    return false
  }

  /// Resolves pattern name references to cached vector indices on each node.
  /// Must be called after all patterns and nodes have been loaded.
  pub fn resolve_pattern_indices(&mut self) {
    for node in self.nodes.iter_mut() {
      match &mut node.node_type {
        NodeType::Junction(junction) => {
          junction.pattern_index = junction.pattern.as_ref()
            .and_then(|id| self.pattern_map.get(id).copied());
        }
        NodeType::Reservoir(reservoir) => {
          reservoir.head_pattern_index = reservoir.head_pattern.as_ref()
            .and_then(|id| self.pattern_map.get(id).copied());
        }
        _ => {}
      }
    }
  }
}

/// Deserialize a network from a JSON or MessagePack file and build the node_map and link_map.
impl<'de> Deserialize<'de> for Network {
  fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
  where
    D: Deserializer<'de>,
  {
    #[derive(Deserialize)]
    struct NetworkData {
      title: Option<Box<str>>,
      options: SimulationOptions,
      nodes: Vec<Node>,
      links: Vec<Link>,
      curves: Vec<Curve>,
      patterns: Vec<Pattern>,
      controls: Vec<Control>,
    }

    let mut data = NetworkData::deserialize(deserializer)?;

    let node_map: HashMap<Box<str>, usize> = data.nodes
      .iter().enumerate().map(|(i, n)| (n.id.clone(), i)).collect();

    let link_map: HashMap<Box<str>, usize> = data.links
      .iter().enumerate().map(|(i, l)| (l.id.clone(), i)).collect();

    let curve_map: HashMap<Box<str>, usize> = data.curves
      .iter().enumerate().map(|(i, c)| (c.id.clone(), i)).collect();

    let pattern_map: HashMap<Box<str>, usize> = data.patterns
      .iter().enumerate().map(|(i, p)| (p.id.clone(), i)).collect();

    let mut contains_pressure_control_valve = false;
    for link in data.links.iter_mut() {
      link.start_node = *node_map.get(&link.start_node_id).unwrap();
      link.end_node = *node_map.get(&link.end_node_id).unwrap();
      if let LinkType::Valve(valve) = &mut link.link_type {
        if valve.valve_type == ValveType::PSV || valve.valve_type == ValveType::PRV {
          contains_pressure_control_valve = true;
        }
      }
    }

    let mut network = Network {
      title: data.title,
      options: data.options,
      nodes: data.nodes,
      links: data.links,
      curves: data.curves,
      patterns: data.patterns,
      controls: data.controls,
      node_map,
      link_map,
      curve_map,
      pattern_map,
      contains_pressure_control_valve,
      network_version: 0,
    };

    network.resolve_pattern_indices();
    network.convert_to_standard(&network.options.clone());

    Ok(network)
  }
}

/// Cargo-private Utility methods to add nodes and links
impl Network {
  pub(crate) fn add_node(&mut self, node: Node) -> Result<(), InputError> {
    if self.node_map.contains_key(&node.id) {
      return Err(InputError::NodeExists { node_id: node.id.clone() });
    }
    self.node_map.insert(node.id.clone(), self.nodes.len());
    self.nodes.push(node);
    Ok(())
  }
  pub(crate) fn add_link(&mut self, link: Link) -> Result<(), InputError> {
    if self.link_map.contains_key(&link.id) {
      return Err(InputError::LinkExists { link_id: link.id.clone() });
    }
    self.link_map.insert(link.id.clone(), self.links.len());
    self.links.push(link);
    Ok(())
  }
}

impl UnitConversion for Network {
  fn convert_to_standard(&mut self, options: &SimulationOptions) {
    // convert the nodes to standard units
    for node in self.nodes.iter_mut() {
      node.convert_to_standard(options);
    }
    // convert the links to standard units
    for link in self.links.iter_mut() {
      link.convert_to_standard(options);
    }

    // convert the options to standard units
    self.options.convert_to_standard();
  }
  fn convert_from_standard(&mut self, options: &SimulationOptions) {
    // convert the nodes from standard units
    for node in self.nodes.iter_mut() {
      node.convert_from_standard(options);
    }
    // convert the links from standard units
    for link in self.links.iter_mut() {
      link.convert_from_standard(options);
    }
    // convert the options from standard units
    self.options.convert_from_standard();
  }
}