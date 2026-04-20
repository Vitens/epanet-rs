use hashbrown::{HashMap, HashSet};

use crate::model::link::{Link, LinkType};
use crate::model::node::{Node, NodeType};
use crate::model::curve::Curve;
use crate::model::valve::ValveType;
use crate::model::junction::Junction;

use crate::model::pattern::Pattern;
use crate::model::options::{SimulationOptions, HeadlossFormula};
use crate::model::control::{Control, ControlCondition};
use crate::model::units::{UnitConversion, FlowUnits};

use crate::constants::*;

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
    
    // Network change tracking properties
    /// Flag to indicate if the network topology has changed, forces complete solver and state re-initialization
    #[serde(skip)]
    pub topology_changed: bool,
    /// indices of nodes that have been updated
    #[serde(skip)]
    pub updated_nodes: HashSet<usize>,
    /// indices of links that have been updated
    #[serde(skip)]
    pub updated_links: HashSet<usize>,
}

impl Network {
  pub fn new(flow_units: FlowUnits, headloss_formula: HeadlossFormula) -> Self {
    let mut network = Self::default();
    network.options = SimulationOptions::new(flow_units, headloss_formula);
    network
  }
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
      topology_changed: false,
      updated_nodes: HashSet::new(),
      updated_links: HashSet::new()
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


// Network addition / modification methods
pub struct JunctionData {
  pub elevation: f64,
  pub basedemand: f64,
  pub emitter_coefficient: f64,
  pub pattern: Option<Box<str>>,
  pub coordinates: Option<(f64, f64)>,
}

pub struct JunctionUpdate {
  pub elevation: Option<f64>,
  pub basedemand: Option<f64>,
  pub emitter_coefficient: Option<f64>,
  pub pattern: Option<Box<str>>,
  pub coordinates: Option<(f64, f64)>,
}

impl Network {
  // Add a new junction to the network and return the index of the new junction
  pub fn add_junction(&mut self, id: &str, data: &JunctionData) -> Result<(), InputError> {

    // resolve the pattern index
    let pattern_index = if let Some(pattern) = &data.pattern {
      Some(self.pattern_map.get(pattern).ok_or(InputError::PatternNotFound { pattern_id: pattern.clone() })?)
    } else {
      None
    };

    let junction = Junction {
      basedemand: data.basedemand,
      emitter_coefficient: data.emitter_coefficient,
      pattern: data.pattern.clone(),
      pattern_index: pattern_index.copied(),
    };
    let mut node = Node {
      id: id.into(),
      elevation: data.elevation,
      node_type: NodeType::Junction(junction),
      coordinates: data.coordinates,
    };

    // convert the node to standard units
    node.convert_to_standard(&self.options);
    // add the node to the network
    self.add_node(node)?;
    // set the topology changed flag to true (force re-initialization of the solver)
    self.topology_changed = true;

    Ok(())
  }

  pub fn update_junction(&mut self, id: &str, update: &JunctionUpdate) -> Result<(), InputError> {
    let node_index = self.node_map.get(id).ok_or(InputError::NodeNotFound { node_id: id.into() })?;

    let node = &mut self.nodes[*node_index];
    // convert to user units
    node.convert_from_standard(&self.options);

    if let NodeType::Junction(junction) = &mut node.node_type {
      node.elevation = update.elevation.unwrap_or(node.elevation);

      if let Some(coordinates) = update.coordinates {
        node.coordinates = Some(coordinates);
      }

      junction.basedemand = update.basedemand.unwrap_or(junction.basedemand);

      junction.emitter_coefficient = update.emitter_coefficient.unwrap_or(junction.emitter_coefficient);

      junction.pattern = update.pattern.clone();

      junction.pattern_index = update.pattern.as_ref().and_then(|pattern| self.pattern_map.get(pattern).copied());
    }
    else {
      return Err(InputError::NodeNotAJunction { node_id: id.into() });
    }
    // convert back to standard units
    node.convert_to_standard(&self.options);

    // add the node to the updated nodes set, forcing state update
    self.updated_nodes.insert(*node_index);

    Ok(())
  }

}


mod tests {
  use super::*;

  #[test]
  fn test_add_and_update_junction() {
    let mut network = Network::default();
    let data = JunctionData {
      elevation: 100.0,
      basedemand: 10.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    };
    network.add_junction("J1", &data).unwrap();

    assert_eq!(network.nodes.len(), 1);
    assert_eq!(network.nodes[0].id, "J1".into());
    assert_eq!(network.nodes[0].elevation, 100.0);
    assert_eq!(network.topology_changed, true);

    let update = JunctionUpdate {
      elevation: Some(200.0),
      basedemand: Some(20.0),
      emitter_coefficient: Some(0.5),
      pattern: None,
      coordinates: None,
    };
    network.update_junction("J1", &update).unwrap();
    assert_eq!(network.nodes[0].elevation, 200.0);
    let NodeType::Junction(junction) = &network.nodes[0].node_type else {
      panic!("Expected Junction node type");
    };
    assert_eq!(junction.basedemand, 20.0);
    assert!((junction.emitter_coefficient - 9.231479).abs() < 1e-6);
  }

  #[test]
  fn test_add_junction_si() {
    let mut network = Network::new(FlowUnits::CMH, HeadlossFormula::HazenWilliams);
    let data = JunctionData {
      elevation: 100.0,
      basedemand: 10.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    };
    network.add_junction("J1", &data).unwrap();
    assert_eq!(network.nodes[0].elevation, 100.0 / MperFT);
  }

}