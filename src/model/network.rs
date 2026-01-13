use hashbrown::HashMap;

use crate::model::link::Link;
use crate::model::node::Node;
use crate::model::curve::Curve;
use crate::model::pattern::Pattern;
use crate::model::options::SimulationOptions;

#[derive(Default)]
pub struct Network {
    pub options: SimulationOptions,
    pub nodes: Vec<Node>,
    pub links: Vec<Link>,

    pub curves: HashMap<Box<str>, Curve>,
    pub patterns: HashMap<Box<str>, Pattern>,

    pub node_map: HashMap<Box<str>, usize>,
    pub link_map: HashMap<Box<str>, usize>,
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