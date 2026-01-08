use hashbrown::HashMap;

pub const UCF_Q: f64 = 0.009809629503517; // Unit conversion factor for flow (m3/h to cfs)
pub const UCF_H: f64 = 3.2808399; // Unit conversion factor for head (m to ft)
pub const UCF_D: f64 = 0.0032808399; // Unit conversion factor for diameter (mm to ft)

#[derive(Default)]
pub struct Network {
    pub nodes: Vec<Node>,
    pub links: Vec<Link>,

    pub node_map: HashMap<Box<str>, usize>,
    pub link_map: HashMap<Box<str>, usize>,
}

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

pub enum NodeType {
    Reservoir,
    Tank,
    Junction { basedemand: f64 },
}

pub enum LinkType {
    Pipe { diameter: f64, length: f64, roughness: f64 },
    Pump,
    Valve,
}

#[derive(Default)]
pub struct NodeResult {
  pub head: f64
}

#[derive(Default)]
pub struct LinkResult {
  pub flow: f64,
}


/// CSC (Compressed Sparse Column) indices for the Jacobian matrix used in the Global Gradient Algorithm
#[derive(Default)]
pub struct CSCIndex {
  pub diag_u: Option<usize>,      // CSC index for J[u,u]
  pub diag_v: Option<usize>,      // CSC index for J[v,v]
  pub off_diag_uv: Option<usize>, // CSC index for J[u,v]
  pub off_diag_vu: Option<usize>, // CSC index for J[v,u]
}

pub struct Node {
    pub id: Box<str>,
    pub node_type: NodeType,
    pub elevation: f64,
    pub demand: f64,

    pub result: NodeResult,
}
impl Node {
  pub fn is_fixed(&self) -> bool {
    matches!(self.node_type, NodeType::Reservoir | NodeType::Tank)
  }
}

pub struct Link {
  pub id: Box<str>,
  pub link_type: LinkType,
  pub start_node: usize,
  pub end_node: usize,

  pub resistance: f64,

  pub result: LinkResult,

  pub csc_index: CSCIndex,
}