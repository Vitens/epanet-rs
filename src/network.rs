use hashbrown::HashMap;

pub const UCF_Q: f64 = 0.009809629503517; // Unit conversion factor for flow (m3/h to cfs)
pub const UCF_H: f64 = 3.2808399; // Unit conversion factor for head (m to ft)
pub const UCF_D: f64 = 0.0032808399; // Unit conversion factor for diameter (mm to ft)
// pub const UCF_Q: f64 = 1.0;
// pub const UCF_H: f64 = 1.0;
// pub const UCF_D: f64 = 1.0 / 12.0;

pub const C_SMALL: f64 = 1e-6;
pub const C_BIG: f64 = 1e8;

pub const H_EXPONENT: f64 = 1.852; // Hazen-Williams exponent

#[derive(Default)]
pub struct Network {
    pub nodes: Vec<Node>,
    pub links: Vec<Link>,

    pub curves: HashMap<Box<str>, Curve>,

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

pub enum ValveType {
  PRV, // Pressure Reducing Valve
  PSV, // Pressure Sensing Valve
  PBV, // Pressure Breaking Valve
  FCV, // Flow Control Valve
  TCV, // Throttle Control Valve
  PCV, // Positional Control Valve
  GPV, // General Purpose Valve
}
pub enum LinkType {
    Pipe { diameter: f64, length: f64, roughness: f64 },
    Pump { speed: f64, head_curve: Box<str>, power: f64 },
    Valve { diameter: f64, setting: f64, curve: Option<Box<str>>, valve_type: ValveType },
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
  pub minor_loss: f64,
  pub start_node: usize,
  pub end_node: usize,

  pub resistance: f64,

  pub result: LinkResult,

  pub csc_index: CSCIndex,

  pub g_inv: f64,
  pub y: f64,
}

impl Link {
  /// Calculate the 1/G_ij and Y_ij coefficients for the link
  pub fn coefficients(&self, ) -> (f64, f64) {
    match &self.link_type {
      LinkType::Pipe { .. } => self.pipe_coefficients(),
      LinkType::Pump { .. } => self.pump_coefficients(),
      LinkType::Valve { diameter, setting, curve, valve_type } => self.valve_coefficients(*diameter, *setting, valve_type, curve),
    }
  }
  /// Calculate the 1/G_ij and Y_ij coefficients for a pipe
  pub fn pipe_coefficients(&self) -> (f64, f64) {
    let q = self.result.flow;
    let q_abs = q.abs().max(1e-8);
    let r = self.resistance;
    let m = self.minor_loss;
    let n = H_EXPONENT;

    let q_pow = q_abs.powf(n - 1.0);


    let g = n * r * q_pow + 2.0 * m * q_abs;
    let g_inv = 1.0 / g;

    let y = (r * q_abs * q_pow + m * q_abs.powi(2)) * q.signum();

    (g_inv, y)
  }
  pub fn pump_coefficients(&self) -> (f64, f64) {
    (0.0, 0.0)
  }
  pub fn valve_coefficients(&self, diameter: f64, setting: f64, valve_type: &ValveType, curve: &Option<Box<str>>) -> (f64, f64) {
    
    let q = self.result.flow;
    let q_abs = q.abs().max(1e-8);


    let m = 0.02517 * setting / diameter.powi(2) / diameter.powi(2);

    if setting == 0.0 {
      // TCF with a setting of 0 is equivalent to a open valve with almost no resistance
      return (1.0 / C_SMALL, q * C_SMALL);
    }

    let g = 2.0 * m * q_abs;
    let g_inv = 1.0 / g;
    let y = m * q_abs.powf(2.0) * q.signum();

    (g_inv, y)

  }
}

#[derive(Debug, Clone)]
pub struct Curve {
  pub id: Box<str>,
  pub x: Vec<f64>,
  pub y: Vec<f64>,
}