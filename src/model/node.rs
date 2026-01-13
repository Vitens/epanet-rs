use crate::model::reservoir::Reservoir;
use crate::model::tank::Tank;
use crate::model::junction::Junction;
/// Node struct
pub struct Node {
    pub id: Box<str>,
    pub node_type: NodeType,
    pub elevation: f64,
}

/// Node types
pub enum NodeType {
    Reservoir(Reservoir),
    Tank(Tank),
    Junction(Junction),
}

// helper methods for nodes to check if they are fixed head
impl Node {
  pub fn is_fixed(&self) -> bool {
    matches!(self.node_type, NodeType::Reservoir(_) | NodeType::Tank(_))
  }
}