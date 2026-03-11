use crate::model::link::{LinkTrait, LinkStatus};
use crate::model::network::Network;
use crate::model::node::NodeType;

/// The solver state is the initial/final state of the solver for a single step
#[derive(Debug, Clone)]
pub struct SolverState {
  pub flows: Vec<f64>,
  pub heads: Vec<f64>,
  pub demands: Vec<f64>,
  pub emitter_flows: Vec<f64>,
  pub demand_flows: Vec<f64>,
  pub statuses: Vec<LinkStatus>,
  pub settings: Vec<f64>,
  pub resistances: Vec<f64>,
}

impl SolverState {
  /// Create a new solver state with the initial values for the flows, heads, demands and statuses and calculate resistances
  pub fn new_with_initial_values(network: &Network) -> Self {

    // get the initial emitter flows
    let initial_emitter_flows = network.nodes.iter().map(|n| if let NodeType::Junction(junction) = &n.node_type {
      if junction.emitter_coefficient > 0.0 { 1.0 } else { 0.0 }
    } else { 0.0 }
    ).collect::<Vec<f64>>();


    Self { flows: network.links.iter().map(|l| l.initial_flow()).collect::<Vec<f64>>(), 
           heads: network.nodes.iter().map(|n| n.initial_head()).collect::<Vec<f64>>(), 
           emitter_flows: initial_emitter_flows,
           demands: vec![0.0; network.nodes.len()], 
           demand_flows: vec![0.0; network.nodes.len()],
           settings: network.links.iter().map(|l| l.initial_setting()).collect::<Vec<f64>>(),
           statuses: network.links.iter().map(|l| l.initial_status).collect::<Vec<LinkStatus>>(),
           resistances: network.links.iter().map(|l| l.resistance()).collect::<Vec<f64>>(),
         }
  }
}