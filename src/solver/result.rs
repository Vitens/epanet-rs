use crate::solver::state::SolverState;
use crate::model::units::{Ft, Cfs, FlowUnits, UnitSystem};
use serde::Serialize;

#[derive(Serialize)]
pub struct SolverResult {
  pub flows: Vec<Vec<Cfs>>,
  pub heads: Vec<Vec<Ft>>,
  pub demands: Vec<Vec<Cfs>>,
}

impl SolverResult {
  pub fn new(n_links: usize, n_nodes: usize, n_steps: usize) -> Self {
    Self { flows: vec![vec![0.0; n_links]; n_steps], heads: vec![vec![0.0; n_nodes]; n_steps], demands: vec![vec![0.0; n_nodes]; n_steps] }
  }

  pub fn append(&mut self, state: &SolverState, step: usize) {
    self.flows[step] = state.flows.clone();
    self.heads[step] = state.heads.clone();
    self.demands[step] = state.demands.clone();
  }

  // convert the solver units back to the original units
  pub fn convert_units(&mut self, to_flow_units: &FlowUnits, to_unit_system: &UnitSystem) {

    let flow_scale = to_flow_units.per_cfs();
    let head_scale = to_unit_system.per_feet();

    self.flows.iter_mut().flatten().for_each(|flow| *flow *= flow_scale);
    self.heads.iter_mut().flatten().for_each(|head| *head *= head_scale);
    self.demands.iter_mut().flatten().for_each(|demand| *demand *= flow_scale);
  }
}