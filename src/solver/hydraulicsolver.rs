//! Hydraulic solver based on the Global Gradient Algorithm (Todini & Pilati, 1987).

use faer::linalg::cholesky::llt::factor::LltError::NonPositivePivot;
use faer::prelude::*;
use faer::sparse::linalg::LltError;
use faer::sparse::linalg::solvers::{Llt, SymbolicLlt};
use faer::sparse::{SparseColMat, SymbolicSparseColMat};
use faer::{Mat, Side};

use crate::solver::matrix::*;
use crate::solver::state::SolverState;

use crate::model::units::Cfs;

use simplelog::{debug, error, warn};

use crate::constants::{BIG_VALUE, PDA_MIN_DIFF, Q_ZERO};
use crate::error::SolverError;
use crate::model::control::ControlCondition;
use crate::model::link::{LinkStatus, LinkTrait, LinkType};
use crate::model::network::Network;
use crate::model::node::NodeType;
use crate::model::options::{DemandModel, SimulationOptions};
use crate::model::valve::ValveType;

/// Flow balance structure containing the total demand, total supply and error
pub struct FlowBalance {
    pub total_demand: Cfs,
    pub total_supply: Cfs,
    pub error: Cfs,
}

/// Iteration statistics structure containing the sum of the flow changes, the sum of the flows, the maximum flow change and the index of the link with the maximum flow change
#[derive(Default)]
pub struct IterationStatistics {
    pub sum_dq: Cfs,
    pub sum_q: Cfs,
    pub max_dq: Cfs,
    pub max_dq_index: usize,
    pub status_changed: bool,
}

impl IterationStatistics {
    pub fn relative_change(&self, options: &SimulationOptions) -> f64 {
        if self.sum_q > options.accuracy {
            self.sum_dq / (self.sum_q + Q_ZERO)
        } else {
            self.sum_dq
        }
    }
    pub fn max_dq_converted(&self, options: &SimulationOptions) -> Cfs {
        self.max_dq * options.flow_units.per_cfs()
    }
}

/// The main hydraulic solver struct
pub struct HydraulicSolver {
    /// global unknown-numbering map: node_to_unknown[node_index] = unknown_index
    pub node_to_unknown: Vec<Option<usize>>,
    /// symbolic sparsity pattern
    pub sparsity_pattern: SymbolicSparseColMat<usize>,
    /// symbolic Cholesky factorization
    pub symbolic_llt: SymbolicLlt<usize>,
    /// AMD fill-reducing permutation: perm_fwd[permuted] = original
    pub perm_fwd: Vec<usize>,
    /// precomputed Jacobian matrix
    pub jac: SparseColMat<usize, f64>,
    /// precomputed CSC indices for the links
    pub csc_indices: Vec<CSCIndex>,
    /// precomputed indices for the rows of the Jacobian matrix for each node
    pub node_rows: Vec<Option<usize>>,
    /// version of the topology of the network on which the solver was created
    pub topology_version: u32,
}

impl HydraulicSolver {
    pub fn new(network: &Network) -> Result<Self, SolverError> {
        // build the sparsity pattern and the global unknown-numbering map
        let node_to_unknown = build_unknown_numbering_map(network);
        let sparsity_pattern = build_sparsity_pattern(network, &node_to_unknown);
        // precompute the CSC indices for the links
        let csc_indices = map_links_to_csc_indices(network, &sparsity_pattern, &node_to_unknown);
        // precompute the indices for the rows of the Jacobian matrix for each node
        let node_rows = map_nodes_to_rows(network, &sparsity_pattern, &node_to_unknown);

        // generate the Jacobian matrix
        let values = vec![0.0; sparsity_pattern.as_ref().row_idx().len()];
        let jac = SparseColMat::new(sparsity_pattern.clone(), values.clone());

        // precompute the AMD fill-reducing permutation for error mapping
        let perm_fwd = compute_amd_permutation(&sparsity_pattern, &node_to_unknown);
        // precompute the symbolic Cholesky factorization
        let symbolic_llt = SymbolicLlt::try_new(jac.symbolic(), Side::Lower)
            .map_err(|e| SolverError::Factorization(e.to_string()))?;

        Ok(Self {
            node_to_unknown,
            sparsity_pattern,
            symbolic_llt,
            perm_fwd,
            jac,
            csc_indices,
            node_rows,
            topology_version: network.topology_version,
        })
    }

    /// Solve the network for a single timestep using the Global Gradient Algorithm (Todini & Pilati, 1987).
    /// Takes a solver state and returns a new state after convergence.
    pub fn solve(
        &self,
        network: &Network,
        state: &SolverState,
    ) -> Result<SolverState, SolverError> {
        // check if there have been changes to the network since the solver was created

        // if the topology version has changed, return an error
        if network.topology_version != self.topology_version {
            return Err(SolverError::StaleTopology);
        }
        // if the properties version has changed, return an error
        if network.properties_version != state.properties_version
            || network.topology_version != self.topology_version
        {
            return Err(SolverError::StaleProperties);
        }

        // clone the solver state to avoid modifying the original
        let mut state = state.clone();

        let unknown_nodes = self.node_to_unknown.iter().filter(|&x| x.is_some()).count();

        let mut values = vec![0.0; self.sparsity_pattern.as_ref().row_idx().len()];
        let mut rhs = vec![0.0; unknown_nodes];

        let mut link_coefficients = ResistanceCoefficients::new(network.links.len());
        let mut jac = self.jac.clone();

        let mut excess_flows = vec![0.0; network.nodes.len()];

        let mut grounded_nodes = vec![false; network.nodes.len()];

        'gga: for iteration in 1..=network.options.max_trials {
            values.fill(0.0);
            rhs.fill(0.0);

            if network.contains_pressure_control_valve {
                self.calculate_excess_flows(network, &state, &mut excess_flows);
            }

            // assemble the Jacobian matrix
            self.assemble_jacobian(
                network,
                &mut state,
                &mut values,
                &mut rhs,
                &mut link_coefficients,
                &excess_flows,
                &grounded_nodes,
            );

            // copy the values to the Jacobian matrix
            jac.val_mut().copy_from_slice(&values);

            // try to factorize the Jacobian matrix
            let llt = match Llt::try_new_with_symbolic(
                self.symbolic_llt.clone(),
                jac.as_ref(),
                Side::Lower,
            ) {
                Ok(llt) => llt,
                Err(LltError::Numeric(NonPositivePivot { index })) => {
                    // get the original unknown index from the AMD permutation and translate it to a node index in the network
                    let original_unknown = self.perm_fwd[index - 1];
                    let node_index = self
                        .node_to_unknown
                        .iter()
                        .position(|&x| x.is_some() && x.unwrap() == original_unknown)
                        .unwrap();

                    // if the factorization failed, attempt to fix the problem by first fixing a possible bad valve
                    if self.fix_bad_valve(network, node_index, &mut state.statuses) {
                        continue 'gga;
                    }
                    // otherwise, ground the node causing the disconnect with a virtual reservoir
                    if !grounded_nodes[node_index] {
                        warn!(
                            "Grounding node '{}' with virtual reservoir (elevation 0) to fix singular matrix",
                            network.nodes[node_index].id
                        );
                        grounded_nodes[node_index] = true;
                        continue 'gga;
                    }
                    let err = SolverError::SingularMatrix {
                        node_id: network.nodes[node_index].id.clone(),
                    };
                    error!("{}", err);
                    return Err(err);
                }
                Err(e) => {
                    return Err(SolverError::Factorization(e.to_string()));
                }
            };

            // solve the system of equations
            let dh = llt.solve(&Mat::from_fn(unknown_nodes, 1, |r, _| rhs[r]));

            // update the heads for the unknown nodes
            for (global, &head_id) in self.node_to_unknown.iter().enumerate() {
                if let Some(i) = head_id {
                    state.heads[global] = dh[(i, 0)];
                }
            }

            // update the links and emitters and gather iteration statistics
            let mut stats = self.update_links(network, &mut state, &link_coefficients);
            self.update_emitter_flows(network, &mut state, &mut stats);

            if network.options.demand_model == DemandModel::PDA {
                self.update_demand_flows(network, &mut state, &mut stats);
            }

            // close/open links connected to tanks based on tank level
            self.update_tank_links(network, &mut state);

            debug!(
                ">> Iteration {}: Relative change: {:.6}, Status changed: {}",
                iteration,
                stats.relative_change(&network.options),
                stats.status_changed
            );
            debug!(
                ">>>> Max flow change: {:.6} for link {}",
                stats.max_dq_converted(&network.options),
                network.links[stats.max_dq_index].id
            );

            let max_flow_change = network.options.max_flow_change.unwrap_or(BIG_VALUE);

            // check for convergence:
            // - not the first iteration
            // - relative change is less than the accuracy
            // - no status changes
            // - maximum flow change is less than the maximum flow change allowed
            if iteration > 1
                && stats.relative_change(&network.options) < network.options.accuracy
                && !stats.status_changed
                && stats.max_dq_converted(&network.options) < max_flow_change
            {
                // if pressure controls are active, apply them and continue the iteration
                if self.apply_pressure_controls(network, &mut state) {
                    continue 'gga;
                }

                // calculate the flow balance, update the node demands
                let flow_balance = self.flow_balance(network, &state.demands, &state.flows);
                if network.options.demand_model == DemandModel::PDA {
                    state.demands = state.demand_flows.clone();
                }
                for i in 0..state.emitter_flows.len() {
                    state.demands[i] += state.emitter_flows[i];
                }
                debug!(
                    "Converged in {} iterations: Error = {:.4}, Supply = {:.4}, Demand = {:.4}",
                    iteration,
                    flow_balance.error,
                    flow_balance.total_supply,
                    flow_balance.total_demand
                );

                return Ok(state);
            }
        }
        Err(SolverError::MaxIterations {
            max_trials: network.options.max_trials,
        })
    }

    fn calculate_excess_flows(
        &self,
        network: &Network,
        state: &SolverState,
        excess_flows: &mut [Cfs],
    ) {
        for (i, emitter_flow) in state.emitter_flows.iter().enumerate() {
            excess_flows[i] = -emitter_flow;
        }
        for (i, demand) in state.demands.iter().enumerate() {
            if network.options.demand_model == DemandModel::PDA {
                excess_flows[i] -= state.demand_flows[i];
            } else {
                excess_flows[i] -= demand;
            }
        }
        for (i, link) in network.links.iter().enumerate() {
            let q = state.flows[i];
            excess_flows[link.start_node] -= q;
            excess_flows[link.end_node] += q;
        }
    }

    /// Update links and gather iteration statistics
    fn update_links(
        &self,
        network: &Network,
        state: &mut SolverState,
        coefficients: &ResistanceCoefficients,
    ) -> IterationStatistics {
        let mut stats = IterationStatistics::default();

        for (i, link) in network.links.iter().enumerate() {
            let dh = state.heads[link.start_node] - state.heads[link.end_node];
            let g_inv = coefficients.g_inv[i];
            let y = coefficients.y[i];

            // calculate the flow change
            let dq = y - g_inv * dh;

            // update the maximum flow change
            if dq.abs() > stats.max_dq {
                stats.max_dq = dq.abs();
                stats.max_dq_index = i;
            }

            // update the link flow
            state.flows[i] -= dq;

            // check if the status of the link has changed
            let new_status = link.update_status(
                state.settings[i],
                state.statuses[i],
                state.flows[i],
                state.heads[link.start_node],
                state.heads[link.end_node],
            );
            if let Some(status) = new_status {
                if state.statuses[i] != LinkStatus::TempClosed
                    && state.statuses[i] != LinkStatus::Xhead
                {
                    stats.status_changed = true;
                }
                debug!(
                    "<yellow>Status changed for link {} from {:?} to {:?}</>",
                    link.id, state.statuses[i], status
                );
                state.statuses[i] = status;
            }

            // update the sum of the flow changes and the sum of the flows
            stats.sum_dq += dq.abs();
            stats.sum_q += state.flows[i].abs();
        }
        stats
    }

    /// Update emitter flows and update iteration statistics
    fn update_emitter_flows(
        &self,
        network: &Network,
        state: &mut SolverState,
        stats: &mut IterationStatistics,
    ) {
        for (i, node) in network.nodes.iter().enumerate() {
            if let NodeType::Junction(junction) = &node.node_type
                && junction.emitter_coefficient > 0.0
            {
                let dh = state.heads[i] - node.elevation;
                let (g_inv, y) = junction
                    .emitter_coefficients(state.emitter_flows[i], network.options.emitter_exponent);
                let dq = (y - dh) * g_inv;
                state.emitter_flows[i] -= dq;
                stats.sum_dq += dq.abs();
                stats.sum_q += state.emitter_flows[i].abs();
                if dq.abs() > stats.max_dq {
                    stats.max_dq = dq.abs();
                }
            }
        }
    }

    /// Update demand flows and update iteration statistics
    fn update_demand_flows(
        &self,
        network: &Network,
        state: &mut SolverState,
        stats: &mut IterationStatistics,
    ) {
        let options = &network.options;

        let dp = (options.required_pressure - options.minimum_pressure).max(PDA_MIN_DIFF);
        let n = 1.0 / options.pressure_exponent;

        for (i, node) in network.nodes.iter().enumerate() {
            if let NodeType::Junction(junction) = &node.node_type
                && state.demands[i] > 0.0
            {
                let (g_inv, y) =
                    junction.demand_coefficients(state.demand_flows[i], state.demands[i], dp, n);

                let dh = state.heads[i] - node.elevation - options.minimum_pressure;
                let dq = (y - dh) * g_inv;

                state.demand_flows[i] -= dq;
                stats.sum_dq += dq.abs();
                stats.sum_q += state.demand_flows[i].abs();
                if dq.abs() > stats.max_dq {
                    stats.max_dq = dq.abs();
                }
            } else {
                // preserve demands for tanks and reservoirs
                state.demand_flows[i] = state.demands[i];
            }
        }
    }

    /// Apply pressure controls to the state
    /// Returns true if any pressure controls were applied
    fn apply_pressure_controls(&self, network: &Network, state: &mut SolverState) -> bool {
        let mut changed = false;

        for control in &network.controls {
            if matches!(
                control.condition,
                ControlCondition::LowPressure { .. } | ControlCondition::HighPressure { .. }
            ) {
                let active = control.is_active(state, network, 0, 0);
                if active {
                    changed = changed || control.activate(state, network);
                }
            }
        }
        changed
    }

    /// Update the links connected to tanks and gather flow balance into/out of tanks
    fn update_tank_links(&self, network: &Network, state: &mut SolverState) {
        for (tank_index, node) in network.nodes.iter().enumerate() {
            if let NodeType::Tank(tank) = &node.node_type {
                state.demands[tank_index] = 0.0;
                let fill_closed =
                    state.heads[tank_index] >= tank.elevation + tank.max_level && !tank.overflow;
                let empty_closed = state.heads[tank_index] <= tank.elevation + tank.min_level;

                for link_index in &tank.links_to {
                    state.demands[tank_index] += state.flows[*link_index];
                    if fill_closed && state.flows[*link_index] > 0.0 {
                        state.statuses[*link_index] = LinkStatus::TempClosed;
                    }
                    if empty_closed && state.flows[*link_index] < 0.0 {
                        state.statuses[*link_index] = LinkStatus::TempClosed;
                    }
                }
                for link_index in &tank.links_from {
                    state.demands[tank_index] -= state.flows[*link_index];
                    if fill_closed && state.flows[*link_index] < 0.0 {
                        state.statuses[*link_index] = LinkStatus::TempClosed;
                    }
                    if empty_closed && state.flows[*link_index] > 0.0 {
                        state.statuses[*link_index] = LinkStatus::TempClosed;
                    }
                }
            }
        }
    }

    /// Fix a bad control (PSV/PRV) valve by setting its status to Closed
    fn fix_bad_valve(
        &self,
        network: &Network,
        node_index: usize,
        statuses: &mut [LinkStatus],
    ) -> bool {
        for (i, link) in network.links.iter().enumerate() {
            if link.start_node != node_index && link.end_node != node_index {
                continue;
            }
            if let LinkType::Valve(valve) = &link.link_type
                && (valve.valve_type == ValveType::PSV || valve.valve_type == ValveType::PRV)
                && statuses[i] == LinkStatus::Active
            {
                debug!("Fixing bad valve for node index: {}", node_index);
                statuses[i] = LinkStatus::XPressure;
                return true;
            }
        }
        false
    }

    /// Calculate the flow balance error
    fn flow_balance(&self, network: &Network, demands: &[Cfs], flows: &[Cfs]) -> FlowBalance {
        let sum_demand: Cfs = demands.iter().sum();

        let mut sum_supply: Cfs = 0.0;
        for (i, link) in network.links.iter().enumerate() {
            if !matches!(
                network.nodes[link.end_node].node_type,
                NodeType::Junction { .. }
            ) {
                sum_supply -= flows[i];
            }
            if !matches!(
                network.nodes[link.start_node].node_type,
                NodeType::Junction { .. }
            ) {
                sum_supply += flows[i];
            }
        }
        let error = sum_demand - sum_supply;
        FlowBalance {
            total_demand: sum_demand,
            total_supply: sum_supply,
            error,
        }
    }
}
