//! Assembly of the Jacobian matrix and right-hand side vector for the GGA system.

use crate::constants::{PDA_MIN_DIFF, SMALL_VALUE};
use crate::model::link::LinkTrait;
use crate::model::network::Network;
use crate::model::node::NodeType;
use crate::model::options::DemandModel;
use crate::solver::hydraulicsolver::HydraulicSolver;
use crate::solver::matrix::ResistanceCoefficients;
use crate::solver::state::SolverState;

impl HydraulicSolver {
    // Assemble the Jacobian matrix and the right-hand side vector for the system of equations
    // First assemble the contributions from the links
    // Then assemble the contributions from the emitters (virtual links)
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_jacobian(
        &self,
        network: &Network,
        state: &mut SolverState,
        values: &mut [f64],
        rhs: &mut [f64],
        link_coefficients: &mut ResistanceCoefficients,
        excess_flows: &[f64],
        grounded_nodes: &[bool],
    ) {
        // assemble the contributions from the links
        self.link_contributions(network, state, values, rhs, link_coefficients, excess_flows);
        // assemble the contributions from the emitters
        self.emitter_contributions(network, state, values, rhs);
        // if Pressure Dependent Analysis, assemble the contributions from the pressure dependent demand
        if network.options.demand_model == DemandModel::PDA {
            self.demand_contributions(network, state, values, rhs);
        }

        // update the contributions from the nodes
        self.node_contributions(network, state, rhs);

        // add virtual reservoir connections for grounded nodes
        self.grounded_node_contributions(values, grounded_nodes);
    }

    /// For nodes grounded to a virtual reservoir at elevation 0, add a small
    /// conductance (SMALL_VALUE) to the diagonal. This effectively pins the head
    /// at 0, equivalent to a direct connection to a reservoir, making the matrix
    /// non-singular for disconnected or ill-conditioned sub-networks.
    fn grounded_node_contributions(&self, values: &mut [f64], grounded_nodes: &[bool]) {
        for (i, &grounded) in grounded_nodes.iter().enumerate() {
            if grounded && let Some(row) = self.node_rows[i] {
                values[row] += SMALL_VALUE;
            }
        }
    }

    fn node_contributions(&self, network: &Network, state: &mut SolverState, rhs: &mut [f64]) {
        for (i, node) in network.nodes.iter().enumerate() {
            if let NodeType::Junction(_) = &node.node_type {
                let idx = self.node_to_unknown[i].unwrap();
                if network.options.demand_model == DemandModel::PDA {
                    rhs[idx] -= state.demand_flows[i];
                } else {
                    rhs[idx] -= state.demands[i];
                }
            }
        }
    }

    fn link_contributions(
        &self,
        network: &Network,
        state: &mut SolverState,
        values: &mut [f64],
        rhs: &mut [f64],
        link_coefficients: &mut ResistanceCoefficients,
        excess_flows: &[f64],
    ) {
        // iterate over the links
        for (i, link) in network.links.iter().enumerate() {
            let q = state.flows[i];
            let csc_index = &self.csc_indices[i];
            let coefficients = link.coefficients(
                q,
                state.resistances[i],
                state.settings[i],
                state.statuses[i],
                excess_flows[link.start_node],
                excess_flows[link.end_node],
            );

            link_coefficients.g_inv[i] = coefficients.g_inv;
            link_coefficients.y[i] = coefficients.y;

            // update the status of the link if it has a new status
            if let Some(status) = coefficients.new_status {
                state.statuses[i] = status;
            }

            // Get the CSC indices for the start and end nodes
            let u = self.node_to_unknown[link.start_node];
            let v = self.node_to_unknown[link.end_node];

            if let Some(i) = u {
                values[csc_index.diag_u.unwrap()] += coefficients.g_inv;
                rhs[i] -= q - coefficients.y;
                if network.nodes[link.end_node].is_fixed() {
                    rhs[i] += coefficients.g_inv * state.heads[link.end_node];
                }
            }
            if let Some(j) = v {
                values[csc_index.diag_v.unwrap()] += coefficients.g_inv;
                rhs[j] += q - coefficients.y;
                if network.nodes[link.start_node].is_fixed() {
                    rhs[j] += coefficients.g_inv * state.heads[link.start_node];
                }
            }
            if let (Some(_i), Some(_j)) = (u, v) {
                values[csc_index.off_diag.unwrap()] -= coefficients.g_inv;
            }

            // apply the upstream/downstream modifications to the nodes (for PSV/PRV valves)
            if let Some(upstream_modification) = coefficients.upstream_modification {
                rhs[u.unwrap()] += upstream_modification.rhs_add;
                values[csc_index.diag_u.unwrap()] += upstream_modification.diagonal_add;
            }

            if let Some(downstream_modification) = coefficients.downstream_modification {
                rhs[v.unwrap()] += downstream_modification.rhs_add;
                values[csc_index.diag_v.unwrap()] += downstream_modification.diagonal_add;
            }
        }
    }

    fn emitter_contributions(
        &self,
        network: &Network,
        state: &mut SolverState,
        values: &mut [f64],
        rhs: &mut [f64],
    ) {
        // iterate over emitters
        for (i, node) in network.nodes.iter().enumerate() {
            if let NodeType::Junction(junction) = &node.node_type
                && junction.emitter_coefficient > 0.0
            {
                // get the index for the diagonal entry in the Jacobian matrix
                let row = self.node_rows[i].unwrap();
                // get the index for the unknown node in the RHS vector
                let idx = self.node_to_unknown[i].unwrap();

                let (g_inv, y) = junction
                    .emitter_coefficients(state.emitter_flows[i], network.options.emitter_exponent);
                // update RHS
                rhs[idx] += (y + node.elevation) * g_inv - state.emitter_flows[i];
                // update matrix diagonal
                values[row] += g_inv;
            }
        }
    }

    fn demand_contributions(
        &self,
        network: &Network,
        state: &mut SolverState,
        values: &mut [f64],
        rhs: &mut [f64],
    ) {
        let options = &network.options;

        // Get demand function parameters
        let dp = (options.required_pressure - options.minimum_pressure).max(PDA_MIN_DIFF);
        let n = 1.0 / options.pressure_exponent;

        // Iterate over all junctions
        for (i, node) in network.nodes.iter().enumerate() {
            if let NodeType::Junction(junction) = &node.node_type {
                // only consider junctions with a positive demand
                if state.demands[i] > 0.0 {
                    // get the index for the diagonal entry in the Jacobian matrix
                    let row = self.node_rows[i].unwrap();
                    // get the index for the unknown node in the RHS vector
                    let idx = self.node_to_unknown[i].unwrap();
                    // get coefficients for the demand function
                    let (g_inv, y) = junction.demand_coefficients(
                        state.demand_flows[i],
                        state.demands[i],
                        dp,
                        n,
                    );

                    if (1.0 / g_inv) > 0.0 {
                        // update RHS
                        rhs[idx] += (y + node.elevation + options.minimum_pressure) * g_inv;
                        // update matrix diagonal
                        values[row] += g_inv;
                    }
                }
            }
        }
    }
}
