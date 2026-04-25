//! Per-step mutable solver state (flows, heads, demands, statuses, settings, resistances).

use crate::model::link::{LinkStatus, LinkTrait};
use crate::model::network::Network;
use crate::model::node::NodeType;

use crate::model::control::ControlCondition;
use crate::model::options::DemandModel;
use crate::model::units::{Cfs, Ft3};

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
    /// version of the topology of the network on which the state was created
    pub topology_version: u32,
    /// version of the properties of the network on which the state was created
    pub properties_version: u32,
}

impl SolverState {
    /// Create a new solver state with the initial values for the flows, heads, demands and statuses and calculate resistances
    pub fn new_with_initial_values(network: &Network) -> Self {
        // get the initial emitter flows
        let initial_emitter_flows = network
            .nodes
            .iter()
            .map(|n| {
                if let NodeType::Junction(junction) = &n.node_type {
                    if junction.emitter_coefficient > 0.0 {
                        1.0
                    } else {
                        0.0
                    }
                } else {
                    0.0
                }
            })
            .collect::<Vec<f64>>();

        Self {
            flows: network
                .links
                .iter()
                .map(|l| l.initial_flow())
                .collect::<Vec<f64>>(),
            heads: network
                .nodes
                .iter()
                .map(|n| n.initial_head())
                .collect::<Vec<f64>>(),
            emitter_flows: initial_emitter_flows,
            demands: vec![0.0; network.nodes.len()],
            demand_flows: vec![0.0; network.nodes.len()],
            settings: network
                .links
                .iter()
                .map(|l| l.initial_setting())
                .collect::<Vec<f64>>(),
            statuses: network
                .links
                .iter()
                .map(|l| l.initial_status)
                .collect::<Vec<LinkStatus>>(),
            resistances: network
                .links
                .iter()
                .map(|l| l.resistance())
                .collect::<Vec<f64>>(),
            topology_version: network.topology_version,
            properties_version: network.properties_version,
        }
    }

    /// Applies demand and head patterns to the state at the given time.
    /// TODO: Add support for multiple patterns
    pub fn apply_patterns(&mut self, network: &Network, time: usize) {
        let time_options = &network.options.time_options;
        let pattern_time = time_options.pattern_start + time;
        let pattern_index = pattern_time / time_options.pattern_timestep;

        for (i, node) in network.nodes.iter().enumerate() {
            let NodeType::Reservoir(reservoir) = &node.node_type else {
                continue;
            };
            let Some(pat_idx) = reservoir.head_pattern_index else {
                continue;
            };
            let pattern = &network.patterns[pat_idx];
            self.heads[i] =
                node.elevation * pattern.multipliers[pattern_index % pattern.multipliers.len()];
        }

        let default_pattern_idx = network
            .options
            .pattern
            .as_ref()
            .or(Some(&Box::from("1")))
            .and_then(|id| network.pattern_map.get(id).copied());

        self.demands = network
            .nodes
            .iter()
            .map(|n| {
                let NodeType::Junction(junction) = &n.node_type else {
                    return 0.0;
                };
                let pat_idx = junction.pattern_index.or(default_pattern_idx);
                let pattern = pat_idx.map(|idx| &network.patterns[idx]);
                let multiplier = match pattern {
                    Some(p) => p.multipliers[pattern_index % p.multipliers.len()],
                    None => 1.0,
                };
                junction.basedemand * multiplier * network.options.demand_multiplier
            })
            .collect::<Vec<Cfs>>();

        if network.options.demand_model == DemandModel::PDA {
            self.demand_flows = self.demands.clone();
        }
    }

    /// Applies controls to the state at the given time.
    pub fn apply_controls(&mut self, network: &Network, time: usize) {
        let clocktime = (time + network.options.time_options.start_clocktime) % (24 * 3600);

        for control in &network.controls {
            // skip controls that are not pressure or level controls
            if matches!(
                control.condition,
                ControlCondition::LowPressure { .. } | ControlCondition::HighPressure { .. }
            ) {
                continue;
            }
            if matches!(
                control.condition,
                ControlCondition::LowLevel { .. } | ControlCondition::HighLevel { .. }
            ) && time == 0
            {
                continue;
            }
            if control.is_active(self, network, time, clocktime) {
                control.activate(self, network);
            }
        }
    }

    /// Updates the tank levels for the given time step.
    pub fn update_tanks(&mut self, network: &Network, timestep: usize) {
        for (tank_index, node) in network.nodes.iter().enumerate() {
            if let NodeType::Tank(tank) = &node.node_type {
                // calculate the delta volume and new head for the tank
                let delta_volume = self.demands[tank_index] * timestep as Ft3;
                let new_head = tank.new_head(delta_volume, self.heads[tank_index]);
                self.heads[tank_index] = new_head;
            }
        }
    }

    /// Updates the state with the changes to the network.
    pub fn update_with_network_changes(&mut self, network: &mut Network) {
        // if the topology has changed, reset the state to the initial values
        if network.topology_version != self.topology_version {
            *self = SolverState::new_with_initial_values(network);
            // clear the updated nodes and links
            network.updated_nodes.clear();
            network.updated_links.clear();
        } else {
            for node_index in network.updated_nodes.iter() {
                let node = &network.nodes[*node_index];

                match &node.node_type {
                    NodeType::Junction(junction) => {
                        // set the emitter flow to 1.0 if the emitter coefficient is greater than 0.0, otherwise set it to 0.0
                        self.emitter_flows[*node_index] = if junction.emitter_coefficient > 0.0 {
                            1.0
                        } else {
                            0.0
                        };
                    }
                    NodeType::Tank(_) | NodeType::Reservoir(_) => {
                        // update the head
                        self.heads[*node_index] = node.initial_head();
                    }
                }
            }

            for link_index in network.updated_links.iter() {
                let link = &network.links[*link_index];
                // update the resistance of the link
                self.resistances[*link_index] = link.resistance();
                // reset the flow of the link to the initial flow
                self.flows[*link_index] = link.initial_flow();
                // reset the status of the link to the initial status
                // TODO: investigate whether this leads to inconsistent results when changing non-status parameters of the link
                self.statuses[*link_index] = link.initial_status;
            }
        }

        // clear the updated nodes and links
        network.updated_nodes.clear();
        network.updated_links.clear();

        self.properties_version = network.properties_version;
    }
}
