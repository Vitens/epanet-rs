#![allow(non_snake_case)]

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

use epanet::ffi::enums;
use epanet::ffi::error_codes::ErrorCode;
use epanet::simulation::Simulation;
use epanet::solver::state::SolverState;

mod native;

fn solved_state(sim: &Simulation) -> Option<&SolverState> {
    if sim.solved { sim.state.as_ref() } else { None }
}

#[pyclass]
#[pyo3(name = "Project")]
struct PyProject {
    sim: Option<Simulation>,
}

#[pymethods]
impl PyProject {
    #[new]
    fn new() -> Self {
        Self { sim: None }
    }

    fn open(&mut self, inp_file: &str) -> PyResult<()> {
        let sim =
            Simulation::from_file(inp_file).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        self.sim = Some(sim);
        Ok(())
    }

    fn close(&mut self) {
        self.sim = None;
    }

    fn saveinpfile(&self, path: &str) -> PyResult<()> {
        let sim = self.sim()?;
        sim.network
            .save_network(path)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    // --- Hydraulic solver ---

    fn openH(&mut self) -> PyResult<()> {
        let sim = self.sim_mut()?;
        sim.initialize_hydraulics()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    #[pyo3(signature = (initflag=0))]
    fn initH(&mut self, initflag: i32) -> PyResult<()> {
        let sim = self.sim_mut()?;
        let _ = initflag;
        if sim.solver.is_none() {
            return Err(PyRuntimeError::new_err("Hydraulic solver not opened"));
        }
        use epanet::solver::state::SolverState;
        sim.state = Some(SolverState::new_with_initial_values(&sim.network));
        Ok(())
    }

    fn runH(&mut self, time: i32) -> PyResult<()> {
        let sim = self.sim_mut()?;
        if sim.solver.is_none() {
            return Err(PyRuntimeError::new_err("Hydraulic solver not opened"));
        }
        if time < 0 {
            return Err(PyRuntimeError::new_err("Illegal numeric value"));
        }
        sim.time = time as usize;
        sim.run_hydraulics()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn nextH(&mut self) -> PyResult<i64> {
        let sim = self.sim_mut()?;
        if sim.solver.is_none() {
            return Err(PyRuntimeError::new_err("Hydraulic solver not opened"));
        }
        let dt = sim.next_hydraulic_timestep();
        Ok(dt as i64)
    }

    #[pyo3(signature = (parallel=false))]
    fn solveH(&mut self, parallel: bool) -> PyResult<PySolverResult> {
        let sim = self.sim_mut()?;
        let result = sim
            .solve_hydraulics(parallel)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(PySolverResult {
            flows: result.flows,
            heads: result.heads,
            demands: result.demands,
        })
    }

    fn closeH(&mut self) -> PyResult<()> {
        let sim = self.sim_mut()?;
        sim.solver = None;
        Ok(())
    }

    // --- Water quality solver (stubs) ---

    fn openQ(&self) -> PyResult<()> {
        Ok(())
    }

    #[pyo3(signature = (_initflag=0))]
    fn initQ(&self, _initflag: i32) -> PyResult<()> {
        Ok(())
    }

    fn closeQ(&self) -> PyResult<()> {
        Ok(())
    }

    // --- Counts ---

    fn getcount(&self, object: i32) -> PyResult<i32> {
        let sim = self.sim()?;
        let count_type = enums::CountType::from_repr(object)
            .ok_or_else(|| PyRuntimeError::new_err("Invalid parameter code"))?;
        use enums::CountType::*;
        use epanet::model::node::NodeType;
        let net = &sim.network;
        let count = match count_type {
            NodeCount => net.nodes.len(),
            TankCount => net
                .nodes
                .iter()
                .filter(|n| matches!(n.node_type, NodeType::Tank(_) | NodeType::Reservoir(_)))
                .count(),
            LinkCount => net.links.len(),
            PatCount => net.patterns.len(),
            CurveCount => net.curves.len(),
            ControlCount => net.controls.len(),
            RuleCount => 0,
        };
        Ok(count as i32)
    }

    // --- Nodes ---

    fn addnode(&mut self, id: &str, node_type: i32) -> PyResult<i32> {
        use enums::NodeType as ENNodeType;
        use epanet::model::demand::Demand;
        use epanet::model::network::modify::{JunctionData, ReservoirData, TankData};

        let sim = self.sim_mut()?;
        let nt = ENNodeType::from_repr(node_type)
            .ok_or_else(|| PyRuntimeError::new_err("Invalid parameter code"))?;

        let result = match nt {
            ENNodeType::Junction => sim.network.add_junction(
                id,
                &JunctionData {
                    elevation: 0.0,
                    demands: vec![Demand {
                        basedemand: 0.0,
                        pattern: None,
                        pattern_index: None,
                        name: None,
                    }],
                    emitter_coefficient: 0.0,
                    coordinates: None,
                },
            ),
            ENNodeType::Reservoir => sim.network.add_reservoir(
                id,
                &ReservoirData {
                    elevation: 0.0,
                    head_pattern: None,
                    coordinates: None,
                },
            ),
            ENNodeType::Tank => sim.network.add_tank(
                id,
                &TankData {
                    elevation: 0.0,
                    initial_level: 0.0,
                    min_level: 0.0,
                    max_level: 0.0,
                    diameter: 0.0,
                    min_volume: 0.0,
                    volume_curve_id: None,
                    overflow: false,
                    coordinates: None,
                },
            ),
        };

        result.map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let index = *sim.network.node_map.get(id).unwrap() + 1;
        Ok(index as i32)
    }

    fn getnodeindex(&self, id: &str) -> PyResult<i32> {
        let sim = self.sim()?;
        let index = sim
            .network
            .node_map
            .get(id)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined node"))?;
        Ok((*index + 1) as i32)
    }

    fn getnodeid(&self, index: i32) -> PyResult<String> {
        let sim = self.sim()?;
        let node = sim
            .network
            .nodes
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined node"))?;
        Ok(node.id.to_string())
    }

    fn setnodeid(&mut self, index: i32, id: &str) -> PyResult<()> {
        let sim = self.sim_mut()?;
        let idx = (index - 1) as usize;

        if sim.network.node_map.contains_key(id) {
            return Err(PyRuntimeError::new_err("Duplicate ID"));
        }

        let node = sim
            .network
            .nodes
            .get_mut(idx)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined node"))?;

        let old_id = node.id.clone();
        sim.network.node_map.remove(&old_id);
        node.id = id.into();
        sim.network.node_map.insert(id.into(), idx);

        for link in sim.network.links.iter_mut() {
            if link.start_node_id == old_id {
                link.start_node_id = id.into();
            }
            if link.end_node_id == old_id {
                link.end_node_id = id.into();
            }
        }
        Ok(())
    }

    fn getnodetype(&self, index: i32) -> PyResult<i32> {
        let sim = self.sim()?;
        let node = sim
            .network
            .nodes
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined node"))?;
        use epanet::model::node::NodeType;
        let t = match &node.node_type {
            NodeType::Junction(_) => 0,
            NodeType::Reservoir(_) => 1,
            NodeType::Tank(_) => 2,
        };
        Ok(t)
    }

    fn getnodevalue(&self, index: i32, property: i32) -> PyResult<f64> {
        let sim = self.sim()?;
        let idx = (index - 1) as usize;
        let node = sim
            .network
            .nodes
            .get(idx)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined node"))?;

        use enums::NodeProperty;
        use epanet::model::node::NodeType;

        let prop = NodeProperty::from_repr(property)
            .ok_or_else(|| PyRuntimeError::new_err("Invalid parameter code"))?;

        let options = &sim.network.options;
        let unit_system = &options.unit_system;
        let flow_units = &options.flow_units;
        let pressure_units = &options.pressure_units;

        let value = match prop {
            NodeProperty::Elevation => node.elevation * unit_system.per_feet(),
            NodeProperty::BaseDemand => match &node.node_type {
                NodeType::Junction(j) => j
                    .demands
                    .first()
                    .map(|d| d.basedemand * flow_units.per_cfs())
                    .unwrap_or(0.0),
                _ => 0.0,
            },
            NodeProperty::Pattern => match &node.node_type {
                NodeType::Junction(j) => j
                    .demands
                    .first()
                    .map(|d| d.pattern_index.map(|i| (i + 1) as f64).unwrap_or(0.0))
                    .unwrap_or(0.0),
                NodeType::Reservoir(r) => {
                    r.head_pattern_index.map(|i| (i + 1) as f64).unwrap_or(0.0)
                }
                _ => 0.0,
            },
            NodeProperty::Emitter => match &node.node_type {
                NodeType::Junction(j) => {
                    let e = j.emitter_coefficient;
                    if e > 0.0 {
                        flow_units.per_cfs()
                            / (pressure_units.per_feet() * e).powf(1.0 / options.emitter_exponent)
                    } else {
                        0.0
                    }
                }
                _ => 0.0,
            },
            NodeProperty::TankLevel => match &node.node_type {
                NodeType::Tank(_) => {
                    solved_state(sim).map_or(0.0, |state| state.heads[idx] - node.elevation)
                        * unit_system.per_feet()
                }
                _ => 0.0,
            },
            NodeProperty::Demand => {
                solved_state(sim).map_or(0.0, |state| state.demands[idx] * flow_units.per_cfs())
            }
            NodeProperty::Head => {
                solved_state(sim).map_or(0.0, |state| state.heads[idx] * unit_system.per_feet())
            }
            NodeProperty::Pressure => solved_state(sim).map_or(0.0, |state| {
                (state.heads[idx] - node.elevation) * pressure_units.per_feet()
            }),
            NodeProperty::TankDiam => match &node.node_type {
                NodeType::Tank(tank) => tank.diameter * unit_system.per_feet(),
                _ => 0.0,
            },
            NodeProperty::MinLevel => match &node.node_type {
                NodeType::Tank(tank) => tank.min_level * unit_system.per_feet(),
                _ => 0.0,
            },
            NodeProperty::MaxLevel => match &node.node_type {
                NodeType::Tank(tank) => tank.max_level * unit_system.per_feet(),
                _ => 0.0,
            },
            NodeProperty::InitVolume => match &node.node_type {
                NodeType::Tank(tank) => {
                    tank.volume_at_level(tank.initial_level) * unit_system.per_cubic_feet()
                }
                _ => 0.0,
            },
            NodeProperty::MinVolume => match &node.node_type {
                NodeType::Tank(tank) => tank.min_volume() * unit_system.per_cubic_feet(),
                _ => 0.0,
            },
            NodeProperty::MaxVolume => match &node.node_type {
                NodeType::Tank(tank) => tank.max_volume() * unit_system.per_cubic_feet(),
                _ => 0.0,
            },
            NodeProperty::TankVolume => match &node.node_type {
                NodeType::Tank(tank) => {
                    if let Some(state) = solved_state(sim) {
                        tank.volume_at_head(state.heads[idx]) * unit_system.per_cubic_feet()
                    } else {
                        tank.volume_at_level(tank.initial_level) * unit_system.per_cubic_feet()
                    }
                }
                _ => 0.0,
            },
            NodeProperty::CanOverflow => match &node.node_type {
                NodeType::Tank(tank) if tank.overflow => 1.0,
                _ => 0.0,
            },
            _ => 0.0,
        };
        Ok(value)
    }

    fn setnodevalue(&mut self, index: i32, property: i32, value: f64) -> PyResult<()> {
        use enums::NodeProperty;
        use epanet::model::network::modify::{JunctionUpdate, NodeUpdate, TankUpdate};

        let sim = self.sim_mut()?;
        let idx = (index - 1) as usize;

        let node_id = sim
            .network
            .nodes
            .get(idx)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined node"))?
            .id
            .clone();

        let prop = NodeProperty::from_repr(property)
            .ok_or_else(|| PyRuntimeError::new_err("Invalid parameter code"))?;

        let result = match prop {
            NodeProperty::Elevation => sim.network.update_node(
                &node_id,
                &NodeUpdate {
                    elevation: Some(value),
                    ..Default::default()
                },
            ),
            NodeProperty::BaseDemand => sim.network.update_junction(
                &node_id,
                &JunctionUpdate {
                    basedemand: Some(value),
                    ..Default::default()
                },
            ),
            NodeProperty::Emitter => sim.network.update_junction(
                &node_id,
                &JunctionUpdate {
                    emitter_coefficient: Some(value),
                    ..Default::default()
                },
            ),
            NodeProperty::TankDiam => sim.network.update_tank(
                &node_id,
                &TankUpdate {
                    diameter: Some(value),
                    ..Default::default()
                },
            ),
            NodeProperty::MinLevel => sim.network.update_tank(
                &node_id,
                &TankUpdate {
                    min_level: Some(value),
                    ..Default::default()
                },
            ),
            NodeProperty::MaxLevel => sim.network.update_tank(
                &node_id,
                &TankUpdate {
                    max_level: Some(value),
                    ..Default::default()
                },
            ),
            NodeProperty::TankLevel => sim.network.update_tank(
                &node_id,
                &TankUpdate {
                    initial_level: Some(value),
                    ..Default::default()
                },
            ),
            NodeProperty::Pattern => {
                let pattern_id = sim
                    .network
                    .patterns
                    .get(value as usize - 1)
                    .ok_or_else(|| PyRuntimeError::new_err("Undefined pattern"))?
                    .id
                    .clone();
                sim.network.update_junction(
                    &node_id,
                    &JunctionUpdate {
                        pattern: Some(Some(pattern_id)),
                        ..Default::default()
                    },
                )
            }
            _ => return Err(PyRuntimeError::new_err("Invalid parameter code")),
        };

        result.map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn getcoord(&self, index: i32) -> PyResult<(f64, f64)> {
        let sim = self.sim()?;
        let node = sim
            .network
            .nodes
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined node"))?;
        node.coordinates
            .ok_or_else(|| PyRuntimeError::new_err("Node has no coordinates"))
    }

    fn setcoord(&mut self, index: i32, x: f64, y: f64) -> PyResult<()> {
        let sim = self.sim_mut()?;
        let node = sim
            .network
            .nodes
            .get_mut((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined node"))?;
        node.coordinates = Some((x, y));
        Ok(())
    }

    #[pyo3(signature = (index, action_code=0))]
    fn deletenode(&mut self, index: i32, action_code: i32) -> PyResult<()> {
        let sim = self.sim_mut()?;
        let node_id = sim
            .network
            .nodes
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined node"))?
            .id
            .clone();
        sim.network
            .remove_node(&node_id, action_code == 0)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    // --- Links ---

    fn addlink(
        &mut self,
        id: &str,
        link_type: i32,
        start_node: &str,
        end_node: &str,
    ) -> PyResult<i32> {
        use enums::LinkType as ENLinkType;
        use epanet::model::link::LinkStatus;
        use epanet::model::network::modify::{PipeData, PumpData, ValveData};
        use epanet::model::options::HeadlossFormula;
        use epanet::model::valve::ValveType;

        let sim = self.sim_mut()?;
        let lt = ENLinkType::from_repr(link_type)
            .ok_or_else(|| PyRuntimeError::new_err("Invalid parameter code"))?;

        let pipe_data = PipeData {
            start_node: start_node.into(),
            end_node: end_node.into(),
            length: 330.0,
            diameter: 10.0 / 12.0 * sim.network.options.unit_system.per_feet(),
            roughness: match sim.network.options.headloss_formula {
                HeadlossFormula::DarcyWeisbach => 0.0015,
                HeadlossFormula::HazenWilliams => 130.0,
                HeadlossFormula::ChezyManning => 0.01,
            },
            minor_loss: 0.0,
            check_valve: false,
            vertices: None,
            initial_status: LinkStatus::Open,
        };

        let valve_data = ValveData {
            start_node: start_node.into(),
            end_node: end_node.into(),
            diameter: 10.0 / 12.0 * sim.network.options.unit_system.per_feet(),
            valve_type: ValveType::PRV,
            setting: 0.0,
            minor_loss: 0.0,
            initial_status: LinkStatus::Active,
            vertices: None,
            curve_id: None,
        };

        let result = match lt {
            ENLinkType::CVPipe => {
                let mut d = pipe_data;
                d.check_valve = true;
                sim.network.add_pipe(id, &d)
            }
            ENLinkType::Pipe => sim.network.add_pipe(id, &pipe_data),
            ENLinkType::Pump => {
                let pump = PumpData {
                    start_node: start_node.into(),
                    end_node: end_node.into(),
                    speed: 1.0,
                    head_curve_id: None,
                    power: 0.0,
                    initial_status: LinkStatus::Open,
                    vertices: None,
                };
                sim.network.add_pump(id, &pump)
            }
            ENLinkType::PRV => {
                let mut v = valve_data;
                v.valve_type = ValveType::PRV;
                sim.network.add_valve(id, &v)
            }
            ENLinkType::PSV => {
                let mut v = valve_data;
                v.valve_type = ValveType::PSV;
                sim.network.add_valve(id, &v)
            }
            ENLinkType::PBV => {
                let mut v = valve_data;
                v.valve_type = ValveType::PBV;
                sim.network.add_valve(id, &v)
            }
            ENLinkType::FCV => {
                let mut v = valve_data;
                v.valve_type = ValveType::FCV;
                sim.network.add_valve(id, &v)
            }
            ENLinkType::TCV => {
                let mut v = valve_data;
                v.valve_type = ValveType::TCV;
                sim.network.add_valve(id, &v)
            }
            ENLinkType::GPV => {
                let mut v = valve_data;
                v.valve_type = ValveType::GPV;
                sim.network.add_valve(id, &v)
            }
            ENLinkType::PCV => {
                let mut v = valve_data;
                v.valve_type = ValveType::PCV;
                sim.network.add_valve(id, &v)
            }
        };

        result.map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let index = *sim.network.link_map.get(id).unwrap() + 1;
        Ok(index as i32)
    }

    #[pyo3(signature = (index, action_code=0))]
    fn deletelink(&mut self, index: i32, action_code: i32) -> PyResult<()> {
        let sim = self.sim_mut()?;
        let link_id = sim
            .network
            .links
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined link"))?
            .id
            .clone();
        let unconditional = action_code == 0;
        sim.network
            .remove_link(&link_id, unconditional)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn getlinkindex(&self, id: &str) -> PyResult<i32> {
        let sim = self.sim()?;
        let index = sim
            .network
            .link_map
            .get(id)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined link"))?;
        Ok((*index + 1) as i32)
    }

    fn getlinkid(&self, index: i32) -> PyResult<String> {
        let sim = self.sim()?;
        let link = sim
            .network
            .links
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined link"))?;
        Ok(link.id.to_string())
    }

    fn setlinkid(&mut self, index: i32, id: &str) -> PyResult<()> {
        let sim = self.sim_mut()?;
        let idx = (index - 1) as usize;

        if sim.network.link_map.contains_key(id) {
            return Err(PyRuntimeError::new_err("Duplicate ID"));
        }

        let link = sim
            .network
            .links
            .get_mut(idx)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined link"))?;

        sim.network.link_map.remove(&link.id);
        link.id = id.into();
        sim.network.link_map.insert(id.into(), idx);
        Ok(())
    }

    fn getlinktype(&self, index: i32) -> PyResult<i32> {
        use epanet::model::link::LinkType;
        use epanet::model::valve::ValveType;

        let sim = self.sim()?;
        let link = sim
            .network
            .links
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined link"))?;

        let t = match &link.link_type {
            LinkType::Pipe(pipe) => {
                if pipe.check_valve {
                    0
                } else {
                    1
                }
            }
            LinkType::Pump(_) => 2,
            LinkType::Valve(valve) => match valve.valve_type {
                ValveType::PRV => 3,
                ValveType::PSV => 4,
                ValveType::PBV => 5,
                ValveType::FCV => 6,
                ValveType::TCV => 7,
                ValveType::GPV => 8,
                ValveType::PCV => 9,
            },
        };
        Ok(t)
    }

    fn getlinknodes(&self, index: i32) -> PyResult<(i32, i32)> {
        let sim = self.sim()?;
        let link = sim
            .network
            .links
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined link"))?;
        Ok(((link.start_node + 1) as i32, (link.end_node + 1) as i32))
    }

    fn setlinknodes(&mut self, index: i32, start_node: i32, end_node: i32) -> PyResult<()> {
        use epanet::model::network::modify::LinkUpdate;

        let sim = self.sim_mut()?;

        if start_node == end_node {
            return Err(PyRuntimeError::new_err(
                "Link assigned same start and end nodes",
            ));
        }

        let link_id = sim
            .network
            .links
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined link"))?
            .id
            .clone();

        let start_node_id = sim
            .network
            .nodes
            .get((start_node - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined node"))?
            .id
            .clone();

        let end_node_id = sim
            .network
            .nodes
            .get((end_node - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined node"))?
            .id
            .clone();

        sim.network
            .update_link(
                &link_id,
                &LinkUpdate {
                    start_node: Some(start_node_id),
                    end_node: Some(end_node_id),
                    ..Default::default()
                },
            )
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn getlinkvalue(&self, index: i32, property: i32) -> PyResult<f64> {
        use enums::{LinkProperty, MISSING_VALUE};
        use epanet::constants::MperFT;
        use epanet::model::link::{LinkStatus, LinkType};
        use epanet::model::options::HeadlossFormula;
        use epanet::model::units::UnitSystem;

        let sim = self.sim()?;
        let idx = (index - 1) as usize;
        let link = sim
            .network
            .links
            .get(idx)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined link"))?;

        let prop = LinkProperty::from_repr(property)
            .ok_or_else(|| PyRuntimeError::new_err("Invalid parameter code"))?;
        let options = &sim.network.options;

        let value = match prop {
            LinkProperty::Diameter => {
                let conversion = match options.unit_system {
                    UnitSystem::US => 12.0,
                    UnitSystem::SI => MperFT * 1e3,
                };
                match &link.link_type {
                    LinkType::Pipe(pipe) => pipe.diameter * conversion,
                    LinkType::Valve(valve) => valve.diameter * conversion,
                    LinkType::Pump(_) => 0.0,
                }
            }
            LinkProperty::Length => match &link.link_type {
                LinkType::Pipe(pipe) => pipe.length * options.unit_system.per_feet(),
                _ => 0.0,
            },
            LinkProperty::Roughness => match &link.link_type {
                LinkType::Pipe(pipe) => {
                    if pipe.headloss_formula == HeadlossFormula::DarcyWeisbach {
                        pipe.roughness * options.unit_system.per_feet()
                    } else {
                        pipe.roughness
                    }
                }
                _ => 0.0,
            },
            LinkProperty::MinorLoss => match &link.link_type {
                LinkType::Pipe(pipe) => pipe.minor_loss * pipe.diameter.powi(4) / 0.02517,
                LinkType::Valve(valve) => valve.minor_loss * valve.diameter.powi(4) / 0.02517,
                _ => 0.0,
            },
            LinkProperty::HeadLoss => {
                if let Some(state) = solved_state(sim) {
                    let headloss = state.heads[link.start_node] - state.heads[link.end_node];
                    match &link.link_type {
                        LinkType::Pipe(_) | LinkType::Valve(_) => {
                            headloss.abs() * options.unit_system.per_feet()
                        }
                        LinkType::Pump(_) => headloss,
                    }
                } else {
                    0.0
                }
            }
            LinkProperty::Status => {
                if let Some(state) = solved_state(sim) {
                    match state.statuses[idx] {
                        LinkStatus::Xhead => 2.0,
                        LinkStatus::Closed | LinkStatus::FixedClosed => 0.0,
                        _ => 1.0,
                    }
                } else {
                    1.0
                }
            }
            LinkProperty::InitStatus => match &link.initial_status {
                LinkStatus::Closed => 0.0,
                _ => 1.0,
            },
            LinkProperty::Setting | LinkProperty::InitSetting => match &link.link_type {
                LinkType::Pipe(_) => 0.0,
                LinkType::Valve(valve) => match valve.valve_type {
                    epanet::model::valve::ValveType::FCV => {
                        valve.setting * options.flow_units.per_cfs()
                    }
                    epanet::model::valve::ValveType::TCV | epanet::model::valve::ValveType::PCV => {
                        valve.setting
                    }
                    _ => valve.setting * options.pressure_units.per_feet(),
                },
                LinkType::Pump(pump) => pump.speed,
            },
            LinkProperty::Flow => solved_state(sim).map_or(MISSING_VALUE, |state| {
                state.flows[idx] * options.flow_units.per_cfs()
            }),
            LinkProperty::Velocity => solved_state(sim).map_or(MISSING_VALUE, |state| {
                let flow = state.flows[idx].abs();
                match &link.link_type {
                    LinkType::Pipe(pipe) => {
                        flow / (pipe.diameter.powi(2) * std::f64::consts::PI / 4.0)
                            * options.unit_system.per_feet()
                    }
                    LinkType::Valve(valve) => {
                        flow / (valve.diameter.powi(2) * std::f64::consts::PI / 4.0)
                            * options.unit_system.per_feet()
                    }
                    LinkType::Pump(_) => 0.0,
                }
            }),
            _ => 0.0,
        };
        Ok(value)
    }

    fn setlinkvalue(&mut self, index: i32, property: i32, value: f64) -> PyResult<()> {
        use enums::LinkProperty;
        use epanet::model::link::{LinkStatus, LinkType};
        use epanet::model::network::modify::{LinkUpdate, PipeUpdate, PumpUpdate, ValveUpdate};

        let sim = self.sim_mut()?;
        let idx = (index - 1) as usize;

        let link = sim
            .network
            .links
            .get(idx)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined link"))?;
        let link_id = link.id.clone();
        let link_type_clone = link.link_type.clone();

        let prop = LinkProperty::from_repr(property)
            .ok_or_else(|| PyRuntimeError::new_err("Invalid parameter code"))?;

        let result = match prop {
            LinkProperty::Diameter => match &link_type_clone {
                LinkType::Pipe(_) => sim.network.update_pipe(
                    &link_id,
                    &PipeUpdate {
                        diameter: Some(value),
                        ..Default::default()
                    },
                ),
                LinkType::Valve(_) => sim.network.update_valve(
                    &link_id,
                    &ValveUpdate {
                        diameter: Some(value),
                        ..Default::default()
                    },
                ),
                _ => return Err(PyRuntimeError::new_err("Invalid parameter code")),
            },
            LinkProperty::Length => sim.network.update_pipe(
                &link_id,
                &PipeUpdate {
                    length: Some(value),
                    ..Default::default()
                },
            ),
            LinkProperty::Roughness => sim.network.update_pipe(
                &link_id,
                &PipeUpdate {
                    roughness: Some(value),
                    ..Default::default()
                },
            ),
            LinkProperty::MinorLoss => match &link_type_clone {
                LinkType::Pipe(_) => sim.network.update_pipe(
                    &link_id,
                    &PipeUpdate {
                        minor_loss: Some(value),
                        ..Default::default()
                    },
                ),
                LinkType::Valve(_) => sim.network.update_valve(
                    &link_id,
                    &ValveUpdate {
                        minor_loss: Some(value),
                        ..Default::default()
                    },
                ),
                _ => return Err(PyRuntimeError::new_err("Invalid parameter code")),
            },
            LinkProperty::InitStatus | LinkProperty::Status => {
                let status = match value as i32 {
                    0 => LinkStatus::Closed,
                    1 => LinkStatus::Open,
                    _ => return Err(PyRuntimeError::new_err("Invalid parameter code")),
                };
                sim.network.update_link(
                    &link_id,
                    &LinkUpdate {
                        initial_status: Some(status),
                        ..Default::default()
                    },
                )
            }
            LinkProperty::InitSetting | LinkProperty::Setting => match &link_type_clone {
                LinkType::Pipe(_) => sim.network.update_pipe(
                    &link_id,
                    &PipeUpdate {
                        roughness: Some(value),
                        ..Default::default()
                    },
                ),
                LinkType::Valve(_) => sim.network.update_valve(
                    &link_id,
                    &ValveUpdate {
                        setting: Some(value),
                        ..Default::default()
                    },
                ),
                LinkType::Pump(_) => sim.network.update_pump(
                    &link_id,
                    &PumpUpdate {
                        speed: Some(value),
                        ..Default::default()
                    },
                ),
            },
            _ => return Err(PyRuntimeError::new_err("Invalid parameter code")),
        };

        result.map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn getheadcurveindex(&self, index: i32) -> PyResult<i32> {
        use epanet::model::link::LinkType;

        let sim = self.sim()?;
        let link = sim
            .network
            .links
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined link"))?;

        match &link.link_type {
            LinkType::Pump(pump) => {
                if let Some(ref curve_id) = pump.head_curve_id {
                    if let Some(&ci) = sim.network.curve_map.get(curve_id) {
                        return Ok((ci + 1) as i32);
                    }
                }
                Ok(0)
            }
            _ => Err(PyRuntimeError::new_err("Undefined pump")),
        }
    }

    fn setheadcurveindex(&mut self, index: i32, head_curve_index: i32) -> PyResult<()> {
        use epanet::model::link::LinkType;
        use epanet::model::network::modify::PumpUpdate;

        let sim = self.sim_mut()?;
        let link = sim
            .network
            .links
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined link"))?;
        let link_id = link.id.clone();

        if !matches!(&link.link_type, LinkType::Pump(_)) {
            return Err(PyRuntimeError::new_err("Undefined pump"));
        }

        let curve_id = if head_curve_index == 0 {
            None
        } else {
            let c = sim
                .network
                .curves
                .get((head_curve_index - 1) as usize)
                .ok_or_else(|| PyRuntimeError::new_err("Undefined curve"))?;
            Some(c.id.clone())
        };

        sim.network
            .update_pump(
                &link_id,
                &PumpUpdate {
                    head_curve_id: Some(curve_id),
                    ..Default::default()
                },
            )
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    // --- Patterns ---

    fn addpattern(&mut self, id: &str) -> PyResult<()> {
        use epanet::model::network::modify::PatternData;

        let sim = self.sim_mut()?;
        sim.network
            .add_pattern(
                id,
                &PatternData {
                    multipliers: vec![1.0],
                },
            )
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn deletepattern(&mut self, index: i32) -> PyResult<()> {
        let sim = self.sim_mut()?;
        let pattern_id = sim
            .network
            .patterns
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined pattern"))?
            .id
            .clone();
        sim.network
            .remove_pattern(&pattern_id)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn getpatternindex(&self, id: &str) -> PyResult<i32> {
        let sim = self.sim()?;
        let index = sim
            .network
            .pattern_map
            .get(id)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined pattern"))?;
        Ok((*index + 1) as i32)
    }

    fn getpatternid(&self, index: i32) -> PyResult<String> {
        let sim = self.sim()?;
        let pattern = sim
            .network
            .patterns
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined pattern"))?;
        Ok(pattern.id.to_string())
    }

    fn setpatternid(&mut self, index: i32, id: &str) -> PyResult<()> {
        use epanet::model::node::NodeType;

        let sim = self.sim_mut()?;
        let idx = (index - 1) as usize;

        if sim.network.pattern_map.contains_key(id) {
            return Err(PyRuntimeError::new_err("Duplicate ID"));
        }

        let pattern = sim
            .network
            .patterns
            .get_mut(idx)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined pattern"))?;

        let old_id = pattern.id.clone();
        sim.network.pattern_map.remove(&pattern.id);
        pattern.id = id.into();
        sim.network.pattern_map.insert(id.into(), idx);

        for node in sim.network.nodes.iter_mut() {
            match &mut node.node_type {
                NodeType::Junction(junction) => {
                    for demand in junction.demands.iter_mut() {
                        if demand.pattern.as_deref() == Some(&old_id) {
                            demand.pattern = Some(id.into());
                        }
                    }
                }
                NodeType::Reservoir(reservoir)
                    if reservoir.head_pattern.as_deref() == Some(&old_id) =>
                {
                    reservoir.head_pattern = Some(id.into());
                }
                _ => {}
            }
        }
        Ok(())
    }

    fn getpatternlen(&self, index: i32) -> PyResult<i32> {
        let sim = self.sim()?;
        let pattern = sim
            .network
            .patterns
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined pattern"))?;
        Ok(pattern.multipliers.len() as i32)
    }

    fn getpatternvalue(&self, index: i32, time: i32) -> PyResult<f64> {
        let sim = self.sim()?;
        let pattern = sim
            .network
            .patterns
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined pattern"))?;
        let t = (time - 1) as usize;
        Ok(pattern.multipliers[t % pattern.multipliers.len()])
    }

    fn setpattern(&mut self, index: i32, multipliers: Vec<f64>) -> PyResult<()> {
        let sim = self.sim_mut()?;
        let pattern = sim
            .network
            .patterns
            .get_mut((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined pattern"))?;
        pattern.multipliers = multipliers;
        Ok(())
    }

    fn getaveragepatternvalue(&self, index: i32) -> PyResult<f64> {
        let sim = self.sim()?;
        let pattern = sim
            .network
            .patterns
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined pattern"))?;
        let avg = pattern.multipliers.iter().sum::<f64>() / pattern.multipliers.len() as f64;
        Ok(avg)
    }

    // --- Curves ---

    fn addcurve(&mut self, id: &str) -> PyResult<()> {
        use epanet::model::network::modify::CurveData;

        let sim = self.sim_mut()?;
        sim.network
            .add_curve(
                id,
                &CurveData {
                    x: vec![1.0],
                    y: vec![1.0],
                },
            )
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn getcurveindex(&self, id: &str) -> PyResult<i32> {
        let sim = self.sim()?;
        let index = sim
            .network
            .curve_map
            .get(id)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined curve"))?;
        Ok((*index + 1) as i32)
    }

    fn getcurveid(&self, index: i32) -> PyResult<String> {
        let sim = self.sim()?;
        let curve = sim
            .network
            .curves
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined curve"))?;
        Ok(curve.id.to_string())
    }

    fn setcurveid(&mut self, index: i32, id: &str) -> PyResult<()> {
        let sim = self.sim_mut()?;
        let idx = (index - 1) as usize;

        if sim.network.curve_map.contains_key(id) {
            return Err(PyRuntimeError::new_err("Duplicate ID"));
        }

        let curve = sim
            .network
            .curves
            .get_mut(idx)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined curve"))?;

        sim.network.curve_map.remove(&curve.id);
        curve.id = id.into();
        sim.network.curve_map.insert(id.into(), idx);
        Ok(())
    }

    fn getcurvelen(&self, index: i32) -> PyResult<i32> {
        let sim = self.sim()?;
        let curve = sim
            .network
            .curves
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined curve"))?;
        Ok(curve.x.len() as i32)
    }

    fn getcurvevalue(&self, index: i32, point_index: i32) -> PyResult<(f64, f64)> {
        let sim = self.sim()?;
        let curve = sim
            .network
            .curves
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined curve"))?;
        let pi = (point_index - 1) as usize;
        if pi >= curve.x.len() {
            return Err(PyRuntimeError::new_err("Undefined curve"));
        }
        Ok((curve.x[pi], curve.y[pi]))
    }

    fn getcurve(&self, index: i32) -> PyResult<(Vec<f64>, Vec<f64>)> {
        let sim = self.sim()?;
        let curve = sim
            .network
            .curves
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined curve"))?;
        Ok((curve.x.clone(), curve.y.clone()))
    }

    fn setcurve(&mut self, index: i32, x: Vec<f64>, y: Vec<f64>) -> PyResult<()> {
        use epanet::model::network::modify::CurveUpdate;

        let sim = self.sim_mut()?;
        let curve_id = sim
            .network
            .curves
            .get((index - 1) as usize)
            .ok_or_else(|| PyRuntimeError::new_err("Undefined curve"))?
            .id
            .clone();

        sim.network
            .update_curve(
                &curve_id,
                &CurveUpdate {
                    x: Some(x),
                    y: Some(y),
                },
            )
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    // --- Time parameters ---

    fn settimeparam(&mut self, param: i32, value: i64) -> PyResult<()> {
        use enums::TimeParameter;

        let sim = self.sim_mut()?;
        let tp = TimeParameter::from_repr(param)
            .ok_or_else(|| PyRuntimeError::new_err("Invalid parameter code"))?;

        let v = value as usize;
        let time_options = &mut sim.network.options.time_options;

        match tp {
            TimeParameter::Duration => {
                time_options.duration = v;
                if time_options.duration > time_options.pattern_start {
                    time_options.pattern_start = time_options.duration;
                }
            }
            TimeParameter::HydStep => {
                if v == 0 {
                    return Err(PyRuntimeError::new_err("Invalid option value"));
                }
                time_options.hydraulic_timestep = v;
                time_options.hydraulic_timestep = time_options
                    .hydraulic_timestep
                    .min(time_options.pattern_timestep);
                time_options.hydraulic_timestep = time_options
                    .hydraulic_timestep
                    .min(time_options.report_timestep);
            }
            TimeParameter::PatternStep => {
                if v == 0 {
                    return Err(PyRuntimeError::new_err("Invalid option value"));
                }
                time_options.pattern_timestep = v;
                if time_options.hydraulic_timestep > time_options.pattern_timestep {
                    time_options.hydraulic_timestep = time_options.pattern_timestep;
                }
            }
            TimeParameter::PatternStart => {
                time_options.pattern_start = v;
            }
            TimeParameter::ReportStep => {
                if v == 0 {
                    return Err(PyRuntimeError::new_err("Invalid option value"));
                }
                time_options.report_timestep = v;
                if time_options.hydraulic_timestep > time_options.report_timestep {
                    time_options.hydraulic_timestep = time_options.report_timestep;
                }
            }
            _ => return Err(PyRuntimeError::new_err("Invalid parameter code")),
        }
        Ok(())
    }

    // --- Simulation options ---

    fn getoption(&self, option: i32) -> PyResult<f64> {
        use enums::SimOption;

        let sim = self.sim()?;
        let opt = SimOption::from_repr(option)
            .ok_or_else(|| PyRuntimeError::new_err("Invalid parameter code"))?;

        let v = match opt {
            SimOption::Trials => sim.network.options.max_trials as f64,
            SimOption::Accuracy => sim.network.options.accuracy,
            SimOption::EmitExpon => {
                if sim.network.options.emitter_exponent > 0.0 {
                    1.0 / sim.network.options.emitter_exponent
                } else {
                    0.0
                }
            }
            SimOption::DemandMult => sim.network.options.demand_multiplier,
            SimOption::FlowChange => sim.network.options.max_flow_change.unwrap_or(0.0),
            SimOption::SpGravity => sim.network.options.specific_gravity,
            SimOption::SpViscos => sim.network.options.viscosity,
            SimOption::CheckFreq => sim.network.options.check_frequency as f64,
            SimOption::MaxCheck => sim.network.options.max_check as f64,
            _ => return Err(PyRuntimeError::new_err("Invalid parameter code")),
        };
        Ok(v)
    }

    fn setoption(&mut self, option: i32, value: f64) -> PyResult<()> {
        use enums::SimOption;

        let sim = self.sim_mut()?;
        let opt = SimOption::from_repr(option)
            .ok_or_else(|| PyRuntimeError::new_err("Invalid parameter code"))?;

        if !value.is_finite() {
            return Err(PyRuntimeError::new_err("Illegal numeric value"));
        }
        if value < 0.0 {
            return Err(PyRuntimeError::new_err("Invalid option value"));
        }

        match opt {
            SimOption::Trials => {
                if value < 1.0 {
                    return Err(PyRuntimeError::new_err("Invalid option value"));
                }
                sim.network.options.max_trials = value as usize;
            }
            SimOption::Accuracy => {
                if !(1.0e-8..=1.0e-1).contains(&value) {
                    return Err(PyRuntimeError::new_err("Invalid option value"));
                }
                sim.network.options.accuracy = value;
            }
            SimOption::EmitExpon => {
                if value <= 0.0 {
                    return Err(PyRuntimeError::new_err("Invalid option value"));
                }
                sim.network.options.emitter_exponent = 1.0 / value;
            }
            SimOption::DemandMult => {
                sim.network.options.demand_multiplier = value;
            }
            SimOption::FlowChange => {
                sim.network.options.max_flow_change = Some(value);
            }
            SimOption::SpGravity => {
                if value <= 0.0 {
                    return Err(PyRuntimeError::new_err("Invalid option value"));
                }
                sim.network.options.specific_gravity = value;
            }
            SimOption::SpViscos => {
                if value <= 0.0 {
                    return Err(PyRuntimeError::new_err("Invalid option value"));
                }
                sim.network.options.viscosity = value;
            }
            SimOption::CheckFreq => {
                sim.network.options.check_frequency = value as usize;
            }
            SimOption::MaxCheck => {
                sim.network.options.max_check = value as usize;
            }
            _ => return Err(PyRuntimeError::new_err("Invalid parameter code")),
        }
        Ok(())
    }

    fn setdemandmodel(
        &mut self,
        demand_model: i32,
        minimum_pressure: f64,
        required_pressure: f64,
        pressure_exponent: f64,
    ) -> PyResult<()> {
        use enums::DemandModel as ENDemandModel;
        use epanet::model::options::DemandModel;

        let sim = self.sim_mut()?;
        let dm = ENDemandModel::from_repr(demand_model)
            .ok_or_else(|| PyRuntimeError::new_err("Invalid parameter code"))?;

        if !minimum_pressure.is_finite()
            || !required_pressure.is_finite()
            || !pressure_exponent.is_finite()
        {
            return Err(PyRuntimeError::new_err("Illegal numeric value"));
        }

        if minimum_pressure < 0.0
            || required_pressure <= minimum_pressure
            || pressure_exponent <= 0.0
        {
            return Err(PyRuntimeError::new_err("Illegal PDA pressure limits"));
        }

        match dm {
            ENDemandModel::Pda => sim.network.options.demand_model = DemandModel::PDA,
            ENDemandModel::Dda => sim.network.options.demand_model = DemandModel::DDA,
        }
        sim.network.options.minimum_pressure = minimum_pressure;
        sim.network.options.required_pressure = required_pressure;
        sim.network.options.pressure_exponent = pressure_exponent;
        Ok(())
    }

    // --- Error ---

    #[staticmethod]
    fn geterror(errcode: i32) -> String {
        match ErrorCode::from_repr(errcode) {
            Some(code) => code.to_string(),
            None => format!("Unknown error code: {errcode}"),
        }
    }
}

impl PyProject {
    fn sim(&self) -> PyResult<&Simulation> {
        self.sim
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("No network data available"))
    }

    fn sim_mut(&mut self) -> PyResult<&mut Simulation> {
        self.sim
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("No network data available"))
    }
}

#[pyclass]
#[pyo3(name = "SolverResult")]
pub(crate) struct PySolverResult {
    #[pyo3(get)]
    pub(crate) flows: Vec<Vec<f64>>,
    #[pyo3(get)]
    pub(crate) heads: Vec<Vec<f64>>,
    #[pyo3(get)]
    pub(crate) demands: Vec<Vec<f64>>,
}

// --- Constants module ---

fn add_constants(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Node types
    m.add("EN_JUNCTION", 0)?;
    m.add("EN_RESERVOIR", 1)?;
    m.add("EN_TANK", 2)?;

    // Link types
    m.add("EN_CVPIPE", 0)?;
    m.add("EN_PIPE", 1)?;
    m.add("EN_PUMP", 2)?;
    m.add("EN_PRV", 3)?;
    m.add("EN_PSV", 4)?;
    m.add("EN_PBV", 5)?;
    m.add("EN_FCV", 6)?;
    m.add("EN_TCV", 7)?;
    m.add("EN_GPV", 8)?;
    m.add("EN_PCV", 9)?;

    // Count types
    m.add("EN_NODECOUNT", 0)?;
    m.add("EN_TANKCOUNT", 1)?;
    m.add("EN_LINKCOUNT", 2)?;
    m.add("EN_PATCOUNT", 3)?;
    m.add("EN_CURVECOUNT", 4)?;
    m.add("EN_CONTROLCOUNT", 5)?;
    m.add("EN_RULECOUNT", 6)?;

    // Node properties
    m.add("EN_ELEVATION", 0)?;
    m.add("EN_BASEDEMAND", 1)?;
    m.add("EN_PATTERN", 2)?;
    m.add("EN_EMITTER", 3)?;
    m.add("EN_INITQUAL", 4)?;
    m.add("EN_SOURCEQUAL", 5)?;
    m.add("EN_SOURCEPAT", 6)?;
    m.add("EN_SOURCETYPE", 7)?;
    m.add("EN_TANKLEVEL", 8)?;
    m.add("EN_DEMAND", 9)?;
    m.add("EN_HEAD", 10)?;
    m.add("EN_PRESSURE", 11)?;
    m.add("EN_QUALITY", 12)?;
    m.add("EN_SOURCEMASS", 13)?;
    m.add("EN_INITVOLUME", 14)?;
    m.add("EN_MIXMODEL", 15)?;
    m.add("EN_MIXZONEVOL", 16)?;
    m.add("EN_TANKDIAM", 17)?;
    m.add("EN_MINVOLUME", 18)?;
    m.add("EN_VOLCURVE", 19)?;
    m.add("EN_MINLEVEL", 20)?;
    m.add("EN_MAXLEVEL", 21)?;
    m.add("EN_MIXFRACTION", 22)?;
    m.add("EN_TANK_KBULK", 23)?;
    m.add("EN_TANKVOLUME", 24)?;
    m.add("EN_MAXVOLUME", 25)?;
    m.add("EN_CANOVERFLOW", 26)?;
    m.add("EN_DEMANDDEFICIT", 27)?;

    // Link properties
    m.add("EN_DIAMETER", 0)?;
    m.add("EN_LENGTH", 1)?;
    m.add("EN_ROUGHNESS", 2)?;
    m.add("EN_MINORLOSS", 3)?;
    m.add("EN_INITSTATUS", 4)?;
    m.add("EN_INITSETTING", 5)?;
    m.add("EN_KBULK", 6)?;
    m.add("EN_KWALL", 7)?;
    m.add("EN_FLOW", 8)?;
    m.add("EN_VELOCITY", 9)?;
    m.add("EN_HEADLOSS", 10)?;
    m.add("EN_STATUS", 11)?;
    m.add("EN_SETTING", 12)?;
    m.add("EN_ENERGY", 13)?;
    m.add("EN_LINKQUAL", 14)?;
    m.add("EN_LINKPATTERN", 15)?;
    m.add("EN_PUMPSTATE", 16)?;
    m.add("EN_PUMPEFFIC", 17)?;
    m.add("EN_PUMPPOWER", 18)?;
    m.add("EN_PUMPHCURVE", 19)?;
    m.add("EN_PUMPECURVE", 20)?;
    m.add("EN_PUMPECOST", 21)?;
    m.add("EN_PUMPEPAT", 22)?;

    // Time parameters
    m.add("EN_DURATION", 0)?;
    m.add("EN_HYDSTEP", 1)?;
    m.add("EN_QUALSTEP", 2)?;
    m.add("EN_PATTERNSTEP", 3)?;
    m.add("EN_PATTERNSTART", 4)?;
    m.add("EN_REPORTSTEP", 5)?;
    m.add("EN_REPORTSTART", 6)?;
    m.add("EN_RULESTEP", 7)?;
    m.add("EN_STATISTIC", 8)?;
    m.add("EN_PERIODS", 9)?;
    m.add("EN_STARTTIME", 10)?;
    m.add("EN_HTIME", 11)?;
    m.add("EN_QTIME", 12)?;
    m.add("EN_HALTFLAG", 13)?;
    m.add("EN_NEXTEVENT", 14)?;
    m.add("EN_NEXTEVENTTANK", 15)?;

    // Simulation options
    m.add("EN_TRIALS", 0)?;
    m.add("EN_ACCURACY", 1)?;
    m.add("EN_TOLERANCE", 2)?;
    m.add("EN_EMITEXPON", 3)?;
    m.add("EN_DEMANDMULT", 4)?;
    m.add("EN_HEADERROR", 5)?;
    m.add("EN_FLOWCHANGE", 6)?;
    m.add("EN_HEADLOSSFORM", 7)?;
    m.add("EN_GLOBALEFFIC", 8)?;
    m.add("EN_GLOBALPRICE", 9)?;
    m.add("EN_GLOBALPATTERN", 10)?;
    m.add("EN_DEMANDCHARGE", 11)?;
    m.add("EN_SP_GRAVITY", 12)?;
    m.add("EN_SP_VISCOS", 13)?;
    m.add("EN_UNBALANCED", 14)?;
    m.add("EN_CHECKFREQ", 15)?;
    m.add("EN_MAXCHECK", 16)?;
    m.add("EN_DAMPLIMIT", 17)?;
    m.add("EN_SP_DIFFUS", 18)?;
    m.add("EN_BULKORDER", 19)?;
    m.add("EN_WALLORDER", 20)?;
    m.add("EN_TANKORDER", 21)?;
    m.add("EN_CONCENLIMIT", 22)?;

    // Demand models
    m.add("EN_DDA", 0)?;
    m.add("EN_PDA", 1)?;

    // Head loss types
    m.add("EN_HW", 0)?;
    m.add("EN_DW", 1)?;
    m.add("EN_CM", 2)?;

    // Flow units
    m.add("EN_CFS", 0)?;
    m.add("EN_GPM", 1)?;
    m.add("EN_MGD", 2)?;
    m.add("EN_IMGD", 3)?;
    m.add("EN_AFD", 4)?;
    m.add("EN_LPS", 5)?;
    m.add("EN_LPM", 6)?;
    m.add("EN_MLD", 7)?;
    m.add("EN_CMH", 8)?;
    m.add("EN_CMD", 9)?;

    // Missing value
    m.add("EN_MISSING", enums::MISSING_VALUE)?;

    Ok(())
}

#[pymodule]
fn epanet_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyProject>()?;
    m.add_class::<PySolverResult>()?;
    add_constants(m)?;
    native::register(m)?;
    Ok(())
}
