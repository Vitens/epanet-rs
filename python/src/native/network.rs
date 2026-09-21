//! `Network`: a Pythonic wrapper around `epanet_rs::model::network::Network`'s
//! `add_*`/`update_*`/`remove_*` methods (see `Network::modify` in the Rust
//! crate) plus read accessors returning the snapshot types in
//! [`crate::native::types`].
//!
//! A `Network` instance does not own its data: it holds a reference to the
//! owning [`crate::native::simulation::Simulation`] and borrows it (through
//! the GIL) for the duration of each call, mirroring how `sim.network` is
//! simply a field access in Rust.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

use epanet::error::InputError;
use epanet::model::demand::Demand as RsDemand;
use epanet::model::network::{
    CurveData, CurveUpdate, JunctionData, JunctionUpdate, LinkUpdate, NodeUpdate, PatternData,
    PatternUpdate, PipeData, PipeUpdate, PumpData, PumpUpdate, ReservoirData, ReservoirUpdate,
    TankData, TankUpdate, ValveData, ValveUpdate,
};

use super::enums::{LinkStatus, ValveType};
use super::input_err;
use super::maps::{CurveMapView, LinkMapView, NodeMapView, PatternMapView};
use super::simulation::Simulation as PySimulation;
use super::types::{self, Curve, Link, Node, Pattern};
use super::views::{ControlView, CurveView, LinkView, NodeView, PatternView};

#[pyclass(module = "epanet_rs.native", name = "Network")]
pub struct Network {
    pub(crate) parent: Py<PySimulation>,
}

#[pymethods]
impl Network {
    // --- junctions / tanks / reservoirs ---

    #[pyo3(signature = (id, elevation, basedemand=0.0, pattern=None, emitter_coefficient=0.0, coordinates=None))]
    #[allow(clippy::too_many_arguments)]
    fn add_junction(
        &self,
        py: Python<'_>,
        id: &str,
        elevation: f64,
        basedemand: f64,
        pattern: Option<String>,
        emitter_coefficient: f64,
        coordinates: Option<(f64, f64)>,
    ) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .add_junction(
                id,
                &JunctionData {
                    elevation,
                    demands: vec![RsDemand {
                        basedemand,
                        pattern: pattern.map(|p| p.into_boxed_str()),
                        pattern_index: None,
                        name: None,
                    }],
                    emitter_coefficient,
                    coordinates,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, elevation, initial_level, min_level, max_level, diameter, min_volume=0.0, volume_curve_id=None, overflow=false, coordinates=None))]
    #[allow(clippy::too_many_arguments)]
    fn add_tank(
        &self,
        py: Python<'_>,
        id: &str,
        elevation: f64,
        initial_level: f64,
        min_level: f64,
        max_level: f64,
        diameter: f64,
        min_volume: f64,
        volume_curve_id: Option<String>,
        overflow: bool,
        coordinates: Option<(f64, f64)>,
    ) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .add_tank(
                id,
                &TankData {
                    elevation,
                    initial_level,
                    min_level,
                    max_level,
                    diameter,
                    min_volume,
                    volume_curve_id: volume_curve_id.map(|s| s.into_boxed_str()),
                    overflow,
                    coordinates,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, elevation, head_pattern=None, coordinates=None))]
    fn add_reservoir(
        &self,
        py: Python<'_>,
        id: &str,
        elevation: f64,
        head_pattern: Option<String>,
        coordinates: Option<(f64, f64)>,
    ) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .add_reservoir(
                id,
                &ReservoirData {
                    elevation,
                    head_pattern: head_pattern.map(|s| s.into_boxed_str()),
                    coordinates,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, elevation=None, coordinates=None))]
    fn update_node(
        &self,
        py: Python<'_>,
        id: &str,
        elevation: Option<f64>,
        coordinates: Option<(f64, f64)>,
    ) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .update_node(
                id,
                &NodeUpdate {
                    elevation,
                    coordinates,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, elevation=None, basedemand=None, emitter_coefficient=None, clear_pattern=false, pattern=None, coordinates=None))]
    #[allow(clippy::too_many_arguments)]
    fn update_junction(
        &self,
        py: Python<'_>,
        id: &str,
        elevation: Option<f64>,
        basedemand: Option<f64>,
        emitter_coefficient: Option<f64>,
        clear_pattern: bool,
        pattern: Option<String>,
        coordinates: Option<(f64, f64)>,
    ) -> PyResult<()> {
        let pattern_update = if clear_pattern {
            Some(None)
        } else {
            pattern.map(|p| Some(p.into_boxed_str()))
        };
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .update_junction(
                id,
                &JunctionUpdate {
                    elevation,
                    demands: None,
                    basedemand,
                    pattern: pattern_update,
                    emitter_coefficient,
                    coordinates,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, elevation=None, coordinates=None, initial_level=None, min_level=None, max_level=None, diameter=None, min_volume=None, overflow=None))]
    #[allow(clippy::too_many_arguments)]
    fn update_tank(
        &self,
        py: Python<'_>,
        id: &str,
        elevation: Option<f64>,
        coordinates: Option<(f64, f64)>,
        initial_level: Option<f64>,
        min_level: Option<f64>,
        max_level: Option<f64>,
        diameter: Option<f64>,
        min_volume: Option<f64>,
        overflow: Option<bool>,
    ) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .update_tank(
                id,
                &TankUpdate {
                    elevation,
                    coordinates,
                    initial_level,
                    min_level,
                    max_level,
                    diameter,
                    min_volume,
                    overflow,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, elevation=None, coordinates=None, clear_head_pattern=false, head_pattern=None))]
    fn update_reservoir(
        &self,
        py: Python<'_>,
        id: &str,
        elevation: Option<f64>,
        coordinates: Option<(f64, f64)>,
        clear_head_pattern: bool,
        head_pattern: Option<String>,
    ) -> PyResult<()> {
        let head_pattern_update = if clear_head_pattern {
            Some(None)
        } else {
            head_pattern.map(|p| Some(p.into_boxed_str()))
        };
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .update_reservoir(
                id,
                &ReservoirUpdate {
                    elevation,
                    coordinates,
                    head_pattern: head_pattern_update,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, unconditional=false))]
    fn remove_node(&self, py: Python<'_>, id: &str, unconditional: bool) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .remove_node(id, unconditional)
            .map_err(input_err)
    }

    // --- pipes / pumps / valves ---

    #[pyo3(signature = (id, start_node, end_node, length, diameter, roughness, minor_loss=0.0, check_valve=false, initial_status=LinkStatus::Open, vertices=None))]
    #[allow(clippy::too_many_arguments)]
    fn add_pipe(
        &self,
        py: Python<'_>,
        id: &str,
        start_node: &str,
        end_node: &str,
        length: f64,
        diameter: f64,
        roughness: f64,
        minor_loss: f64,
        check_valve: bool,
        initial_status: LinkStatus,
        vertices: Option<Vec<(f64, f64)>>,
    ) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .add_pipe(
                id,
                &PipeData {
                    start_node: start_node.into(),
                    end_node: end_node.into(),
                    length,
                    diameter,
                    roughness,
                    minor_loss,
                    check_valve,
                    initial_status: initial_status.into(),
                    vertices,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, start_node, end_node, speed=1.0, head_curve_id=None, power=0.0, initial_status=LinkStatus::Open, vertices=None))]
    #[allow(clippy::too_many_arguments)]
    fn add_pump(
        &self,
        py: Python<'_>,
        id: &str,
        start_node: &str,
        end_node: &str,
        speed: f64,
        head_curve_id: Option<String>,
        power: f64,
        initial_status: LinkStatus,
        vertices: Option<Vec<(f64, f64)>>,
    ) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .add_pump(
                id,
                &PumpData {
                    start_node: start_node.into(),
                    end_node: end_node.into(),
                    speed,
                    head_curve_id: head_curve_id.map(|s| s.into_boxed_str()),
                    power,
                    initial_status: initial_status.into(),
                    vertices,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, start_node, end_node, diameter, valve_type, setting=0.0, curve_id=None, minor_loss=0.0, initial_status=LinkStatus::Active, vertices=None))]
    #[allow(clippy::too_many_arguments)]
    fn add_valve(
        &self,
        py: Python<'_>,
        id: &str,
        start_node: &str,
        end_node: &str,
        diameter: f64,
        valve_type: ValveType,
        setting: f64,
        curve_id: Option<String>,
        minor_loss: f64,
        initial_status: LinkStatus,
        vertices: Option<Vec<(f64, f64)>>,
    ) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .add_valve(
                id,
                &ValveData {
                    start_node: start_node.into(),
                    end_node: end_node.into(),
                    diameter,
                    valve_type: valve_type.into(),
                    setting,
                    curve_id: curve_id.map(|s| s.into_boxed_str()),
                    minor_loss,
                    initial_status: initial_status.into(),
                    vertices,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, length=None, diameter=None, roughness=None, minor_loss=None, check_valve=None))]
    #[allow(clippy::too_many_arguments)]
    fn update_pipe(
        &self,
        py: Python<'_>,
        id: &str,
        length: Option<f64>,
        diameter: Option<f64>,
        roughness: Option<f64>,
        minor_loss: Option<f64>,
        check_valve: Option<bool>,
    ) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .update_pipe(
                id,
                &PipeUpdate {
                    length,
                    diameter,
                    roughness,
                    minor_loss,
                    check_valve,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, speed=None, power=None, clear_head_curve=false, head_curve_id=None))]
    fn update_pump(
        &self,
        py: Python<'_>,
        id: &str,
        speed: Option<f64>,
        power: Option<f64>,
        clear_head_curve: bool,
        head_curve_id: Option<String>,
    ) -> PyResult<()> {
        let head_curve_update = if clear_head_curve {
            Some(None)
        } else {
            head_curve_id.map(|c| Some(c.into_boxed_str()))
        };
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .update_pump(
                id,
                &PumpUpdate {
                    speed,
                    power,
                    head_curve_id: head_curve_update,
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, diameter=None, setting=None, minor_loss=None, clear_curve=false, curve_id=None, valve_type=None))]
    #[allow(clippy::too_many_arguments)]
    fn update_valve(
        &self,
        py: Python<'_>,
        id: &str,
        diameter: Option<f64>,
        setting: Option<f64>,
        minor_loss: Option<f64>,
        clear_curve: bool,
        curve_id: Option<String>,
        valve_type: Option<ValveType>,
    ) -> PyResult<()> {
        let curve_update = if clear_curve {
            Some(None)
        } else {
            curve_id.map(|c| Some(c.into_boxed_str()))
        };
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .update_valve(
                id,
                &ValveUpdate {
                    diameter,
                    setting,
                    minor_loss,
                    curve_id: curve_update,
                    valve_type: valve_type.map(Into::into),
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, start_node=None, end_node=None, vertices=None, initial_status=None))]
    fn update_link(
        &self,
        py: Python<'_>,
        id: &str,
        start_node: Option<String>,
        end_node: Option<String>,
        vertices: Option<Vec<(f64, f64)>>,
        initial_status: Option<LinkStatus>,
    ) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .update_link(
                id,
                &LinkUpdate {
                    start_node: start_node.map(|s| s.into_boxed_str()),
                    end_node: end_node.map(|s| s.into_boxed_str()),
                    vertices,
                    initial_status: initial_status.map(Into::into),
                },
            )
            .map_err(input_err)
    }

    #[pyo3(signature = (id, unconditional=false))]
    fn remove_link(&self, py: Python<'_>, id: &str, unconditional: bool) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .remove_link(id, unconditional)
            .map_err(input_err)
    }

    // --- patterns ---

    fn add_pattern(&self, py: Python<'_>, id: &str, multipliers: Vec<f64>) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .add_pattern(id, &PatternData { multipliers })
            .map_err(input_err)
    }

    fn update_pattern(&self, py: Python<'_>, id: &str, multipliers: Vec<f64>) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .update_pattern(
                id,
                &PatternUpdate {
                    multipliers: Some(multipliers),
                },
            )
            .map_err(input_err)
    }

    fn remove_pattern(&self, py: Python<'_>, id: &str) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner.network.remove_pattern(id).map_err(input_err)
    }

    // --- curves ---

    fn add_curve(&self, py: Python<'_>, id: &str, x: Vec<f64>, y: Vec<f64>) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .add_curve(id, &CurveData { x, y })
            .map_err(input_err)
    }

    #[pyo3(signature = (id, x=None, y=None))]
    fn update_curve(
        &self,
        py: Python<'_>,
        id: &str,
        x: Option<Vec<f64>>,
        y: Option<Vec<f64>>,
    ) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner
            .network
            .update_curve(id, &CurveUpdate { x, y })
            .map_err(input_err)
    }

    fn remove_curve(&self, py: Python<'_>, id: &str) -> PyResult<()> {
        let mut sim = self.parent.borrow_mut(py);
        sim.inner.network.remove_curve(id).map_err(input_err)
    }

    // --- reads ---

    fn get_node(&self, py: Python<'_>, id: &str) -> PyResult<Node> {
        let sim = self.parent.borrow(py);
        let network = &sim.inner.network;
        let idx = *network
            .node_map
            .get(id)
            .ok_or_else(|| InputError::NodeNotFound { node_id: id.into() })
            .map_err(input_err)?;
        Ok(types::node_snapshot(&network.nodes[idx], &network.options))
    }

    fn get_link(&self, py: Python<'_>, id: &str) -> PyResult<Link> {
        let sim = self.parent.borrow(py);
        let network = &sim.inner.network;
        let idx = *network
            .link_map
            .get(id)
            .ok_or_else(|| InputError::LinkNotFound { link_id: id.into() })
            .map_err(input_err)?;
        Ok(types::link_snapshot(&network.links[idx], &network.options))
    }

    fn get_pattern(&self, py: Python<'_>, id: &str) -> PyResult<Pattern> {
        let sim = self.parent.borrow(py);
        let network = &sim.inner.network;
        let idx = *network
            .pattern_map
            .get(id)
            .ok_or_else(|| InputError::PatternNotFound {
                pattern_id: id.into(),
            })
            .map_err(input_err)?;
        Ok(types::pattern_snapshot(&network.patterns[idx]))
    }

    fn get_curve(&self, py: Python<'_>, id: &str) -> PyResult<Curve> {
        let sim = self.parent.borrow(py);
        let network = &sim.inner.network;
        let idx = *network
            .curve_map
            .get(id)
            .ok_or_else(|| InputError::CurveNotFound {
                curve_id: id.into(),
            })
            .map_err(input_err)?;
        Ok(types::curve_snapshot(&network.curves[idx]))
    }

    /// Lazily view every node without materializing all of them up front
    /// (see [`crate::native::views`]); use [`Network::get_node`] for a
    /// single lookup by id, or [`Network::node_map`] when only ids/indices are
    /// needed.
    fn nodes(&self, py: Python<'_>) -> NodeView {
        NodeView {
            parent: self.parent.clone_ref(py),
        }
    }

    /// Lazily view every link. See [`Network::nodes`].
    fn links(&self, py: Python<'_>) -> LinkView {
        LinkView {
            parent: self.parent.clone_ref(py),
        }
    }

    /// Lazily view every pattern. See [`Network::nodes`].
    fn patterns(&self, py: Python<'_>) -> PatternView {
        PatternView {
            parent: self.parent.clone_ref(py),
        }
    }

    /// Lazily view every curve. See [`Network::nodes`].
    fn curves(&self, py: Python<'_>) -> CurveView {
        CurveView {
            parent: self.parent.clone_ref(py),
        }
    }

    /// Lazily view every control. See [`Network::nodes`].
    fn controls(&self, py: Python<'_>) -> ControlView {
        ControlView {
            parent: self.parent.clone_ref(py),
        }
    }

    /// Live view mapping every node id to its position in [`Network::nodes`]
    /// — mirrors the real `Network::node_map` field. Cheaper than `nodes()`
    /// for existence checks (`"J1" in net.node_map()`) and is the way to
    /// correlate a node id with the positional arrays in `SolverState`/
    /// `SolverResult` (`state.heads[net.node_map()["J1"]]`), which are
    /// indexed the same way as `Network.nodes`.
    ///
    /// Like [`Network::nodes`]'s `NodeView`, this is live, not a snapshot:
    /// `in`/indexing/`len` always reflect the network's current state, so
    /// they keep working correctly across `add_*`/`remove_*` calls made
    /// after this method returns. Only `keys()`/`values()`/`items()`/
    /// `__iter__` capture a point-in-time snapshot (matching how Python's
    /// own `dict` iterators behave once iteration has started). Call
    /// `dict(net.node_map())` if you want an ordinary, fully-independent
    /// `dict` snapshot instead. One caveat regardless: indices shift after
    /// `remove_node`/`remove_link` (they use swap-remove), so an index
    /// looked up before either call can point at the wrong node afterwards.
    fn node_map(&self, py: Python<'_>) -> NodeMapView {
        NodeMapView {
            parent: self.parent.clone_ref(py),
        }
    }

    /// Live view mapping every link id to its position in [`Network::links`]. See [`Network::node_map`].
    fn link_map(&self, py: Python<'_>) -> LinkMapView {
        LinkMapView {
            parent: self.parent.clone_ref(py),
        }
    }

    /// Live view mapping every pattern id to its position in [`Network::patterns`]. See [`Network::node_map`].
    fn pattern_map(&self, py: Python<'_>) -> PatternMapView {
        PatternMapView {
            parent: self.parent.clone_ref(py),
        }
    }

    /// Live view mapping every curve id to its position in [`Network::curves`]. See [`Network::node_map`].
    fn curve_map(&self, py: Python<'_>) -> CurveMapView {
        CurveMapView {
            parent: self.parent.clone_ref(py),
        }
    }

    fn node_count(&self, py: Python<'_>) -> usize {
        self.parent.borrow(py).inner.network.nodes.len()
    }

    fn link_count(&self, py: Python<'_>) -> usize {
        self.parent.borrow(py).inner.network.links.len()
    }

    fn pattern_count(&self, py: Python<'_>) -> usize {
        self.parent.borrow(py).inner.network.patterns.len()
    }

    fn curve_count(&self, py: Python<'_>) -> usize {
        self.parent.borrow(py).inner.network.curves.len()
    }

    fn title(&self, py: Python<'_>) -> Option<String> {
        self.parent
            .borrow(py)
            .inner
            .network
            .title
            .as_deref()
            .map(String::from)
    }

    fn has_tanks(&self, py: Python<'_>) -> bool {
        self.parent.borrow(py).inner.network.has_tanks()
    }

    fn has_pressure_controls(&self, py: Python<'_>) -> bool {
        self.parent.borrow(py).inner.network.has_pressure_controls()
    }

    /// Save the network back out to an INP file. Mirrors `Network::save_network`.
    fn save_network(&self, py: Python<'_>, path: &str) -> PyResult<()> {
        self.parent
            .borrow(py)
            .inner
            .network
            .save_network(path)
            .map_err(PyRuntimeError::new_err)
    }
}
