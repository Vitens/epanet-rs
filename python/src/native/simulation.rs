//! `Simulation`: a thin, idiomatic wrapper around `epanet_rs::simulation::Simulation`.
//!
//! Method names and behaviour mirror the Rust API directly (see
//! `epanet_rs::simulation::Simulation`'s doc comments for the exact
//! semantics of each step).

use pyo3::prelude::*;

use epanet::simulation::Simulation as RsSimulation;

use super::enums::{FlowUnits, HeadlossFormula, LinkStatus};
use super::network::Network;
use super::{input_err, solved_state, solver_err};
use crate::PySolverResult as SolverResult;

/// Snapshot of the hydraulic solver state at a single point in time
/// (converted to the network's configured units), available after
/// `run_hydraulics()`. Mirrors `epanet_rs::solver::state::SolverState`.
#[pyclass(frozen, get_all, module = "epanet_rs.native")]
#[derive(Clone, Debug)]
pub struct SolverState {
    pub flows: Vec<f64>,
    pub heads: Vec<f64>,
    pub demands: Vec<f64>,
    pub statuses: Vec<LinkStatus>,
}

#[pyclass(module = "epanet_rs.native", name = "Simulation")]
pub struct Simulation {
    pub(crate) inner: RsSimulation,
}

#[pymethods]
impl Simulation {
    /// Load a network from an INP file. Mirrors `Simulation::from_file`.
    #[staticmethod]
    fn from_file(path: &str) -> PyResult<Self> {
        let inner = RsSimulation::from_file(path).map_err(input_err)?;
        Ok(Self { inner })
    }

    /// Create an empty simulation with the given units. Mirrors `Simulation::init`.
    #[staticmethod]
    fn new(flow_units: FlowUnits, headloss_formula: HeadlossFormula) -> Self {
        Self {
            inner: RsSimulation::init(flow_units.into(), headloss_formula.into()),
        }
    }

    /// Returns a `Network` handle for reading/modifying this simulation's network.
    fn network(slf: &Bound<'_, Self>) -> Network {
        Network {
            parent: slf.clone().unbind(),
        }
    }

    /// Resets the simulation to initial conditions. Mirrors `Simulation::initialize_hydraulics`
    /// (equivalent to `EN_openH` followed by `EN_initH`).
    fn initialize_hydraulics(&mut self) -> PyResult<()> {
        self.inner.initialize_hydraulics().map_err(solver_err)
    }

    /// Applies patterns/controls and solves hydraulics at the current time.
    /// Mirrors `Simulation::run_hydraulics` (equivalent to `EN_runH`).
    fn run_hydraulics(&mut self) -> PyResult<usize> {
        self.inner.run_hydraulics().map_err(solver_err)
    }

    /// Advances to the next hydraulic timestep, returning its length in
    /// seconds (0 when the simulation is complete). Mirrors
    /// `Simulation::next_hydraulic_timestep` (equivalent to `EN_nextH`).
    fn next_hydraulic_timestep(&mut self) -> usize {
        self.inner.next_hydraulic_timestep()
    }

    /// Runs a complete extended-period simulation. Mirrors
    /// `Simulation::solve_hydraulics` (equivalent to `EN_solveH`); pass
    /// `parallel=True` to solve report steps in parallel (only supported for
    /// networks without tanks or pressure controls).
    #[pyo3(signature = (parallel=false))]
    fn solve_hydraulics(&mut self, parallel: bool) -> PyResult<SolverResult> {
        let result = self.inner.solve_hydraulics(parallel).map_err(solver_err)?;
        Ok(SolverResult {
            flows: result.flows,
            heads: result.heads,
            demands: result.demands,
        })
    }

    /// The current simulation time, in seconds.
    #[getter]
    fn time(&self) -> usize {
        self.inner.time
    }

    /// Whether the simulation has a valid solved state.
    #[getter]
    fn solved(&self) -> bool {
        self.inner.solved
    }

    /// The solved hydraulic state at the current time (converted to the
    /// network's configured units), or `None` if nothing has been solved yet.
    fn state(&self) -> Option<SolverState> {
        let options = &self.inner.network.options;
        let flow_scale = options.flow_units.per_cfs();
        let head_scale = options.unit_system.per_feet();

        solved_state(&self.inner).map(|state| SolverState {
            flows: state.flows.iter().map(|f| f * flow_scale).collect(),
            heads: state.heads.iter().map(|h| h * head_scale).collect(),
            demands: state.demands.iter().map(|d| d * flow_scale).collect(),
            statuses: state.statuses.iter().map(|s| (*s).into()).collect(),
        })
    }
}
