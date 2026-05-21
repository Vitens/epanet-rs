//! Python bindings exposing hydraulic simulation results as Polars DataFrames.

use polars::prelude::*;
use pyo3::exceptions::{PyIOError, PyRuntimeError};
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

use crate::error::{InputError, SolverError};
use crate::simulation::Simulation;
use crate::solver::result::SolverResult;

fn input_error(err: InputError) -> PyErr {
    PyIOError::new_err(err.to_string())
}

fn solver_error(err: SolverError) -> PyErr {
    PyRuntimeError::new_err(err.to_string())
}

fn polars_error(err: PolarsError) -> PyErr {
    PyRuntimeError::new_err(err.to_string())
}

fn node_dataframe(
    values: &[Vec<f64>],
    node_ids: &[String],
    report_timestep: usize,
    value_name: &str,
) -> PyResult<PyDataFrame> {
    let n_steps = values.len();
    let n_nodes = node_ids.len();
    let len = n_steps * n_nodes;

    let mut seconds = Vec::with_capacity(len);
    let mut ids = Vec::with_capacity(len);
    let mut flattened = Vec::with_capacity(len);

    for (step, row) in values.iter().enumerate() {
        let time = (step * report_timestep) as u32;
        for (node_idx, node_id) in node_ids.iter().enumerate() {
            seconds.push(time);
            ids.push(node_id.as_str());
            flattened.push(row[node_idx]);
        }
    }

    let df = DataFrame::new(
        len,
        vec![
            Column::new("time".into(), seconds),
            Column::new("node_id".into(), ids),
            Column::new(value_name.into(), flattened),
        ],
    )
    .map_err(polars_error)?;

    Ok(PyDataFrame(df))
}

fn link_dataframe(
    values: &[Vec<f64>],
    link_ids: &[String],
    report_timestep: usize,
    value_name: &str,
) -> PyResult<PyDataFrame> {
    let n_steps = values.len();
    let n_links = link_ids.len();
    let len = n_steps * n_links;

    let mut seconds = Vec::with_capacity(len);
    let mut ids = Vec::with_capacity(len);
    let mut flattened = Vec::with_capacity(len);

    for (step, row) in values.iter().enumerate() {
        let time = (step * report_timestep) as u32;
        for (link_idx, link_id) in link_ids.iter().enumerate() {
            seconds.push(time);
            ids.push(link_id.as_str());
            flattened.push(row[link_idx]);
        }
    }

    let df = DataFrame::new(
        len,
        vec![
            Column::new("time".into(), seconds),
            Column::new("link_id".into(), ids),
            Column::new(value_name.into(), flattened),
        ],
    )
    .map_err(polars_error)?;

    Ok(PyDataFrame(df))
}

/// Hydraulic simulation results exposed as Polars DataFrames.
///
/// Numeric columns are transferred to Python without copying the underlying buffers.
#[pyclass(name = "HydraulicResults")]
pub struct PyHydraulicResults {
    results: SolverResult,
    node_ids: Vec<String>,
    link_ids: Vec<String>,
    report_timestep: usize,
}

#[pymethods]
impl PyHydraulicResults {
    /// Node heads at each report step.
    ///
    /// Columns: `time`, `node_id`, `head`.
    fn heads(&self) -> PyResult<PyDataFrame> {
        node_dataframe(
            &self.results.heads,
            &self.node_ids,
            self.report_timestep,
            "head",
        )
    }

    /// Link flows at each report step.
    ///
    /// Columns: `time`, `link_id`, `flow`.
    fn flows(&self) -> PyResult<PyDataFrame> {
        link_dataframe(
            &self.results.flows,
            &self.link_ids,
            self.report_timestep,
            "flow",
        )
    }

    /// Node demands at each report step.
    ///
    /// Columns: `time`, `node_id`, `demand`.
    fn demands(&self) -> PyResult<PyDataFrame> {
        node_dataframe(
            &self.results.demands,
            &self.node_ids,
            self.report_timestep,
            "demand",
        )
    }

    /// All hydraulic results in one long-format DataFrame.
    ///
    /// Columns: `time`, `element_id`, `element_type`, `head`, `demand`, `flow`.
    /// Node rows have `head` and `demand`; link rows have `flow` (other value columns are null).
    fn dataframe(&self) -> PyResult<PyDataFrame> {
        let n_steps = self.results.heads.len();
        let n_nodes = self.node_ids.len();
        let n_links = self.link_ids.len();
        let len = n_steps * (n_nodes + n_links);

        let mut seconds = Vec::with_capacity(len);
        let mut element_ids = Vec::with_capacity(len);
        let mut element_types = Vec::with_capacity(len);
        let mut heads = Vec::with_capacity(len);
        let mut demands = Vec::with_capacity(len);
        let mut flows = Vec::with_capacity(len);

        for (step, head_row) in self.results.heads.iter().enumerate() {
            let time = (step * self.report_timestep) as u32;
            let demand_row = &self.results.demands[step];
            for (node_idx, node_id) in self.node_ids.iter().enumerate() {
                seconds.push(time);
                element_ids.push(node_id.as_str());
                element_types.push("node");
                heads.push(Some(head_row[node_idx]));
                demands.push(Some(demand_row[node_idx]));
                flows.push(None);
            }
            let flow_row = &self.results.flows[step];
            for (link_idx, link_id) in self.link_ids.iter().enumerate() {
                seconds.push(time);
                element_ids.push(link_id.as_str());
                element_types.push("link");
                heads.push(None);
                demands.push(None);
                flows.push(Some(flow_row[link_idx]));
            }
        }

        let df = DataFrame::new(
            len,
            vec![
                Column::new("time".into(), seconds),
                Column::new("element_id".into(), element_ids),
                Column::new("element_type".into(), element_types),
                Column::new("head".into(), heads),
                Column::new("demand".into(), demands),
                Column::new("flow".into(), flows),
            ],
        )
        .map_err(polars_error)?;

        Ok(PyDataFrame(df))
    }

    #[getter]
    fn report_steps(&self) -> usize {
        self.results.heads.len()
    }
}

/// EPANET hydraulic simulation driver.
#[pyclass(name = "Simulation")]
pub struct PySimulation {
    inner: Simulation,
}

#[pymethods]
impl PySimulation {
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let inner = Simulation::from_file(path).map_err(input_error)?;
        Ok(Self { inner })
    }

    /// Run a full hydraulic simulation and return results as Polars DataFrames.
    #[pyo3(signature = (parallel = false))]
    fn solve(&mut self, parallel: bool) -> PyResult<PyHydraulicResults> {
        let report_timestep = self.inner.network.options.time_options.report_timestep;
        let node_ids = self
            .inner
            .network
            .nodes
            .iter()
            .map(|node| node.id.to_string())
            .collect();
        let link_ids = self
            .inner
            .network
            .links
            .iter()
            .map(|link| link.id.to_string())
            .collect();

        let results = self
            .inner
            .solve_hydraulics(parallel)
            .map_err(solver_error)?;

        Ok(PyHydraulicResults {
            results,
            node_ids,
            link_ids,
            report_timestep,
        })
    }
}
