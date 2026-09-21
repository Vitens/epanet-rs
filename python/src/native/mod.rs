//! Idiomatic Python bindings for the native `epanet-rs` Rust API, exposed as
//! the `epanet_rs.native` submodule — as opposed to [`crate::PyProject`],
//! which mirrors the EPANET 2.3 C toolkit API (numeric property codes,
//! 1-based indices).
//!
//! This module exposes the same functions the Rust API exposes: methods on
//! [`simulation::Simulation`] and [`network::Network`] map directly onto
//! `epanet_rs::simulation::Simulation` and
//! `epanet_rs::model::network::Network`'s `add_*`/`update_*`/`remove_*`
//! methods, using Python keyword arguments in place of the `*Data`/`*Update`
//! structs. Reads return the snapshot types in [`types`] instead of
//! numeric property codes; [`network::Network`]'s collection accessors
//! (`nodes()`, `links()`, ...) return the lazy, live, index-able views in
//! [`views`] rather than eagerly building a snapshot per element, and its
//! `node_map()`/`link_map()`/... return the live, dict-like views in
//! [`maps`].
//!
//! # Stability
//!
//! This is public API intended for external use, but the package is still
//! pre-1.0 (`0.x`); breaking changes here are called out explicitly in
//! release notes rather than assumed safe for any `0.x` bump. See the
//! "Stability" section of `python/README.md` for the specific design
//! choices (frozen snapshot reads, method- rather than property-style
//! collection accessors, explicit `clear_*` flags for tri-state updates)
//! most likely to be revisited before `1.0`.

pub mod enums;
pub mod maps;
pub mod network;
pub mod simulation;
pub mod types;
pub mod views;

use pyo3::create_exception;
use pyo3::prelude::*;

use epanet::error::{InputError as RsInputError, SolverError as RsSolverError};
use epanet::simulation::Simulation as RsSimulation;
use epanet::solver::state::SolverState as RsSolverState;

create_exception!(epanet_rs.native, InputError, pyo3::exceptions::PyException);
create_exception!(epanet_rs.native, SolverError, pyo3::exceptions::PyException);

pub(crate) fn input_err(e: RsInputError) -> PyErr {
    InputError::new_err(e.to_string())
}

pub(crate) fn solver_err(e: RsSolverError) -> PyErr {
    SolverError::new_err(e.to_string())
}

/// Returns the solved state only if the simulation is currently in a solved
/// state (mirrors the crate-private `Simulation::solved_state` helper, which
/// is not visible outside the `epanet-rs` crate).
pub(crate) fn solved_state(sim: &RsSimulation) -> Option<&RsSolverState> {
    if sim.solved { sim.state.as_ref() } else { None }
}

/// Registers the `epanet_rs.native` submodule on the parent module and makes
/// it importable as `epanet_rs.native` / `from epanet_rs.native import ...`
/// (plain `add_submodule` alone does not register submodules in
/// `sys.modules`, see https://github.com/PyO3/pyo3/issues/759).
pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = parent.py();
    let m = PyModule::new(py, "native")?;
    m.setattr(
        "__doc__",
        "Idiomatic bindings for the native epanet-rs Rust API (Simulation, Network, ...).\n\n\
         Public API intended for external use, but pre-1.0 (0.x): breaking changes are \
         called out explicitly in release notes. See the \"Stability\" section of \
         python/README.md.",
    )?;

    m.add("InputError", py.get_type::<InputError>())?;
    m.add("SolverError", py.get_type::<SolverError>())?;

    enums::register(&m)?;
    types::register(&m)?;
    views::register(&m)?;
    maps::register(&m)?;
    m.add_class::<network::Network>()?;
    m.add_class::<simulation::Simulation>()?;
    m.add_class::<simulation::SolverState>()?;
    m.add_class::<crate::PySolverResult>()?;

    parent.add_submodule(&m)?;

    py.import("sys")?
        .getattr("modules")?
        .set_item("epanet_rs.native", &m)?;

    Ok(())
}
