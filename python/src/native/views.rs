//! Lazy, index-able views over `Network` collections (`nodes()`, `links()`,
//! `patterns()`, `curves()`, `controls()`).
//!
//! `Network::nodes()` and friends used to eagerly build a `Vec<Node>`
//! (and friends), which meant every call fully materialized *every* element
//! of the network into nested Python objects before returning — expensive
//! and wasteful for large networks or for callers that only wanted to peek
//! at the first few elements, index a single one, or hold a reference to
//! iterate later. The types here instead return a lightweight view
//! (`len()`/indexing/`repr` only touch the underlying `Network`, no
//! snapshots are built) whose `__iter__` yields snapshots one at a time on
//! demand, so `for x in net.nodes(): ...`, `net.nodes()[0]`,
//! `next(iter(net.nodes()))`, and `itertools.islice(net.nodes(), 5)` only
//! pay for the elements actually touched.
//!
//! For cases that only need identifiers or O(1) existence/index lookups
//! (e.g. to correlate an id with the positional arrays in `SolverState`/
//! `SolverResult`, which are indexed the same way as these collections),
//! `Network::node_map()`/`link_map()`/`pattern_map()`/`curve_map()` mirror
//! the real `Network::node_map` field and avoid snapshot construction
//! entirely.

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

use epanet::model::control::Control as RsControl;
use epanet::model::curve::Curve as RsCurve;
use epanet::model::link::Link as RsLink;
use epanet::model::network::Network as RsNetwork;
use epanet::model::node::Node as RsNode;
use epanet::model::pattern::Pattern as RsPattern;

use super::simulation::Simulation as PySimulation;
use super::types;
use super::types::{Control, Curve, Link, Node, Pattern};

/// Generates a `{Name}View` (len/index/iter over `network.$field`, no
/// eager snapshot construction) and matching `{Name}Iter` (yields one
/// snapshot per `__next__`) pair. `$snapshot` is `|item, network| -> $item`.
macro_rules! define_view {
    ($view:ident, $iter:ident, $item:ty, $field:ident, $snapshot:expr) => {
        #[doc = concat!(
                            "Lazy, index-able view over `Network::",
                            stringify!($field),
                            "`. See [`crate::native::views`] for why this isn't a plain list."
                        )]
        #[pyclass(module = "epanet_rs.native")]
        pub struct $view {
            pub(crate) parent: Py<PySimulation>,
        }

        #[pymethods]
        impl $view {
            fn __len__(&self, py: Python<'_>) -> usize {
                self.parent.borrow(py).inner.network.$field.len()
            }

            fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<$item> {
                let sim = self.parent.borrow(py);
                let network = &sim.inner.network;
                let len = network.$field.len() as isize;
                let resolved = if index < 0 { index + len } else { index };
                if resolved < 0 || resolved >= len {
                    return Err(PyIndexError::new_err("index out of range"));
                }
                let item = &network.$field[resolved as usize];
                Ok(($snapshot)(item, network))
            }

            fn __iter__(&self, py: Python<'_>) -> $iter {
                $iter {
                    parent: self.parent.clone_ref(py),
                    index: 0,
                }
            }

            fn __repr__(&self, py: Python<'_>) -> String {
                format!(concat!(stringify!($view), "(len={})"), self.__len__(py))
            }
        }

        #[doc = concat!("Lazy iterator produced by `", stringify!($view), ".__iter__`.")]
        #[pyclass(module = "epanet_rs.native")]
        pub struct $iter {
            parent: Py<PySimulation>,
            index: usize,
        }

        #[pymethods]
        impl $iter {
            fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
                slf
            }

            fn __next__(&mut self, py: Python<'_>) -> Option<$item> {
                let sim = self.parent.borrow(py);
                let network = &sim.inner.network;
                let item = network.$field.get(self.index)?;
                let snapshot = ($snapshot)(item, network);
                self.index += 1;
                Some(snapshot)
            }
        }
    };
}

define_view!(
    NodeView,
    NodeIter,
    Node,
    nodes,
    |n: &RsNode, network: &RsNetwork| { types::node_snapshot(n, &network.options) }
);
define_view!(
    LinkView,
    LinkIter,
    Link,
    links,
    |l: &RsLink, network: &RsNetwork| { types::link_snapshot(l, &network.options) }
);
define_view!(
    PatternView,
    PatternIter,
    Pattern,
    patterns,
    |p: &RsPattern, _network: &RsNetwork| { types::pattern_snapshot(p) }
);
define_view!(
    CurveView,
    CurveIter,
    Curve,
    curves,
    |c: &RsCurve, _network: &RsNetwork| { types::curve_snapshot(c) }
);
define_view!(
    ControlView,
    ControlIter,
    Control,
    controls,
    |c: &RsControl, network: &RsNetwork| { types::control_snapshot(c, network, &network.options) }
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<NodeView>()?;
    m.add_class::<NodeIter>()?;
    m.add_class::<LinkView>()?;
    m.add_class::<LinkIter>()?;
    m.add_class::<PatternView>()?;
    m.add_class::<PatternIter>()?;
    m.add_class::<CurveView>()?;
    m.add_class::<CurveIter>()?;
    m.add_class::<ControlView>()?;
    m.add_class::<ControlIter>()?;
    Ok(())
}
