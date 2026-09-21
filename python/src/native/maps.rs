//! Live, dict-like views over `Network`'s id -> index maps (`node_map`,
//! `link_map`, `pattern_map`, `curve_map`).
//!
//! Unlike a plain Python `dict` snapshot, `__len__`/`__contains__`/
//! `__getitem__`/`get`/`keys`/`values`/`items` on these types always read
//! the network's *current* state — the same "no snapshot" model as
//! [`crate::native::views`]'s `NodeView`/`LinkView`/etc, and for the same
//! reason: `in`/indexing/`len` are all direct hash lookups against whatever
//! the network looks like *right now*, so there's no stale-snapshot problem
//! to begin with. `__iter__` snapshots the key list at the moment it's
//! called (matching how Python's own `dict` iterators behave once
//! iteration has started), but that only affects iteration order/set, not
//! `in`/indexing/`len`.
//!
//! Call `dict(net.node_map())` to materialize an independent, ordinary
//! `dict` snapshot when that's what you want instead.

use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;

use super::simulation::Simulation as PySimulation;

/// Iterator over a pre-captured snapshot of a map's keys, shared by every
/// `*MapView::__iter__` (enumeration order/set is fixed at the point
/// `__iter__` is called, matching Python's own dict iterators).
#[pyclass(module = "epanet_rs.native")]
pub struct MapKeyIter {
    keys: std::vec::IntoIter<String>,
}

#[pymethods]
impl MapKeyIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<String> {
        self.keys.next()
    }
}

/// Generates a live, dict-like `{Name}` view over `network.$field`
/// (a `HashMap<Box<str>, usize>`). See the module docs for why this can be
/// live where a plain returned `dict` could not.
macro_rules! define_map_view {
    ($view:ident, $field:ident) => {
        #[doc = concat!(
                            "Live, dict-like view over `Network::",
                            stringify!($field),
                            "`. See [`crate::native::maps`] for why this isn't a plain `dict`."
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

            fn __contains__(&self, py: Python<'_>, key: &str) -> bool {
                self.parent
                    .borrow(py)
                    .inner
                    .network
                    .$field
                    .contains_key(key)
            }

            fn __getitem__(&self, py: Python<'_>, key: &str) -> PyResult<usize> {
                self.parent
                    .borrow(py)
                    .inner
                    .network
                    .$field
                    .get(key)
                    .copied()
                    .ok_or_else(|| PyKeyError::new_err(key.to_string()))
            }

            /// Like `dict.get`. `default` is `i64` (not `usize`) so common
            /// sentinels like `-1` work; a found index is always
            /// non-negative and fits `i64` for any realistic network size.
            #[pyo3(signature = (key, default=None))]
            fn get(&self, py: Python<'_>, key: &str, default: Option<i64>) -> Option<i64> {
                self.parent
                    .borrow(py)
                    .inner
                    .network
                    .$field
                    .get(key)
                    .map(|idx| *idx as i64)
                    .or(default)
            }

            /// Current keys, snapshotted (also usable as the source for
            /// `dict(view)`, since `dict()` looks for a `.keys()` method).
            fn keys(&self, py: Python<'_>) -> Vec<String> {
                self.parent
                    .borrow(py)
                    .inner
                    .network
                    .$field
                    .keys()
                    .map(|k| k.to_string())
                    .collect()
            }

            /// Current values, snapshotted.
            fn values(&self, py: Python<'_>) -> Vec<usize> {
                self.parent
                    .borrow(py)
                    .inner
                    .network
                    .$field
                    .values()
                    .copied()
                    .collect()
            }

            /// Current (key, value) pairs, snapshotted.
            fn items(&self, py: Python<'_>) -> Vec<(String, usize)> {
                self.parent
                    .borrow(py)
                    .inner
                    .network
                    .$field
                    .iter()
                    .map(|(k, v)| (k.to_string(), *v))
                    .collect()
            }

            fn __iter__(&self, py: Python<'_>) -> MapKeyIter {
                MapKeyIter {
                    keys: self.keys(py).into_iter(),
                }
            }

            fn __repr__(&self, py: Python<'_>) -> String {
                format!(concat!(stringify!($view), "(len={})"), self.__len__(py))
            }
        }
    };
}

define_map_view!(NodeMapView, node_map);
define_map_view!(LinkMapView, link_map);
define_map_view!(PatternMapView, pattern_map);
define_map_view!(CurveMapView, curve_map);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<MapKeyIter>()?;
    m.add_class::<NodeMapView>()?;
    m.add_class::<LinkMapView>()?;
    m.add_class::<PatternMapView>()?;
    m.add_class::<CurveMapView>()?;
    Ok(())
}
