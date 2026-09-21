# epanet-rs

A fast, modern re-implementation of the EPANET2 hydraulic solver, written in Rust with Python bindings.

## Installation

```bash
pip install epanet-rs
```

## Usage

The API mirrors the EPANET 2.3 toolkit, so existing EPANET code translates directly:

```python
from epanet_rs import Project, EN_NODECOUNT, EN_LINKCOUNT, EN_HEAD, EN_FLOW

# Create and open a project
p = Project()
p.open("network.inp")

# Run a full simulation
result = p.solveH()
print(f"Time steps: {len(result.heads)}")
print(f"Head at node 1, step 0: {result.heads[0][0]}")

# Or use the step-by-step API
p.openH()
p.initH()
while True:
    t = p.runH(p.time if hasattr(p, 'time') else 0)
    dt = p.nextH()
    if dt == 0:
        break
p.closeH()

# Query network
n_nodes = p.getcount(EN_NODECOUNT)
n_links = p.getcount(EN_LINKCOUNT)

# Get results after solving
head = p.getnodevalue(1, EN_HEAD)
flow = p.getlinkvalue(1, EN_FLOW)

# Parallel solving (networks without tanks/pressure controls)
result = p.solveH(parallel=True)

p.close()
```

## API Compatibility

This package exposes the same functions as the EPANET 2.3 C toolkit, using Pythonic conventions:

- `Project()` replaces `EN_createproject` / `EN_deleteproject`
- Methods like `open`, `openH`, `runH`, `nextH`, `solveH`, `closeH` map directly
- Node/link/pattern/curve accessors: `getnodeindex`, `getnodevalue`, `setnodevalue`, etc.
- All EPANET constants are available: `EN_ELEVATION`, `EN_HEAD`, `EN_FLOW`, etc.

## Native Rust API

`epanet_rs.native` exposes the underlying Rust API directly, instead of the
numeric property codes and 1-based indices used by the EPANET-2.3-compatible
`Project`. Values are read and written in the network's configured units, and
`add_*`/`update_*`/`remove_*` methods use keyword arguments in place of the
Rust `*Data`/`*Update` structs:

```python
from epanet_rs.native import Simulation, FlowUnits, HeadlossFormula, ValveType

# Load an existing network
sim = Simulation.from_file("network.inp")
net = sim.network()

# Build up a network from scratch instead
sim = Simulation.new(FlowUnits.GPM, HeadlossFormula.HazenWilliams)
net = sim.network()
net.add_reservoir("R1", elevation=700.0)
net.add_junction("J1", elevation=650.0, basedemand=100.0)
net.add_pipe("P1", start_node="R1", end_node="J1", length=1000.0, diameter=12.0, roughness=120.0)

# Modify the network
net.update_pipe("P1", roughness=140.0)

# Solve (sequentially or in parallel)
result = sim.solve_hydraulics(parallel=True)
print(result.heads[-1])

# Inspect the network with typed, named-field snapshots instead of property codes
node = net.get_node("J1")
print(node.elevation, node.junction.demands[0].basedemand)

for link in net.links():
    if link.pipe is not None:
        print(link.id, link.pipe.diameter, link.pipe.roughness)
```

`nodes()`/`links()`/`patterns()`/`curves()`/`controls()` return a lazy,
index-able view rather than a plain list: `len(net.nodes())` and
`net.nodes()[0]` don't build a snapshot for every node, and iterating only
builds one snapshot at a time, so `next(iter(net.nodes()))` or
`itertools.islice(net.links(), 5)` doesn't pay for the rest of a large
network. `node_map()`/`link_map()`/`pattern_map()`/`curve_map()` mirror the
Rust `Network::node_map` field (id -> position): a dict-like, *live* view
that gives O(1) existence checks and index lookups without building any
node/link snapshot at all — and they're the way to correlate an id with the
positional arrays in `SolverState`/`SolverResult`, which are indexed the
same way as `Network.nodes`/`links`:

```python
first_node = net.nodes()[0]          # builds exactly one Node
n_nodes = len(net.nodes())           # no Node objects built at all
node_map = net.node_map()            # dict-like live view, no Node objects built
"J1" in node_map                     # O(1) existence check
node_map["J1"]                       # O(1) index lookup, raises KeyError if missing
node_map.get("J1")                   # like dict.get: None (or a given default) if missing

state = sim.state()
head_at_j1 = state.heads[node_map["J1"]]  # correlate id -> flat solver array
```

Like the views above, `node_map()` is live rather than a snapshot:
`in`/indexing/`len`/`get` on a `node_map()` you're holding onto keep
reflecting the network even after a later `add_junction`/`remove_node`/etc.
Only `keys()`/`values()`/`items()`/iteration capture a point-in-time
snapshot (the same way Python's own `dict` iterators behave once iteration
has started) — call `dict(net.node_map())` if you want an ordinary, fully
independent `dict` snapshot instead. One caveat regardless of live-ness:
indices shift after a `remove_node`/`remove_link` call specifically (they
use swap-remove), so an index read before either call can silently point at
a *different* node afterwards — re-read the index after mutating.

`InputError` and `SolverError` (`epanet_rs.native.InputError` /
`epanet_rs.native.SolverError`) mirror the corresponding Rust error enums and
are raised in place of generic runtime errors.

### Stability

`epanet_rs.native` is public API and intended for external use, but the
package is still pre-1.0 (`0.x`). Until a `1.0` release:

- Breaking changes to `epanet_rs.native` will be called out explicitly in
  release notes rather than assumed safe for any `0.x` bump.
- The following are deliberate, documented design choices rather than
  incidental details, and are the parts most likely to be revisited based on
  real-world usage before `1.0`: individual reads (`get_node`, indexing a
  view, iterating a view) return frozen snapshot objects — mutation always
  goes through `Network`'s `add_*`/`update_*`/`remove_*` methods — but the
  *containers* that produce them (`nodes()`/`links()`/`patterns()`/
  `curves()`/`controls()`'s views, and `node_map()`/`link_map()`/
  `pattern_map()`/`curve_map()`) are live rather than snapshots, and reflect
  later mutations without needing to be re-fetched; collection accessors
  are plain methods rather than properties, both because `network()`
  allocates a new wrapper and because the lazy/live nature of the views and
  maps (see above) would be obscured by property syntax; single-integer
  indexing on the sequence views does not support slicing; and tri-state
  "clear vs. leave untouched vs. set" updates (e.g. a junction's demand
  pattern) are expressed as an explicit `clear_pattern=True` flag alongside
  the value argument, rather than a sentinel or nested
  `Optional[Optional[...]]`.
- The EPANET-2.3-toolkit-mirroring `Project`/`SolverResult`/`EN_*` API is
  unaffected by any of the above and follows its own, longer-established
  compatibility expectations.

## Parallel Solving

For networks without tanks or pressure controls, `solveH(parallel=True)` runs the
extended-period simulation in parallel using all available CPU cores, providing up to
5x speedup on large networks.
