/*!
EPANET-RS: A fast, modern re-implementation of the EPANET2 hydraulic solver

## Background

The EPANET2 solver has been the industry standard for hydraulic network simulation for decades due
to its robustness, numerical stability, and extensive validation. Its algorithms and results are
widely trusted in both academic and industrial applications.

The core EPANET2 codebase, however, is several decades old and written in C, making it difficult
to maintain, extend, and optimize using modern software engineering practices. In addition, the
original implementation predates modern CPU architectures and therefore does not fully exploit
multi-core processors or SIMD (Single Instruction, Multiple Data) capabilities.

Modern applications of the EPANET solver, such as monte-carlo simulations, leak detection
algorithms and real-time digital twins of huge networks require a more modern implementation,
with better performance and maintainability.

`epanet-rs` is a modern reimplementation of the EPANET2 hydraulic solver written in Rust,
designed to preserve the original algorithms and numerical behavior while enabling safer memory
management, improved maintainability, and performance optimizations through multi-threading and
SIMD acceleration.

`epanet-rs` runs about as fast as the original EPANET2_3 solver in sequential mode, and up to 5
times faster in parallel mode for extended period simulations for supported networks (no
tanks/controls).

## Design Goals

- **Numerical Parity** with EPANET2_3 solver (`hyd_2_3` branch), ensuring identical or equivalent
  results for the same input
- **High Performance** through multi-threading, SIMD acceleration and a modern,
  [faer](https://crates.io/crates/faer) based sparse solver
- **Parallelization** of the solver loop by rewriting the solver to return vectors of heads and
  flows instead of in-place assignment
- **Memory Safety** through Rust's ownership and borrowing system
- **Modern API** with a focus on ease of use and parallelization
- **Backwards Compatibility** with EPANET2_3 network API methods

## Crate layout

The public API is organised into a few top-level modules:

- [`model`] — data structures describing a hydraulic network (nodes, links, patterns, curves,
  controls, options). Networks can be loaded from EPANET `.inp` files or constructed
  programmatically.
- [`io`] — parsers and serializers for EPANET `.inp` files and JSON/MessagePack snapshots of a
  [`model::network::Network`].
- [`solver`] — the low-level hydraulic solver. [`solver::hydraulicsolver::HydraulicSolver`]
  performs a single-step pressure/flow solve on a pre-built sparsity pattern, and
  [`solver::state::SolverState`] holds the per-step mutable state (flows, heads, demands, link
  statuses, ...).
- [`simulation`] — the high-level [`simulation::Simulation`] driver that mirrors the EPANET
  `EN_openH` / `EN_initH` / `EN_runH` / `EN_nextH` / `EN_solveH` workflow, including tanks,
  controls and reporting.
- [`ffi`] — a subset of the EPANET C API exposed as `extern "C"` functions, for use from other
  languages.
- [`error`] — error types returned by the parser ([`error::InputError`]) and solver
  ([`error::SolverError`]).

## Examples

### Loading and solving a network

Use the high-level [`simulation::Simulation`] API to run a full extended-period simulation and
collect flows, heads and demands at every report step:

```no_run
use epanet_rs::simulation::Simulation;

# fn main() -> Result<(), Box<dyn std::error::Error>> {
let mut simulation = Simulation::from_file("tests/pump.inp")?;

// `parallel = false` runs the EPANET-style sequential loop, which supports
// tanks, controls, and quality. Pass `true` for parallel solving when the
// network has no tanks or pressure controls.
let results = simulation.solve_hydraulics(false)?;

// `results.flows[step][link_index]` / `results.heads[step][node_index]` /
// `results.demands[step][node_index]` are indexed by report step.
println!("solved {} report steps", results.heads.len());
# Ok(())
# }
```

### Running a parallel simulation using the low-level API

The following example shows how to drive the solver directly. A leak-detection scenario is
simulated for every junction in parallel with [rayon](https://crates.io/crates/rayon), starting
from a shared, pre-solved initial state. For each junction an extra demand of 5 is added to
mimic a leak, the network is re-solved, and the resulting head at node `"1"` is recorded.
This pattern generalises to monte-carlo and fire-flow simulations on networks without tanks or
pressure controls.

```no_run
use epanet_rs::model::network::Network;
use epanet_rs::model::node::NodeType;
use epanet_rs::solver::hydraulicsolver::HydraulicSolver;
use epanet_rs::solver::state::SolverState;
use rayon::prelude::*;

# fn main() -> Result<(), Box<dyn std::error::Error>> {
let network = Network::from_file("tests/pump.inp")?;
let solver = HydraulicSolver::new(&network)?;

// Build the initial state and solve once. The resulting state is used as a
// warm start for every parallel scenario, which accelerates convergence.
let mut initial = SolverState::new_with_initial_values(&network);
initial.apply_patterns(&network, 0);
let initial = solver.solve(&network, &initial)?;

// Indices of all junctions (only junctions carry a demand that can leak).
let junctions: Vec<usize> = network.nodes.iter().enumerate()
    .filter(|(_, n)| matches!(n.node_type, NodeType::Junction(_)))
    .map(|(i, _)| i)
    .collect();

// Resolve the monitoring node once. Pressure = head - elevation.
let monitor = network.node_map["1"];
let monitor_elevation = network.nodes[monitor].elevation;

// For each junction, add a leak demand of 5, re-solve in parallel, and
// collect the resulting pressure at node "1".
let pressures: Vec<f64> = junctions.par_iter()
    .map(|&j| {
        let mut state = initial.clone();
        state.demands[j] += 5.0;
        let solved = solver.solve(&network, &state)?;
        Ok(solved.heads[monitor] - monitor_elevation)
    })
    .collect::<Result<_, epanet_rs::error::SolverError>>()?;

println!("simulated {} leak scenarios", pressures.len());
# Ok(())
# }
```

### Modifying and constructing networks

Networks can also be built and mutated programmatically through the strongly-typed `*Data` /
`*Update` structs re-exported from [`model::network`]. Topology changes (adding or removing
nodes/links, changing a link's endpoints) invalidate any cached solver; property-only changes
(e.g. pipe roughness, pump speed) are tracked incrementally and applied to the solver state when the network is solved.

```no_run
use epanet_rs::simulation::Simulation;
use epanet_rs::model::network::{Network, JunctionData, PipeData, ReservoirData, PipeUpdate};
use epanet_rs::model::link::LinkStatus;
use epanet_rs::model::options::HeadlossFormula;
use epanet_rs::model::units::FlowUnits;

# fn main() -> Result<(), Box<dyn std::error::Error>> {
let mut network = Network::new(FlowUnits::LPS, HeadlossFormula::DarcyWeisbach);

network.add_reservoir("R1", &ReservoirData {
    elevation: 100.0,
    ..Default::default()
})?;

network.add_junction("J1", &JunctionData {
    elevation: 50.0,
    basedemand: 1.0,
    ..Default::default()
})?;

network.add_pipe("P1", &PipeData {
    start_node: "R1".into(),
    end_node: "J1".into(),
    length: 1000.0,
    diameter: 200.0,
    roughness: 0.1,
    minor_loss: 0.0,
    check_valve: false,
    initial_status: LinkStatus::Open,
    vertices: None,
})?;

let mut simulation = Simulation::new(network);

// solve the network for a single time step
simulation.run_hydraulics()?;

// update the pipe roughness  
simulation.network.update_pipe("P1", &PipeUpdate {
    roughness: Some(0.2),
    ..Default::default()
})?;

// solve the network again
simulation.run_hydraulics()?;

# Ok(())
# }
```
*/
pub mod io;
pub mod model;
pub mod solver;
pub mod constants;
pub mod utils;
pub mod simulation;
pub mod error;
pub mod ffi;
