# EPANET-RS
A fast, modern and safe re-implementation of the EPANET2 hydraulic solver, written in Rust.

[![Rust](https://github.com/vitens/epanet-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/AbelHeinsbroek/epanet-rs/actions/workflows/rust.yml)
[![crate](https://img.shields.io/crates/v/epanet-rs.svg)](https://crates.io/crates/epanet-rs)

Dual-licensed under [Apache 2.0](LICENSE-APACHE) or [MIT](LICENSE-MIT).

## Background

The EPANET2 solver has been the industry standard for hydraulic network simulation for decades due to its robustness, numerical stability, and extensive validation. Its algorithms and results are widely trusted in both academic and industrial applications.
The core EPANET2 codebase, however, is several decades old and written in C, making it difficult to maintain, extend, and optimize using modern software engineering practices. In addition, the original implementation predates modern CPU architectures and therefore does not fully exploit multi-core processors or SIMD (Single Instruction, Multiple Data) capabilities.

Modern applications of the EPANET solver, such as monte-carlo simulations, leak detection algorithms and real-time digital twins of huge networks require a more modern implementation, with better performance and maintainability.

`epanet-rs` is a modern reimplementation of the EPANET2 hydraulic solver written in Rust, designed to preserve the original algorithms and numerical behavior while enabling safer memory management, improved maintainability, and performance optimizations through multi-threading and SIMD acceleration.  

`epanet-rs` runs about as fast as the original EPANET2_3 solver in sequential mode, and up to 5 times faster in parallel mode for extended period simulations for supported networks (no tanks/controls)!

## Design Goals

- **Numerical Parity** with EPANET2_3 solver (hyd_2_3 branch), ensuring identical or equivalent results for the same input
- **High Performance** through multi-threading, SIMD acceleration and a modern, [faer](https://crates.io/crates/faer) based sparse solver
- **Parallelization** of the solver loop by rewriting the solver to return vectors of heads and flows instead of in-place assignment
- **Memory Safety** through Rust's ownership and borrowing system
- **Modern API** with a focus on ease of use and parallelization
- **Backwards Compatibility** with EPANET2_3 network API methods

## Usage

```bash
# Run simulation
epanet-rs run <network_file.inp>

# Run with output file (.json format)
epanet-rs run <network_file.inp> output.json

# Run with parallel solving (for extended period simulations without tanks/controls)
epanet-rs run <network_file.inp> < --parallel

# Convert network to different format (JSON or MessagePack)
epanet-rs convert <network_file.inp> output.json
epanet-rs convert <network_file.inp> output.msgpack

# Validate results against EPANET
epanet-rs validate <network_file.inp>
```

## Building

```bash
# Debug build
cargo build

# Release build (optimized)
cargo build --release
```

## Testing

```bash
# Run all tests
cargo test

# Run solver tests
cargo test --test solver_test
```

## Supported Features and To-Do

- [x] INP file support
- [x] Parallel Hydraulic Solver for extended period simulations
- [ ] Quality simulations
- [ ] Pressure dependent demand simulation
- [ ] CONTROLS and RULES
- [ ] EMITTERS and LEAKAGE
- [ ] ENERGY
- [ ] Backwards compatibility with EPANET2_3 network API methods

## Dependencies

- [faer](https://crates.io/crates/faer) - Sparse linear algebra
- [rayon](https://crates.io/crates/rayon) - Parallel iteration
- [hashbrown](https://crates.io/crates/hashbrown) - Fast hash maps
- [serde](https://crates.io/crates/serde) - Serialization
- [clap](https://crates.io/crates/clap) - Command line parsing


## License

EPANET-RS is licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Acknowledgments

EPANET-RS is heavily inspired by [EPANET 2.3](https://github.com/OpenWaterAnalytics/EPANET) and implements the same hydraulic algorithms and models, with slight modifications to improve performance and memory safety. We acknowledge the original EPANET developers, the [OpenWaterAnalytics](https://github.com/OpenWaterAnalytics) team and the broader water distribution modeling community for their decades of work establishing these foundations.

## Author

Abel Heinsbroek (Vitens N.V.)
