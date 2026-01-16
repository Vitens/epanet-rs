# EPANET-RS

An extremely fast, modern and safe re-implementation of the EPANET2 hydraulic solver, written in Rust.

## Background

The EPANET2 solver has been the industry standard for hydraulic network simulation for decades due to its robustness, numerical stability, and extensive validation. Its algorithms and results are widely trusted in both academic and industrial applications.
The core EPANET2 codebase, however, is several decades old and written in C, making it difficult to maintain, extend, and optimize using modern software engineering practices. In addition, the original implementation predates modern CPU architectures and therefore does not fully exploit multi-core processors or SIMD (Single Instruction, Multiple Data) capabilities.

Modern applications of the EPANET solver, such as monte-carlo simulations, leak detection algorithms and real-time digital twins of huge networks require a more modern implementation, with better performance and maintainability.

`epanet-rs` is a modern reimplementation of the EPANET2 hydraulic solver written in Rust, designed to preserve the original algorithms and numerical behavior while enabling safer memory management, improved maintainability, and performance optimizations through multi-threading and SIMD acceleration.  

`epanet-rs` runs about 3 times faster on a modern CPU then the original EPANET2_3 solver in sequential mode, and 8 times faster in parallel mode for extended period simulations for supported networks (no tanks/controls)!

## Design Goals

- **Numerical Parity** with EPANET2_3 solver, ensuring identical or equivalent results for the same input
- **High Performance** through multi-threading, SIMD acceleration and a modern, [faer](https://crates.io/crates/faer) based sparse solver
- **Parallelization** of the solver loop by rewriting the solver to return vectors of heads and flows instead of in-place assignment
- **Memory Safety** through Rust's ownership and borrowing system

## Usage

```bash
# Run simulation
epanet-rs <network_file.inp>

# Run with output file (JSON or MessagePack)
epanet-rs <network_file.inp> output.json

# Run with verbose output
epanet-rs <network_file.inp> --verbose

# Run with parallel solving (for extended period simulations)
epanet-rs <network_file.inp> --parallel

# Print results to console
epanet-rs <network_file.inp> --print-results
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
- [-] Parallel Hydraulic Solver for extended period simulations
- [ ] Quality simulations
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

See [LICENSE](LICENSE) file.

## Author

Abel Heinsbroek (Vitens N.V.)
