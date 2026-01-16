# EPANET-RS

An extremely fast, modern and safe re-implementation of the EPANET2 hydraulic solver, written in Rust.

## Why?

The EPANET2 solver has been the industry standard for hydraulic network simulation for decades due to its robustness, numerical stability, and extensive validation. Its algorithms and results are widely trusted in both academic and industrial applications.
The core EPANET2 codebase, however, is several decades old and written in C, making it difficult to maintain, extend, and optimize using modern software engineering practices. In addition, the original implementation predates modern CPU architectures and therefore does not fully exploit multi-core processors or SIMD (Single Instruction, Multiple Data) capabilities.
`epanet-rs` is a modern reimplementation of the EPANET2 hydraulic solver written in Rust, designed to preserve the original algorithms and numerical behavior while enabling safer memory management, improved maintainability, and performance optimizations through multi-threading and SIMD acceleration.

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

## Supported Network Elements

| Element | Status |
|---------|--------|
| Junctions | Supported |
| Reservoirs | Supported |
| Tanks | Not yet implemented |
| Pipes | Supported |
| Pumps | Supported (single-point curves) |
| Valves | Partial |

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
