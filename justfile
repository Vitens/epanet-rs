# List available recipes
default:
    @just --list

# Run all tests
test:
    cargo test --all-features --workspace

# Generate and display coverage summary
coverage:
    #!/usr/bin/env bash
    if ! command -v cargo-llvm-cov &> /dev/null; then
        echo "cargo-llvm-cov not found. Installing..."
        cargo install cargo-llvm-cov
    fi
    cargo llvm-cov --all-features --workspace

# Generate HTML coverage report
coverage-html:
    #!/usr/bin/env bash
    if ! command -v cargo-llvm-cov &> /dev/null; then
        echo "cargo-llvm-cov not found. Installing..."
        cargo install cargo-llvm-cov
    fi
    cargo llvm-cov --all-features --workspace --html
    echo "Coverage report generated at target/llvm-cov/html/index.html"

# Clean coverage artifacts
clean:
    cargo llvm-cov clean --workspace
    rm -f lcov.info
