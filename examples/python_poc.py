#!/usr/bin/env python3
"""Proof-of-concept: load an INP file, solve hydraulics, inspect Polars results."""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

from epanet_rs import Simulation


def main() -> None:
    inp_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("tests/pump.inp")

    simulation = Simulation(str(inp_path))
    results = simulation.solve()

    df = results.dataframe()
    last_time = df["time"].max()

    print(f"Loaded and solved: {inp_path}")
    print(f"Report steps: {results.report_steps}")
    print("\nNodes (last report step):")
    print(
        df.filter(pl.col("time") == last_time, pl.col("element_type") == "node")
        .sort("element_id")
        .head(10)
    )
    print("\nLinks (last report step):")
    print(
        df.filter(pl.col("time") == last_time, pl.col("element_type") == "link")
        .sort("element_id")
        .head(10)
    )


if __name__ == "__main__":
    main()
