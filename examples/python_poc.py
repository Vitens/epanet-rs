#!/usr/bin/env python3
"""Proof-of-concept: load an INP file, solve hydraulics, inspect Polars results."""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl
import matplotlib.pyplot as plt

from epanet_rs import Simulation


def main() -> None:

    simulation = Simulation("benchmarks/L-TOWN.inp")
    results = simulation.solve()

    df = results.dataframe()

    d = df.filter(pl.col("element_id") == "T1")

    plt.plot(d["time"], d["head"])
    plt.show()

if __name__ == "__main__":
    main()
