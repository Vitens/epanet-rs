"""Shared pytest fixtures for the epanet-rs Python binding tests."""

from __future__ import annotations

from pathlib import Path

import pytest

# Repo-level `tests/` directory (shared .inp fixtures used by the Rust test
# suite too), two levels up from `python/tests/`.
REPO_TESTS_DIR = Path(__file__).resolve().parents[2] / "tests"


@pytest.fixture
def inp_path():
    """Returns a function that resolves a fixture `.inp` filename to a path
    in the repo-level `tests/` directory."""

    def _resolve(name: str) -> str:
        path = REPO_TESTS_DIR / name
        assert path.exists(), f"missing fixture INP file: {path}"
        return str(path)

    return _resolve
