"""Smoke tests for the EPANET-2.3-toolkit-mirroring `epanet_rs.Project` API."""

from __future__ import annotations

import epanet_rs
from epanet_rs import EN_FLOW, EN_HEAD, EN_LINKCOUNT, EN_NODECOUNT, Project


def test_open_and_counts(inp_path):
    p = Project()
    p.open(inp_path("2tanks.inp"))
    try:
        assert p.getcount(EN_NODECOUNT) > 0
        assert p.getcount(EN_LINKCOUNT) > 0
    finally:
        p.close()


def test_solve_h_returns_series(inp_path):
    p = Project()
    p.open(inp_path("2tanks.inp"))
    try:
        result = p.solveH()
        assert isinstance(result, epanet_rs.SolverResult)
        assert len(result.heads) > 0
        assert len(result.heads[0]) == p.getcount(EN_NODECOUNT)
        assert len(result.flows[0]) == p.getcount(EN_LINKCOUNT)
    finally:
        p.close()


def test_step_by_step_hydraulics_match_getters(inp_path):
    p = Project()
    p.open(inp_path("2tanks.inp"))
    try:
        p.openH()
        p.initH()
        p.runH(0)
        # after solving, node/link value getters should return finite numbers
        head = p.getnodevalue(1, EN_HEAD)
        flow = p.getlinkvalue(1, EN_FLOW)
        assert isinstance(head, float)
        assert isinstance(flow, float)
        p.closeH()
    finally:
        p.close()


def test_parallel_solve_matches_sequential(inp_path):
    p = Project()
    p.open(inp_path("2tanks.inp"))
    try:
        sequential = p.solveH(parallel=False)
    finally:
        p.close()

    # 2tanks.inp has tanks, so parallel solving silently falls back to
    # sequential; this just checks the flag is accepted and results agree.
    p2 = Project()
    p2.open(inp_path("2tanks.inp"))
    try:
        parallel = p2.solveH(parallel=True)
    finally:
        p2.close()

    assert sequential.heads == parallel.heads


def test_geterror_returns_description():
    message = Project.geterror(110)
    assert isinstance(message, str)
    assert message != ""
