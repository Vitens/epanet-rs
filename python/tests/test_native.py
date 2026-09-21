"""Tests for the native `epanet_rs.native` Rust API bindings."""

from __future__ import annotations

import math

import pytest
from epanet_rs import native

# --- Simulation / loading ---


def test_from_file_loads_network(inp_path):
    sim = native.Simulation.from_file(inp_path("2tanks.inp"))
    net = sim.network()
    assert net.node_count() > 0
    assert net.link_count() > 0


def test_from_file_missing_file_raises_input_error():
    with pytest.raises(native.InputError):
        native.Simulation.from_file("does-not-exist.inp")


def test_import_from_submodule_path(inp_path):
    # `epanet_rs.native` must be importable both as an attribute of
    # `epanet_rs` and as `from epanet_rs.native import ...` / `import
    # epanet_rs.native` (plain pyo3 `add_submodule` does not register
    # submodules in `sys.modules` by default).
    from epanet_rs.native import Simulation

    sim = Simulation.from_file(inp_path("2tanks.inp"))
    assert sim.network().node_count() > 0


# --- building a network from scratch ---


def test_build_network_from_scratch():
    sim = native.Simulation.new(native.FlowUnits.GPM, native.HeadlossFormula.HazenWilliams)
    net = sim.network()

    net.add_reservoir("R1", elevation=700.0)
    net.add_junction("J1", elevation=650.0, basedemand=100.0)
    net.add_junction("J2", elevation=640.0, basedemand=50.0)
    net.add_pipe("P1", start_node="R1", end_node="J1", length=1000.0, diameter=12.0, roughness=120.0)
    net.add_pipe("P2", start_node="J1", end_node="J2", length=500.0, diameter=8.0, roughness=120.0)

    assert net.node_count() == 3
    assert net.link_count() == 2

    result = sim.solve_hydraulics()
    assert len(result.heads) == 1
    assert len(result.heads[0]) == 3
    # heads should decrease downstream of the reservoir under positive demand
    heads = result.heads[0]
    assert heads[0] == pytest.approx(700.0)
    assert heads[1] < heads[0]


def test_add_tank_and_pattern_and_curve():
    sim = native.Simulation.new(native.FlowUnits.CFS, native.HeadlossFormula.HazenWilliams)
    net = sim.network()

    net.add_pattern("PAT1", multipliers=[1.0, 1.2, 0.8])
    net.add_junction("J1", elevation=100.0, basedemand=10.0, pattern="PAT1")
    net.add_reservoir("R1", elevation=150.0)
    net.add_tank(
        "T1",
        elevation=120.0,
        initial_level=5.0,
        min_level=0.0,
        max_level=20.0,
        diameter=50.0,
        coordinates=(1.0, 2.0),
    )
    net.add_curve("C1", x=[0.0, 500.0, 1000.0], y=[100.0, 80.0, 0.0])
    net.add_pipe("P1", start_node="R1", end_node="J1", length=1000.0, diameter=12.0, roughness=100.0)
    net.add_pump("PU1", start_node="R1", end_node="T1", head_curve_id="C1")

    assert net.has_tanks()

    node = net.get_node("J1")
    assert node.node_type == native.NodeType.Junction
    assert node.junction.demands[0].basedemand == pytest.approx(10.0)
    assert node.junction.demands[0].pattern == "PAT1"

    tank_node = net.get_node("T1")
    assert tank_node.node_type == native.NodeType.Tank
    assert tank_node.coordinates == (1.0, 2.0)
    assert tank_node.tank.max_level == pytest.approx(20.0)

    pump_link = net.get_link("PU1")
    assert pump_link.link_type == native.LinkType.Pump
    assert pump_link.pump.head_curve_id == "C1"

    curve = net.get_curve("C1")
    assert curve.x == [0.0, 500.0, 1000.0]
    assert curve.y == [100.0, 80.0, 0.0]

    pattern = net.get_pattern("PAT1")
    assert pattern.multipliers == [1.0, 1.2, 0.8]


def test_add_valve_and_update():
    sim = native.Simulation.new(native.FlowUnits.CFS, native.HeadlossFormula.HazenWilliams)
    net = sim.network()
    net.add_junction("J1", elevation=100.0)
    net.add_junction("J2", elevation=90.0)
    net.add_valve(
        "V1",
        start_node="J1",
        end_node="J2",
        diameter=6.0,
        valve_type=native.ValveType.PRV,
        setting=20.0,
    )

    valve = net.get_link("V1").valve
    assert valve.valve_type == native.ValveType.PRV
    assert valve.setting == pytest.approx(20.0)

    net.update_valve("V1", setting=30.0)
    assert net.get_link("V1").valve.setting == pytest.approx(30.0)


# --- update semantics (explicit clear_* flags instead of tri-state Optionals) ---


def test_update_junction_pattern_clear_and_set():
    sim = native.Simulation.new(native.FlowUnits.CFS, native.HeadlossFormula.HazenWilliams)
    net = sim.network()
    net.add_pattern("PAT1", multipliers=[1.0])
    net.add_junction("J1", elevation=100.0, basedemand=10.0, pattern="PAT1")

    # basedemand-only update must not disturb the existing pattern
    net.update_junction("J1", basedemand=20.0)
    node = net.get_node("J1")
    assert node.junction.demands[0].basedemand == pytest.approx(20.0)
    assert node.junction.demands[0].pattern == "PAT1"

    # explicit clear
    net.update_junction("J1", clear_pattern=True)
    assert net.get_node("J1").junction.demands[0].pattern is None


def test_update_pattern_unknown_id_raises_input_error():
    sim = native.Simulation.new(native.FlowUnits.CFS, native.HeadlossFormula.HazenWilliams)
    net = sim.network()
    with pytest.raises(native.InputError):
        net.update_pattern("missing", multipliers=[1.0])


# --- removal ---


def test_remove_node_cascades_links():
    sim = native.Simulation.new(native.FlowUnits.CFS, native.HeadlossFormula.HazenWilliams)
    net = sim.network()
    net.add_junction("J1", elevation=100.0)
    net.add_junction("J2", elevation=90.0)
    net.add_pipe("P1", start_node="J1", end_node="J2", length=100.0, diameter=12.0, roughness=100.0)

    net.remove_node("J2", unconditional=True)
    assert net.node_count() == 1
    assert net.link_count() == 0
    with pytest.raises(native.InputError):
        net.get_node("J2")


# --- unit conversion round-trips ---


def test_si_units_round_trip():
    sim = native.Simulation.new(native.FlowUnits.CMH, native.HeadlossFormula.HazenWilliams)
    net = sim.network()
    net.add_junction("J1", elevation=100.0)
    net.add_junction("J2", elevation=90.0)
    net.add_pipe(
        "P1",
        start_node="J1",
        end_node="J2",
        length=1000.0,
        diameter=300.0,  # mm
        roughness=120.0,
        minor_loss=2.0,
    )

    link = net.get_link("P1")
    assert link.pipe.diameter == pytest.approx(300.0)
    assert link.pipe.length == pytest.approx(1000.0)
    assert link.pipe.minor_loss == pytest.approx(2.0)

    node = net.get_node("J1")
    assert node.elevation == pytest.approx(100.0)


def test_us_units_pipe_diameter_in_inches():
    sim = native.Simulation.new(native.FlowUnits.CFS, native.HeadlossFormula.HazenWilliams)
    net = sim.network()
    net.add_junction("J1", elevation=0.0)
    net.add_junction("J2", elevation=0.0)
    net.add_pipe("P1", start_node="J1", end_node="J2", length=1000.0, diameter=12.0, roughness=100.0)

    assert net.get_link("P1").pipe.diameter == pytest.approx(12.0)


# --- solving / state ---


def test_step_by_step_hydraulics_and_state(inp_path):
    sim = native.Simulation.from_file(inp_path("2tanks.inp"))
    assert sim.state() is None
    assert not sim.solved

    sim.initialize_hydraulics()
    sim.run_hydraulics()

    assert sim.solved
    state = sim.state()
    assert state is not None
    assert all(math.isfinite(h) for h in state.heads)
    assert len(state.statuses) == sim.network().link_count()


def test_run_hydraulics_without_init_raises_solver_error():
    sim = native.Simulation.new(native.FlowUnits.CFS, native.HeadlossFormula.HazenWilliams)
    with pytest.raises(native.SolverError):
        sim.run_hydraulics()


def test_solve_hydraulics_parallel_matches_sequential(inp_path):
    sequential = native.Simulation.from_file(inp_path("2tanks.inp")).solve_hydraulics(parallel=False)
    # 2tanks.inp has tanks, so parallel mode silently falls back to sequential.
    parallel = native.Simulation.from_file(inp_path("2tanks.inp")).solve_hydraulics(parallel=True)
    assert sequential.heads == parallel.heads


# --- enums ---


def test_enum_equality_and_int_conversion():
    assert native.LinkStatus.Open == native.LinkStatus.Open
    assert native.LinkStatus.Open != native.LinkStatus.Closed
    assert isinstance(int(native.FlowUnits.CFS), int)


# --- lazy views over nodes/links/patterns/curves/controls ---


def _sample_network():
    sim = native.Simulation.new(native.FlowUnits.CFS, native.HeadlossFormula.HazenWilliams)
    net = sim.network()
    net.add_pattern("PAT1", multipliers=[1.0])
    net.add_curve("C1", x=[0.0, 1.0], y=[1.0, 0.0])
    net.add_junction("J1", elevation=100.0)
    net.add_junction("J2", elevation=90.0)
    net.add_junction("J3", elevation=80.0)
    net.add_pipe("P1", start_node="J1", end_node="J2", length=100.0, diameter=12.0, roughness=100.0)
    net.add_pipe("P2", start_node="J2", end_node="J3", length=100.0, diameter=12.0, roughness=100.0)
    return net


def test_node_view_len_and_getitem_do_not_require_iteration():
    net = _sample_network()
    view = net.nodes()

    assert len(view) == 3
    assert isinstance(view[0], native.Node)
    # negative indexing, like a normal Python sequence
    assert view[-1].id == view[2].id
    with pytest.raises(IndexError):
        view[3]
    with pytest.raises(IndexError):
        view[-4]


def test_node_view_iteration_matches_node_map_keys():
    # node_map()/link_map()/... are hash maps: their *key set* matches the
    # view's ids, but iteration order is unspecified and need not match
    # nodes()'s storage order — use the view (or `net.get_node(id)`) when
    # order matters, and the map only for lookups/existence checks.
    net = _sample_network()
    assert {n.id for n in net.nodes()} == set(net.node_map())
    assert {l.id for l in net.links()} == set(net.link_map())
    assert {p.id for p in net.patterns()} == set(net.pattern_map())
    assert {c.id for c in net.curves()} == set(net.curve_map())


def test_node_view_iterator_is_lazy_and_resumable():
    net = _sample_network()
    it = iter(net.nodes())

    first = next(it)
    assert first.id == net.nodes()[0].id
    # the iterator holds its own position independent of a fresh view/call
    second = next(it)
    assert second.id == net.nodes()[1].id

    # a *second*, independent iterator starts from the beginning again
    fresh = list(net.nodes())
    assert [n.id for n in fresh] == [n.id for n in net.nodes()]


def test_node_view_supports_islice_without_building_everything():
    import itertools

    net = _sample_network()
    first_two = [n.id for n in itertools.islice(net.nodes(), 2)]
    assert first_two == [net.nodes()[0].id, net.nodes()[1].id]


def test_link_view_and_control_view():
    net = _sample_network()
    link_view = net.links()
    assert len(link_view) == 2
    assert link_view[0].pipe is not None

    control_view = net.controls()
    assert len(control_view) == 0
    assert list(control_view) == []


def test_views_reflect_mutation_after_they_were_obtained():
    # views hold a reference to the simulation, not a snapshot, so a view
    # obtained before a mutation still reflects it afterwards.
    net = _sample_network()
    view = net.nodes()
    assert len(view) == 3
    net.add_junction("J4", elevation=70.0)
    assert len(view) == 4


# --- id -> index maps (mirror Network::node_map/link_map/...) ---


def test_node_map_gives_index_into_the_view():
    net = _sample_network()
    node_map = net.node_map()

    assert "J1" in node_map  # O(1) existence check, no Node built
    assert "missing" not in node_map

    view = net.nodes()
    for node_id, index in node_map.items():
        assert view[index].id == node_id


def test_map_correlates_ids_with_solver_state_arrays(inp_path):
    # SolverState/SolverResult arrays are indexed the same way as
    # Network.nodes/links, which is exactly what node_map()/link_map() are
    # for: going from an id to a position in those flat arrays without
    # building any Node/Link snapshot.
    sim = native.Simulation.from_file(inp_path("2tanks.inp"))
    net = sim.network()
    node_map = net.node_map()
    link_map = net.link_map()

    sim.initialize_hydraulics()
    sim.run_hydraulics()
    state = sim.state()

    assert len(state.heads) == len(node_map)
    assert len(state.flows) == len(link_map)
    # the index for an id must match the same id's position in the view
    nodes = net.nodes()
    for node_id, index in node_map.items():
        assert nodes[index].id == node_id
        assert math.isfinite(state.heads[index])


def test_node_map_getitem_and_get():
    net = _sample_network()
    node_map = net.node_map()

    assert node_map["J1"] == 0
    with pytest.raises(KeyError):
        node_map["missing"]

    assert node_map.get("J1") == 0
    assert node_map.get("missing") is None
    assert node_map.get("missing", -1) == -1  # negative sentinels work


def test_node_map_converts_to_a_plain_dict():
    net = _sample_network()
    snapshot = dict(net.node_map())
    assert snapshot == {"J1": 0, "J2": 1, "J3": 2}
    # an ordinary dict, fully independent of the network from here on
    net.remove_node("J1", unconditional=True)
    assert snapshot == {"J1": 0, "J2": 1, "J3": 2}


def test_node_map_is_live_like_the_views():
    # a *held* node_map()/nodes() both keep reflecting the network, unlike
    # a `dict(net.node_map())` snapshot (see test above).
    net = _sample_network()
    held_map = net.node_map()
    held_view = net.nodes()
    assert len(held_map) == 3
    assert len(held_view) == 3
    assert "J4" not in held_map

    net.add_junction("J4", elevation=70.0)

    assert len(held_map) == 4  # live: reflects the new junction
    assert len(held_view) == 4  # live: reflects the new junction
    assert "J4" in held_map
    assert held_map["J4"] == 3


def test_node_map_indices_shift_after_removal():
    # indices shift the moment they change (remove_node/remove_link use
    # swap-remove internally) — even a *held* node_map() reflects this
    # immediately, since it's live, but an index read *before* the removal
    # can point at the wrong node afterwards.
    net = _sample_network()
    held_map = net.node_map()
    assert held_map["J1"] == 0

    net.remove_node("J1", unconditional=True)

    assert "J1" not in held_map  # live: the removal is visible immediately
    assert set(held_map.keys()) == {"J2", "J3"}
