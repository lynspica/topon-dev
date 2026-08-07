"""Unit tests for the periodic cell recorded by the topology generators.

Background
----------
``infer_dims_from_graph`` used to derive the box from the node positions
alone, as ``max - min + 1``. That is exact for simple cubic (sites on
integer coordinates, unit spacing) but wrong for every lattice with
fractional basis sites: BCC and FCC body/face sites sit at +0.5 and
Diamond sites at quarter-cell offsets, so the coordinates stop short of
the cell edge and the estimate overshoots.

The overshoot is not cosmetic. The box feeds every minimum-image
calculation in the pipeline, so an inflated box makes a large fraction of
edges resolve to the wrong periodic replica and get built at twice their
true bond length (measured on 4x4x4: 169/512 BCC edges, 360/1536 FCC).

Generators now record the exact cell in ``G.graph["box"]`` and the
positional estimate remains only as a fallback for graphs written before
they did.
"""

import networkx as nx
import numpy as np
import pytest

from topon.topology.generator_python import PythonTopologyGenerator
from topon.topology.generator_python_diamond import create_diamond_lattice
from topon.topology.loader import (
    format_box_header,
    infer_dims_from_graph,
    load_graph,
    read_box_header,
    remove_vacancies,
    save_nodes_edges,
)


class _Cfg:
    """Minimal stand-in for a GeneratorConfig."""

    def __init__(self, lattice_type, lattice_size=(4, 4, 4)):
        self.lattice_type = lattice_type
        self.lattice_size = lattice_size
        self.max_functionality = 6
        self.degree_distribution = ""


def _build(lattice_type, dims=(4, 4, 4)):
    gen = PythonTopologyGenerator(_Cfg(lattice_type, dims))
    return gen._create_lattice(dims, lattice_type)


def _edge_lengths(G, box):
    """Minimum-image edge lengths under ``box``."""
    box = np.asarray(box, dtype=float)
    pos = {n: np.asarray(d["pos"], dtype=float) for n, d in G.nodes(data=True)}
    out = []
    for u, v in G.edges():
        d = pos[u] - pos[v]
        out.append(np.linalg.norm(d - box * np.round(d / box)))
    return np.asarray(out)


# ---------------------------------------------------------------------------
# The recorded cell
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("lattice_type", ["SC", "BCC", "FCC"])
def test_generator_records_true_cell(lattice_type):
    """The cell is the lattice repeat, not the extent of the coordinates."""
    G = _build(lattice_type, (4, 4, 4))
    assert G.graph["box"] == (4.0, 4.0, 4.0)
    assert np.allclose(infer_dims_from_graph(G), [4.0, 4.0, 4.0])


def test_diamond_records_true_cell():
    G = create_diamond_lattice(3, 3, 3)
    assert G.graph["box"] == (3.0, 3.0, 3.0)
    assert np.allclose(infer_dims_from_graph(G), [3.0, 3.0, 3.0])


def test_non_cubic_cell_is_recorded_per_axis():
    """A 3x4x5 lattice must not collapse to a single box length."""
    G = _build("BCC", (3, 4, 5))
    assert np.allclose(infer_dims_from_graph(G), [3.0, 4.0, 5.0])


# ---------------------------------------------------------------------------
# The bug this prevents: wrong periodic replica, doubled bond lengths
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "lattice_type,expected_nn",
    [("SC", 1.0), ("BCC", np.sqrt(3) / 2), ("FCC", np.sqrt(2) / 2)],
)
def test_all_edges_are_nearest_neighbours(lattice_type, expected_nn):
    """Every edge of a full lattice spans exactly one nearest-neighbour bond.

    This is the assertion that fails under the old heuristic: an inflated
    box sends some edges to the wrong periodic image, where they measure
    twice the true bond length.
    """
    G = _build(lattice_type, (4, 4, 4))
    lengths = _edge_lengths(G, infer_dims_from_graph(G))
    assert np.allclose(lengths, expected_nn), (
        f"{lattice_type}: {np.sum(~np.isclose(lengths, expected_nn))} of "
        f"{len(lengths)} edges are not nearest-neighbour bonds "
        f"(max {lengths.max():.4f}, expected {expected_nn:.4f})"
    )


@pytest.mark.parametrize("lattice_type", ["BCC", "FCC"])
def test_positional_fallback_would_be_wrong(lattice_type):
    """Pin the old behaviour so the regression is detectable, not silent.

    Without a recorded box the estimate inflates by half a cell and edges
    start landing on the wrong replica. If a future change makes the
    fallback correct for these lattices, this test fails and should be
    deleted rather than worked around.
    """
    G = _build(lattice_type, (4, 4, 4))
    del G.graph["box"]

    fallback = infer_dims_from_graph(G)
    assert np.allclose(fallback, [4.5, 4.5, 4.5])

    lengths = _edge_lengths(G, fallback)
    assert not np.allclose(lengths, lengths.min())


def test_sc_is_unaffected_by_the_change():
    """SC must be byte-for-byte what it was: the fallback was exact there.

    Simple cubic is the lattice the frozen regression outputs were built
    on, so the recorded box has to agree with what the estimate produced.
    """
    G = _build("SC", (5, 5, 5))
    with_box = infer_dims_from_graph(G)

    del G.graph["box"]
    fallback = infer_dims_from_graph(G)

    assert np.allclose(with_box, fallback)
    assert np.allclose(with_box, [5.0, 5.0, 5.0])


# ---------------------------------------------------------------------------
# The cell survives the transforms the pipeline applies
# ---------------------------------------------------------------------------

def test_cell_survives_multigraph_and_vacancy_removal():
    """pipeline._generate_topology does both before reading dims."""
    G = nx.MultiGraph(_build("FCC", (4, 4, 4)))
    G.add_node(99_999, pos=(0.0, 0.0, 0.0))  # vacancy: degree 0
    assert remove_vacancies(G) == 1
    assert np.allclose(infer_dims_from_graph(G), [4.0, 4.0, 4.0])


def test_cell_survives_sculpting():
    """run_single_trial copies the base graph; graph attrs must come along."""
    gen = PythonTopologyGenerator(_Cfg("BCC", (3, 3, 3)))
    graphs = gen.generate(trials=5, max_saves=1)
    assert graphs, "sculpting produced no graph"
    assert np.allclose(infer_dims_from_graph(graphs[0]), [3.0, 3.0, 3.0])


# ---------------------------------------------------------------------------
# .nodes round-trip
# ---------------------------------------------------------------------------

def test_nodes_round_trip_preserves_cell(tmp_path):
    G = _build("BCC", (4, 4, 4))
    nodes_p, edges_p = tmp_path / "n.nodes", tmp_path / "n.edges"
    save_nodes_edges(G, nodes_p, edges_p)

    assert read_box_header(nodes_p) == (4.0, 4.0, 4.0)

    loaded, dims = load_graph(nodes_path=nodes_p, edges_path=edges_p)
    assert np.allclose(dims, [4.0, 4.0, 4.0])
    assert np.allclose(_edge_lengths(loaded, dims), np.sqrt(3) / 2)


def test_nodes_without_box_header_still_load(tmp_path):
    """Backwards compatibility: existing .nodes files have no BOX line.

    The frozen regression inputs are exactly this shape, so the fallback
    has to keep working unchanged.
    """
    nodes_p, edges_p = tmp_path / "n.nodes", tmp_path / "n.edges"
    nodes_p.write_text(
        "# NodeID X Y Z Degree\n"
        "0 0.000000 0.000000 0.000000 1\n"
        "1 1.000000 0.000000 0.000000 1\n"
    )
    edges_p.write_text("0 1\n")

    assert read_box_header(nodes_p) is None

    G, dims = load_graph(nodes_path=nodes_p, edges_path=edges_p)
    assert "box" not in G.graph
    assert np.allclose(dims, [2.0, 1.0, 1.0])  # the max-min+1 fallback


def test_malformed_box_header_falls_back(tmp_path):
    nodes_p, edges_p = tmp_path / "n.nodes", tmp_path / "n.edges"
    nodes_p.write_text(
        "# BOX not a number\n"
        "# NodeID X Y Z Degree\n"
        "0 0.000000 0.000000 0.000000 1\n"
        "1 1.000000 0.000000 0.000000 1\n"
    )
    edges_p.write_text("0 1\n")

    assert read_box_header(nodes_p) is None
    _, dims = load_graph(nodes_path=nodes_p, edges_path=edges_p)
    assert np.allclose(dims, [2.0, 1.0, 1.0])


def test_non_positive_box_header_rejected(tmp_path):
    nodes_p = tmp_path / "n.nodes"
    nodes_p.write_text("# BOX 4 0 4\n# NodeID X Y Z Degree\n0 0.0 0.0 0.0 0\n")
    assert read_box_header(nodes_p) is None


def test_box_header_format_is_stable():
    """The C generator has to emit this exact spelling to stay in lockstep."""
    assert format_box_header((6, 6, 6)) == "# BOX 6 6 6"
    assert format_box_header((3.0, 4.0, 5.5)) == "# BOX 3 4 5.5"


# ---------------------------------------------------------------------------
# gpickle path
# ---------------------------------------------------------------------------

def test_recorded_cell_overrides_stale_stored_dims(tmp_path):
    """A gpickle saved with a fallback-derived dims must not win.

    ``save_graph`` stores ``(graph, dims)``. If that dims came from the
    old estimate it is half a cell too large, and the graph's own box is
    the trustworthy value.
    """
    import pickle

    G = nx.MultiGraph(_build("BCC", (4, 4, 4)))
    path = tmp_path / "g.gpickle"
    with open(path, "wb") as f:
        pickle.dump((G, np.array([4.5, 4.5, 4.5])), f)  # stale dims

    _, dims = load_graph(gpickle_path=path)
    assert np.allclose(dims, [4.0, 4.0, 4.0])


def test_gpickle_without_box_keeps_stored_dims(tmp_path):
    import pickle

    G = nx.MultiGraph(_build("SC", (4, 4, 4)))
    del G.graph["box"]
    path = tmp_path / "g.gpickle"
    with open(path, "wb") as f:
        pickle.dump((G, np.array([4.0, 4.0, 4.0])), f)

    _, dims = load_graph(gpickle_path=path)
    assert np.allclose(dims, [4.0, 4.0, 4.0])
