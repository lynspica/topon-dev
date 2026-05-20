"""Regression guard: npz_writer edge_index must be 0-based row positions.

The bug: ``write_npz`` stored *original simulation node IDs* in
``edge_index`` -- chain dual-ids offset above ``max_xlink_id`` and
sparse crosslink node-ids -- instead of 0-based positions into
``node_features``. With vacancy removal leaving ID gaps, the raw-ID
range exceeds N, so PyTorch Geometric throws CUDA "index out of bounds".

These tests build a graph with deliberately SPARSE, non-contiguous node
IDs (contiguous IDs would mask the bug because id == position by
coincidence) and assert:
  1. every edge_index value is a valid 0-based row position [0, N)
  2. each chemical edge connects a chain row (type 0) to a crosslink
     row (type 1) -- i.e. the remap points at the *correct* rows
  3. write_npz -> load_npz round-trips chains + entanglements
"""
from __future__ import annotations

import networkx as nx
import numpy as np

from topon.topology.loader import load_npz
from topon.writers.npz_writer import write_npz


def _sparse_graph() -> nx.MultiGraph:
    """MultiGraph with non-contiguous node IDs (0, 5, 10, 42).

    In the buggy writer this yields node_ids = [43,44,45,46, 0,5,10,42]
    (N=8) but edge_index values up to 46 -- well past N-1.
    """
    G = nx.MultiGraph()
    for nid, pos in [(0, (0.0, 0.0, 0.0)), (5, (1.0, 0.0, 0.0)),
                     (10, (0.0, 1.0, 0.0)), (42, (1.0, 1.0, 1.0))]:
        G.add_node(nid, pos=pos)
    G.add_edge(0, 5, dp=10)
    G.add_edge(5, 10, dp=12)
    G.add_edge(10, 42, dp=8)
    G.add_edge(0, 42, dp=15)
    # one entanglement pair: chain (0,5,0) <-> chain (10,42,0)
    G[0][5][0]["entangled_with"] = (10, 42, 0)
    G[0][5][0]["entanglement_count"] = 1
    G[10][42][0]["entangled_with"] = (0, 5, 0)
    G[10][42][0]["entanglement_count"] = 1
    return G


def test_edge_index_is_zero_based_positions(tmp_path):
    """edge_index values must all be valid 0-based rows into node_features."""
    G = _sparse_graph()
    out = tmp_path / "net.npz"
    write_npz(G, str(out), dp=10, dims=np.array([2.0, 2.0, 2.0]))

    data = np.load(out)
    node_features = data["node_features"]
    edge_index = data["edge_index"]
    n_nodes = node_features.shape[0]

    assert edge_index.shape[1] > 0, "expected edges in the test graph"
    assert edge_index.min() >= 0
    assert edge_index.max() < n_nodes, (
        f"edge_index max {int(edge_index.max())} >= N {n_nodes} -- "
        f"edge_index is still in original-ID space (the bug)"
    )


def test_edge_index_points_at_correct_rows(tmp_path):
    """Each chemical edge connects a chain row (type 0) to a crosslink
    row (type 1) -- proves the remap targets the RIGHT rows, not just
    in-range ones."""
    G = _sparse_graph()
    out = tmp_path / "net.npz"
    write_npz(G, str(out), dp=10, dims=np.array([2.0, 2.0, 2.0]))

    data = np.load(out)
    types = data["node_features"][:, 0]   # column 0 = type (0 chain / 1 xlink)
    edge_index = data["edge_index"]
    edge_type = data["edge_type"]

    for k in range(edge_index.shape[1]):
        if int(edge_type[k]) != 0:        # chemical edges only
            continue
        a, b = int(edge_index[0, k]), int(edge_index[1, k])
        pair = {types[a], types[b]}
        assert pair == {0.0, 1.0}, (
            f"chemical edge rows {a},{b} have types {pair}; expected one "
            f"chain (0) + one crosslink (1)"
        )


def test_npz_round_trip_through_loader(tmp_path):
    """write_npz -> load_npz reconstructs the same chains + entanglement."""
    G = _sparse_graph()
    out = tmp_path / "net.npz"
    write_npz(G, str(out), dp=10, dims=np.array([2.0, 2.0, 2.0]))

    G2, _ = load_npz(out)

    assert set(G2.nodes()) == set(G.nodes())
    assert G2.number_of_edges() == G.number_of_edges()

    chains_in = sorted((min(u, v), max(u, v), d.get("dp"))
                       for u, v, _, d in G.edges(keys=True, data=True))
    chains_out = sorted((min(u, v), max(u, v), d.get("dp"))
                        for u, v, _, d in G2.edges(keys=True, data=True))
    assert chains_in == chains_out

    n_ent_out = sum(1 for _, _, _, d in G2.edges(keys=True, data=True)
                    if d.get("entangled_with")) // 2
    assert n_ent_out == 1
