"""The entangled-edge realization helper, both methods.

The pipeline and both canonical workflow modules delegate here, so what
these tests pin down is the contract they all rely on: one path per
entangled edge, one bead per chain atom, junctions excluded, and the
waypoint pair actually wound rather than merely bent.
"""
import networkx as nx
import numpy as np
import pytest

from topon.config.schema import EntanglementsConfig
from topon.conformation.entanglement.realize import entangled_backbone_paths

DIMS = np.array([4.0, 4.0, 4.0])


def _pair_graph():
    """Two disjoint parallel chords one lattice unit apart, entangled."""
    G = nx.MultiGraph()
    G.add_node("a0", pos=(1.0, 1.0, 1.0))
    G.add_node("a1", pos=(3.0, 1.0, 1.0))
    G.add_node("b0", pos=(1.0, 2.0, 1.0))
    G.add_node("b1", pos=(3.0, 2.0, 1.0))
    ka = G.add_edge("a0", "a1")
    kb = G.add_edge("b0", "b1")
    ea, eb = ("a0", "a1", ka), ("b0", "b1", kb)
    G.edges[ea]["entangled_with"] = eb
    G.edges[ea]["entanglement_count"] = 2
    G.edges[eb]["entangled_with"] = ea
    G.edges[eb]["entanglement_count"] = 2
    return G, ea, eb


def _length(path):
    p = np.asarray(path, float)
    return float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum())


@pytest.mark.parametrize("method", ["waypoint", "kink"])
def test_counts_and_interior(method):
    G, ea, eb = _pair_graph()
    atoms = {ea: list(range(20)), eb: list(range(30))}
    out = entangled_backbone_paths(G, DIMS, atoms, method=method,
                                   kink_params={"overshoot": 0.2,
                                                "z_amp": 0.5,
                                                "sigma": 0.15})
    assert set(out) == {ea, eb}
    assert len(out[ea]) == 20
    assert len(out[eb]) == 30
    # Junctions are not in the path: the first bead sits off the node.
    assert np.linalg.norm(np.asarray(out[ea][0])
                          - np.array([1.0, 1.0, 1.0])) > 1e-3


def test_waypoint_pair_is_wound():
    G, ea, eb = _pair_graph()
    atoms = {ea: list(range(40)), eb: list(range(40))}
    out = entangled_backbone_paths(G, DIMS, atoms, method="waypoint")

    # Wound, not straight: each path is longer than its 2.0 chord.
    assert _length(out[ea]) > 2.2
    assert _length(out[eb]) > 2.2

    # Wound about each other: the two paths interleave across the midline
    # between the chords (y = 1.5), which a kink pointed at the partner
    # never does from both sides.
    ya = np.asarray(out[ea], float)[:, 1]
    yb = np.asarray(out[eb], float)[:, 1]
    assert ya.max() > 1.5 and yb.min() < 1.5


def test_untangled_edges_untouched():
    G, ea, eb = _pair_graph()
    kc = G.add_edge("a0", "b0")
    ec = ("a0", "b0", kc)
    atoms = {ea: list(range(10)), eb: list(range(10)),
             ec: list(range(10))}
    out = entangled_backbone_paths(G, DIMS, atoms, method="waypoint")
    assert ec not in out


def test_methods_differ():
    G, ea, eb = _pair_graph()
    atoms = {ea: list(range(24)), eb: list(range(24))}
    wp = entangled_backbone_paths(G, DIMS, atoms, method="waypoint")
    kk = entangled_backbone_paths(G, DIMS, atoms, method="kink")
    assert not np.allclose(np.asarray(wp[ea]), np.asarray(kk[ea]))


def test_config_default_is_waypoint():
    assert EntanglementsConfig().method == "waypoint"
    assert EntanglementsConfig(method="kink").method == "kink"
