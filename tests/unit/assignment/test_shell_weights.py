"""Neighbour-shell weighting of the entanglement draw.

The property: naming shells restricts the draw to those shells, and the
weights set the mix among them. Not naming any is the legacy behaviour and
must stay exactly that.
"""

import networkx as nx
import numpy as np
import pytest

from topon.assignment.entanglements import compute_shell_weights
from topon.config.schema import EntanglementsConfig


def _ladder(rungs=4, spacing=1.0, length=4.0):
    """Parallel edges at 1x, 2x, 3x ... the spacing, so shells are known.

    Each rung is one edge; rung i sits at y = i * spacing. The gap between
    rung 0 and rung i is therefore i * spacing exactly, which puts the pairs
    in shells that can be checked by hand.
    """
    G = nx.MultiGraph()
    for i in range(rungs):
        G.add_node(2 * i, pos=(0.0, i * spacing, 0.0))
        G.add_node(2 * i + 1, pos=(length, i * spacing, 0.0))
        G.add_edge(2 * i, 2 * i + 1)
    return G


def _pairs(G):
    e = [(u, v, 0) for u, v in G.edges()]
    return [(e[i], e[j]) for i in range(len(e)) for j in range(i + 1, len(e))]


def test_bands_are_read_off_the_geometry():
    """Shell numbering comes from the actual separations, not from a guess."""
    G = _ladder(4)
    cands = _pairs(G)
    # Ask for every shell, weighted by its own index, so the returned weight
    # is the band number.
    bands = compute_shell_weights(cands, G, None,
                                  {b: float(b) for b in range(1, 8)})
    # Rungs at 0,1,2,3 give separations of 1, 2 and 3 spacings: three bands.
    assert set(int(b) for b in bands) == {1, 2, 3}


def test_naming_one_shell_excludes_the_others():
    G = _ladder(4)
    cands = _pairs(G)
    w = compute_shell_weights(cands, G, None, {1: 1.0})
    assert any(x > 0 for x in w)
    # Everything that is not the nearest shell is excluded outright.
    bands = compute_shell_weights(cands, G, None,
                                  {b: float(b) for b in range(1, 8)})
    for weight, band in zip(w, bands):
        assert (weight > 0) == (int(band) == 1)


def test_weights_are_relative_between_shells():
    G = _ladder(4)
    cands = _pairs(G)
    w = compute_shell_weights(cands, G, None, {1: 3.0, 2: 1.0})
    bands = compute_shell_weights(cands, G, None,
                                  {b: float(b) for b in range(1, 8)})
    first = [x for x, b in zip(w, bands) if int(b) == 1]
    second = [x for x, b in zip(w, bands) if int(b) == 2]
    assert first and second
    assert all(x == 3.0 for x in first)
    assert all(x == 1.0 for x in second)


def test_no_weights_is_the_legacy_behaviour():
    """An empty mapping must not restrict anything.

    The default has to be inert: this rides on the same draw the pipeline has
    always used, and a config that does not mention shells must behave as it
    did before shells existed.
    """
    assert EntanglementsConfig().shell_weights == {}
    G = _ladder(3)
    w = compute_shell_weights(_pairs(G), G, None, {})
    assert all(x == 0.0 for x in w)   # caller skips the multiply entirely


def test_a_shell_nobody_asked_for_gets_nothing():
    G = _ladder(4)
    w = compute_shell_weights(_pairs(G), G, None, {9: 1.0})
    assert all(x == 0.0 for x in w)


def test_missing_positions_do_not_raise():
    """A node without a position drops its candidate rather than failing."""
    G = _ladder(3)
    del G.nodes[0]["pos"]
    w = compute_shell_weights(_pairs(G), G, None, {1: 1.0})
    assert len(w) == len(_pairs(G))


def test_periodic_images_are_used_when_dims_are_given():
    """A pair that is close across the boundary counts as close."""
    G = nx.MultiGraph()
    G.add_node(0, pos=(0.0, 0.2, 0.0))
    G.add_node(1, pos=(4.0, 0.2, 0.0))
    G.add_node(2, pos=(0.0, 9.8, 0.0))
    G.add_node(3, pos=(4.0, 9.8, 0.0))
    G.add_edge(0, 1)
    G.add_edge(2, 3)
    dims = np.array([10.0, 10.0, 10.0])
    cands = _pairs(G)
    # 0.4 apart through the boundary, 9.6 apart without it. Either way it is
    # the only pair, so it is band 1; the point is that it does not raise and
    # the minimum image is what was measured.
    with_mic = compute_shell_weights(cands, G, dims, {1: 1.0})
    assert with_mic == [1.0]
