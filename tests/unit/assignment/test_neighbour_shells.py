"""Neighbour shells and shell-resolved selection.

`find_crossing_candidates` gives each chain its single nearest disjoint
neighbour, so the pool it returns is entirely first-shell -- measured on an SC
4x4x4 network, 89 of 89 pairs. No reweighting of that pool can reach a second
or third neighbour that is not in it, which is what these functions add.
"""
import networkx as nx
import numpy as np
import pytest

from topon.assignment.entanglements import (
    chain_distances,
    neighbour_shells,
    select_by_shells,
)


def grid(n=4, spacing=1.0):
    """A simple cubic network, periodic, with every node degree 6."""
    G = nx.MultiGraph()
    for i in range(n):
        for j in range(n):
            for k in range(n):
                G.add_node((i, j, k),
                           pos=np.array([i, j, k], float) * spacing)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    G.add_edge((i, j, k),
                               ((i + d[0]) % n, (j + d[1]) % n,
                                (k + d[2]) % n))
    return G, np.array([n, n, n], float) * spacing


def test_shells_land_on_the_lattice_distances():
    """Shell k is the k-th distinct distance, not the k-th sorted neighbour.

    On a simple cubic lattice those are 1, sqrt(2), sqrt(3), 2, and each is
    heavily degenerate, so a sorted-list reading of "second neighbour" would
    return another chain at distance 1.
    """
    G, dims = grid()
    d = chain_distances(G, dims)
    sh = neighbour_shells(G, dims, max_shell=4, distances=d)
    seen = {}
    for chain, by in sh.items():
        for s, others in by.items():
            for o in others:
                r = d.get((chain, o), d.get((o, chain)))
                seen.setdefault(s, []).append(r)
    for s, want in ((1, 1.0), (2, np.sqrt(2)), (3, np.sqrt(3))):
        assert np.mean(seen[s]) == pytest.approx(want, abs=0.02)


def test_a_shell_holds_every_neighbour_at_that_distance():
    G, dims = grid()
    sh = neighbour_shells(G, dims, max_shell=3)
    counts = [len(by[1]) for by in sh.values() if 1 in by]
    # Degenerate by construction: one neighbour would mean the shells were
    # being read off a sorted list rather than off the geometry.
    assert min(counts) > 1


def test_chains_sharing_a_junction_are_not_neighbours():
    G, dims = grid()
    d = chain_distances(G, dims)
    for (a, b) in d:
        assert not ({a[0], a[1]} & {b[0], b[1]})


def test_density_follows_the_existing_convention():
    """total = per_chain * 0.5 * num_chains, as distribution mode uses."""
    G, dims = grid()
    sh = neighbour_shells(G, dims, max_shell=3)
    rng = np.random.default_rng(0)
    sel = select_by_shells(G, 2.0, {1: 1.0}, dims, shells=sh, rng=rng)
    assert sum(c for _a, _b, c in sel) == int(
        round(2.0 * 0.5 * G.number_of_edges()))


def test_naming_one_shell_draws_only_from_it():
    G, dims = grid()
    sh = neighbour_shells(G, dims, max_shell=4)
    rng = np.random.default_rng(1)
    sel = select_by_shells(G, 1.0, {2: 1.0}, dims, shells=sh, rng=rng)
    in_two = {(min(c, o), max(c, o))
              for c, by in sh.items() for o in by.get(2, ())}
    for a, b, _c in sel:
        assert (min(a, b), max(a, b)) in in_two


def test_the_mix_is_the_closest_whole_draw_approximation():
    G, dims = grid()
    sh = neighbour_shells(G, dims, max_shell=4)
    rng = np.random.default_rng(2)
    want = {1: 0.2, 2: 0.5, 3: 0.25, 4: 0.05}
    sel = select_by_shells(G, 2.0, want, dims, shells=sh, rng=rng)
    of_shell = {s: {(min(c, o), max(c, o))
                    for c, by in sh.items() for o in by.get(s, ())}
                for s in want}
    got = {s: 0 for s in want}
    for a, b, c in sel:
        for s, pool in of_shell.items():
            if (min(a, b), max(a, b)) in pool:
                got[s] += c
                break
    total = sum(got.values())
    for s, f in want.items():
        assert got[s] / total == pytest.approx(f, abs=0.06)


def test_zero_density_selects_nothing():
    G, dims = grid()
    assert select_by_shells(G, 0.0, {1: 1.0}, dims) == []


def test_an_empty_mix_selects_nothing():
    G, dims = grid()
    assert select_by_shells(G, 2.0, {}, dims) == []


def test_selection_is_reproducible():
    G, dims = grid()
    sh = neighbour_shells(G, dims, max_shell=3)
    a = select_by_shells(G, 1.5, {1: 0.5, 2: 0.5}, dims, shells=sh,
                         rng=np.random.default_rng(7))
    b = select_by_shells(G, 1.5, {1: 0.5, 2: 0.5}, dims, shells=sh,
                         rng=np.random.default_rng(7))
    assert a == b


def test_max_per_pair_caps_repeat_draws():
    G, dims = grid()
    sh = neighbour_shells(G, dims, max_shell=2)
    sel = select_by_shells(G, 4.0, {1: 1.0}, dims, shells=sh,
                           rng=np.random.default_rng(3), max_per_pair=2)
    assert sel and max(c for _a, _b, c in sel) <= 2
