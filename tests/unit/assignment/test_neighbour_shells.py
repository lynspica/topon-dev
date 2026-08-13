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

    Shell 1 is the chains that share a crosslink, at distance zero: they are
    the closest strands in the network and they can wind around each other,
    since fixing the junction constrains where the winding sits rather than
    whether there is one. After that come the simple cubic distances 1,
    sqrt(2), sqrt(3), each heavily degenerate, so a sorted-list reading of
    "second neighbour" would return another chain at the same distance.
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
    for s, want in ((1, 0.0), (2, 1.0), (3, np.sqrt(2)), (4, np.sqrt(3))):
        assert np.mean(seen[s]) == pytest.approx(want, abs=0.02)


def test_a_shell_holds_every_neighbour_at_that_distance():
    G, dims = grid()
    sh = neighbour_shells(G, dims, max_shell=3)
    counts = [len(by[1]) for by in sh.values() if 1 in by]
    # Degenerate by construction: one neighbour would mean the shells were
    # being read off a sorted list rather than off the geometry.
    assert min(counts) > 1


def test_chains_sharing_a_junction_are_the_closest_shell():
    """They are included, and they are shell 1.

    Excluding them was a code limitation mistaken for physics. On SC 4x4x4 it
    discarded 263 pairs, five per chain, each strand carrying 77 sigma of
    contour on a 5.4 sigma chord.
    """
    G, dims = grid()
    d = chain_distances(G, dims)
    shared = [(a, b) for (a, b) in d
              if len({a[0], a[1]} & {b[0], b[1]}) == 1]
    assert shared, "junction-sharing pairs must be ranked, not discarded"
    assert all(d[q] == pytest.approx(0.0, abs=1e-9) for q in shared)

    sh = neighbour_shells(G, dims, max_shell=2)
    a, b = shared[0]
    assert b in sh[a].get(1, ()), "a shared junction is the closest shell"


def test_parallel_edges_are_still_excluded():
    """Two chains joining the same two junctions are one contact, not a pair."""
    G, dims = grid()
    d = chain_distances(G, dims)
    for (a, b) in d:
        assert {a[0], a[1]} != {b[0], b[1]}


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


def test_yield_weighting_asks_for_more_where_a_pair_is_worth_less():
    """The request is a mix of entanglements, not of pairs.

    A pair in an outer shell delivers more, because the routed chain travels
    further and picks up more on the way. Measured over 62 designed pairs on
    SC, asked 0.20/0.50/0.25/0.05 and delivered 0.11/0.57/0.32/0.00.
    """
    G, dims = grid()
    sh = neighbour_shells(G, dims, max_shell=3)
    want = {1: 0.5, 2: 0.5}
    # Shell 2 pairs are worth four times shell 1 pairs here.
    yields = {1: 0.05, 2: 0.20}
    plain = select_by_shells(G, 2.0, want, dims, shells=sh,
                             rng=np.random.default_rng(5))
    tuned = select_by_shells(G, 2.0, want, dims, shells=sh,
                             rng=np.random.default_rng(5),
                             yield_by_shell=yields)

    of_shell = {s: {(min(c, o), max(c, o))
                    for c, by in sh.items() for o in by.get(s, ())}
                for s in want}

    def share(sel, s):
        n = sum(c for a, b, c in sel
                if (min(a, b), max(a, b)) in of_shell[s])
        tot = sum(c for _a, _b, c in sel)
        return n / tot if tot else 0.0

    # Weighting must move draws toward the shell whose pairs are worth less.
    assert share(tuned, 1) > share(plain, 1)
    assert share(tuned, 2) < share(plain, 2)


def test_no_yield_given_is_the_old_behaviour():
    G, dims = grid()
    sh = neighbour_shells(G, dims, max_shell=3)
    a = select_by_shells(G, 1.5, {1: 0.5, 2: 0.5}, dims, shells=sh,
                         rng=np.random.default_rng(9))
    b = select_by_shells(G, 1.5, {1: 0.5, 2: 0.5}, dims, shells=sh,
                         rng=np.random.default_rng(9), yield_by_shell=None)
    assert a == b


def test_equal_yields_change_nothing():
    G, dims = grid()
    sh = neighbour_shells(G, dims, max_shell=3)
    a = select_by_shells(G, 1.5, {1: 0.5, 2: 0.5}, dims, shells=sh,
                         rng=np.random.default_rng(11))
    b = select_by_shells(G, 1.5, {1: 0.5, 2: 0.5}, dims, shells=sh,
                         rng=np.random.default_rng(11),
                         yield_by_shell={1: 0.3, 2: 0.3})
    assert a == b


# ---------------------------------------------------------------------------
# Spatial placement
# ---------------------------------------------------------------------------

def pair_positions(G, sel, dims):
    """Where each selected pair sits, as fractional coordinates."""
    out = []
    for a, b, c in sel:
        mids = []
        for u, v, _k in (a, b):
            pu = np.asarray(G.nodes[u]["pos"], float)
            pv = np.asarray(G.nodes[v]["pos"], float)
            d = pv - pu
            d -= dims * np.round(d / dims)
            mids.append(pu + 0.5 * d)
        d = mids[1] - mids[0]
        d -= dims * np.round(d / dims)
        p = mids[0] + 0.5 * d
        out += [(p - dims * np.floor(p / dims)) / dims] * c
    return np.array(out)


def test_a_region_concentrates_the_entanglements_inside_it():
    """The stretch goal: put them where you want them, not everywhere."""
    G, dims = grid(6)
    sh = neighbour_shells(G, dims, max_shell=2)
    plain = select_by_shells(G, 2.0, {1: 1.0}, dims, shells=sh,
                             rng=np.random.default_rng(0))
    biased = select_by_shells(G, 2.0, {1: 1.0}, dims, shells=sh,
                              rng=np.random.default_rng(0),
                              bias_kind="region",
                              bias_params={"center": [0.5, 0.5, 0.5],
                                           "radius": 0.25,
                                           "strength": 20.0})

    def inside(sel):
        p = pair_positions(G, sel, dims)
        r = np.linalg.norm(p - 0.5, axis=1)
        return float((r < 0.25).mean())

    assert inside(biased) > inside(plain) + 0.2


def test_a_gradient_puts_more_at_one_end_than_the_other():
    G, dims = grid(6)
    sh = neighbour_shells(G, dims, max_shell=2)
    sel = select_by_shells(G, 3.0, {1: 1.0}, dims, shells=sh,
                           rng=np.random.default_rng(1),
                           bias_kind="gradient",
                           bias_params={"axis": "z", "strength": 3.0})
    z = pair_positions(G, sel, dims)[:, 2]
    assert (z > 0.5).mean() > 0.65


def test_anti_region_depletes_instead_of_concentrating():
    G, dims = grid(6)
    sh = neighbour_shells(G, dims, max_shell=2)
    plain = select_by_shells(G, 2.0, {1: 1.0}, dims, shells=sh,
                             rng=np.random.default_rng(2))
    hole = select_by_shells(G, 2.0, {1: 1.0}, dims, shells=sh,
                            rng=np.random.default_rng(2),
                            bias_kind="anti_region",
                            bias_params={"center": [0.5, 0.5, 0.5],
                                         "radius": 0.3,
                                         "strength": 20.0})

    def inside(sel):
        p = pair_positions(G, sel, dims)
        return float((np.linalg.norm(p - 0.5, axis=1) < 0.3).mean())

    assert inside(hole) < inside(plain)


def test_uniform_is_the_behaviour_without_a_bias():
    G, dims = grid()
    sh = neighbour_shells(G, dims, max_shell=2)
    a = select_by_shells(G, 1.5, {1: 1.0}, dims, shells=sh,
                         rng=np.random.default_rng(4))
    b = select_by_shells(G, 1.5, {1: 1.0}, dims, shells=sh,
                         rng=np.random.default_rng(4), bias_kind="uniform")
    assert a == b


def test_a_spatial_bias_and_a_shell_mix_compose():
    """Both requests are honoured at once, not one instead of the other."""
    G, dims = grid(6)
    sh = neighbour_shells(G, dims, max_shell=3)
    sel = select_by_shells(G, 2.0, {1: 0.5, 2: 0.5}, dims, shells=sh,
                           rng=np.random.default_rng(6),
                           bias_kind="region",
                           bias_params={"center": [0.5, 0.5, 0.5],
                                        "radius": 0.3, "strength": 20.0})
    of_shell = {s: {(min(c, o), max(c, o))
                    for c, by in sh.items() for o in by.get(s, ())}
                for s in (1, 2)}
    got = {s: sum(c for a, b, c in sel
                  if (min(a, b), max(a, b)) in of_shell[s]) for s in (1, 2)}
    tot = sum(got.values())
    # The shell mix survives the spatial bias.
    assert 0.35 < got[1] / tot < 0.65
    # And the spatial bias is still doing something.
    p = pair_positions(G, sel, dims)
    assert (np.linalg.norm(p - 0.5, axis=1) < 0.3).mean() > 0.35
