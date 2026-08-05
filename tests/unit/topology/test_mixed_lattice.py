"""Unit tests for the mixed SC/BCC/FCC lattice (``lattice_type="MIX"``).

Semantics under test
--------------------
All three lattices share the cubic cell corner and each adds sites on
top of it: BCC one body centre, FCC three face centres. So ``MIX``
places the corner in every cell, the body centre with probability
``mix_fractions["BCC"]`` and each face centre with probability
``mix_fractions["FCC"]``. The ``"SC"`` entry is the remainder and
contributes no site of its own, which is what makes the three fractions
a partition summing to 1.

Edges join every pair within ``mix_cutoff`` under the minimum image,
because a mixed point set has no single neighbour shell to enumerate.

The load-bearing guarantee is that ``MIX`` at ``(1, 0, 0)`` reproduces
``SC`` exactly. It deliberately does *not* hold at the other two
corners, and `test_mix_bcc_corner_is_not_canonical_bcc` pins that so the
difference stays visible rather than becoming a surprise.
"""

import random

import networkx as nx
import numpy as np
import pytest

from topon.config.schema import GeneratorConfig
from topon.topology.generator_python import PythonTopologyGenerator
from topon.topology.loader import infer_dims_from_graph

DIMS = (4, 4, 4)
N_CELLS = DIMS[0] * DIMS[1] * DIMS[2]


class _Cfg:
    """Minimal stand-in for a GeneratorConfig."""

    def __init__(self, lattice_type, mix=None, cutoff=1.0,
                 lattice_size=DIMS, degree_distribution="", max_func=6):
        self.lattice_type = lattice_type
        self.lattice_size = lattice_size
        self.max_functionality = max_func
        self.degree_distribution = degree_distribution
        self.mix_fractions = mix or {"SC": 1.0, "BCC": 0.0, "FCC": 0.0}
        self.mix_cutoff = cutoff


def build(lattice_type, mix=None, cutoff=1.0, dims=DIMS, seed=42):
    random.seed(seed)
    gen = PythonTopologyGenerator(_Cfg(lattice_type, mix, cutoff, dims))
    return gen._create_lattice(dims, lattice_type)


def edge_lengths(G):
    box = np.asarray(G.graph["box"], dtype=float)
    pos = {n: np.asarray(d["pos"], dtype=float) for n, d in G.nodes(data=True)}
    return np.asarray([
        np.linalg.norm((pos[u] - pos[v]) - box * np.round((pos[u] - pos[v]) / box))
        for u, v in G.edges()
    ])


# ---------------------------------------------------------------------------
# The anchor: the SC corner is exactly the existing SC lattice
# ---------------------------------------------------------------------------

def test_mix_sc_corner_reproduces_sc_exactly():
    """MIX at (1, 0, 0) must be indistinguishable from lattice_type="SC".

    This is what makes MIX additive: an existing SC study re-expressed as
    a mix is the same network, down to node ids.
    """
    sc = build("SC")
    mix = build("MIX", {"SC": 1.0, "BCC": 0.0, "FCC": 0.0})

    assert sorted(mix.nodes()) == sorted(sc.nodes())
    for n in sc.nodes():
        assert mix.nodes[n]["pos"] == sc.nodes[n]["pos"]
    assert {frozenset(e) for e in mix.edges()} == {frozenset(e) for e in sc.edges()}
    assert mix.graph["box"] == sc.graph["box"]


def test_mix_bcc_corner_is_not_canonical_bcc():
    """Pin the documented discontinuity at the BCC and FCC corners.

    MIX uses a distance cutoff, so at (0, 1, 0) it admits the
    corner-corner shell at 1.0 on top of the corner-body shell at 0.866.
    That gives 14 neighbours per node where canonical BCC has 8. Same
    story at (0, 0, 1): 18 against FCC's 12. The site sets do match; only
    the connectivity differs. If this ever starts passing, MIX changed
    meaning and the docstrings need revisiting.
    """
    for frac, pure, mix_deg, pure_deg in (
        ({"SC": 0.0, "BCC": 1.0, "FCC": 0.0}, "BCC", 14, 8),
        ({"SC": 0.0, "BCC": 0.0, "FCC": 1.0}, "FCC", 18, 12),
    ):
        mix = build("MIX", frac)
        canonical = build(pure)

        # Same sites ...
        assert mix.number_of_nodes() == canonical.number_of_nodes()
        assert (sorted(tuple(d["pos"]) for _, d in mix.nodes(data=True))
                == sorted(tuple(d["pos"]) for _, d in canonical.nodes(data=True)))
        # ... different connectivity.
        assert {d for _, d in mix.degree()} == {mix_deg}
        assert {d for _, d in canonical.degree()} == {pure_deg}


# ---------------------------------------------------------------------------
# Site counts
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "mix,expected",
    [
        ({"SC": 1.0, "BCC": 0.0, "FCC": 0.0}, 1 * N_CELLS),   # SC
        ({"SC": 0.0, "BCC": 1.0, "FCC": 0.0}, 2 * N_CELLS),   # BCC
        ({"SC": 0.0, "BCC": 0.0, "FCC": 1.0}, 4 * N_CELLS),   # FCC
    ],
)
def test_pure_corners_have_exact_site_counts(mix, expected):
    """N, 2N and 4N sites, matching the pure lattices they stand in for."""
    assert build("MIX", mix).number_of_nodes() == expected


def test_mixture_site_count_matches_the_formula_on_average():
    """Expected count is N*(1 + f_bcc + 3*f_fcc); draws are binomial.

    Averaged over seeds so the test is about the formula rather than one
    lucky draw. The tolerance is loose enough that it cannot flake, but
    tight enough to catch a wrong coefficient (dropping the factor 3 on
    the face sites would land ~26% low).
    """
    mix = {"SC": 0.2, "BCC": 0.4, "FCC": 0.4}
    expected = N_CELLS * (1 + mix["BCC"] + 3 * mix["FCC"])

    counts = [build("MIX", mix, seed=s).number_of_nodes() for s in range(12)]
    assert abs(np.mean(counts) - expected) / expected < 0.05


def test_zero_fraction_places_no_sites_of_that_kind():
    """f_fcc = 0 must leave no half-integer face site behind."""
    G = build("MIX", {"SC": 0.5, "BCC": 0.5, "FCC": 0.0})
    for _, data in G.nodes(data=True):
        x, y, z = data["pos"]
        fracs = {v % 1.0 for v in (x, y, z)}
        # Corner sites are all-integer; body sites are all-half.
        assert fracs in ({0.0}, {0.5}), f"unexpected site at {data['pos']}"


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def test_box_is_the_cell_count_not_the_coordinate_extent():
    """MIX inherits the recorded-cell fix; half-integer sites must not inflate it."""
    G = build("MIX", {"SC": 0.2, "BCC": 0.4, "FCC": 0.4})
    assert G.graph["box"] == (4.0, 4.0, 4.0)
    assert np.allclose(infer_dims_from_graph(G), [4.0, 4.0, 4.0])


def test_non_cubic_cell():
    G = build("MIX", {"SC": 0.5, "BCC": 0.5, "FCC": 0.0}, dims=(3, 4, 5))
    assert np.allclose(infer_dims_from_graph(G), [3.0, 4.0, 5.0])


def test_every_edge_is_within_the_cutoff():
    for cutoff in (0.9, 1.0, 1.2):
        G = build("MIX", {"SC": 0.2, "BCC": 0.4, "FCC": 0.4}, cutoff=cutoff)
        if G.number_of_edges() == 0:
            continue
        assert edge_lengths(G).max() <= cutoff + 1e-9


def test_cutoff_below_the_corner_shell_drops_corner_corner_bonds():
    """A cutoff under 1.0 excludes the simple-cubic shell.

    Documents why 1.0 is the default: below it the always-present corner
    sublattice stops being bonded to itself.
    """
    frac = {"SC": 0.0, "BCC": 1.0, "FCC": 0.0}
    tight = build("MIX", frac, cutoff=0.9)
    assert np.isclose(edge_lengths(tight).max(), np.sqrt(3) / 2)
    assert {d for _, d in tight.degree()} == {8}      # canonical BCC after all


def test_mixing_adds_distance_shells():
    """The actual point of mixing: more neighbour distances, smoother P(Ree).

    SC alone offers one edge length. Adding BCC and FCC sites brings in
    the 0.866 and 0.707 shells, plus the 0.5 body-to-face contact.
    """
    def shells(mix):
        return set(np.round(edge_lengths(build("MIX", mix)), 4))

    sc_only = shells({"SC": 1.0, "BCC": 0.0, "FCC": 0.0})
    all_three = shells({"SC": 0.2, "BCC": 0.4, "FCC": 0.4})

    assert sc_only == {1.0}
    assert all_three == {0.5, 0.7071, 0.866, 1.0}
    assert len(all_three) > len(sc_only)


def test_body_and_face_sites_can_sit_half_a_cell_apart():
    """Pin the closest contact a mixture can produce.

    A BCC body centre and an FCC face centre on the shared cell face are
    0.5 apart, tighter than any pure lattice's nearest neighbour (SC 1.0,
    BCC 0.866, FCC 0.707). DP is assigned independently of edge length,
    so this sets the bond-length spread at fixed DP and is a real
    consequence of mixing rather than a defect.
    """
    L = edge_lengths(build("MIX", {"SC": 0.2, "BCC": 0.4, "FCC": 0.4}))
    assert np.isclose(L.min(), 0.5)
    assert np.isclose(L.max() / L.min(), 2.0)


# ---------------------------------------------------------------------------
# Reproducibility and sculpting
# ---------------------------------------------------------------------------

def test_same_seed_gives_the_same_lattice():
    mix = {"SC": 0.2, "BCC": 0.4, "FCC": 0.4}
    a = build("MIX", mix, seed=7)
    b = build("MIX", mix, seed=7)
    assert a.number_of_nodes() == b.number_of_nodes()
    assert {frozenset(e) for e in a.edges()} == {frozenset(e) for e in b.edges()}


def test_different_seeds_give_different_lattices():
    mix = {"SC": 0.2, "BCC": 0.4, "FCC": 0.4}
    a = build("MIX", mix, seed=1)
    b = build("MIX", mix, seed=2)
    assert a.number_of_nodes() != b.number_of_nodes() or (
        {frozenset(e) for e in a.edges()} != {frozenset(e) for e in b.edges()}
    )


def test_sculpting_runs_on_a_mixed_lattice():
    """The strict-sculpting stages must work unchanged on a mixed base graph."""
    random.seed(42)
    cfg = _Cfg("MIX", {"SC": 0.34, "BCC": 0.33, "FCC": 0.33},
               degree_distribution="0:0,1:0", max_func=4)
    graphs = PythonTopologyGenerator(cfg).generate(trials=4000, max_saves=1)

    assert graphs, "sculpting produced no graph from a mixed lattice"
    G = graphs[0]
    assert max(d for _, d in G.degree()) <= 4
    assert np.allclose(infer_dims_from_graph(G), [4.0, 4.0, 4.0])


def test_unreachable_target_still_fails_fast_on_a_mixed_lattice():
    """The V40 guard reads bounds off the built graph, so it adapts to MIX."""
    cfg = _Cfg("MIX", {"SC": 1.0, "BCC": 0.0, "FCC": 0.0},
               degree_distribution="e:100000")
    gen = PythonTopologyGenerator(cfg)
    with pytest.raises(ValueError, match=r"e:100000 exceeds"):
        gen.generate(trials=1_000_000, time_limit=5)


def test_error_message_names_the_fractions():
    """An unreachable target should say which mixture it was unreachable on.

    Degree 99 is well past what any cutoff-built mixture reaches (the
    densest, all-FCC case tops out at 18), so this exercises the
    max-degree branch of the V40 guard.
    """
    cfg = _Cfg("MIX", {"SC": 0.2, "BCC": 0.4, "FCC": 0.4},
               degree_distribution="99:5")
    gen = PythonTopologyGenerator(cfg)
    with pytest.raises(ValueError, match=r"MIX \(SC:0\.2,BCC:0\.4,FCC:0\.4\)"):
        gen.generate(trials=10, time_limit=5)


# ---------------------------------------------------------------------------
# Config validation (fails at load time, before a long run starts)
# ---------------------------------------------------------------------------

def test_config_accepts_mix():
    cfg = GeneratorConfig(lattice_type="MIX",
                          mix_fractions={"SC": 0.2, "BCC": 0.4, "FCC": 0.4})
    assert cfg.lattice_type == "MIX"
    assert cfg.mix_cutoff == 1.0


def test_config_default_is_unchanged_sc():
    cfg = GeneratorConfig()
    assert cfg.lattice_type == "SC"
    assert cfg.mix_fractions == {"SC": 1.0, "BCC": 0.0, "FCC": 0.0}


@pytest.mark.parametrize(
    "bad,match",
    [
        ({"SC": 0.5, "BCC": 0.2, "FCC": 0.0}, "must sum to 1"),
        ({"SC": 0.5, "BCC": 0.5, "FCC": 0.5}, "must sum to 1"),
        ({"SC": 1.5, "BCC": -0.5, "FCC": 0.0}, "non-negative"),
        ({"SC": 0.5, "Diamond": 0.5}, "unknown key"),
    ],
)
def test_config_rejects_bad_fractions(bad, match):
    with pytest.raises(ValueError, match=match):
        GeneratorConfig(mix_fractions=bad)


def test_config_rejects_non_positive_cutoff():
    with pytest.raises(ValueError):
        GeneratorConfig(mix_cutoff=0.0)


def test_unknown_lattice_type_is_rejected_by_the_generator():
    gen = PythonTopologyGenerator(_Cfg("HCP"))
    with pytest.raises(NotImplementedError, match="SC, BCC, FCC, or MIX"):
        gen._create_lattice(DIMS, "HCP")
