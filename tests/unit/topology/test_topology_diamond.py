"""
Unit tests for the Diamond-lattice topology generator.

Focus: the fail-fast guard against structurally-unreachable
``degree_distribution`` targets (sculpting only removes edges, so the
freshly-built lattice is a hard ceiling), mirroring the guard in
``PythonTopologyGenerator``.
"""

import pytest
from collections import namedtuple

from topon.topology.generator_python_diamond import (
    DiamondTopologyGenerator,
    create_diamond_lattice,
)

# Duck-typed config: DiamondTopologyGenerator reads these three attributes.
DiamondConfig = namedtuple(
    "DiamondConfig", ["lattice_size", "max_functionality", "degree_distribution"]
)

# A 2x2x2 diamond lattice: 64 nodes, 128 edges, 4-regular.
BASE_NODES = 64
BASE_EDGES = 128


def _cfg(degree_distribution=""):
    return DiamondConfig(
        lattice_size=(2, 2, 2),
        max_functionality=4,
        degree_distribution=degree_distribution,
    )


def test_base_lattice_shape():
    """Sanity: the 2x2x2 diamond lattice matches the assumed counts."""
    g = create_diamond_lattice(2, 2, 2)
    assert g.number_of_nodes() == BASE_NODES
    assert g.number_of_edges() == BASE_EDGES
    assert max(d for _, d in g.degree()) == 4


def test_over_target_edge_count_raises_fast():
    """e:N above the lattice edge count must fail fast, not hang.

    The time_limit is a belt-and-suspenders guard: the raise happens before
    the trial loop, so a future regression fails this test in seconds rather
    than hanging CI.
    """
    gen = DiamondTopologyGenerator(_cfg("e:200"))
    with pytest.raises(ValueError, match=r"e:200 exceeds the 128 edges"):
        gen.generate(trials=1_000_000, time_limit=5)


def test_reachable_edge_count_target_ok():
    """e:N below the full lattice sculpts to exactly N edges."""
    gen = DiamondTopologyGenerator(_cfg("e:124"))
    graphs = gen.generate(trials=50)
    assert len(graphs) > 0
    assert graphs[0].number_of_edges() == 124


def test_empty_distribution_returns_full_lattice():
    """No constraints -> the raw 4-regular lattice, untouched (no false raise)."""
    gen = DiamondTopologyGenerator(_cfg(""))
    graphs = gen.generate(trials=1)
    assert len(graphs) == 1
    assert graphs[0].number_of_edges() == BASE_EDGES


def test_over_target_per_degree_count_raises():
    """A per-degree target asking for more nodes than exist fails fast."""
    gen = DiamondTopologyGenerator(_cfg("3:100"))
    with pytest.raises(ValueError, match=r"3:100 exceeds the 64 nodes"):
        gen.generate(trials=1_000_000, time_limit=5)


def test_over_target_per_degree_degree_raises():
    """A per-degree target above the lattice's max degree (4) fails fast."""
    gen = DiamondTopologyGenerator(_cfg("9:5"))
    with pytest.raises(ValueError, match=r"maximum degree in a 2x2x2 Diamond lattice is 4"):
        gen.generate(trials=1_000_000, time_limit=5)
