"""
Unit tests for Topology Module.
Tests Python-based topology generation and graph loading.
"""

import pytest
import networkx as nx
from collections import namedtuple
from topon.topology.generator_python import PythonTopologyGenerator

# Mock configuration object
TopologyConfig = namedtuple('TopologyConfig', [
    'lattice_source', 'lattice_size', 'periodicity', 
    'degree_distribution', 'max_functionality'
])

@pytest.fixture
def small_lattice_config():
    """Config for a 3x3x3 Simple Cubic lattice."""
    return TopologyConfig(
        lattice_source="SC",
        lattice_size="3x3x3",
        periodicity="111",
        degree_distribution="",
        max_functionality=6
    )

def test_lattice_creation(small_lattice_config):
    """Test initial SC lattice creation."""
    gen = PythonTopologyGenerator(small_lattice_config)
    g = gen._create_lattice((3,3,3), "SC")
    
    assert g.number_of_nodes() == 27
    # SC lattice has 6 neighbors per node (periodic)
    # Total edges = 27 * 6 / 2 = 81
    assert g.number_of_edges() == 81
    
    # Verify degree 6 for all
    degrees = [d for n, d in g.degree()]
    assert all(d == 6 for d in degrees)

def test_strict_sculpting_d0(small_lattice_config):
    """Test Stage 1: Creating explicit degree-0 nodes."""
    # Request 2 nodes of degree 0
    config = small_lattice_config._replace(degree_distribution="0:2")
    gen = PythonTopologyGenerator(config)
    
    graphs = gen.generate(trials=10)
    assert len(graphs) > 0
    g = graphs[0]
    
    degrees = [d for n, d in g.degree()]
    d0_count = degrees.count(0)
    assert d0_count == 2
    assert nx.number_connected_components(g) > 1 # isolated nodes

def test_reduced_functionality(small_lattice_config):
    """Test Stage 3: Enforcing max functionality."""
    # Max func 4
    config = small_lattice_config._replace(max_functionality=4)
    gen = PythonTopologyGenerator(config)
    
    graphs = gen.generate(trials=5)
    assert len(graphs) > 0
    g = graphs[0]
    
    degrees = [d for n, d in g.degree()]
    assert max(degrees) <= 4

def test_connectivity_check():
    """Test _is_subgraph_connected Logic."""
    gen = PythonTopologyGenerator(TopologyConfig("SC", "2x2x2", "111", "", 6))
    g = nx.Graph()
    g.add_edges_from([(0,1), (1,2), (2,3)]) # 0-1-2-3 line
    
    # All active -> Connected
    status = {0:"ACTIVE", 1:"ACTIVE", 2:"ACTIVE", 3:"ACTIVE"}
    assert gen._is_subgraph_connected(g, status) == True
    
    # 3 active, 0 removed (isolated)
    g.remove_edge(0,1)
    # 0 isolated
    status[0] = "IS_DEGREE_0"
    # 1-2-3 still connected
    assert gen._is_subgraph_connected(g, status) == True
    
    # Break 1-2 -> 1 isolated, 2-3 connected
    g.remove_edge(1,2)
    # 1 is ACTIVE, 2-3 ACTIVE
    # Active nodes: 1, 2, 3. Edges: 2-3.
    # Disconnected (1 is separate from 2-3)
    assert gen._is_subgraph_connected(g, status) == False


# ---------------------------------------------------------------------------
# Unreachable-target validation (fail fast instead of churning through trials)
# ---------------------------------------------------------------------------

def test_over_target_edge_count_raises_fast(small_lattice_config):
    """e:N above the base lattice's edge count must fail fast, not hang.

    A 3x3x3 SC lattice has only 81 edges. Sculpting only ever removes edges,
    so e:128 is structurally unreachable and used to make generate() grind
    through hundreds of thousands of doomed trials. The time_limit is a
    belt-and-suspenders guard: the raise happens before the trial loop, so a
    future regression fails this test in seconds rather than hanging CI.
    """
    config = small_lattice_config._replace(degree_distribution="e:128")
    gen = PythonTopologyGenerator(config)
    with pytest.raises(ValueError, match=r"e:128 exceeds the 81 edges"):
        gen.generate(trials=1_000_000, time_limit=5)


def test_edge_count_target_at_capacity_ok(small_lattice_config):
    """e:N equal to the full lattice edge count is reachable (boundary)."""
    config = small_lattice_config._replace(degree_distribution="e:81")
    gen = PythonTopologyGenerator(config)
    graphs = gen.generate(trials=5)
    assert len(graphs) > 0
    assert graphs[0].number_of_edges() == 81


def test_reachable_edge_count_target_ok(small_lattice_config):
    """e:N below the full lattice edge count still sculpts to exactly N edges."""
    config = small_lattice_config._replace(degree_distribution="e:70")
    gen = PythonTopologyGenerator(config)
    graphs = gen.generate(trials=50)
    assert len(graphs) > 0
    assert graphs[0].number_of_edges() == 70


def test_over_target_per_degree_count_raises(small_lattice_config):
    """A per-degree target asking for more nodes than exist fails fast."""
    # 3x3x3 SC has only 27 nodes; 100 nodes of degree 3 is impossible.
    config = small_lattice_config._replace(degree_distribution="3:100")
    gen = PythonTopologyGenerator(config)
    with pytest.raises(ValueError, match=r"3:100 exceeds the 27 nodes"):
        gen.generate(trials=1_000_000, time_limit=5)


def test_over_target_per_degree_degree_raises(small_lattice_config):
    """A per-degree target above the lattice's maximum degree fails fast."""
    # SC max degree is 6; requesting degree-7 nodes is unreachable by sculpting.
    config = small_lattice_config._replace(degree_distribution="7:5")
    gen = PythonTopologyGenerator(config)
    with pytest.raises(ValueError, match=r"maximum degree in a 3x3x3 SC lattice is 6"):
        gen.generate(trials=1_000_000, time_limit=5)


# ---------------------------------------------------------------------------
# Connectivity check: hand-rolled traversal must match NetworkX exactly
# ---------------------------------------------------------------------------

def _networkx_reference(g, node_status):
    """What `_is_subgraph_connected` used to do, kept as the oracle."""
    active = [n for n in g.nodes() if node_status[n] == "ACTIVE"]
    if not active:
        return True
    return nx.is_connected(g.subgraph(active))


@pytest.mark.parametrize("lattice_type,dims", [("SC", (4, 4, 4)),
                                               ("BCC", (3, 3, 3)),
                                               ("FCC", (3, 3, 3))])
def test_connectivity_matches_networkx_on_random_states(lattice_type, dims):
    """The direct traversal must agree with NetworkX on every state.

    `_is_subgraph_connected` walks `g._adj` itself instead of building a
    `g.subgraph(...)` view, because the view re-evaluates its node filter
    on every neighbour access and that call is ~99% of the generator's
    runtime. The speedup is only worth having if the answers are
    identical, so this fuzzes random edge removals against random status
    assignments and compares.
    """
    import random

    gen = PythonTopologyGenerator(
        TopologyConfig(lattice_type, "x".join(map(str, dims)), "111", "", 6)
    )
    base = gen._create_lattice(dims, lattice_type)

    random.seed(20260805)
    for _ in range(120):
        g = base.copy()
        edges = list(g.edges())
        for e in random.sample(edges, random.randint(0, len(edges) // 2)):
            g.remove_edge(*e)
        status = {}
        for n in g.nodes():
            r = random.random()
            status[n] = ("ACTIVE" if r < 0.7
                         else "IS_DEGREE_0" if r < 0.85 else "IS_DEGREE_1")
        assert gen._is_subgraph_connected(g, status) == _networkx_reference(g, status)


def test_connectivity_handles_degenerate_states():
    """No ACTIVE nodes at all, and a single ACTIVE node, both count as connected."""
    gen = PythonTopologyGenerator(TopologyConfig("SC", "2x2x2", "111", "", 6))
    g = nx.Graph()
    g.add_edges_from([(0, 1), (1, 2)])

    none_active = {n: "IS_DEGREE_0" for n in g.nodes()}
    assert gen._is_subgraph_connected(g, none_active) is True
    assert _networkx_reference(g, none_active) is True

    one_active = {0: "ACTIVE", 1: "IS_DEGREE_1", 2: "IS_DEGREE_0"}
    assert gen._is_subgraph_connected(g, one_active) == _networkx_reference(g, one_active)

