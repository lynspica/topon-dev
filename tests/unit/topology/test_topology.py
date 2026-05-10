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

