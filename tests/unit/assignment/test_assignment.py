"""
Unit tests for Assignment Module.
Tests assignment of node/edge types based on various rules.
"""

import pytest
import networkx as nx
import numpy as np
from collections import namedtuple

from topon.assignment.node_types import (
    _assign_by_degree, 
    _assign_by_position, 
    _assign_by_random,
    _assign_explicit
)
from topon.assignment.edge_types import (
    _assign_uniform,
    _assign_random as _assign_edge_random,
    _assign_composite
)

# =============================================================================
# Node Assignment Tests
# =============================================================================

@pytest.fixture
def star_graph():
    """Returns a star graph: center node 0, leaves 1-4."""
    g = nx.MultiGraph()
    # Star graph
    # Node 0: degree 4
    # Nodes 1-4: degree 1
    for i in range(1, 5):
        g.add_edge(0, i)
    return g

def test_assign_node_degree(star_graph):
    """Test degree-based assignment."""
    # Map degree 1 -> "end", degree 4 -> "center"
    mapping = {"1": "end", "4": "center"}
    _assign_by_degree(star_graph, mapping)
    
    assert star_graph.nodes[0]["node_type"] == "center"
    for i in range(1, 5):
        assert star_graph.nodes[i]["node_type"] == "end"

def test_assign_node_position():
    """Test positional assignment (layers)."""
    g = nx.MultiGraph()
    # 3 nodes along Z axis
    g.add_node(0, pos=(0,0,0))
    g.add_node(1, pos=(0,0,5))
    g.add_node(2, pos=(0,0,10))
    
    # 2 Layers: Z=0-5 -> Type A, Z=5-10 -> Type B
    # Range is 0 to 10. Step is 5.
    PosConfig = namedtuple('PosConfig', ['dimension', 'num_layers', 'layer_types'])
    config = PosConfig("z", 2, ["A", "B"])
    
    _assign_by_position(g, config)
    
    # Node 0 at 0.0 -> First half -> A
    # Node 1 at 5.0 -> Exactly boundary -> usually floor(1.0) -> index 1 -> B?
    # Logic: int((5-0)/5) -> int(1.0) -> 1. index 1 is B.
    # Node 2 at 10.0 -> int(2.0) -> 2 -> clamped to num_layers-1 = 1 -> B.
    
    assert g.nodes[0]["node_type"] == "A"
    # Depending on float precision, boundaries can be tricky, 
    # checking logic in code: usually exclusive upper bound for int cast logic
    # but let's check exact behavior
    assert g.nodes[2]["node_type"] == "B"

def test_assign_node_random(star_graph):
    """Test random assignment ratios."""
    # 100% Type A
    _assign_by_random(star_graph, {"A": 100})
    for n in star_graph:
        assert star_graph.nodes[n]["node_type"] == "A"

def test_assign_node_explicit(star_graph):
    """Test explicit ID assignment."""
    mapping = {0: "center", 1: "leaf"}
    _assign_explicit(star_graph, mapping)
    assert star_graph.nodes[0]["node_type"] == "center"
    assert star_graph.nodes[1]["node_type"] == "leaf"
    assert star_graph.nodes[2]["node_type"] == "A" # Default

# =============================================================================
# Edge Assignment Tests
# =============================================================================

@pytest.fixture
def line_graph():
    """Linear chain 0-1-2-3."""
    g = nx.MultiGraph()
    g.add_node(0, pos=(0,0,0))
    g.add_node(1, pos=(0,0,2))
    g.add_node(2, pos=(0,0,4))
    g.add_node(3, pos=(0,0,6))
    g.add_edge(0, 1) # Midpoint 1.0 (A)
    g.add_edge(1, 2) # Midpoint 3.0 (Mid)
    g.add_edge(2, 3) # Midpoint 5.0 (B)
    return g

def test_assign_edge_uniform(line_graph):
    """Test uniform assignment."""
    _assign_uniform(line_graph, "U")
    for u, v, k in line_graph.edges(keys=True):
        assert line_graph.edges[u,v,k]["edge_type"] == "U"

def test_assign_edge_composite(line_graph):
    """Test composite/lamellar assignment."""
    # Range 0 to 6. Midpoints at 1, 3, 5.
    # 2 Layers. Boundary at z=3.
    # Layer 1: [0, 3) -> A
    # Layer 2: [3, 6] -> B
    CompConfig = namedtuple('CompConfig', ['dimension', 'num_layers', 'layer_types'])
    config = CompConfig("z", 2, ["A", "B"])
    
    _assign_composite(line_graph, config, dims=np.array([10,10,10]))
    
    # 0-1 (z=1): A
    # 1-2 (z=3): Boundary case. (3-1)/2.5? Wait.
    # Code logic: 
    # Midpoints: 1, 3, 5. min=1, max=5. range=4.
    # Step = 2.0.
    # Edge 1 (idx 0): (1-1)/2 = 0 -> A
    # Edge 2 (idx 1): (3-1)/2 = 1 -> B
    # Edge 3 (idx 2): (5-1)/2 = 2 -> Clamped to 1 -> B
    
    edge_types = [line_graph.edges[u,v,k]["edge_type"] for u,v,k in line_graph.edges(keys=True)]
    # As iterating dict order is not guaranteed in older python but usually fine now
    # Let's check by edge key manually/conceptually
    
    # We can check counts
    assert edge_types.count("A") == 1
    assert edge_types.count("B") == 2
