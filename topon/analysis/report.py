"""
Graph analysis for Topon.

Analyzes the graph to report:
- Node/edge counts
- Degree distribution
- Max possible defects (primary/secondary loops)
- Max possible entanglements
"""

import networkx as nx
import numpy as np
from typing import Optional


def analyze_graph(
    G: nx.MultiGraph, 
    dims: Optional[np.ndarray] = None,
    verbose: bool = True
) -> dict:
    """
    Analyze the graph and return a report.
    
    Args:
        G: Graph to analyze.
        dims: Box dimensions.
        verbose: Print report to stdout.
        
    Returns:
        Analysis dictionary.
    """
    analysis = {}
    
    # Basic counts
    analysis["num_nodes"] = G.number_of_nodes()
    analysis["num_edges"] = G.number_of_edges()
    
    # Degree distribution
    degrees = [d for _, d in G.degree()]
    degree_counts = {}
    for d in range(max(degrees) + 1 if degrees else 1):
        count = degrees.count(d)
        if count > 0:
            degree_counts[d] = count
    analysis["degree_distribution"] = degree_counts
    
    # Node type counts after assignment
    node_types = {}
    for _, data in G.nodes(data=True):
        nt = data.get("node_type")
        if nt:
            node_types[nt] = node_types.get(nt, 0) + 1
    analysis["node_types"] = node_types
    
    # Edge type counts after assignment
    edge_types = {}
    for _, _, _, data in G.edges(keys=True, data=True):
        et = data.get("edge_type")
        if et:
            edge_types[et] = edge_types.get(et, 0) + 1
    analysis["edge_types"] = edge_types
    
    # Max possible primary loops (edges where endpoints share another edge)
    analysis["max_primary_loops"] = count_primary_loop_candidates(G)
    
    # Max possible secondary loops (parallel edges)
    analysis["max_secondary_loops"] = count_secondary_loop_candidates(G)
    
    # Max possible entanglements (disjoint parallel edge pairs)
    analysis["max_entanglements"] = count_entanglement_candidates(G, dims)
    
    if verbose:
        print_analysis(analysis)
    
    return analysis


def count_primary_loop_candidates(G: nx.MultiGraph) -> int:
    """
    Count edges that could become primary loops.
    
    A primary loop is when an edge's endpoints share another edge,
    creating a triangle-like defect.
    """
    count = 0
    processed = set()
    
    for u, v, key in G.edges(keys=True):
        edge_id = (min(u, v), max(u, v), key)
        if edge_id in processed:
            continue
        processed.add(edge_id)
        
        # Check if u and v share any other neighbors
        neighbors_u = set(G.neighbors(u)) - {v}
        neighbors_v = set(G.neighbors(v)) - {u}
        shared = neighbors_u & neighbors_v
        
        if shared:
            count += 1
    
    return count


def count_secondary_loop_candidates(G: nx.MultiGraph) -> int:
    """
    Count potential secondary loop locations.
    
    A secondary loop would be parallel edges between same node pair.
    For regular graphs, this counts potential locations, not existing ones.
    """
    # For a MultiGraph, just count existing multi-edges
    count = 0
    processed = set()
    
    for u, v in G.edges():
        pair = (min(u, v), max(u, v))
        if pair in processed:
            continue
        processed.add(pair)
        
        # Count edges between u and v
        num_edges = G.number_of_edges(u, v)
        if num_edges > 1:
            count += num_edges - 1  # Additional edges are secondary loops
    
    return count


def count_entanglement_candidates(
    G: nx.MultiGraph, 
    dims: Optional[np.ndarray] = None
) -> int:
    """
    Estimate max possible entanglement pairs.
    
    Counts disjoint edge pairs where both endpoints have degree > 1.
    This is an upper bound estimate.
    """
    # Get valid edges (both endpoints have degree > 1)
    valid_edges = []
    for u, v, key in G.edges(keys=True):
        if G.degree(u) > 1 and G.degree(v) > 1:
            valid_edges.append((u, v, key))
    
    # Count disjoint pairs (no shared nodes)
    count = 0
    for i, (u1, v1, k1) in enumerate(valid_edges):
        nodes1 = {u1, v1}
        for j in range(i + 1, len(valid_edges)):
            u2, v2, k2 = valid_edges[j]
            nodes2 = {u2, v2}
            if nodes1.isdisjoint(nodes2):
                count += 1
    
    # This can be very large, so we return actual candidate count
    # from the find_crossing_candidates algorithm instead
    # For now, return a reasonable estimate
    return min(count, len(valid_edges) // 2)


def print_analysis(analysis: dict) -> None:
    """Print formatted analysis report."""
    print("=" * 50)
    print("GRAPH ANALYSIS REPORT")
    print("=" * 50)
    print(f"Nodes: {analysis['num_nodes']}")
    print(f"Edges: {analysis['num_edges']}")
    print()
    print("Degree Distribution:")
    for d, count in sorted(analysis["degree_distribution"].items()):
        print(f"  d={d}: {count}")
    print()
    if analysis["node_types"]:
        print("Node Types:")
        for t, count in sorted(analysis["node_types"].items()):
            print(f"  {t}: {count}")
        print()
    if analysis["edge_types"]:
        print("Edge Types:")
        for t, count in sorted(analysis["edge_types"].items()):
            print(f"  {t}: {count}")
        print()
    print("Defect/Entanglement Potential:")
    print(f"  Max primary loops: {analysis['max_primary_loops']}")
    print(f"  Max secondary loops: {analysis['max_secondary_loops']}")
    print(f"  Max entanglement pairs: {analysis['max_entanglements']}")
    print("=" * 50)
