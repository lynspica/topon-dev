"""
Defect Injection Module for Topon.

Implements primary loop injection (parallel edges between connected node pairs).
"""

import random
from typing import Tuple, List
import networkx as nx


def count_primary_loops(G: nx.MultiGraph) -> int:
    """
    Count existing primary loops (parallel edges between same node pairs).
    
    Returns:
        Number of node pairs with more than one edge between them.
    """
    edge_pairs = {}
    for u, v in G.edges():
        key = (min(u, v), max(u, v))
        edge_pairs[key] = edge_pairs.get(key, 0) + 1
    
    return sum(1 for count in edge_pairs.values() if count > 1)


def get_eligible_pairs(G: nx.MultiGraph, max_degree: int = None) -> List[Tuple[int, int]]:
    """
    Get node pairs that can receive a primary loop.
    
    Conditions:
    1. Pair has exactly one edge currently.
    2. If max_degree is set, both nodes must have degree < max_degree.
    
    Returns:
        List of (u, v) tuples of eligible node pairs.
    """
    edge_pairs = {}
    for u, v in G.edges():
        key = (min(u, v), max(u, v))
        edge_pairs[key] = edge_pairs.get(key, 0) + 1
    
    candidates = []
    degrees = dict(G.degree())
    
    for (u, v), count in edge_pairs.items():
        if count != 1:
            continue
            
        # Check degree limits if specified
        if max_degree is not None:
            # We are adding an edge, so degree will increase by 1
            if degrees[u] >= max_degree or degrees[v] >= max_degree:
                continue
        
        candidates.append((u, v))
        
    return candidates


def inject_primary_loops(G: nx.MultiGraph, target: int, target_type: str = "count", 
                         inherit_dp: bool = True, max_degree: int = None) -> int:
    """
    Inject primary loops by adding parallel edges to existing connected node pairs.
    
    Args:
        G: NetworkX MultiGraph to modify IN PLACE.
        target: Number/percentage of loops to add.
        target_type: "count" or "percentage".
        inherit_dp: If True, new edges inherit 'dp' from the existing edge.
        max_degree: If Set, do not add edges to nodes already at/above this degree.
                    Useful for avoiding chemical valence violations (e.g. 4 for Si).
    
    Returns:
        Number of primary loops actually injected.
    """
    eligible = get_eligible_pairs(G, max_degree=max_degree)
    
    if not eligible:
        print(f"    Warning: No eligible pairs for primary loop injection (max_degree={max_degree}).")
        return 0
    
    # Calculate actual number to inject
    if target_type == "percentage":
        num_to_inject = max(1, int(len(eligible) * target / 100))
    else:
        num_to_inject = min(target, len(eligible))
    
    # Random selection
    random.shuffle(eligible)
    selected = eligible[:num_to_inject]
    
    injected = 0
    for u, v in selected:
        # Get attributes from existing edge
        existing_edge_data = G.get_edge_data(u, v)
        if existing_edge_data:
            first_key = list(existing_edge_data.keys())[0]
            attrs = existing_edge_data[first_key].copy()
            if not inherit_dp:
                attrs.pop('dp', None)
        else:
            attrs = {}
        
        attrs['is_primary_loop'] = True
        G.add_edge(u, v, **attrs)
        injected += 1
    
    print(f"    Injected {injected} primary loops (limit max_degree={max_degree}).")
    return injected


def analyze_primary_loop_potential(G: nx.MultiGraph, max_degree: int = None) -> dict:
    """
    Analyze the graph for primary loop injection potential.
    """
    eligible = get_eligible_pairs(G, max_degree=max_degree)
    existing = count_primary_loops(G)
    
    return {
        "max_possible_primary_loops": len(eligible),
        "existing_primary_loops": existing,
        "eligible_pairs": len(eligible),
        "constraints": {"max_degree": max_degree}
    }
