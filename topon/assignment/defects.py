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


def get_eligible_pairs(
    G: nx.MultiGraph,
    max_degree: int = None,
    exclude_node_types: tuple = ("end",),
) -> List[Tuple[int, int]]:
    """
    Get node pairs that can receive a primary loop.

    Conditions:
    1. Pair has exactly one edge currently.
    2. If max_degree is set, both nodes must have degree < max_degree.
    3. Neither endpoint's `node_type` attribute is in `exclude_node_types`.
       Defaults to excluding end-cap nodes, which have only 1 free chemistry
       valence (e.g. a Si chain-cap with 3 methyls leaves 1 bond for the
       chain); adding a parallel edge there would over-valence at the
       chemistry stage. Primary loops physically represent a chain that
       returns to the same *junction* pair anyway, so end caps are also
       semantically wrong as endpoints.

    Returns:
        List of (u, v) tuples of eligible node pairs.
    """
    edge_pairs = {}
    for u, v in G.edges():
        key = (min(u, v), max(u, v))
        edge_pairs[key] = edge_pairs.get(key, 0) + 1

    candidates = []
    degrees = dict(G.degree())
    excluded = set(exclude_node_types or ())

    def node_type(n):
        attrs = G.nodes[n]
        return attrs.get("node_type", attrs.get("type", "A"))

    for (u, v), count in edge_pairs.items():
        if count != 1:
            continue

        if max_degree is not None:
            if degrees[u] >= max_degree or degrees[v] >= max_degree:
                continue

        if excluded and (node_type(u) in excluded or node_type(v) in excluded):
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
    
    # Random selection — shuffle the eligibility list but re-check the
    # degree constraint at injection time. A single node can appear in
    # multiple eligible pairs (e.g. u-v1 and u-v2); injecting both bumps
    # u's degree twice. Without a per-injection re-check we silently
    # over-valence at the chemistry stage. See defect-demo P1-K notes.
    random.shuffle(eligible)
    candidates = eligible  # try the full shuffled list, stop when budget hit

    injected = 0
    for u, v in candidates:
        if injected >= num_to_inject:
            break
        if max_degree is not None:
            if G.degree(u) >= max_degree or G.degree(v) >= max_degree:
                continue
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
