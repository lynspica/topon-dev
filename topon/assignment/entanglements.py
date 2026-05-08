"""
Entanglement selection for Topon.

Selects pairs of edges for parametric entanglement based on:
- Nearest disjoint neighbor algorithm (parallel edges closest to each other)
- Target count or percentage
"""

import random
from typing import Optional
import networkx as nx
import numpy as np

from topon.config.schema import EntanglementsConfig


def select_entanglements(
    G: nx.MultiGraph, 
    config: EntanglementsConfig,
    dims: Optional[np.ndarray] = None,
    max_possible: Optional[int] = None,
    candidates: Optional[list] = None,
    num_chains: Optional[int] = None
) -> list[tuple[tuple, tuple, int]]:
    """
    Select entanglement pairs from the graph.
    
    Args:
        G: Graph with node positions.
        config: Entanglement configuration.
        dims: Box dimensions for MIC.
        max_possible: Max possible entanglements (from analysis).
        candidates: Optional pre-computed candidates.
        num_chains: Number of chains (for distribution mode).
        
    Returns:
        List of ((u1, v1, key1), (u2, v2, key2), count) tuples.
        In strict mode, count is always 1.
    """
    if not config.enabled:
        return []
    
    # Find candidate pairs
    if candidates is None:
        candidates = find_crossing_candidates(G, dims)
    
    if not candidates:
        print("    No entanglement candidates found")
        return []
    
    # Helper to calculate center
    def get_center(c):
        e1, e2 = c
        u1, v1, _ = e1
        u2, v2, _ = e2
        p1u = np.array(G.nodes[u1]['pos'])
        p1v = np.array(G.nodes[v1]['pos'])
        p2u = np.array(G.nodes[u2]['pos'])
        p2v = np.array(G.nodes[v2]['pos'])
        
        # Midpoint 1
        vec1 = p1v - p1u
        if dims is not None: vec1 = vec1 - dims * np.round(vec1 / dims)
        m1 = p1u + 0.5 * vec1
        if dims is not None: m1 = m1 - dims * np.floor(m1 / dims)
        
        # Midpoint 2
        vec2 = p2v - p2u
        if dims is not None: vec2 = vec2 - dims * np.round(vec2 / dims)
        m2 = p2u + 0.5 * vec2
        if dims is not None: m2 = m2 - dims * np.floor(m2 / dims)
        
        # Center
        diff = m2 - m1
        if dims is not None: diff = diff - dims * np.round(diff / dims)
        center = m1 + 0.5 * diff
        if dims is not None: center = center - dims * np.floor(center / dims)
        
        return center

    # ============ DISTRIBUTION MODE ============
    if config.avg_crosslinks_per_chain is not None:
        if num_chains is None:
            num_chains = G.number_of_edges()
        
        total_draws = int(config.avg_crosslinks_per_chain * 0.5 * num_chains)
        print(f"    Distribution mode: {config.avg_crosslinks_per_chain} avg crosslinks/chain")
        print(f"    Total draws to distribute: {total_draws}")
        
        # Track kink data
        kink_candidates = {}  # kink_idx -> (e1, e2)
        kink_counts = {}      # kink_idx -> count
        kink_centers = {}     # kink_idx -> center
        
        # Track which edges are locked to which kink
        edge_to_kink = {}
        
        # Valid candidate indices (initially all)
        valid_candidates = set(range(len(candidates)))
        
        min_dist_sq = 1e-4
        
        for draw in range(total_draws):
            # Build draw pool: valid candidates + existing kink indices (as negative numbers)
            # We use negative indices for existing kinks to distinguish from candidate indices
            draw_pool = list(valid_candidates) + [-k-1 for k in kink_candidates.keys()]
            
            if not draw_pool:
                print(f"    Warning: Empty draw pool at draw {draw}")
                break
            
            # Pick randomly from pool
            pick = random.choice(draw_pool)
            
            if pick >= 0:
                # Picked a valid candidate -> create new kink
                cand_idx = pick
                cand = candidates[cand_idx]
                e1, e2 = cand
                
                # Calculate center and check spatial exclusivity
                center = get_center(cand)
                collision = False
                for ec in kink_centers.values():
                    diff = center - ec
                    if dims is not None:
                        diff = diff - dims * np.round(diff / dims)
                    if np.dot(diff, diff) < min_dist_sq:
                        collision = True
                        break
                
                if collision:
                    # Remove from valid and retry
                    valid_candidates.discard(cand_idx)
                    continue
                
                # Create new kink
                kink_idx = len(kink_candidates)
                kink_candidates[kink_idx] = cand
                kink_counts[kink_idx] = 1
                kink_centers[kink_idx] = center
                edge_to_kink[e1] = kink_idx
                edge_to_kink[e2] = kink_idx
                
                # Remove this candidate from valid
                valid_candidates.discard(cand_idx)
                
                # Remove all candidates that share edges with this kink
                to_remove = set()
                for other_idx in valid_candidates:
                    oe1, oe2 = candidates[other_idx]
                    if oe1 == e1 or oe1 == e2 or oe2 == e1 or oe2 == e2:
                        to_remove.add(other_idx)
                valid_candidates -= to_remove
                
            else:
                # Picked an existing kink -> increment count
                kink_idx = -pick - 1
                kink_counts[kink_idx] += 1
        
        # Build result
        selected = []
        for kink_idx, cand in kink_candidates.items():
            count = kink_counts[kink_idx]
            selected.append((cand[0], cand[1], count))
        
        # Store in graph edges
        for (e1, e2, count) in selected:
            G.edges[e1]["entangled_with"] = e2
            G.edges[e1]["entanglement_count"] = count
            G.edges[e2]["entangled_with"] = e1
            G.edges[e2]["entanglement_count"] = count
        
        total_count = sum(kink_counts.values())
        print(f"    Created {len(selected)} unique kinks with total count {total_count}")
        return selected
    
    # ============ STRICT MODE (Legacy) ============
    # Determine target count
    if config.target_type == "percentage":
        max_val = max_possible if max_possible else len(candidates)
        target = int(config.target * max_val / 100)
    else:
        target = config.target
    
    # Randomly shuffle candidates for unbiased selection
    random.shuffle(candidates)
    
    selected = []
    used_edges = set()
    existing_centers = []
    min_dist_sq = 1e-4
    
    for cand in candidates:
        if len(selected) >= target:
            break
            
        e1, e2 = cand
        
        # 1. Edge Exclusivity Check
        if e1 in used_edges or e2 in used_edges:
            continue
            
        # 2. Location Exclusivity Check
        center = get_center(cand)
        collision = False
        for ec in existing_centers:
            diff = center - ec
            if dims is not None:
                diff = diff - dims * np.round(diff / dims)
            dist_sq = np.dot(diff, diff)
            
            if dist_sq < min_dist_sq:
                collision = True
                break
        
        if collision:
            continue
            
        # Select!
        selected.append((cand[0], cand[1], 1))  # count=1 in strict mode
        used_edges.add(e1)
        used_edges.add(e2)
        existing_centers.append(center)
    
    # Store in graph edges
    for (e1, e2, count) in selected:
        G.edges[e1]["entangled_with"] = e2
        G.edges[e2]["entangled_with"] = e1
    
    print(f"    Selected {len(selected)} entanglement pairs (strict mode)")
    return selected


def find_crossing_candidates(
    G: nx.MultiGraph, 
    dims: Optional[np.ndarray] = None
) -> list[tuple[tuple, tuple]]:
    """
    Find edge pairs suitable for entanglement using nearest disjoint neighbor.
    
    Criteria:
    - Both endpoints have degree > 1 (not dangling ends)
    - Edges don't share any nodes (disjoint)
    - Closest midpoint distance
    
    Args:
        G: Graph with node positions.
        dims: Box dimensions for MIC.
        
    Returns:
        List of (edge1, edge2) tuples where each edge is (u, v, key).
    """
    print("    Finding crossing candidates...")
    
    # Get edges with valid positions and degree > 1 endpoints
    edges = []
    midpoints = []
    edge_nodes = {}
    
    for u, v, key in G.edges(keys=True):
        # Check degree
        if G.degree(u) <= 1 or G.degree(v) <= 1:
            continue
        
        # Get positions
        pos_u = G.nodes[u].get("pos")
        pos_v = G.nodes[v].get("pos")
        
        if pos_u is None or pos_v is None:
            continue
        
        pos_u = np.array(pos_u)
        pos_v = np.array(pos_v)
        
        # Calculate midpoint with MIC
        if dims is not None:
            vec = pos_v - pos_u
            vec = vec - dims * np.round(vec / dims)
            midpoint = pos_u + 0.5 * vec
            # Wrap to box
            midpoint = midpoint - dims * np.floor(midpoint / dims)
        else:
            midpoint = 0.5 * (pos_u + pos_v)
        
        edge_key = (u, v, key)
        edges.append(edge_key)
        midpoints.append(midpoint)
        edge_nodes[edge_key] = {u, v}
    
    if len(edges) < 2:
        return []
    
    midpoints = np.array(midpoints)
    candidates = []
    processed = set()
    unique_geometries = set()
    
    # For each edge, find nearest disjoint neighbor
    for i, edge_a in enumerate(edges):
        nodes_a = edge_nodes[edge_a]
        mid_a = midpoints[i]
        
        # Calculate distances to all other edges
        best_dist = float('inf')
        best_match = None
        
        for j, edge_b in enumerate(edges):
            if i == j:
                continue
            
            nodes_b = edge_nodes[edge_b]
            
            # Must be disjoint (no shared nodes)
            if not nodes_a.isdisjoint(nodes_b):
                continue
            
            # Calculate distance with MIC
            mid_b = midpoints[j]
            if dims is not None:
                vec = mid_b - mid_a
                vec = vec - dims * np.round(vec / dims)
                dist = np.linalg.norm(vec)
            else:
                dist = np.linalg.norm(mid_b - mid_a)
            
            if dist < best_dist:
                best_dist = dist
                best_match = edge_b
        
        if best_match is not None:
            # Create a unique key based on the node sets of the two edges
            # (u1, v1) and (u2, v2)
            nodes_pair = frozenset([
                frozenset(nodes_a),
                frozenset(edge_nodes[best_match])
            ])
            
            pair = tuple(sorted([edge_a, best_match]))
            
            if pair not in processed and nodes_pair not in unique_geometries:
                candidates.append((edge_a, best_match))
                processed.add(pair)
                unique_geometries.add(nodes_pair)
    
    print(f"    Found {len(candidates)} candidate pairs")
    return candidates


def get_kink_params(config: EntanglementsConfig) -> dict:
    """Get kink parameters for entanglement geometry."""
    return {
        "overshoot": config.kink_params.overshoot,
        "z_amp": config.kink_params.z_amp,
        "sigma": config.kink_params.sigma,
    }
