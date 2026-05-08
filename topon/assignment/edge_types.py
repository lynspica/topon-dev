"""
Edge type assignment for Topon.

Assigns abstract types (e.g., 'A', 'B') to edges based on:
- Uniform (all same type)
- Random assignment
- Composite/lamellar (positional)
"""

import random
from typing import Optional
import networkx as nx
import numpy as np

from topon.config.schema import EdgeTypesConfig


def assign_edge_types(
    G: nx.MultiGraph, 
    config: EdgeTypesConfig, 
    dims: Optional[np.ndarray] = None
) -> None:
    """
    Assign edge types to all edges in the graph.
    
    Args:
        G: Graph to modify (in place).
        config: Edge types configuration.
        dims: Box dimensions for MIC calculations.
    """
    method = config.method
    
    if method == "uniform":
        _assign_uniform(G, config.uniform.type)
    elif method == "random":
        _assign_random(G, config.random.type_ratios)
    elif method == "composite":
        _assign_composite(G, config.composite, dims)
    else:
        raise ValueError(f"Unknown edge type assignment method: {method}")
    
    # Count results
    type_counts = {}
    for u, v, key, data in G.edges(keys=True, data=True):
        edge_type = data.get("edge_type", "unassigned")
        type_counts[edge_type] = type_counts.get(edge_type, 0) + 1
    
    print(f"    Method: {method}")
    for t, count in sorted(type_counts.items()):
        print(f"      {t}: {count}")


def _assign_uniform(G: nx.MultiGraph, edge_type: str) -> None:
    """
    Assign the same type to all edges.
    
    Args:
        G: Graph to modify.
        edge_type: Type to assign.
    """
    for u, v, key in G.edges(keys=True):
        G.edges[u, v, key]["edge_type"] = edge_type


def _assign_random(G: nx.MultiGraph, type_ratios: dict[str, float]) -> None:
    """
    Assign edge types randomly based on ratios.
    
    Args:
        G: Graph to modify.
        type_ratios: Dict mapping type to weight.
    """
    types = list(type_ratios.keys())
    weights = list(type_ratios.values())
    total = sum(weights)
    probabilities = [w / total for w in weights]
    
    for u, v, key in G.edges(keys=True):
        edge_type = random.choices(types, weights=probabilities, k=1)[0]
        G.edges[u, v, key]["edge_type"] = edge_type


def _assign_composite(
    G: nx.MultiGraph, 
    config, 
    dims: Optional[np.ndarray]
) -> None:
    """
    Assign edge types based on edge midpoint position (composite/lamellar).
    
    Args:
        G: Graph to modify.
        config: Composite config with dimension, num_layers, layer_types.
        dims: Box dimensions for MIC.
    """
    dim_map = {"x": 0, "y": 1, "z": 2}
    dim_idx = dim_map[config.dimension]
    
    # Calculate edge midpoints
    midpoints = []
    edges = list(G.edges(keys=True))
    
    for u, v, key in edges:
        pos_u = G.nodes[u].get("pos")
        pos_v = G.nodes[v].get("pos")
        
        if pos_u is None or pos_v is None:
            midpoints.append(0)
            continue
        
        pos_u = np.array(pos_u)
        pos_v = np.array(pos_v)
        
        # Apply MIC for midpoint calculation
        if dims is not None:
            vec = pos_v - pos_u
            vec = vec - dims * np.round(vec / dims)
            midpoint = pos_u + 0.5 * vec
        else:
            midpoint = 0.5 * (pos_u + pos_v)
        
        midpoints.append(midpoint[dim_idx])
    
    if not midpoints or all(m == 0 for m in midpoints):
        print("    Warning: No edge positions found, using uniform fallback")
        _assign_uniform(G, config.layer_types[0])
        return
    
    # Determine layer boundaries
    min_pos = min(midpoints)
    max_pos = max(midpoints)
    range_size = (max_pos - min_pos) / config.num_layers if max_pos > min_pos else 1
    
    # Assign types
    for (u, v, key), midpoint in zip(edges, midpoints):
        layer_idx = min(int((midpoint - min_pos) / range_size), config.num_layers - 1)
        edge_type = config.layer_types[layer_idx % len(config.layer_types)]
        G.edges[u, v, key]["edge_type"] = edge_type
