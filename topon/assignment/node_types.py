"""
Node type assignment for Topon.

Assigns abstract types (e.g., 'A', 'B', 'end') to nodes based on:
- Degree mapping
- Positional (layer-based)
- Random assignment
- Explicit per-node assignment
"""

import random
from typing import Optional
import networkx as nx
import numpy as np

from topon.config.schema import NodeTypesConfig


def assign_node_types(G: nx.MultiGraph, config: NodeTypesConfig) -> None:
    """
    Assign node types to all nodes in the graph.
    
    Args:
        G: Graph to modify (in place).
        config: Node types configuration.
    """
    method = config.method
    
    if method == "degree":
        _assign_by_degree(G, config.degree.mapping)
    elif method == "positional":
        _assign_by_position(G, config.positional)
    elif method == "random":
        _assign_by_random(G, config.random.type_ratios)
    elif method == "explicit":
        _assign_explicit(G, config.explicit)
    else:
        raise ValueError(f"Unknown node type assignment method: {method}")
    
    # Count results
    type_counts = {}
    for node, data in G.nodes(data=True):
        node_type = data.get("node_type", "unassigned")
        type_counts[node_type] = type_counts.get(node_type, 0) + 1
    
    print(f"    Method: {method}")
    for t, count in sorted(type_counts.items()):
        print(f"      {t}: {count}")


def _assign_by_degree(G: nx.MultiGraph, mapping: dict[str, str]) -> None:
    """
    Assign node types based on node degree.
    
    Args:
        G: Graph to modify.
        mapping: Dict mapping degree (as string) to node type.
    """
    for node in G.nodes():
        degree = G.degree(node)
        degree_str = str(degree)
        
        if degree_str in mapping:
            node_type = mapping[degree_str]
        else:
            # Default: use 'A' for any unmapped degree
            node_type = "A"
        
        G.nodes[node]["node_type"] = node_type


def _assign_by_position(G: nx.MultiGraph, config) -> None:
    """
    Assign node types based on spatial position (layers).
    
    Args:
        G: Graph to modify.
        config: Positional config with dimension, num_layers, layer_types.
    """
    dim_map = {"x": 0, "y": 1, "z": 2}
    dim_idx = dim_map[config.dimension]
    
    # Get position range
    positions = []
    for node, data in G.nodes(data=True):
        if "pos" in data:
            positions.append(data["pos"][dim_idx])
    
    if not positions:
        print("    Warning: No node positions found, using degree method fallback")
        _assign_by_degree(G, {"1": "end"})
        return
    
    min_pos = min(positions)
    max_pos = max(positions)
    range_size = (max_pos - min_pos) / config.num_layers if max_pos > min_pos else 1
    
    for node, data in G.nodes(data=True):
        if "pos" in data:
            pos = data["pos"][dim_idx]
            layer_idx = min(int((pos - min_pos) / range_size), config.num_layers - 1)
            node_type = config.layer_types[layer_idx % len(config.layer_types)]
        else:
            node_type = config.layer_types[0]
        
        G.nodes[node]["node_type"] = node_type


def _assign_by_random(G: nx.MultiGraph, type_ratios: dict[str, float]) -> None:
    """
    Assign node types randomly based on ratios.
    
    Args:
        G: Graph to modify.
        type_ratios: Dict mapping type to weight.
    """
    types = list(type_ratios.keys())
    weights = list(type_ratios.values())
    total = sum(weights)
    probabilities = [w / total for w in weights]
    
    for node in G.nodes():
        node_type = random.choices(types, weights=probabilities, k=1)[0]
        G.nodes[node]["node_type"] = node_type


def _assign_explicit(G: nx.MultiGraph, explicit_map: dict[int, str]) -> None:
    """
    Assign node types from explicit per-node mapping.
    
    Args:
        G: Graph to modify.
        explicit_map: Dict mapping node ID to type.
    """
    for node in G.nodes():
        if node in explicit_map:
            node_type = explicit_map[node]
        else:
            # Default for unmapped nodes
            node_type = "A"
        
        G.nodes[node]["node_type"] = node_type
