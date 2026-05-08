"""
Graph loader for Topon.

Handles loading topology from various file formats.
"""

import pickle
from pathlib import Path
from typing import Optional, Union

import networkx as nx
import numpy as np
import pandas as pd


def load_graph(
    gpickle_path: Optional[Union[str, Path]] = None,
    nodes_path: Optional[Union[str, Path]] = None,
    edges_path: Optional[Union[str, Path]] = None,
) -> tuple[nx.MultiGraph, Optional[np.ndarray]]:
    """
    Load a topology graph from file(s).
    
    Args:
        gpickle_path: Path to a .gpickle file (takes precedence).
        nodes_path: Path to a .nodes file.
        edges_path: Path to a .edges file.
        
    Returns:
        Tuple of (NetworkX MultiGraph, box dimensions array or None).
        
    Raises:
        ValueError: If no valid file paths provided.
        FileNotFoundError: If specified files don't exist.
    """
    if gpickle_path:
        return _load_from_gpickle(gpickle_path)
    elif nodes_path and edges_path:
        return _load_from_nodes_edges(nodes_path, edges_path)
    else:
        raise ValueError(
            "Must provide either gpickle_path or both nodes_path and edges_path"
        )


def _load_from_gpickle(path: Union[str, Path]) -> tuple[nx.MultiGraph, Optional[np.ndarray]]:
    """
    Load graph from a .gpickle file.
    
    Args:
        path: Path to .gpickle file.
        
    Returns:
        Tuple of (graph, dims).
    """
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"Gpickle file not found: {path}")
    
    with open(path, "rb") as f:
        data = pickle.load(f)
    
    # Handle different gpickle formats
    if isinstance(data, tuple) and len(data) == 2:
        # Format: (graph, dims)
        G, dims = data
    elif isinstance(data, nx.Graph):
        # Format: just the graph
        G = data
        dims = _infer_dims_from_graph(G)
    elif isinstance(data, dict) and "graph" in data:
        # Format: dict with 'graph' and optionally 'dims'
        G = data["graph"]
        dims = data.get("dims")
    else:
        raise ValueError(f"Unrecognized gpickle format in {path}")
    
    # Ensure it's a MultiGraph
    if not isinstance(G, nx.MultiGraph):
        G = nx.MultiGraph(G)
    
    # Remove vacancies (degree-0 nodes)
    n_removed = _remove_vacancies(G)
    
    print(f"Loaded graph from {path.name}")
    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    if n_removed > 0:
        print(f"  Removed {n_removed} vacancies (degree-0 nodes)")
    
    return G, dims


def _load_from_nodes_edges(
    nodes_path: Union[str, Path],
    edges_path: Union[str, Path]
) -> tuple[nx.MultiGraph, Optional[np.ndarray]]:
    """
    Load graph from .nodes and .edges files.
    
    File formats:
    - .nodes: NodeID X Y Z Degree (whitespace-separated, # comments)
    - .edges: Node1 Node2 (whitespace-separated, # comments)
    
    Args:
        nodes_path: Path to .nodes file.
        edges_path: Path to .edges file.
        
    Returns:
        Tuple of (graph, dims).
    """
    nodes_path = Path(nodes_path)
    edges_path = Path(edges_path)
    
    if not nodes_path.exists():
        raise FileNotFoundError(f"Nodes file not found: {nodes_path}")
    if not edges_path.exists():
        raise FileNotFoundError(f"Edges file not found: {edges_path}")
    
    # Load nodes
    nodes_df = pd.read_csv(
        nodes_path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["id", "x", "y", "z", "degree"]
    )
    
    # Load edges
    edges_df = pd.read_csv(
        edges_path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["u", "v"]
    )
    
    # Build graph
    G = nx.MultiGraph()
    
    for _, row in nodes_df.iterrows():
        G.add_node(
            int(row["id"]),
            pos=(float(row["x"]), float(row["y"]), float(row["z"]))
        )
    
    for _, row in edges_df.iterrows():
        u, v = int(row["u"]), int(row["v"])
        if G.has_node(u) and G.has_node(v):
            G.add_edge(u, v)
    
    # Infer dimensions from positions
    dims = _infer_dims_from_graph(G)
    
    # Remove vacancies (degree-0 nodes)
    n_removed = _remove_vacancies(G)
    
    print(f"Loaded graph from {nodes_path.name} + {edges_path.name}")
    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    if n_removed > 0:
        print(f"  Removed {n_removed} vacancies (degree-0 nodes)")
    
    return G, dims


def _remove_vacancies(G: nx.Graph) -> int:
    """
    Remove degree-0 nodes (vacancies) from graph.
    
    Vacancies are lattice positions with no edges - they should not
    become atoms in the simulation.
    
    Args:
        G: Graph to modify in-place.
        
    Returns:
        Number of nodes removed.
    """
    vacancies = [n for n in G.nodes() if G.degree(n) == 0]
    if vacancies:
        G.remove_nodes_from(vacancies)
    return len(vacancies)


def _infer_dims_from_graph(G: nx.Graph) -> Optional[np.ndarray]:
    """
    Infer box dimensions from node positions.
    
    Args:
        G: Graph with 'pos' node attributes.
        
    Returns:
        Box dimensions as numpy array, or None if no positions.
    """
    positions = []
    for node, data in G.nodes(data=True):
        if "pos" in data:
            positions.append(data["pos"])
    
    if not positions:
        return None
    
    positions = np.array(positions)
    
    # Assume box starts at 0, dimensions are max + 1 (for lattice spacing)
    # This is a heuristic - actual dims should be stored in gpickle
    max_pos = positions.max(axis=0)
    min_pos = positions.min(axis=0)
    
    # For integer lattice positions, dims = max - min + 1
    dims = max_pos - min_pos + 1
    
    return dims


def save_graph(
    G: nx.Graph,
    output_path: Union[str, Path],
    dims: Optional[np.ndarray] = None
) -> None:
    """
    Save graph to a .gpickle file.
    
    Args:
        G: NetworkX graph to save.
        output_path: Path for output file.
        dims: Optional box dimensions to save with graph.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    data = (G, dims) if dims is not None else G
    
    with open(output_path, "wb") as f:
        pickle.dump(data, f)
    
    print(f"Saved graph to {output_path}")


def get_node_positions(G: nx.Graph) -> dict[int, np.ndarray]:
    """
    Extract node positions as a dictionary.
    
    Args:
        G: Graph with 'pos' node attributes.
        
    Returns:
        Dict mapping node ID to position array.
    """
    positions = {}
    for node, data in G.nodes(data=True):
        if "pos" in data:
            positions[node] = np.array(data["pos"])
    return positions
