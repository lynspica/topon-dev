"""
topon.topology.network
======================
Stage 1 backend for polymer network graph generation.

Wraps the existing topology loader so it returns a NetworkGraph,
which is the canonical inter-stage type defined in topon.core.

The underlying generation logic lives in:
  topon/_legacy/topology/generator.py   (C-extension wrapper)
  topon/_legacy/topology/generator_python.py  (pure Python port)

Usage
-----
    from topon.topology.network import load

    ng = load(
        nodes_path="path/to/graph.nodes",
        edges_path="path/to/graph.edges",
    )
    # ng is a NetworkGraph — pass it to a chemistry backend
"""

from pathlib import Path
from typing import Optional, Union

from topon.core import NetworkGraph
from topon.topology.loader import load_graph   # existing loader, unchanged


def load(
    nodes_path: Optional[Union[str, Path]] = None,
    edges_path: Optional[Union[str, Path]] = None,
    gpickle_path: Optional[Union[str, Path]] = None,
    **metadata,
) -> NetworkGraph:
    """
    Load a polymer network topology from file and return a NetworkGraph.

    Accepts the same file formats as the legacy loader:
      - .nodes / .edges file pair
      - .gpickle file

    Parameters
    ----------
    nodes_path, edges_path : path to .nodes and .edges files
    gpickle_path : path to .gpickle file (takes precedence)
    **metadata : arbitrary provenance info stored in NetworkGraph.metadata

    Returns
    -------
    NetworkGraph
        Stage 1 output ready for a chemistry backend.
    """
    G, dims = load_graph(
        gpickle_path=gpickle_path,
        nodes_path=nodes_path,
        edges_path=edges_path,
    )
    return NetworkGraph.from_legacy(G, dims, **metadata)
