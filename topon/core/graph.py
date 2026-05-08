"""
topon.core.graph
================
Canonical data types shared across all pipeline stages.

NetworkGraph is the output of Stage 1 (topology) and the input to Stage 2
(chemistry). It wraps a NetworkX MultiGraph and carries box dimensions.

Rules:
- Topology engines must return a NetworkGraph.
- Chemistry backends must accept a NetworkGraph.
- Nothing outside core/ should define new inter-stage types.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import networkx as nx
import numpy as np


@dataclass
class NetworkGraph:
    """
    Polymer network topology: a graph of crosslink nodes connected by chains.

    Attributes
    ----------
    graph : nx.MultiGraph
        Nodes are crosslink junctions. Node attribute ``pos`` is a 3-tuple of
        lattice coordinates (integers or floats).
        Edges are polymer chains. Edge attribute ``dp`` is the degree of
        polymerization (number of repeat units on that chain).
    dims : np.ndarray, shape (3,)
        Simulation box dimensions in lattice units.  Multiply by a scale
        factor (computed during chemistry embedding) to get Angstroms.
    metadata : dict
        Arbitrary key/value info (source file, generation parameters, etc.).
        Not used by pipeline stages — for logging/provenance only.

    Notes
    -----
    MultiGraph is required (not Graph) because two nodes can share more than
    one chain (parallel edges).
    """

    graph: nx.MultiGraph
    dims: np.ndarray
    metadata: dict = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Convenience accessors — keep workflow code readable
    # ------------------------------------------------------------------

    @property
    def n_nodes(self) -> int:
        return self.graph.number_of_nodes()

    @property
    def n_edges(self) -> int:
        return self.graph.number_of_edges()

    def node_positions(self) -> dict[int, np.ndarray]:
        """Return {node_id: np.array([x, y, z])} for all nodes with 'pos'."""
        return {
            n: np.asarray(data["pos"])
            for n, data in self.graph.nodes(data=True)
            if "pos" in data
        }

    def edges_with_dp(self) -> list[tuple[int, int, int]]:
        """Return [(u, v, dp), ...] for all edges."""
        return [
            (u, v, data.get("dp", 1))
            for u, v, data in self.graph.edges(data=True)
        ]

    # ------------------------------------------------------------------
    # Interop with legacy code that works directly on (G, dims) tuples
    # ------------------------------------------------------------------

    @classmethod
    def from_legacy(
        cls,
        G: nx.MultiGraph,
        dims: np.ndarray,
        **metadata,
    ) -> "NetworkGraph":
        """Wrap the (G, dims) tuple returned by topology.loader.load_graph."""
        return cls(graph=G, dims=np.asarray(dims), metadata=metadata)

    def to_legacy(self) -> tuple[nx.MultiGraph, np.ndarray]:
        """Return (G, dims) for code that hasn't been ported yet."""
        return self.graph, self.dims

    def __repr__(self) -> str:
        return (
            f"NetworkGraph(nodes={self.n_nodes}, edges={self.n_edges}, "
            f"dims={self.dims})"
        )
