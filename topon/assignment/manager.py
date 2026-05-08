"""
Assignment Manager for Topon.

Orchestrates all graph attribute assignment operations:
- Node types
- Edge types
- DP distribution
- Defects
- Entanglements
- Grafts
- Copolymers
"""

from typing import Optional
import networkx as nx
import numpy as np

from topon.config.schema import AssignmentConfig
from topon.assignment import node_types, edge_types, dp_distribution


class AssignmentManager:
    """
    Manages assignment of attributes to graph nodes and edges.
    
    Assignment order:
    1. Node types (based on degree/position/random)
    2. Edge types (uniform/random/composite)
    3. DP distribution (per edge type)
    4. Defects (after types assigned)
    5. Entanglements (after defects)
    6. Grafts (per edge type)
    7. Copolymers (per edge type)
    """
    
    def __init__(self, G: nx.MultiGraph, dims: Optional[np.ndarray], config: AssignmentConfig):
        """
        Initialize the assignment manager.
        
        Args:
            G: NetworkX MultiGraph with node positions.
            dims: Box dimensions for periodic boundary calculations.
            config: Assignment configuration.
        """
        self.G = G
        self.dims = dims
        self.config = config
        
        # Analysis results (populated by analyze())
        self.analysis = {}
    
    def analyze(self) -> dict:
        """
        Analyze the graph to determine max possible modifications.
        
        Returns:
            Dict with analysis results.
        """
        print("Analyzing graph...")
        
        # Degree distribution
        degrees = [d for _, d in self.G.degree()]
        degree_counts = {}
        for d in range(max(degrees) + 1):
            degree_counts[d] = degrees.count(d)
        
        # Calculate max possible primary loops
        from topon.assignment import defects
        defect_analysis = defects.analyze_primary_loop_potential(self.G)
        
        self.analysis = {
            "num_nodes": self.G.number_of_nodes(),
            "num_edges": self.G.number_of_edges(),
            "degree_distribution": degree_counts,
            "max_primary_loops": defect_analysis["max_possible_primary_loops"],
            "existing_primary_loops": defect_analysis["existing_primary_loops"],
            "max_entanglements": None,  # spatial estimate not yet implemented; falls back to len(candidates)
        }
        
        self._print_analysis()
        return self.analysis
    
    def _print_analysis(self) -> None:
        """Print analysis report."""
        print(f"  Nodes: {self.analysis['num_nodes']}")
        print(f"  Edges: {self.analysis['num_edges']}")
        print(f"  Degree distribution:")
        for d, count in sorted(self.analysis['degree_distribution'].items()):
            if count > 0:
                print(f"    d={d}: {count}")
    
    def run(self) -> nx.MultiGraph:
        """
        Run all assignment operations.
        
        Returns:
            The modified graph with all attributes assigned.
        """
        print("Running assignment...")
        
        # 1. Assign node types
        self._assign_node_types()
        
        # 2. Assign edge types
        self._assign_edge_types()
        
        # 3. Assign DP distribution
        self._assign_dp()
        
        # 4. Apply defects (if enabled)
        if self.config.defects.primary_loops.enabled or self.config.defects.secondary_loops.enabled:
            self._apply_defects()
        
        # 5. Select entanglements (if enabled)
        if self.config.entanglements.enabled:
            self._select_entanglements()
        
        # 6. Assign grafts (if enabled)
        if self.config.grafts.enabled:
            self._assign_grafts()
        
        # 7. Assign copolymers (if enabled)
        if self.config.copolymer.enabled:
            self._assign_copolymers()
        
        return self.G
    
    def _assign_node_types(self) -> None:
        """Assign node types based on configuration."""
        print("  Assigning node types...")
        node_types.assign_node_types(self.G, self.config.node_types)
    
    def _assign_edge_types(self) -> None:
        """Assign edge types based on configuration."""
        print("  Assigning edge types...")
        edge_types.assign_edge_types(self.G, self.config.edge_types, self.dims)
    
    def _assign_dp(self) -> None:
        """Assign DP values to edges."""
        print("  Assigning DP values...")
        dp_distribution.assign_dp(self.G, self.config.dp_distribution)
    
    def _apply_defects(self) -> None:
        """Apply defect modifications."""
        print("  Applying defects...")
        from topon.assignment import defects
        
        if self.config.defects.primary_loops.enabled:
            target = self.config.defects.primary_loops.target
            target_type = self.config.defects.primary_loops.target_type
            defects.inject_primary_loops(self.G, target, target_type)
    
    def _select_entanglements(self) -> None:
        """Select entanglement pairs."""
        print("  Selecting entanglements...")
        from topon.assignment import entanglements
        self.entangled_pairs = entanglements.select_entanglements(
            self.G, 
            self.config.entanglements, 
            self.dims,
            self.analysis.get("max_entanglements")
        )
    
    def _assign_grafts(self) -> None:
        """Assign graft side-chain information to edges.

        For each edge whose ``edge_type`` has a graft config, randomly
        selects backbone positions for side-chain attachment based on
        ``graft_density`` and writes three attributes to the edge:

        * ``graft_positions`` — list of backbone bead indices (0-based)
        * ``graft_dp``        — number of beads per side chain
        * ``graft_monomer``   — monomer/bead-type name for side-chain beads
        """
        import random

        graft_cfg = self.config.grafts.per_edge_type
        if not graft_cfg:
            print("  Assigning grafts... (no per_edge_type config, skipping)")
            return

        print("  Assigning grafts...")
        total = 0
        for u, v, key, data in self.G.edges(keys=True, data=True):
            edge_type = data.get("edge_type", "A")
            conf = graft_cfg.get(edge_type)
            if conf is None:
                continue
            dp = data.get("dp", 1)
            density = conf.graft_density
            positions = [k for k in range(dp) if random.random() < density]
            self.G[u][v][key]["graft_positions"] = positions
            self.G[u][v][key]["graft_dp"] = conf.side_chain_dp
            self.G[u][v][key]["graft_monomer"] = conf.side_chain_monomer
            total += len(positions)

        print(f"    Grafts assigned: {total} side chains across "
              f"{self.G.number_of_edges()} edges")

    def _assign_copolymers(self) -> None:
        """Assign per-position monomer sequences to edges.

        For each edge whose ``edge_type`` has a copolymer config, generates
        a monomer sequence of length ``dp`` using
        :func:`topon.chemistry.sequences.generate_monomer_sequence` and
        writes it as the ``monomer_sequence`` edge attribute.
        """
        from topon.chemistry.sequences import generate_monomer_sequence

        cop_cfg = self.config.copolymer.per_edge_type
        if not cop_cfg:
            print("  Assigning copolymers... (no per_edge_type config, skipping)")
            return

        print("  Assigning copolymers...")
        assigned = 0
        for u, v, key, data in self.G.edges(keys=True, data=True):
            edge_type = data.get("edge_type", "A")
            conf = cop_cfg.get(edge_type)
            if conf is None:
                continue
            dp = data.get("dp", 1)
            seq_cfg = {
                "arrangement": conf.arrangement,
                "composition": [
                    {"monomer": c.monomer, "fraction": c.fraction}
                    for c in conf.composition
                ],
            }
            self.G[u][v][key]["monomer_sequence"] = generate_monomer_sequence(
                dp, seq_cfg, default_monomer=edge_type
            )
            assigned += 1

        print(f"    Copolymer sequences assigned to {assigned} edges")
