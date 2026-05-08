"""
Graph Attributor
=================
Applies node and edge type assignments to a raw graph,
producing an attributed graph ready for chemistry/simulation.
"""

import json
import networkx as nx
from pathlib import Path

from topon.assignment.entanglements import select_entanglements
from topon.config.schema import EntanglementsConfig


class GraphAttributor:
    """
    Reads assignment rules from JSON files and applies them to a graph.
    """
    
    def __init__(self, graph, dims=None):
        """
        Args:
            graph: NetworkX graph to attribute
            dims: Optional lattice dimensions [x, y, z]
        """
        self.G = graph.copy()
        self.dims = dims if dims is not None else [1.0, 1.0, 1.0]
        
    def apply_node_assignment(self, assignment_path):
        """
        Apply node type assignment rules from JSON file.
        
        Supported methods:
        - "degree": Assign based on node degree
        
        Args:
            assignment_path: Path to node assignment JSON
        """
        with open(assignment_path, 'r') as f:
            rules = json.load(f)
            
        method = rules.get('method', 'degree')
        
        if method == 'degree':
            mapping = rules.get('mapping', {})
            for node in self.G.nodes():
                degree = self.G.degree(node)
                node_type = mapping.get(str(degree), 'A')  # Default to 'A'
                self.G.nodes[node]['node_type'] = node_type
                
        print(f"  Applied node assignment ({method}): {self._count_node_types()}")
        return self
        
    def apply_edge_assignment(self, assignment_path):
        """
        Apply edge type assignment rules from JSON file.
        
        Supported methods:
        - "uniform": All edges same type
        
        Args:
            assignment_path: Path to edge assignment JSON
        """
        with open(assignment_path, 'r') as f:
            rules = json.load(f)
            
        method = rules.get('method', 'uniform')
        
        if method == 'uniform':
            edge_type = rules.get('type', 'A')
            for u, v, key, data in self.G.edges(keys=True, data=True):
                data['edge_type'] = edge_type
                
        print(f"  Applied edge assignment ({method}): {self._count_edge_types()}")
        return self
        
    def apply_dp(self, dp=20):
        """
        Apply degree of polymerization to all edges.
        """
        for u, v, key, data in self.G.edges(keys=True, data=True):
            data['dp'] = dp
        print(f"  Applied DP={dp} to {self.G.number_of_edges()} edges")
        return self
        
    def apply_entanglements(self, config_dict: dict):
        """
        Select and mark entangled edge pairs.
        
        Args:
            config_dict: Dictionary matching EntanglementsConfig structure
        """
        config = EntanglementsConfig(**config_dict)
        if not config.enabled:
            return self
        
        # Distribution mode uses avg_crosslinks_per_chain instead of target
        if config.avg_crosslinks_per_chain is not None:
            print(f"  Applying entanglements (distribution mode: {config.avg_crosslinks_per_chain} avg/chain)...")
        else:
            print(f"  Applying entanglements (target: {config.target} {config.target_type})...")
        
        # Pass num_chains for distribution mode
        num_chains = self.G.number_of_edges()
        select_entanglements(self.G, config, self.dims, num_chains=num_chains)
        return self
        
    def save(self, output_path):
        """
        Save the attributed graph.
        """
        output_path = Path(output_path)
        
        if output_path.suffix == '.gpickle':
            nx.write_gpickle(self.G, output_path)
        else:
            # Save as nodes + edges files
            nodes_path = output_path.with_suffix('.nodes')
            edges_path = output_path.with_suffix('.edges')
            
            with open(nodes_path, 'w') as f:
                for node in sorted(self.G.nodes()):
                    pos = self.G.nodes[node].get('pos', (0, 0, 0))
                    node_type = self.G.nodes[node].get('node_type', 'A')
                    f.write(f"{node} {pos[0]} {pos[1]} {pos[2]} {node_type}\n")
                    
            with open(edges_path, 'w') as f:
                for u, v, data in self.G.edges(data=True):
                    edge_type = data.get('edge_type', 'A')
                    dp = data.get('dp', 20)
                    f.write(f"{u} {v} {edge_type} {dp}\n")
                    
            print(f"  Saved: {nodes_path.name}, {edges_path.name}")
            
        return output_path
        
    def get_graph(self):
        """Return the attributed graph."""
        return self.G
        
    def _count_node_types(self):
        """Count nodes by type."""
        counts = {}
        for node in self.G.nodes():
            t = self.G.nodes[node].get('node_type', '?')
            counts[t] = counts.get(t, 0) + 1
        return counts
        
    def _count_edge_types(self):
        """Count edges by type."""
        counts = {}
        for u, v, data in self.G.edges(data=True):
            t = data.get('edge_type', '?')
            counts[t] = counts.get(t, 0) + 1
        return counts
