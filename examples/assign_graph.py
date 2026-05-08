"""
Assign Graph Types Example
==========================

This script demonstrates the two-stage workflow:
1. Load raw graph + assignment rules → Create attributed graph
2. Use attributed graph for chemistry/simulation (in generate_cg.py)

Usage:
    python assign_graph.py
"""

import sys
from pathlib import Path

# Add package root
pkg_dir = Path(__file__).parent.parent
sys.path.insert(0, str(pkg_dir))

from topon.topology.loader import load_graph
from topon.assignment.attributor import GraphAttributor

# === CONFIGURATION ===
GRAPH_PATH = pkg_dir / "tests/sample_graphs/network_N5x5x5_trial3.nodes"
EDGES_PATH = pkg_dir / "tests/sample_graphs/network_N5x5x5_trial3.edges"

NODE_ASSIGNMENT = Path(__file__).parent / "defaults/node_degree.json"
EDGE_ASSIGNMENT = Path(__file__).parent / "defaults/edge_uniform.json"

OUTPUT_DIR = Path(__file__).parent / "output_attributed"
DP = 20

def main():
    print("=" * 60)
    print("GRAPH ATTRIBUTION WORKFLOW")
    print("=" * 60)
    
    # 1. Load raw graph
    print("\n1. Loading raw graph...")
    G, dims = load_graph(nodes_path=str(GRAPH_PATH), edges_path=str(EDGES_PATH))
    print(f"   Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    
    # 2. Create attributor and apply rules
    print("\n2. Applying assignment rules...")
    attributor = GraphAttributor(G, dims)
    attributor.apply_node_assignment(NODE_ASSIGNMENT)
    attributor.apply_edge_assignment(EDGE_ASSIGNMENT)
    attributor.apply_dp(DP)
    
    # 3. Save attributed graph
    print("\n3. Saving attributed graph...")
    OUTPUT_DIR.mkdir(exist_ok=True)
    output_path = OUTPUT_DIR / "graph_attributed"
    attributor.save(output_path)
    
    # 4. Show summary
    print("\n" + "=" * 60)
    print("ATTRIBUTION COMPLETE")
    print("=" * 60)
    print(f"Output: {OUTPUT_DIR.resolve()}")
    print(f"  - graph_attributed.nodes")
    print(f"  - graph_attributed.edges")
    print("\nNext: Use these files in generate_cg.py or generate_atomistic.py")

if __name__ == "__main__":
    main()
