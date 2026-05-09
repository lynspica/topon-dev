"""End-linking topology generation via the pure-Python generator.

Generates a 6x6x6 SC-lattice polymer-network topology graph using the
Python port of the C generator (`topon.topology.generator_python`).
The same generator is invoked when ``topology.generator.exe_path`` is
``null`` in a ``topon generate`` config.

Usage:
    python run.py
"""
from __future__ import annotations

import time
from pathlib import Path

from topon.config.schema import GeneratorConfig
from topon.topology.generator_python import PythonTopologyGenerator
from topon.topology.loader import load_graph

OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


def main() -> None:
    print("--- End-linking topology generation (Python generator) ---")
    cfg = GeneratorConfig(
        lattice_size="6x6x6",
        lattice_type="SC",
        periodicity="111",
        max_functionality=4,
        degree_distribution="0:13,1:25",
        max_trials=1_000_000,
        max_saves=1,
    )
    t0 = time.perf_counter()
    gen = PythonTopologyGenerator(cfg)
    graphs = gen.generate(trials=20)
    dt = time.perf_counter() - t0
    print(f"  generated in {dt:.3f} s")

    if not graphs:
        print("  no graphs produced (constraint may be too strict)")
        return

    G = graphs[0]
    print(f"  nodes: {G.number_of_nodes()}, edges: {G.number_of_edges()}")

    # Save to .nodes / .edges so other demos can load it via topology.source = "load"
    nodes_path = OUTPUT_DIR / "network.nodes"
    edges_path = OUTPUT_DIR / "network.edges"
    with open(nodes_path, "w") as f:
        for n, attrs in G.nodes(data=True):
            f.write(f"{n} {attrs.get('x', 0)} {attrs.get('y', 0)} {attrs.get('z', 0)}\n")
    with open(edges_path, "w") as f:
        for u, v in G.edges():
            f.write(f"{u} {v}\n")
    print(f"  saved: {nodes_path}, {edges_path}")
    print()
    print("Use the generated .nodes/.edges files as the topology source for an")
    print("atomistic or coarse-grained chemistry build by setting")
    print(f"  topology.source = \"load\"")
    print(f"  topology.existing_files.nodes_file = \"{nodes_file}\"")
    print(f"  topology.existing_files.edges_file = \"{edges_file}\"")
    print("in any of the demos under examples/demos/polymer/.")


if __name__ == "__main__":
    main()
