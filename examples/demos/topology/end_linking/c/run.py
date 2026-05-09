"""End-linking topology generation via the C generator.

Generates the same 6x6x6 SC-lattice polymer-network topology as the Python
demo next door, but invokes the compiled C executable
(``generator.exe`` / ``generator``) as a subprocess. Faster than the
Python port for hard cases.

Set ``GENERATOR_EXE`` below (or the ``TOPON_GENERATOR_EXE`` env var) to the
absolute path of your compiled binary before running.

Usage:
    python run.py
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

from topon.config.schema import GeneratorConfig
from topon.topology.generator import run_generator
from topon.topology.loader import load_graph

# Configure the path to your compiled C binary here, or via env var.
GENERATOR_EXE = os.environ.get("TOPON_GENERATOR_EXE", "")

OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


def main() -> None:
    if not GENERATOR_EXE or not Path(GENERATOR_EXE).exists():
        print("ERROR: C generator binary not found.")
        print(f"  GENERATOR_EXE = {GENERATOR_EXE!r}")
        print("Set TOPON_GENERATOR_EXE to the absolute path of generator.exe,")
        print("or edit GENERATOR_EXE at the top of this script. The Python")
        print("generator demo next door is a good alternative for most cases")
        print("(see ../python/run.py).")
        sys.exit(1)

    print(f"--- End-linking topology generation (C generator: {GENERATOR_EXE}) ---")
    cfg = GeneratorConfig(
        lattice_size="6x6x6",
        lattice_type="SC",
        periodicity="111",
        max_functionality=4,
        degree_distribution="0:13,1:25",
        max_trials=1_000_000,
        max_saves=1,
        exe_path=GENERATOR_EXE,
    )
    t0 = time.perf_counter()
    nodes_path, edges_path = run_generator(cfg, OUTPUT_DIR, exe_path=GENERATOR_EXE)
    dt = time.perf_counter() - t0
    print(f"  generated in {dt:.3f} s")

    G, dims = load_graph(nodes_path=str(nodes_path), edges_path=str(edges_path))
    print(f"  nodes: {G.number_of_nodes()}, edges: {G.number_of_edges()}, box: {dims}")
    print(f"  saved: {nodes_path} + {edges_path}")


if __name__ == "__main__":
    main()
