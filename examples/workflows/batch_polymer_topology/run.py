"""Generate a batch of 25 polymer-network topologies and export each
in three formats (.nodes/.edges, GraphML, NPZ), plus a single CSV of
per-graph properties.

This is the kind of dataset you'd feed into a graph-neural-network
training pipeline or a structure-property survey — "give me N graphs
with the same lattice + degree distribution, with their summary stats
in one place."

The generator is deterministic per seed; SEED_BASE + i gives 25
independent realisations. Tune the dataset by editing the knobs at the
top — everything below them is mechanical.

Outputs (under runs/batch_polymer_topology/):
    summary.csv                  one row per graph, columns below
    graph_000/
        network.nodes            x/y/z coordinates, one line per node
        network.edges            u v pairs, one line per edge
        network.graphml          dual-graph form for GNN libraries
        network.npz              dual-graph dense arrays for PyTorch
    graph_001/ ...

CSV columns: seed, n_nodes, n_edges, avg_degree, min_degree, max_degree,
n_chain_ends (degree-1 nodes), n_interior (degree>=2 nodes).

Usage:
    python examples/workflows/batch_polymer_topology/run.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import networkx as nx
import numpy as np

from topon.config.schema import GeneratorConfig
from topon.topology.generator_python import PythonTopologyGenerator
from topon.writers.graphml_writer import write_graphml
from topon.writers.npz_writer import write_npz


# ----- knobs you'd typically change ----------------------------------------
N_GRAPHS = 25
OUTPUT_ROOT = Path(__file__).parent / "output"
LATTICE_SIZE = "5x5x5"
LATTICE_TYPE = "SC"
PERIODICITY = "111"
MAX_FUNCTIONALITY = 4
# "0:13,1:25" = 13% degree-0 nodes (vacancies), 25% degree-1 (chain caps),
# remainder distributed across degree-2..max_functionality. See
# topon.topology.generator_python::PythonTopologyGenerator for the format.
DEGREE_DISTRIBUTION = "0:13,1:25"
SEED_BASE = 42
DP_DEFAULT = 50         # nominal degree of polymerization, written into the NPZ/GraphML
MAX_TRIALS = 1_000_000  # generator retries per graph

# ----- helpers --------------------------------------------------------------

def write_nodes_edges(G: nx.Graph, nodes_path: Path, edges_path: Path) -> None:
    """Persist the raw lattice form so other topon configs can `load` this graph."""
    with nodes_path.open("w") as f:
        for n, attrs in G.nodes(data=True):
            f.write(f"{n} {attrs.get('x', 0)} {attrs.get('y', 0)} {attrs.get('z', 0)}\n")
    with edges_path.open("w") as f:
        for u, v in G.edges():
            f.write(f"{u} {v}\n")


def stats(G: nx.MultiGraph) -> dict:
    degs = [d for _, d in G.degree()]
    return {
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
        "avg_degree": round(float(np.mean(degs)), 4),
        "min_degree": int(min(degs)),
        "max_degree": int(max(degs)),
        "n_chain_ends": sum(1 for d in degs if d == 1),
        "n_interior": sum(1 for d in degs if d >= 2),
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []

    print(f"--- batch_polymer_topology: {N_GRAPHS} graphs, {LATTICE_TYPE} {LATTICE_SIZE} ---")
    for i in range(N_GRAPHS):
        seed = SEED_BASE + i
        cfg = GeneratorConfig(
            lattice_size=LATTICE_SIZE,
            lattice_type=LATTICE_TYPE,
            periodicity=PERIODICITY,
            max_functionality=MAX_FUNCTIONALITY,
            degree_distribution=DEGREE_DISTRIBUTION,
            max_trials=MAX_TRIALS,
            max_saves=1,
            seed=seed,
        )
        gen = PythonTopologyGenerator(cfg)
        graphs = gen.generate(trials=20)
        if not graphs:
            print(f"  graph_{i:03d}: no graph produced (constraint too strict)")
            continue
        G = graphs[0]
        if not isinstance(G, nx.MultiGraph):
            G = nx.MultiGraph(G)

        sub = OUTPUT_ROOT / f"graph_{i:03d}"
        sub.mkdir(exist_ok=True)
        write_nodes_edges(G, sub / "network.nodes", sub / "network.edges")
        write_graphml(G, str(sub / "network.graphml"), dp=DP_DEFAULT)
        write_npz(G, str(sub / "network.npz"), dp=DP_DEFAULT)

        row = {"graph_id": f"graph_{i:03d}", "seed": seed, **stats(G)}
        rows.append(row)
        print(f"  graph_{i:03d}: {row['n_nodes']} nodes, {row['n_edges']} edges, "
              f"avg_deg={row['avg_degree']}")

    if not rows:
        print("No graphs were produced. Check generator config.")
        return

    csv_path = OUTPUT_ROOT / "summary.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {len(rows)} graphs and {csv_path.relative_to(Path.cwd())}")


if __name__ == "__main__":
    main()
