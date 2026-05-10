# Showcase — small reference data files

Small reference files demonstrating topon's input/output formats. Useful for:
- Seeing what a generated topology looks like before running a demo end-to-end.
- Loading as a starting topology in any of the polymer demos via `topology.source = "load"`.
- Quickly inspecting expected file layouts.

## Contents

### `network_5x5x5/`

A 5×5×5 SC-lattice network produced by the topology generator on trial 3.
Tiny by topon standards — 125 nodes, 210 edges; the three files together are well under 10 KB.

| File | Purpose |
|---|---|
| `network.nodes` | Node ID + (x, y, z) coordinates per line — the topology generator's primary output. |
| `network.edges` | Edge list: pair of node IDs per line. |
| `generation.log` | Trial-3 edge-removal log (the sequence of `remove edge N-M` operations the sculpting algorithm took to satisfy the degree-distribution constraint). |

To use this network as the topology for any polymer demo:

```json
"topology": {
    "source": "load",
    "existing_files": {
        "nodes_file": "examples/showcase/network_5x5x5/network.nodes",
        "edges_file": "examples/showcase/network_5x5x5/network.edges"
    }
}
```

Then run:

```bash
topon generate <your-config>.json
```

(Subject to the schema gap noted in `internal/DEVELOPMENT_INTERNAL.md` §1 P0-A — for polymer demos the workflow-script path may be needed instead of the bare CLI.)

## Larger reference data

For paper-scale topologies (6×6×6 and beyond, with mechanics measurements), see [`examples/npjcompmat/data/mechanics/`](../npjcompmat/data/mechanics/) — ~900 networks indexed by their degree-distribution constraint string.
