# `examples/workflows/` — scriptable recipes for topon

These are **standalone Python scripts** — not config-driven demos — that show how to drive topon programmatically for batch generation, parameter sweeps, and pipeline-chaining use cases. Each script is meant to be read top-to-bottom; you'll typically copy it, change the knob block at the top, and run.

| Workflow | What it does | Output |
|---|---|---|
| [`batch_polymer_topology/run.py`](batch_polymer_topology/run.py) | Generates 25 lattice-network graphs (seed-controlled), exports each as `.nodes/.edges` + GraphML + NPZ, and writes a single per-graph property CSV. | `output/graph_NNN/{network.nodes,network.edges,network.graphml,network.npz}` + `output/summary.csv` |
| [`bfm_gel_point_sweep/run.py`](bfm_gel_point_sweep/run.py) | Sweeps BFM parameters (chains, repeats, segs-per-block, intra-chain separation, equilibration) and records the gel-point conversion for each. | `output/{label}.json` per config + `output/summary.csv` |
| [`bfm_to_martini/run.py`](bfm_to_martini/run.py) | Builds a MARTINI 3 coarse-grained LAMMPS system from a BFM topology (resilin sequence by default). | `output/protein_network.{data,in.settings,in.groups}` + `output/relaxation/protein_network_stage{1,2,3}.in` |
| [`bfm_to_charmm/run.py`](bfm_to_charmm/run.py) | Builds a CHARMM36m atomistic LAMMPS system at several water contents from a BFM topology. Uses bundled FF data. | `output/sys/wXX/protein_network.{data,in.settings,in.groups}` + `relaxation/...` per water content |

## How they fit together

```
       PythonTopologyGenerator                       generate_topology
              |                                            |
              v                                            v
   batch_polymer_topology                       bfm_gel_point_sweep
        (lattice networks)                         (BFM topologies + gel)
              |                                            |
              |                                            +--> bfm_to_martini  (MARTINI 3 CG)
              |                                            |
              |                                            +--> bfm_to_charmm   (CHARMM36m AA)
              v
  GraphML / NPZ datasets for downstream
  graph-NN training or structure-property surveys
```

## Common knob-blocks

Every script has a "knobs you'd typically change" comment block right after imports. Edit those values, run the script, read the output. The intent is that each script is short enough to read in one sitting and obvious enough that you can fork it for your own variant.

## Stand-alone vs. demos

The scripts in `examples/demos/` are config-driven reference cases — one JSON per knob combination — and run via `topon generate <config.json>` or `python examples/run_via_api.py <config.json>`. These workflows are for cases where the demo config form is too constraining: you want a loop, a sweep, a CSV, or to chain stages together programmatically.

For the demo configs that *are* directly runnable, see [`examples/demos/`](../demos/) and the top-level [`examples/README.md`](../README.md).
