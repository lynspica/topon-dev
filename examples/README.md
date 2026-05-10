# topon — Examples

Curated demos showcasing what topon can build.

```bash
topon generate examples/demos/<path>/config.json
```

> **Heads-up — schema gap.** The current `topon generate` strict-validates the JSON config against `topon.config.schema.ToponConfig`. Some demo configs include `conformation`, `simulation`, or `execution` sections that are not yet covered by the Pydantic schema; those configs run through the workflow scripts under `tests/workflows/` (which read the raw JSON) but will fail `topon generate`'s strict validation pass. Tracked in [internal notes](../internal/DEVELOPMENT_INTERNAL.md). For the demos that don't use those sections (the topology and POSS demos here), `topon generate` works as-is.

For the full CLI reference and the JSON-config schema, see [`docs/USAGE.md`](../docs/USAGE.md).

## Layout

```
examples/
├── templates/              minimal.json + full.json starter configs (also produced by `topon init`)
├── defaults/               default node-type / edge-type assignment fragments
├── showcase/               small reference data files (input format examples)
├── run_via_api.py          generic Python-API runner (equivalent to `topon generate`)
├── demos/
│   ├── polymer/
│   │   ├── atomistic/      DREIDING — basic, entanglement, copolymer, graft, defect, combined
│   │   └── coarse_grained/ Kremer-Grest — basic, entanglement, copolymer, graft, defect, combined
│   ├── protein/
│   │   ├── charmm/         CHARMM atomistic protein networks (legacy topro path)
│   │   └── martini/        MARTINI 3 protein networks (current — `topon.protein_network`)
│   ├── topology/
│   │   ├── end_linking/    same network from two paths (C generator vs Python port) + speed logs
│   │   └── bfm/            BFM lattice topology used by topro
│   └── poss/               POSS junctions in atomistic networks
└── npjcompmat/             pre-generated dataset + figure notebook from the npj submission
```

## Per-category overviews

- [demos/polymer/](demos/polymer/README.md) — polymer networks, atomistic vs CG with each architecture knob.
- [demos/protein/](demos/protein/README.md) — protein networks (the user-facing name is **topro**), CHARMM legacy + MARTINI current.
- [demos/topology/](demos/topology/README.md) — topology generation only (no chemistry), C vs Python paths plus BFM.
- [demos/poss/](demos/poss/README.md) — POSS-junction node chemistry in atomistic networks.
- [showcase/](showcase/README.md) — small reference data files (input format examples; not generated demos).

## npjcompmat

`examples/npjcompmat/` is a copy of the dataset + figure pipeline that accompanied the npj Computational Materials submission. It is not a topon-generated demo — it ships pre-generated network topologies (`data/mechanics/`), the analysis dataset (`data/dataset.pkl` and `data/csv/`), and the figure-generation notebook (`generate_figures.ipynb`).

To regenerate the figures:

```bash
cd examples/npjcompmat/
jupyter notebook generate_figures.ipynb
```
