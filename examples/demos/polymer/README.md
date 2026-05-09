# Polymer-network demos

Each demo below ships a `config.json`. Run with:

```bash
topon generate examples/demos/polymer/<atomistic|coarse_grained>/<demo>/config.json
```

The two resolutions share the same architecture knobs; pick whichever your downstream LAMMPS workflow expects.

## Atomistic (DREIDING)

| Demo | What it shows |
|---|---|
| [`atomistic/basic/`](atomistic/basic/) | Minimal atomistic network — uniform edge/node types, no architecture knobs. |
| [`atomistic/entanglement/`](atomistic/entanglement/) | Five Gaussian-Kink entanglements injected into the network. |
| [`atomistic/copolymer/`](atomistic/copolymer/) | Block / random / alternating monomer sequences along chains. |
| [`atomistic/graft/`](atomistic/graft/) | Side-chain attachment with dynamic DP scaling (V20). |
| [`atomistic/defect/`](atomistic/defect/) | Primary-loop defects (parallel edges) with valence protection. |
| [`atomistic/combined/`](atomistic/combined/) | Entanglements + grafts together (V21 baseline). |

## Coarse-grained (Kremer-Grest)

The same six demos, with `chemistry.model_type = "coarse_grained"` and KG-specific simulation knobs (`include_angles`, `pair_style`).

| Demo | What it shows |
|---|---|
| [`coarse_grained/basic/`](coarse_grained/basic/) | Bead-spring network with attractive LJ (`pair_style: attractive`). |
| [`coarse_grained/entanglement/`](coarse_grained/entanglement/) | Five Gaussian-Kink entanglements at CG scale. |
| [`coarse_grained/copolymer/`](coarse_grained/copolymer/) | Per-bead monomer sequence assignment. |
| [`coarse_grained/graft/`](coarse_grained/graft/) | Side-chain CG bead attachment. |
| [`coarse_grained/defect/`](coarse_grained/defect/) | Primary loops in CG. |
| [`coarse_grained/combined/`](coarse_grained/combined/) | Entanglements + grafts in CG. |

## Architecture-knob reference

| Knob | Config section | Notes |
|---|---|---|
| Entanglements | `assignment.entanglements` | Gaussian-Kink geometry, count or distribution mode |
| Copolymers | `assignment.copolymer.per_edge_type` | block / random / alternating / gradient |
| Grafts | `assignment.grafts.per_edge_type` | per-edge graft density + side-chain DP/monomer |
| Defects | `assignment.defects.primary_loops` | parallel-edge injection with valence cap |

For the full schema, see [`docs/USAGE.md`](../../../docs/USAGE.md) Appendix A.

## Topology source

These demos generate their topology in-place via the Python generator. To use a pre-generated topology instead, swap `topology.source` to `"load"` and point at any `.nodes`+`.edges` pair (for example, the network produced by [`../topology/end_linking/python/run.py`](../topology/end_linking/python/run.py), or one of the npj-paper topologies under `examples/npjcompmat/data/mechanics/`).
