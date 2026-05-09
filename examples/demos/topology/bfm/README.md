# BFM topology (topro)

The BFM (Bond-Fluctuation Model) topology generator that drives topro's protein-network workflow. Self-avoiding-walk chain placement on a cubic lattice with end / kink-crankshaft / reptation MC moves and stochastic Y-Y crosslinking detected via Union-Find gel-point.

This demo emits *just* the topology JSON snapshot — no chemistry, no LAMMPS files. Useful for inspecting the gel-point structure, sweeping topology parameters before committing to a full MARTINI run, or feeding a custom snapshot into a different chemistry builder.

## Run

```bash
python -m topon.protein_network topology \
    --n-chains 16 --n-repeats 6 \
    --segs-per-block 2 \
    --target-packing 0.45 \
    --equil-steps 5000 \
    --seed 42 \
    --output bfm_topology.json
```

Output is a single JSON file containing snapshots at multiple convergence points (`gel_point`, `post_gel_1`, …). The format is byte-compatible with the legacy topro `topo_*.json`.

## Inspecting

```python
from topon.protein_network.topology_io import load_topology, get_snapshot

topo = load_topology("bfm_topology.json")
snap = get_snapshot(topo, "gel_point")
print(snap.summary())
```

## Feeding into a chemistry build

Pass the snapshot to `build_protein_system` to attach MARTINI 3 chemistry:

```python
from topon.protein_network import bfm, builder, topology_io
from topon.protein_network.martini_ff import MartiniLibrary

topo = topology_io.load_topology("bfm_topology.json")
snap = topology_io.get_snapshot(topo, "gel_point")
sys_ = builder.build_protein_system(
    snap,
    sequence_3letter="GLY-GLY-ARG-PRO-SER-ASP-SER-TYR-GLY-ALA-PRO-GLY-GLY-GLY-ASN" * 6,
    library=MartiniLibrary.from_package_data(),
)
```

This is the lower-level entry point that `python -m topon.protein_network generate` wires up automatically. Use it when you want to swap the chemistry stage (e.g. for a future CHARMM build — see the [`charmm/` README](../../protein/charmm/README.md)).

## Reference

Full BFM and topology-IO API: [`docs/USAGE.md`](../../../../docs/USAGE.md) §4.1.
