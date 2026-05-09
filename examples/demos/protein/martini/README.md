# Topro — MARTINI 3 protein networks

The current protein-network workflow — coarse-grained MARTINI 3, BFM lattice topology with stochastic dityrosine crosslinks. Implemented under `topon.protein_network`; user-facing name **topro**.

## Quick run — dry resilin reference

```bash
python -m topon.protein_network generate \
    --block-seq GGRPSDSYGAPGGGN \
    --n-repeats 6 --n-chains 4 \
    --equil-steps 5000 \
    --water-density 0 \
    --output runs/resilin_dry --seed 42
```

This generates a 4-chain × 6-repeat resilin network, no water, dry. Output:

```
runs/resilin_dry/
  protein_network_topology.json   # BFM snapshots (gel_point, post_gel_1, ...)
  protein_network.data            # LAMMPS data file
  protein_network.in.settings     # pair_coeff / bond_coeff / etc.
  protein_network.in.groups       # protein vs water groups
  relaxation/
    protein_network_stage1.in     # soft-push overlap removal
    protein_network_stage2.in     # LJ-epsilon ramp via nve/limit
    protein_network_stage3.in     # tight CG min + brief NVT/NPT @ 310 K
```

Run the relaxation stages in LAMMPS (serial or MPI), in order. After stage 3, `system_equilibrated.data` sits in the parent directory ready for production.

## Variants

```bash
# Hydrated
python -m topon.protein_network generate \
    --block-seq GGRPSDSYGAPGGGN --n-repeats 6 --n-chains 16 \
    --water-density 4 \
    --output runs/resilin_w4 --seed 42

# Sweep across water densities
python -m topon.protein_network sweep \
    --block-seq GGRPSDSYGAPGGGN --n-repeats 6 --n-chains 16 \
    --water-densities 0,1,4 \
    --output runs/resilin_sweep --seed 42

# Topology only (BFM snapshot JSON, no chemistry yet)
python -m topon.protein_network topology \
    --n-chains 16 --n-repeats 6 \
    --output topo.json
```

## Sequence coverage

The vendored MARTINI 3 protein FF currently covers 8 amino acids (the resilin reference set): `GLY`, `ALA`, `ARG`, `PRO`, `SER`, `ASP`, `TYR`, `ASN`. To extend to other residues, drop a polyply-generated ITP into `tests/_martini_extracted/Martini_Ahmet/itp_files/` and re-run `tools/extract_residues_from_itp.py` — the residue table regenerates without code edits.

`TYR` positions are stochastically linked into dityrosine crosslinks (SC4–SC4 → TN6 bead) by the BFM gel-point stage.

## Full reference

For all CLI flags (water beads, ion packing, hierarchical stage 1, etc.), the Python API, the known approximations, and the full output description, see [`docs/USAGE.md`](../../../../docs/USAGE.md) §4.1.
