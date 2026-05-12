# Topro — CHARMM36m atomistic protein networks

CHARMM36m all-atom protein networks: BFM lattice topology with stochastic dityrosine crosslinks, atomistic placement via the CHARMM RTF/PRM tables, TIP3P solvation + NaCl background, three-stage soft → LJ-ramp → tight relaxation.

Lives at [`topon.protein_network.charmm`](../../../../topon/protein_network/charmm/) — forked verbatim from the legacy `topro` package after the MARTINI port (`topon.protein_network`) inherited the BFM topology stage but not the atomistic chemistry stage.

The bundled `data/` ships CHARMM36m PRM/RTF/CMAP files (the `updated_charges` IDP-optimized variants) so the demo runs out of the box without external FF downloads.

## Quick run — small dry resilin

```bash
# 1. Generate a BFM topology (resilin GGRPSDSYGAPGGGN, 8 chains, 8 repeats)
python - <<'PY'
from pathlib import Path
from topon.protein_network.bfm import generate_topology
from topon.protein_network.charmm.topology_io import save_topology
out = Path("runs/charmm_test"); out.mkdir(parents=True, exist_ok=True)
topo = generate_topology(
    n_chains=8, n_repeats=8, segs_per_block=2,
    equil_steps=10000, n_extra_snapshots=2, snapshot_delta_conv=0.05,
    seed=42, verbose=False,
)
save_topology(topo, str(out / "topo.json"))
PY

# 2. Build atomistic LAMMPS systems for w0 + w35 water contents
python -m topon.protein_network.charmm.build_systems \
    --topology runs/charmm_test/topo.json \
    --snapshot 0 \
    --output runs/charmm_test/sys \
    --water_contents 0,35 \
    --n_repeats 8

# 3. Run the 3-stage relaxation in LAMMPS (per water content)
cd runs/charmm_test/sys/w0/relaxation
lmp -in protein_network_stage1.in   # ~5s   soft overlap removal
lmp -in protein_network_stage2.in   # ~3min epsilon ramp, nve/limit
lmp -in protein_network_stage3.in   # ~30s  CG min + NVT + NPT
```

Output layout per water content:

```
runs/charmm_test/sys/
  w0/
    protein_network.data          # LAMMPS atomistic data (~11k atoms dry)
    protein_network.in.settings   # pair/bond/angle/dihedral coeffs
    protein_network.in.groups     # protein, chain01..chainN, water, ions
    charmm36m.cmap                # backbone CMAP correction grid
    relaxation/
      protein_network_stage1.in   # soft-push overlap removal
      protein_network_stage2.in   # LJ epsilon ramp
      protein_network_stage3.in   # tight CG min -> NVT -> NPT
  w35/  ...                        # 35 wt% water + NaCl background
```

After stage 3, `system_equilibrated.data` lives in the parent `wXX/` folder, ready for production.

## CLI reference

Run `python -m topon.protein_network.charmm.build_systems --help` for the full flag list. Key flags:

| Flag | Default | Description |
|---|---|---|
| `--topology` | (required) | Path to BFM topology JSON |
| `--snapshot` | `gel_point` | Snapshot label or integer index |
| `--block_seq` | `GGRPSDSYGAPGGGN` | One-letter repeat block (resilin default) |
| `--water_contents` | `0,35,55,65,75` | wt% water values to build |
| `--salt_conc` | `0.15` | NaCl background concentration in mol/L |
| `--charmm_prm` | bundled CHARMM36m | Override PRM file |
| `--charmm_rtf` | bundled CHARMM36m | Override RTF file |
| `--charmm_cmap` | bundled CHARMM36m | Override CMAP file |
| `--target_density` | `0.85` g/cm³ | Initial box density for auto lattice scale |

## Differences from the MARTINI port

| | atomistic CHARMM | MARTINI 3 (sibling folder) |
|---|---|---|
| FF | CHARMM36m + CMAP, TIP3P | MARTINI 3 (vendored ITPs) |
| Atom count, 4 chains × 8 repeats | ~6,000 atoms (dry) | ~400 beads (dry) |
| Solvation | Per-water explicit | MARTINI W beads |
| Sequence coverage | All 20 amino acids (via RTF) | 8 residues (resilin set) |
| Implementation | `topon.protein_network.charmm` | `topon.protein_network` (top-level) |

The two share the BFM topology stage (`topon.protein_network.bfm`), so you can generate one topology JSON and feed it into either chemistry path.

## Known limitations

- The default `--snapshot` is `gel_point`; if your topology generation didn't reach gel (small system + few `equil_steps`), pass `--snapshot 0` (first snapshot, post-gel or otherwise).
- Stage 2 is the long pole — for the 8×8 demo it's ~3 minutes serial, scaling roughly linearly with atom count. Run with MPI (`mpirun -np 8 lmp -in protein_network_stage2.in`) for production.
- The legacy topro `scripts/build_systems.py` had a Windows-absolute fallback path; the integrated CLI uses bundled `data/` PRM/RTF/CMAP and resolves them per-platform.

## Where things live

```
topon/protein_network/
  bfm.py               # BFM topology generator (shared with MARTINI)
  sequence.py          # Block-seq -> BFM-node residue mapping
  charmm/              # this demo's chemistry stage
    charmm_ff.py       # RTF/PRM parser
    builder.py         # Atomistic placement, water + ion packer
    lammps_writer.py   # LAMMPS data/settings/groups + 3 stage scripts
    topology_io.py     # Topology JSON load/save
    build_systems.py   # CLI entry point (this demo)
    data/              # CHARMM36m PRM/RTF/CMAP (bundled)
```
