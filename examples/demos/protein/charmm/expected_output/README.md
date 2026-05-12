# Expected output — CHARMM atomistic protein-network demo

Generated 2026-05-10 by `examples/demos/protein/charmm/run.py` against the
config in [`../config.json`](../config.json). Reproducing this directory
exactly requires the bundled CHARMM36m PRM/RTF/CMAP files at
`topon/protein_network/charmm/data/` and a CHARMM-compatible LAMMPS build
(version `2 Apr 2025` here).

## Files in this folder

| File | Stage | Notes |
|---|---|---|
| `topo.json` | BFM | 8 chains × 8 repeats; seed=42, 10 000 MC steps; converges to conv=0.156 (no gel — small system). |
| `protein_network.data` | atomistic build | LAMMPS data file for the **dry** (w=0) snapshot. ~11 086 atoms. |
| `protein_network.in.settings` | atomistic build | All bonded + non-bonded coeffs from CHARMM36m. |
| `protein_network.in.groups` | atomistic build | `protein`, `chain01..chain08`, water, ions groups. |
| `protein_network_stage1.in` | LAMMPS | Soft overlap removal. Runs in ~5 s on a single core. |
| `protein_network_stage2.in` | LAMMPS | LJ epsilon ramp. ~3 min serial for this size. |
| `protein_network_stage3.in` | LAMMPS | Tight CG min + NVT + NPT. ~30 s. |
| `log.stage1.lammps` | LAMMPS | Reference stage-1 log; final `Total wall time: 0:00:05`. |

## Reproducing

```bash
python examples/demos/protein/charmm/run.py --output runs/charmm_demo
cd runs/charmm_demo/sys/w0/relaxation
lmp -in protein_network_stage1.in
lmp -in protein_network_stage2.in
lmp -in protein_network_stage3.in
```

Stage 2 and stage 3 outputs aren't included here — they take
non-trivial wall time and the user is expected to run them themselves.

## What changed vs. the legacy `topro` package

This demo runs through `topon.protein_network.charmm.build_systems`,
which is a thin re-export of the legacy `topro/scripts/build_systems.py`
with two adjustments:

1. PRM/RTF/CMAP files default to `topon/protein_network/charmm/data/`
   instead of an absolute Windows path under `tests/protein_network/`.
2. The `group protein all` output for dry systems was switched to
   `group protein union all` (the original was invalid LAMMPS syntax —
   the legacy code never exercised the dry path through the integrated
   tree's lj/charmm/coul/long stage 1 script).
