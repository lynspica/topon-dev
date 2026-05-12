# Expected output — polymer/coarse_grained/graft

Generated 2026-05-11 by re-running through `topon.pipeline.Pipeline`
after the CG-graft conformation fix (P1-J FIXED).
Topology: `tests/sample_graphs/network_N5x5x5_trial3` (5x5x5 SC).

## Stage results

| Stage | Status |
|---|---|
| Pipeline.run() | ok (8589 atoms, 2 atom types) |
| LAMMPS stage1 | ok in 2.2 s |
| LAMMPS stage2 | ok in 3 min 03 s |
| LAMMPS stage3 | ok in 0 min 57 s |

## Notes

- Pre-fix: Pipeline's CG branch never consulted `_builder.graft_atom_map`,
  so graft beads landed at (0, 0, 0) and LAMMPS stage 1 errored with
  `Neighbor list overflow` (P1-J). The chemistry stage was building the
  side chains correctly; only the placement code was wrong.
- 2026-05-11 fix: unified the CG branch of `_run_chemistry_stage` with
  the atomistic branch. CG now emits `system_backbone.displace` +
  `system_grafts.displace` (replacing the legacy combined
  `system_beads.displace`), with the same entanglement-aware kink
  placement and 3-way graft-length cap as atomistic. Stage 1 now passes
  in 2.2 s; full equilibration runs in under 4 minutes.

## Reproducing

```bash
python - <<'PY'
from pathlib import Path
from topon.config import load_config_full
from topon.config.schema import TopologyConfig, ExistingFilesConfig
from topon.pipeline import Pipeline
cfg, raw = load_config_full(Path("examples/demos/polymer/coarse_grained/graft/config.json"))
cfg.topology = TopologyConfig(source="load", existing_files=ExistingFilesConfig(
    nodes_file="tests/sample_graphs/network_N5x5x5_trial3.nodes",
    edges_file="tests/sample_graphs/network_N5x5x5_trial3.edges"))
cfg.study.output_dir = "runs/cg_graft"; cfg.study.name = "run"
Pipeline(cfg, raw_config=raw).run()
PY
cd runs/cg_graft/run/04_Simulation
lmp -in minimize_1_serial.in
lmp -in minimize_2_parallel.in
lmp -in minimize_3_parallel.in
```
