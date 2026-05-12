# Expected output — polymer/atomistic/combined

Generated 2026-05-10 by re-running through `topon.pipeline.Pipeline`
after the P1-K fix (Option C: Sanitize -> AddHs -> Gasteiger + 5-displace tail).
Topology: `tests/sample_graphs/network_N5x5x5_trial3` (5x5x5 SC, injected).

## Stage results

| Stage | Status |
|---|---|
| Pipeline.run() | ok (42449 atoms, 5 displace files) |
| LAMMPS stage1 | ok in 78.3 s |
| LAMMPS stage2 | ok in 3937.4 s |
| LAMMPS stage3 | partial: CG minimize completed (`system_minimized_final.data` present); NVT/NPT equilibration cut off at 1800 s budget |

## Notes

- Largest atomistic demo (42 k atoms = entanglements + grafts at DP=20). Wall time scales super-linearly under serial PPPM compared to the 21 k single-feature demos.
- Stage 3 minimize ran to completion (production-ready coordinates in `system_minimized_final.data`); only the post-minimize NVT/NPT short equilibration was cut by the 30 min budget. Re-run `lmp -in minimize_3_parallel.in` with a ~60 min budget to also get `system_equilibrated.data`.

## Reproducing

```bash
python - <<'PY'
from pathlib import Path
from topon.config import load_config_full
from topon.config.schema import TopologyConfig, ExistingFilesConfig
from topon.pipeline import Pipeline
cfg, raw = load_config_full(Path("examples/demos/polymer/atomistic/combined/config.json"))
cfg.topology = TopologyConfig(source="load", existing_files=ExistingFilesConfig(
    nodes_file="tests/sample_graphs/network_N5x5x5_trial3.nodes",
    edges_file="tests/sample_graphs/network_N5x5x5_trial3.edges"))
cfg.study.output_dir = "runs/combined"; cfg.study.name = "run"
Pipeline(cfg, raw_config=raw).run()
PY
cd runs/combined/run/04_Simulation
lmp -in minimize_1_serial.in
lmp -in minimize_2_parallel.in
lmp -in minimize_3_parallel.in
```