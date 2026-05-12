# Expected output — polymer/atomistic/defect

Generated 2026-05-10 by re-running through `topon.pipeline.Pipeline`
after the P1-K chemistry fixes (over-valent Si + atomistic graft).
Topology: `tests/sample_graphs/network_N5x5x5_trial3` (5x5x5 SC).

## Stage results

| Stage | Status |
|---|---|
| Pipeline.run() | ok (21944 atoms, neutral net charge) |
| LAMMPS stage1 | ok in 31.2 s |
| LAMMPS stage2 | ok in 1994.1 s |
| LAMMPS stage3 | ok in 790.0 s |

## Notes

- Fixed 2026-05-10 (P1-K): primary-loop injection now respects per-node-type valence (max_degree=4) AND skips end-cap nodes, so the chemistry build never produces over-valent Si. AddHs + Gasteiger return clean neutral charges naturally — no NaN scrub, no charge neutralisation needed.

## Reproducing

```bash
python - <<'PY'
from pathlib import Path
from topon.config import load_config_full
from topon.config.schema import TopologyConfig, ExistingFilesConfig
from topon.pipeline import Pipeline
cfg, raw = load_config_full(Path("examples/demos/polymer/atomistic/defect/config.json"))
cfg.topology = TopologyConfig(source="load", existing_files=ExistingFilesConfig(
    nodes_file="tests/sample_graphs/network_N5x5x5_trial3.nodes",
    edges_file="tests/sample_graphs/network_N5x5x5_trial3.edges"))
cfg.study.output_dir = "runs/defect"; cfg.study.name = "run"
Pipeline(cfg, raw_config=raw).run()
PY
cd runs/defect/run/04_Simulation
lmp -in minimize_1_serial.in
lmp -in minimize_2_parallel.in
lmp -in minimize_3_parallel.in
```