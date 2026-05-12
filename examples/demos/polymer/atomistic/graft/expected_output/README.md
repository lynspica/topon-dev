# Expected output — polymer/atomistic/graft

Generated 2026-05-10 by re-running through `topon.pipeline.Pipeline`
after the P1-K chemistry fixes (over-valent Si + atomistic graft).
Topology: `tests/sample_graphs/network_N5x5x5_trial3` (5x5x5 SC).

## Stage results

| Stage | Status |
|---|---|
| Pipeline.run() | ok (41565 atoms, neutral net charge) |
| LAMMPS stage1 | ok in 56.3 s |
| LAMMPS stage2 | ok in 3807.7 s |
| LAMMPS stage3 | ok in 1982.7 s |

## Notes

- Fixed 2026-05-10: atomistic graft chemistry now uses a per-repeat PDMS builder (`ChemistryBuilder._build_pdms_chain_with_grafts`) that emits real side chains at backbone positions chosen by `graft_density`. Each graft is placed perpendicular to the backbone vector with effective length min(0.5 * edge_len, graft_dp/backbone_dp * edge_len, 0.5 * lattice_spacing) — the third term keeps graft tips inside their own lattice cell. Pre-fix graft demo built only the backbone (no side-chain atoms); chemistry now produces ~2x more atoms.

## Reproducing

```bash
python - <<'PY'
from pathlib import Path
from topon.config import load_config_full
from topon.config.schema import TopologyConfig, ExistingFilesConfig
from topon.pipeline import Pipeline
cfg, raw = load_config_full(Path("examples/demos/polymer/atomistic/graft/config.json"))
cfg.topology = TopologyConfig(source="load", existing_files=ExistingFilesConfig(
    nodes_file="tests/sample_graphs/network_N5x5x5_trial3.nodes",
    edges_file="tests/sample_graphs/network_N5x5x5_trial3.edges"))
cfg.study.output_dir = "runs/graft"; cfg.study.name = "run"
Pipeline(cfg, raw_config=raw).run()
PY
cd runs/graft/run/04_Simulation
lmp -in minimize_1_serial.in
lmp -in minimize_2_parallel.in
lmp -in minimize_3_parallel.in
```