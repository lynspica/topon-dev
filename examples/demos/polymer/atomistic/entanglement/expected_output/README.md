# Expected output — polymer/atomistic/entanglement

Generated 2026-05-10 by re-running through `topon.pipeline.Pipeline`
after the P1-K fix (Option C: Sanitize -> AddHs -> Gasteiger + 5-displace tail).
Topology: `tests/sample_graphs/network_N5x5x5_trial3` (5x5x5 SC, injected).

## Stage results

| Stage | Status |
|---|---|
| Pipeline.run() | ok (21449 atoms, 5 displace files where applicable) |
| LAMMPS stage1 | ok in 35 s |
| LAMMPS stage2 | ok in 26 min 21 s |
| LAMMPS stage3 | ok in 6 min 05 s |

## Notes

- Full end-to-end success: `system_equilibrated.data` and `after_nvt_real.data` both present.

## File inventory

Chemistry stage: `system.data`, `system.in.settings`, `system.groups`,
and the displacement files `system_{nodes,backbone,grafts,pendant,hydrogens}.displace`
(empty `*_grafts.displace` is expected for atomistic — see internal P1-K notes).

Conformation stage: `system_conformed.data` (after `apply_displacements`),
`system_relaxed.data` (after overlap resolution).

Simulation stage: the three calibrated LAMMPS scripts, plus whichever
of `system_after_soft.data` (stage 1), `system_ramped.data` (stage 2),
`system_minimized_final.data` + `system_equilibrated.data` (stage 3),
and `log.stage{1,2,3}.lammps` actually completed in this run.

## Reproducing

```bash
python - <<'PY'
from pathlib import Path
from topon.config import load_config_full
from topon.config.schema import TopologyConfig, ExistingFilesConfig
from topon.pipeline import Pipeline
cfg, raw = load_config_full(Path("examples/demos/polymer/atomistic/entanglement/config.json"))
cfg.topology = TopologyConfig(source="load", existing_files=ExistingFilesConfig(
    nodes_file="tests/sample_graphs/network_N5x5x5_trial3.nodes",
    edges_file="tests/sample_graphs/network_N5x5x5_trial3.edges"))
cfg.study.output_dir = "runs/entanglement"; cfg.study.name = "run"
Pipeline(cfg, raw_config=raw).run()
PY
cd runs/entanglement/run/04_Simulation
lmp -in minimize_1_serial.in
lmp -in minimize_2_parallel.in
lmp -in minimize_3_parallel.in
```