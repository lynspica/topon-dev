# Expected output — polymer/atomistic/basic

Generated 2026-05-10 by re-running the demo end-to-end through `topon.pipeline.Pipeline`
after the P1-K fix (Option C: Sanitize → AddHs → Gasteiger + 5-displace tail).
Topology source: `tests/sample_graphs/network_N5x5x5_trial3` (5×5×5 SC, injected — the
demo config doesn't ship its own topology).

## Stage results

| Stage | Status | Wall time | Final notes |
|---|---|---|---|
| Pipeline.run() | ok | ~2 s | 10 949 atoms, 4 DREIDING types (Si3, C_3, O_3, H_), net charge 0.0000 e |
| LAMMPS stage 1 (soft min) | ok | 4 s | `system_after_soft.data`; TotEng 39 761 |
| LAMMPS stage 2 (LJ ramp) | ok | 5 min 05 s | `system_ramped.data`; final E_pair −163 033; Temp 668 K (cooled from 752) |
| LAMMPS stage 3 (CG min + NVT + NPT) | ok | 1 min 31 s | `system_equilibrated.data`; "All Minimization Stages Complete." |

Total Pipeline + LAMMPS wall time: ~6 m 40 s on one CPU core.

## Files in this folder

**Chemistry (`02_Chemistry/`-equivalent)**

| File | Notes |
|---|---|
| `system.data` | LAMMPS data, 10 949 atoms, atom_style full, charged |
| `system.in.settings` | DREIDING pair/bond/angle/dihedral/improper coeffs |
| `system.groups` | `nodes` (junction atoms) + `beads` (everything else) |
| `system_nodes.displace` | Junction atom target coords (graph node positions) |
| `system_backbone.displace` | Chain backbone Si atoms |
| `system_grafts.displace` | (empty — atomistic graft path doesn't populate `graft_atom_map` yet) |
| `system_pendant.displace` | C/O side-chain atoms (Gaussian propagation through `mol_h`) |
| `system_hydrogens.displace` | H atoms from `AddHs` |

**Conformation (`03_Conformation/`-equivalent)**

| File | Notes |
|---|---|
| `system_conformed.data` | After `apply_displacements` |
| `system_relaxed.data` | After overlap resolution |

**Simulation (`04_Simulation/`-equivalent)**

| File | Notes |
|---|---|
| `minimize_1_serial.in`, `minimize_2_parallel.in`, `minimize_3_parallel.in` | The three calibrated LAMMPS scripts |
| `system_after_soft.data` | After stage 1 |
| `system_ramped.data` | After stage 2 |
| `system_minimized_final.data` | After stage 3 CG min |
| `system_equilibrated.data` | After stage 3 NVT/NPT — production-ready |
| `log.stage1.lammps`, `log.stage2.lammps`, `log.stage3.lammps` | Per-stage LAMMPS logs |

## Reproducing

```bash
python - <<'PY'
from pathlib import Path
from topon.config import load_config_full
from topon.config.schema import TopologyConfig, ExistingFilesConfig
from topon.pipeline import Pipeline
cfg, raw = load_config_full(Path("examples/demos/polymer/atomistic/basic/config.json"))
cfg.topology = TopologyConfig(source="load", existing_files=ExistingFilesConfig(
    nodes_file="tests/sample_graphs/network_N5x5x5_trial3.nodes",
    edges_file="tests/sample_graphs/network_N5x5x5_trial3.edges"))
cfg.study.output_dir = "runs/basic_demo"; cfg.study.name = "run"
Pipeline(cfg, raw_config=raw).run()
PY

cd runs/basic_demo/run/04_Simulation
lmp -in minimize_1_serial.in   # ~4 s
lmp -in minimize_2_parallel.in # ~5 min
lmp -in minimize_3_parallel.in # ~1.5 min
```
