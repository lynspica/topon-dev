# Expected output — examples/demos/poss

Generated 2026-05-10 by `_tmp_run_all_demos.py` (one-off runner).
Topology source: `tests/sample_graphs/network_N5x5x5_trial3` (injected — the demo config doesn't ship its own topology).

## Stage status

- **load_config_full**: ok
- **topology**: as-configured (generate)
- **Pipeline.run**: ok
- **stage1**: ok in 10.0s
- **stage2**: exit 1, 4.4s; ERROR: Must use kspace_modify gewald for uncharged system (src/KSPACE/pppm.cpp:992)

## Files in this folder

Whatever the pipeline + LAMMPS produced before any failure:
`system.data`, `system.in.settings`, `system.groups`, minimize scripts, `log.*.lammps`, intermediate `system_*.data`.