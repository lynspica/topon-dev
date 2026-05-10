# Expected output — atomistic-basic

Snapshot of running the sibling `config.json` (or its programmatic
equivalent) through `Pipeline.run()` plus LAMMPS stage-1 minimize
on the topon project's bundled `tests/sample_graphs/network_N5x5x5_trial3`.

## Files

- `system.data` — Stage-4 chemistry output (LAMMPS data file)
- `system.in.settings` — pair / bond / angle / dihedral coefficients
- `system.groups` — group definitions (`nodes`, `beads`)
- `minimize_1_serial.in` — Stage-6 LAMMPS input script
- `log.lammps` — what LAMMPS prints during the stage-1 run (energies, etc.)
- `system_after_soft.data` — system state after stage-1 minimize

Reproduce with `topon generate config.json` (when this demo
ships a JSON config compatible with `Pipeline.run()`) or via
`examples/run_via_api.py`. Generated 2026-05-10.
