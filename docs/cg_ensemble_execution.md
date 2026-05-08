# Generating Your Coarse-Grained Ensemble

You now have a fully operational 1000-graph Kremer-Grest CG generation suite powered by a multiprocessing Python orchestration script! 

This orchestrator parses the highly detailed topologies previously saved into GraphML, translates them backwards into a pure repulsive `0.01 Å` minimized atomistic system without destroying the tagged entanglements, and drops out all resulting `.data` structures and LAMMPS parallel `.in` simulations automatically!

## How To Run It (Topon root directoy)

All 1000 source topological networks were generated natively from the `generator.exe` inside your environment.
Run the orchestration script directly:

```powershell
# Default (Utilizes 8 CPU workers concurrently to maximize speed)
python tests/workflows/run_cg_ensemble.py
```

### Advanced Usage

You can change the target input locations (where the source GraphML histories are loaded from) or modify the CPU count if you wish to run more/less operations at exactly the same time:

```powershell
python tests/workflows/run_cg_ensemble.py --workers 4 --input_dir path/to/study --output_dir target/output/
```

## Where Do Files Go?
All 1000 generated datasets will perfectly mirror `system_0001`!
These subdirectories will be produced per `system_XXXX` mapping out:
- `02_Chemistry`: Contains the un-conformed `.data` file with 80k beads and `system.in.settings` / `system.groups` files parameterizing the Kremer-Grest forces natively.
- `03_Conformation`: Contains `system_conformed.data` and `system_relaxed.data` mapped precisely backwards against the constraints with `verify_bonds.py` mathematics.
- `04_Simulation`: Contains the output LAMMPS `.in` files `minimize_1_serial.in` and `minimize_3_parallel.in` scripts fully initialized and ready to simulate!
