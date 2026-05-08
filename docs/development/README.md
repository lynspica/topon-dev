# Topon Development Documentation

This folder contains all development-related documentation for the Topon project.

## Contents

| Document | Description |
|----------|-------------|
| [changelog.md](changelog.md) | Version-by-version changes (v1–v32+) |
| [implementation_plan.md](implementation_plan.md) | Architecture, module reference, data flow |
| [walkthrough.md](walkthrough.md) | Step-by-step guide: CG, Atomistic, simbox |
| [tasks.md](tasks.md) | Task tracking (accomplished and planned) |

## Quick Reference

### Config Files
- `examples/config_cg.json` - CG workflow configuration
- `examples/config_atomistic.json` - Atomistic workflow configuration
- `examples/experimental.json` - Production LAMMPS parameters
- `examples/defaults/node_degree.json` - Node assignment rules
- `examples/defaults/edge_uniform.json` - Edge assignment rules

### Running Workflows
```bash
# CG baseline
python tests/workflows/generate_cg.py

# Atomistic baseline
python tests/workflows/generate_atomistic.py

# CG with entanglements
python tests/workflows/generate_cg_entangled.py --output output_dir

# Atomistic with grafts
python tests/workflows/generate_atomistic_combined.py --output output_dir

# Run ensemble (1000 systems, 8 workers)
python tests/workflows/run_cg_ensemble.py

# Run unit tests
pytest tests/unit/

# Verify both workflows end-to-end
python examples/verify_workflows.py
```

### Key Architecture Note
Both the `topon generate` CLI / `Pipeline` class **and** the `tests/workflows/generate_*.py`
scripts are fully functional. The CLI is the recommended entry point for new workflows.
See [tasks.md](tasks.md) for status.
