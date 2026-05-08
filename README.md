# Topon: Topological Polymer Network Generator

**Topon** is a modular Python package for generating complex polymer network structures for Molecular Dynamics (MD) simulations. It bridges the gap between graph-theoretical topology and physical chemical structures, supporting both Coarse-Grained (KG) and Atomistic (All-Atom) models.

## Key Features

- **Graph-Based Topology**: Define network connectivity first, then map to chemistry.
- **Dual Resolution**:
  - **Coarse-Grained (CG)**: Kremer-Grest model (FENE/Harmonic bonds, LJ potentials).
  - **Atomistic**: Full chemistry (e.g., PDMS, Silica, POSS) with DREIDING force field.
- **Advanced Architecture**:
  - **Entanglements**: Physical knots (Gaussian Kinks) preserving topology; single and multi-entanglement supported.
  - **Copolymers**: Block, Random, and Alternating sequences.
  - **Grafts**: Side-chain functionalization with dynamic scaling.
  - **Defects**: Configurable injection of primary loops (parallel edges).
  - **POSS Junctions**: Polyhedral Oligomeric Silsesquioxane node chemistry.
  - **Custom Topology**: Load any graph from `.nodes`/`.edges` or `.graphml` files.
- **Workflow Automation**: Integrated `SimulationRunner` for LAMMPS execution.
- **Ensemble Generation**: Multiprocessing orchestrator for batch production runs.
- **Molecule Packing** (`simbox`): General-purpose box packer for arbitrary molecular systems (Epoxy-PDMS, Amino-PDMS, POSS) with crosslinking simulation templates.

## Installation

```bash
git clone https://github.com/lynspica/topon.git
cd topon
pip install -e .
```

## Quick Start

### 1. Configuration
Topon uses JSON configuration files. See `examples/` for templates.

```json
{
    "chemistry": {
        "model_type": "coarse_grained",
        "degree_of_polymerization": 20
    },
    "assignment": {
        "entanglements": {"enabled": true, "target": 5},
        "grafts": {"enabled": true, "per_edge_type": {"A": {"graft_density": 0.05}}}
    }
}
```

### 2. Running a Workflow
Use the pre-built workflows in `tests/workflows/` to generate systems.

**Coarse-Grained Network with Entanglements**
```bash
python tests/workflows/generate_cg_entangled.py --output output_dir
```

**Atomistic Network with Grafts**
```bash
python tests/workflows/generate_atomistic_combined.py --output output_dir
```

**CG with POSS Junctions**
```bash
python tests/workflows/generate_atomistic_poss.py --output output_dir
```

**Ensemble (1000 systems, parallel)**
```bash
python tests/workflows/run_cg_ensemble.py --workers 8
```

**Custom Topology**
```bash
python tests/workflows/generate_v27_custom_topology.py --output output_dir
```

## Documentation

| Document | Description |
|---|---|
| [docs/cli.md](docs/cli.md) | **CLI reference** — all `topon` commands and options |
| [docs/config_reference.md](docs/config_reference.md) | **Config reference** — every JSON key with examples |
| [docs/development/walkthrough.md](docs/development/walkthrough.md) | Step-by-step usage guide (CG, Atomistic, simbox) |
| [docs/development/implementation_plan.md](docs/development/implementation_plan.md) | Module architecture and data flow |
| [docs/simbox.md](docs/simbox.md) | `simbox` sub-system: molecule packing for crosslinking |
| [docs/cg_ensemble_execution.md](docs/cg_ensemble_execution.md) | Ensemble / batch production workflow |
| [docs/development/changelog.md](docs/development/changelog.md) | Version history (v1–v32+) |
| [tests/VERSION_HISTORY.md](tests/VERSION_HISTORY.md) | Output directory reference map |

## License
Proprietary / Internal Use.
