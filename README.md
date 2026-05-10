# Topon — Topological Polymer Network Generator

**Topon** is a modular Python package for generating complex polymer (and, via **topro**, protein) network structures for Molecular Dynamics simulations. It separates graph-theoretical topology from chemistry, supporting both Coarse-Grained (Kremer-Grest) and Atomistic (DREIDING) models, with a parallel MARTINI 3 protein-network sub-system.

## Key features

- **Graph-based topology** — define network connectivity first, then map to chemistry. The same graph can produce a CG or an atomistic system.
- **Dual resolution**:
  - **Coarse-Grained** — Kremer-Grest model (FENE / Harmonic bonds, LJ potentials).
  - **Atomistic** — DREIDING force field with PDMS, fluorinated, phenyl, and POSS chemistries.
- **Architecture knobs**:
  - **Entanglements** via Gaussian Kink geometry (single + multi-entanglement).
  - **Copolymers** with Block, Random, Alternating, Gradient sequences.
  - **Grafts / side chains** with dynamic DP scaling.
  - **Defects** (primary loops / parallel edges) with valence protection.
  - **POSS junctions** (Si₈O₁₂ cage chemistry).
  - **Custom topology** loaded from `.nodes`/`.edges`/`.gpickle`/`.graphml`.
- **Sub-systems**:
  - **simbox** — independent molecule packer for crosslinking studies (Epoxy-PDMS, Amino-PDMS, AM0270-POSS).
  - **topro** — `topon.protein_network`, MARTINI 3 protein-network generator from a residue-sequence string.
  - **singlechain** — single chain in solvent for solubility studies.
- **Workflow automation** — `SimulationRunner` for LAMMPS subprocess execution.

## Installation

```bash
git clone https://github.com/lynspica/topon.git
cd topon
pip install -e .
```

LAMMPS (`lmp` on `PATH`) is required only to *run* generated systems, not to generate them.

## Quick start

```bash
# Generate a starter config
topon init --output my_run.json

# Run the full six-stage pipeline
topon generate my_run.json --output ./runs
```

For ready-to-use configs, see `examples/config_*.json`. For full CLI options, recipes, and config schema, see [docs/USAGE.md](docs/USAGE.md).

## Documentation

| Document | Description |
|---|---|
| [AGENTS.md](AGENTS.md) | **AI-agent onboarding** — read this first if you (or your AI assistant) are new to the project |
| [docs/USAGE.md](docs/USAGE.md) | **How to run topon** — CLI reference, sub-system APIs, recipes, JSON-config schema (appendix) |
| [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | **How topon is organised** — six-stage pipeline, module map, design principles |
| [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) | Version-by-version history + Q→R→R→R methodology |
| [docs/JOURNAL.md](docs/JOURNAL.md) | Engineering journal — dated entries on changes, issues, fixes |
| [CLAUDE.md](CLAUDE.md) | Rules every change must follow |
| [tests/VERSION_HISTORY.md](tests/VERSION_HISTORY.md) | Test-output directory reference map |

## License

Proprietary / Internal Use.
