# Topon — Implementation Plan & Architecture

**Last updated:** 2026-03-23

This document describes the architecture of the `topon` package, how its modules interact, and the design intent behind each layer.

---

## 1. Design Philosophy

Topon separates **topology** (graph structure) from **chemistry** (atom types, coordinates, force-field) and **I/O** (LAMMPS files). This separation means the same graph can be used to produce a CG Kremer-Grest system *or* a fully atomistic DREIDING system without changing the topology logic.

---

## 2. Four-Stage Pipeline

```
Stage 1: Topology
    PythonTopologyGenerator / generator.exe
    → NetworkX MultiGraph (nodes = junctions, edges = chains)

Stage 2: Assignment
    AssignmentManager
    → node types, edge types, DP distribution, grafts,
      copolymer sequences, entanglements, defects

Stage 3: Chemistry
    ChemistryBuilder (CG or Atomistic)
    → RDKit Mol (all atoms, bonds, 3D coordinates via conformers)
    → overlap resolution (OverlapResolver)

Stage 4: Writers + I/O
    LammpsWriter (CG) / AtomisticLammpsWriter (Atomistic)
    → system.data, minimize.in, nvt.in, npt.in
    SimulationRunner → subprocess LAMMPS
```

The four stages are implemented as separate classes that pass data through well-defined interfaces. The primary entry point for each stage is documented below.

---

## 3. Module Reference

### 3.1 `topon/topology/`

| Module | Class/Function | Purpose |
|--------|---------------|---------|
| `generator.py` | `TopologyGenerator` | Wraps `generator.exe` C binary (subprocess call) |
| `generator_python.py` | `PythonTopologyGenerator` | Pure-Python port; preferred when C binary unavailable |
| `loader.py` | `TopologyLoader` | Loads `.nodes`/`.edges`/`.gpickle`/`.graphml` files |

**Graph format:** `networkx.MultiGraph` with:
- Node attributes: `degree` (int), `node_type` (str)
- Edge attributes: `edge_type` (str), `dp` (int), `sequence` (list), `entangled` (bool)

### 3.2 `topon/assignment/`

| Module | Class | Purpose |
|--------|-------|---------|
| `manager.py` | `AssignmentManager` | Orchestrates all assignment steps from config |
| `attributor.py` | `GraphAttributor` | Re-attributes an existing graph (for round-trip workflows) |
| `node_types.py` | `NodeTypeAssigner` | Degree-based, positional, random, explicit assignment |
| `edge_types.py` | `EdgeTypeAssigner` | Uniform, random, composite assignment |
| `dp_distribution.py` | `DPDistributor` | Schulz-Zimm PDI distribution |
| `defects.py` | `inject_primary_loops()` | Parallel edges (loops) as defects; valence protection |
| `entanglements.py` | `EntanglementAssigner` | Gaussian Kink geometry for topological knots |
| `sequences.py` | `SequenceBuilder` | Block, Random, Alternating, Gradient copolymers |

`AssignmentManager.assign(graph, config)` is the single entry point for Stage 2.

### 3.3 `topon/chemistry/`

| Module | Class | Purpose |
|--------|-------|---------|
| `builder_cg.py` | `CGChemistryBuilder` | Kremer-Grest bead placement (FENE bonds, LJ) |
| `builder_atomistic.py` | `AtomisticChemistryBuilder` | DREIDING atom placement; PDMS/Si/POSS chemistries |
| `sequences.py` | — | Shared sequence helpers |

Both builders accept an annotated NetworkX graph and return an RDKit `Mol` object with 3D coordinates.

### 3.4 `topon/conformation/`

| Module | Class | Purpose |
|--------|-------|---------|
| `conformer.py` | `ConformerGenerator` | Per-chain ETKDGv3 conformer generation |
| `overlap.py` | `OverlapResolver` | Iterative overlap correction (gradient descent) |

### 3.5 `topon/writers/`

| Module | Class | Purpose |
|--------|-------|---------|
| `lammps_cg.py` | `LammpsWriter` | CG LAMMPS data file + input scripts |
| `lammps_atomistic.py` | `AtomisticLammpsWriter` | Atomistic LAMMPS data file + input scripts |
| `graphml.py` | `GraphMLWriter` | GraphML export (dual-graph with entanglement edges) |

### 3.6 `topon/forcefield/`

| Module | Purpose |
|--------|---------|
| `dreiding.py` | DREIDING parameter file parser; bond/angle/dihedral assignment |
| `kg.py` | Kremer-Grest parameters (FENE, LJ) |

### 3.7 `topon/simbox/`

Independent sub-system for molecule packing. See [../simbox.md](../simbox.md) for full documentation.

| Module | Class | Purpose |
|--------|-------|---------|
| `molecule.py` | `Molecule` | RDKit mol + reactive-site annotations |
| `library.py` | `MoleculeLibrary` | Pre-built Epoxy-PDMS, Amino-PDMS, POSS |
| `packer.py` | `BoxPacker` | Grid-based random packing with overlap detection |
| `system.py` | `assemble()` | Merges packed molecules into unified system |
| `writer.py` | `write_lammps()` | DREIDING LAMMPS data + ff_coeffs.in |
| `inputs.py` | `write_inputs()` | 4-stage LAMMPS input scripts |

### 3.8 `topon/config/`

Pydantic schema for JSON configuration files. `ConfigLoader.load(path)` returns a validated `ToponConfig` object used by `AssignmentManager`.

### 3.9 `topon/simulation/`

`SimulationRunner` wraps a LAMMPS subprocess. Configured via `config.execution` (auto_run, executable, n_procs).

---

## 4. Data Flow Diagram

```
config.json
    │
    ▼
ConfigLoader.load()
    │
    ├─► TopologyGenerator.generate()  ──►  NetworkX MultiGraph
    │                                            │
    ├─► AssignmentManager.assign()    ◄──────────┘
    │       │  (annotates graph in-place)
    │       ▼
    │   Annotated MultiGraph
    │       │
    ├─► ChemistryBuilder.build()      ──►  RDKit Mol (3D coords)
    │                                            │
    ├─► OverlapResolver.resolve()     ◄──────────┘
    │       │
    │       ▼
    │   Resolved Mol
    │       │
    └─► LammpsWriter.write()          ──►  system.data, *.in scripts
            │
            ▼
        SimulationRunner.run()        ──►  LAMMPS subprocess
```

---

## 5. Configuration Schema

Key config fields (see `examples/config_cg.json` for full example):

```json
{
  "topology": {
    "source": "generate",           // or "load"
    "lattice_source": "SC",
    "lattice_size": "5x5x5",
    "periodicity": "111"
  },
  "chemistry": {
    "model_type": "coarse_grained", // or "atomistic"
    "degree_of_polymerization": 20
  },
  "assignment": {
    "entanglements": {"enabled": true, "target": 5},
    "grafts": {"enabled": false},
    "defects": {"primary_loops": 0}
  },
  "simulation": {
    "pair_style": "repulsive",      // or "attractive"
    "include_angles": true,
    "remove_cg_angles": false
  },
  "execution": {
    "auto_run": false,
    "executable": "lmp",
    "n_procs": 4
  }
}
```

---

## 6. Workflow Scripts

The `Pipeline` class (`topon/pipeline.py`) is a **stub** — all stage methods raise `NotImplementedError`. The real entry points are the workflow scripts in `tests/workflows/`:

| Script | What it does |
|--------|-------------|
| `generate_cg.py` | Baseline CG network |
| `generate_cg_entangled.py` | CG + Gaussian Kink entanglements |
| `generate_atomistic.py` | Baseline DREIDING atomistic |
| `generate_atomistic_combined.py` | Atomistic + grafts + entanglements |
| `generate_atomistic_poss.py` | Atomistic + POSS junction nodes |
| `generate_simbox_crosslink.py` | simbox: pack Epoxy/Amino/POSS → crosslink setup |
| `run_cg_ensemble.py` | Parallel ensemble (multiprocessing) |

---

## 7. Testing Strategy

| Layer | Location | Coverage |
|-------|----------|---------|
| Unit tests | `tests/unit/` | topology generator, assignment, chemistry, config, simbox |
| Regression tests | `tests/regression/` | CG header/workflow, Atomistic header/workflow, simbox crosslink |
| Run: | `pytest tests/unit/` | Fast (~5s) |
| Run: | `pytest tests/regression/` | Slow (~5 min, runs full workflows) |

---

## 8. Known Limitations

- `Pipeline` class and `topon generate` CLI are not functional (stubs).
- `gui/app.py` does not exist.
- POSS crosslink regression reference only covers `n_poss=0`; POSS fraction sweep has no automated regression test.
- `topon/external/auto_martini/` and `topon/data/martini3/` are vendored but unused (Martini v24.1 experiment was archived).
