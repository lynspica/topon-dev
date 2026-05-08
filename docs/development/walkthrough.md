# Topon — Step-by-Step Walkthrough

This guide walks through generating a polymer network system from scratch using the two main pipelines: **Coarse-Grained (CG)** and **Atomistic (DREIDING)**. A third section covers the **simbox** crosslink packing workflow.

---

## Prerequisites

```bash
pip install -e .           # install topon in editable mode
# LAMMPS must be on PATH as 'lmp' for simulation runs (optional for generation)
```

Dependencies: `numpy`, `networkx`, `rdkit`, `pydantic`, `scipy`.

---

## Part 1 — Coarse-Grained (Kremer-Grest) Network

### Step 1: Generate a topology

The topology is a cubic lattice that is "sculpted" to a target degree distribution. The Python generator is the default:

```python
from topon.topology.generator_python import PythonTopologyGenerator
from collections import namedtuple

TopologyConfig = namedtuple(
    'TopologyConfig',
    ['lattice_source', 'lattice_size', 'periodicity',
     'degree_distribution', 'max_functionality']
)
config = TopologyConfig(
    lattice_source="SC",
    lattice_size="5x5x5",
    periodicity="111",
    degree_distribution="",  # default: no explicit degree targets
    max_functionality=6
)
gen = PythonTopologyGenerator(config)
graphs = gen.generate(trials=20)
graph = graphs[0]   # NetworkX MultiGraph: nodes=junctions, edges=chains
```

### Step 2: Assign types and DP

```python
from topon.assignment.manager import AssignmentManager
from topon.config.loader import ConfigLoader

cfg = ConfigLoader.load("examples/config_cg.json")
manager = AssignmentManager(cfg)
manager.assign(graph)
# graph now has: node_type ('J'/'E'), edge_type ('A'), dp (int), sequence (list)
```

### Step 3: Build CG chemistry

```python
from topon.chemistry.builder_cg import CGChemistryBuilder

builder = CGChemistryBuilder(cfg)
mol = builder.build(graph)
# mol: RDKit Mol with explicit H and 3D coordinates (one bead per repeat unit)
```

### Step 4: Write LAMMPS files

```python
from topon.writers.lammps_cg import LammpsWriter
from pathlib import Path

writer = LammpsWriter(cfg)
output = Path("output/my_cg_system")
output.mkdir(parents=True, exist_ok=True)
paths = writer.write(mol, graph, output_dir=output)
# writes: system.data, minimize_1_parallel.in, nvt.in, npt.in
```

### Shortcut: use a workflow script

```bash
python tests/workflows/generate_cg.py --output output/cg_baseline
```

For entanglements:
```bash
python tests/workflows/generate_cg_entangled.py \
    --output output/cg_entangled \
    --n_entanglements 5
```

For ensemble production (1000 systems, 8 workers):
```bash
python tests/workflows/run_cg_ensemble.py --workers 8 --output output/ensemble
```

---

## Part 2 — Atomistic (DREIDING) Network

The atomistic workflow follows the same four stages but uses DREIDING atom types (Si3, O_3, C_3, H_, N_3) and generates a full-chemistry RDKit Mol.

### Step 1–2: Topology and assignment (same as CG)

Use `examples/config_atomistic.json` instead of `config_cg.json`.

### Step 3: Build atomistic chemistry

```python
from topon.chemistry.builder_atomistic import AtomisticChemistryBuilder

builder = AtomisticChemistryBuilder(cfg)
mol = builder.build(graph)
# mol: RDKit Mol with Si-O-Si backbone, explicit H, 3D conformers
```

Atom counts scale with degree of polymerization. A 125-node network with DP=20 yields ~10,000 atoms.

### Step 4: Write LAMMPS files

```python
from topon.writers.lammps_atomistic import AtomisticLammpsWriter

writer = AtomisticLammpsWriter(cfg)
paths = writer.write(mol, graph, output_dir=output)
# writes: system.data (atom_style full), settings.in, minimize.in, nvt.in, npt.in
```

### Shortcut: use a workflow script

```bash
# Baseline
python tests/workflows/generate_atomistic.py --output output/atomistic_baseline

# With grafts + entanglements
python tests/workflows/generate_atomistic_combined.py \
    --output output/atomistic_combined

# With POSS junction nodes
python tests/workflows/generate_atomistic_poss.py \
    --output output/atomistic_poss
```

---

## Part 3 — simbox Crosslink Workflow

The `simbox` sub-system is independent of the network pipeline. It packs arbitrary molecules (Epoxy-PDMS, Amino-PDMS, AM0270-POSS) into a periodic box and generates LAMMPS input for crosslinking studies.

### Step 1: Build molecules from the library

```python
from topon.simbox.library import MoleculeLibrary

lib = MoleculeLibrary()
epoxy = lib.epoxy_pdms(n_dms=2)   # ~71 atoms, 2 DMS repeat units
amino = lib.amino_pdms(n_dms=8)   # ~117 atoms, 8 DMS repeat units
poss  = lib.am0270_poss()          # ~207 atoms, Si8O12 cage
```

Each `Molecule` carries reactive-site annotations auto-detected by SMARTS:
- `"epoxide"` — oxirane rings
- `"primary_amine"` — `-NH2`
- `"secondary_amine"` — `>NH`

### Step 2: Pack into a box

```python
from topon.simbox.packer import BoxPacker

packer = BoxPacker(
    density=0.85,       # g/cm³
    min_dist=2.0,       # minimum interatomic distance (Å)
    seed=42,            # for reproducibility
)
packed = packer.pack([
    (epoxy, 600),       # 600 Epoxy-PDMS molecules
    (amino, 300),       # 300 Amino-PDMS molecules
    (poss,    0),       # 0 POSS (pure epoxy/amino system)
])
# packed.total_atoms == 77700, packed.box_lengths ~= [168, 168, 168] Å
```

The packer uses grid-based spatial hashing (cell size = `min_dist`) for O(N) overlap detection. If a molecule cannot be placed after `max_attempts`, the box expands by `growth_factor` and placement retries.

### Step 3: Assemble the system

```python
from topon.simbox.system import assemble

system = assemble(packed)
# system.mol  — merged RDKit Mol (77700 atoms)
# system.reactive_sites  — list of ReactiveSiteEntry with global atom indices
```

### Step 4: Write LAMMPS files

```python
from topon.simbox.writer import write_lammps
from topon.simbox.inputs import write_inputs

write_lammps(system, output_dir="output/simbox")
write_inputs(system, output_dir="output/simbox", temperature=300.0)
```

Output files:
```
output/simbox/
  system.data        ← DREIDING data file (atom_style full)
  settings.in        ← pair/bond/angle/dihedral coefficients
  ff_coeffs.in       ← all *_coeff commands (for fix/bond/react templates)
  groups.txt         ← group definitions by reactive-group type
  1_minimize.in      ← soft push-off + CG minimisation
  2_nvt.in           ← NVT thermalisation
  3_npt.in           ← NPT density equilibration
  4b_crosslink.in    ← crosslink script template
```

### Step 5: Run crosslinking in LAMMPS

Edit `4b_crosslink.in` to choose the reaction strategy:
- **Option A** (`fix bond/react`): template-based reactions using `.mol` files
- **Option B** (`fix bond/create`): distance-based bond creation (simpler)

Then run the sequence:
```bash
lmp -in 1_minimize.in
lmp -in 2_nvt.in
lmp -in 3_npt.in
lmp -in 4b_crosslink.in
```

### Shortcut: CLI (recommended)

```bash
topon simbox --output output/simbox_crosslink \
    --n-epoxy 600 --n-amino 300 --n-poss 0 \
    --seed 42
```

### Shortcut: Python API

```python
from topon.simbox.workflow import run_workflow
run_workflow("output/simbox_crosslink", n_epoxy=600, n_amino=300, n_poss=0, seed=42)
```

---

## Common Options

### Entanglements

Add topological knots (Gaussian Kinks) via config:

```json
"assignment": {
    "entanglements": {"enabled": true, "target": 5}
}
```

Or via the entangled workflow script (`generate_cg_entangled.py`, `generate_atomistic_combined.py`).

### Copolymers

```json
"chemistry": {
    "copolymer": {"enabled": true, "sequence": "block", "fraction_B": 0.3}
}
```

### Grafts (side chains)

```json
"assignment": {
    "grafts": {"enabled": true, "per_edge_type": {"A": {"graft_density": 0.05}}}
}
```

### Pair style

```json
"simulation": {"pair_style": "repulsive"}   // WCA, ~4× faster equilibration
"simulation": {"pair_style": "attractive"}  // LJ 2.5σ cutoff
```

---

## Running Tests

```bash
# Unit tests (fast, ~5s)
pytest tests/unit/ -v

# Regression tests (slow, ~5 min each)
pytest tests/regression/test_cg_network.py -v
pytest tests/regression/test_atomistic_network.py -v
pytest tests/regression/test_simbox_crosslink.py -v

# All regression tests
pytest tests/regression/ -v
```

Regression tests regenerate systems from scratch and compare byte-level against versioned reference outputs in `tests/output/`.
