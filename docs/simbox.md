# simbox — Molecule Packing Sub-System

`topon.simbox` is a general-purpose sub-system for packing individual molecules into a periodic simulation box and generating LAMMPS input scripts for crosslinking studies. It operates independently of the polymer network pipeline and is oriented toward Epoxy-PDMS / Amino-PDMS / POSS hybrid systems.

---

## Overview

The simbox workflow:

```
MoleculeLibrary         →  Molecule objects (RDKit mol + reactive-site annotations)
       ↓
BoxPacker.pack()        →  PackedBox (placed molecules with 3D coordinates)
       ↓
assemble(packed)        →  AssembledSystem (merged RDKit mol, reactive-site registry)
       ↓
write_lammps(system)    →  system.data, settings.in, groups.txt, ff_coeffs.in
       ↓
write_inputs(system)    →  1_minimize.in, 2_nvt.in, 3_npt.in, 4b_crosslink.in
```

The quickest way to run the full workflow is the CLI:

```bash
topon simbox --output output/simbox_run --n-epoxy 600 --n-amino 300
```

Or call the Python API directly via `topon.simbox.workflow.run_workflow` (see [Workflow module](#topon-simbox-workflow--full-workflow) below).

The legacy script `tests/workflows/generate_simbox_crosslink.py` delegates to the same `run_workflow` function.

---

## Modules

### `topon.simbox.molecule` — Molecule

```python
from topon.simbox.molecule import Molecule

mol = Molecule.from_smiles("EpoxyPDMS", "C1OC1COCCC[Si](C)(C)O...")
mol = Molecule.from_pdb("MyMol", "path/to/file.pdb")
mol = Molecule.from_mol("MyMol", rdkit_mol_object)
```

Each `Molecule` stores:
- `mol` — RDKit Mol with explicit H and a 3D conformer (ETKDGv3 + MMFF optimised)
- `mw` — exact molecular weight (g/mol)
- `reactive_sites` — dict mapping group name → list of heavy-atom indices, auto-detected via SMARTS:
  - `"epoxide"` — oxirane ring (`[C]1[O][C]1`)
  - `"primary_amine"` — `-NH2` (`[NX3;H2;!$([NH2]C=O)]`)
  - `"secondary_amine"` — `>NH` (`[NX3;H1]([#6])[#6]`)

---

### `topon.simbox.library` — MoleculeLibrary

Pre-built molecules for siloxane crosslinking experiments:

```python
from topon.simbox.library import MoleculeLibrary

lib = MoleculeLibrary()

epoxy  = lib.epoxy_pdms(n_dms=2)   # Glycidoxypropyl-PDMS, ~500 g/mol
amino  = lib.amino_pdms(n_dms=8)   # Aminopropyl-PDMS, ~850 g/mol
poss   = lib.am0270_poss()          # AminopropylIsooctyl POSS, ~1267 g/mol
custom = lib.custom("C1OC1", name="MyEpoxide")
```

**Epoxy-PDMS** (`n_dms` = number of internal DMS repeat units):
```
Epoxide-CH2-O-CH2CH2CH2-Si(Me)-[O-Si(Me)2]_n-O-Si(Me)-CH2CH2CH2-O-CH2-Epoxide
```

**Amino-PDMS**:
```
H2N-CH2CH2CH2-Si(Me)-[O-Si(Me)2]_n-O-Si(Me)-CH2CH2CH2-NH2
```

**AM0270 POSS** (AminopropylIsooctyl POSS):
- Si8O12 cube cage (12 bridging oxygens)
- Corner 0: `-CH2CH2CH2-NH2` (reactive amine)
- Corners 1–7: 2,4,4-trimethylpentyl (isooctyl, inert)

---

### `topon.simbox.packer` — BoxPacker

Packs molecules into a periodic box using grid-based spatial hashing for overlap detection.

```python
from topon.simbox.packer import BoxPacker

packer = BoxPacker(
    density=0.85,        # target density in g/cm³
    min_dist=2.0,        # minimum interatomic distance in Å
    seed=42,             # random seed for reproducibility
    max_attempts=1000,   # placement attempts per molecule
    growth_factor=1.05,  # box expansion factor when packing fails
)

packed = packer.pack([
    (epoxy_mol, 100),    # 100 epoxy molecules
    (amino_mol, 50),     # 50 amino molecules
    (poss_mol, 10),      # 10 POSS cages
])

print(f"Box: {packed.box_lengths} Å")
print(f"Total atoms: {packed.total_atoms}")
```

**Algorithm:**
1. Compute initial box size from total molecular mass and target density.
2. Shuffle insertion order.
3. For each molecule: random rotation (Shoemake quaternion) + random translation; check overlap with minimum-image convention.
4. If placement fails after `max_attempts`, grow the box by `growth_factor` and retry (up to 20 growth rounds).

---

### `topon.simbox.system` — AssembledSystem

Merges all placed molecules into a single RDKit Mol with global bookkeeping:

```python
from topon.simbox.system import assemble

system = assemble(packed)

print(f"Total atoms: {system.mol.GetNumAtoms()}")
print(f"Total molecules: {system.num_molecules}")
print(f"Reactive sites: {len(system.reactive_sites)}")

for site in system.reactive_sites:
    print(f"  mol {site.molecule_id}: {site.group_name} at atom {site.global_atom_idx}")
```

`AssembledSystem` fields:
- `mol` — merged RDKit Mol (all atoms, 3D coordinates)
- `box_lengths` — `ndarray([Lx, Ly, Lz])` in Å
- `molecule_ids` — per-atom LAMMPS molecule ID (1-based)
- `species_names` — per-molecule species name
- `reactive_sites` — list of `ReactiveSiteEntry` (global atom index + group name)

---

### `topon.simbox.writer` — LAMMPS Data File

Writes a DREIDING-compatible `atom_style full` data file:

```python
from topon.simbox.writer import write_lammps

paths = write_lammps(system, output_dir="run/02_Chemistry")
# writes: system.data, settings.in, groups.txt, ff_coeffs.in
```

---

### `topon.simbox.inputs` — LAMMPS Input Scripts

Generates a 4-stage simulation protocol:

```python
from topon.simbox.inputs import write_inputs

paths = write_inputs(
    system,
    output_dir="run/04_Simulation",
    temperature=300.0,
    pressure=1.0,
)
# writes:
#   1_minimize.in   — soft push-off (pair/soft) → CG minimisation
#   2_nvt.in        — NVT thermalisation
#   3_npt.in        — NPT density equilibration
#   4b_crosslink.in — crosslink reaction template (commented, two options)
```

**Stage 1 — Soft push-off + minimisation:**
- Phase A: `pair_style soft` with ramped prefactor (0→60) + brief NVT to resolve overlaps
- Phase B: Switch to `lj/cut` DREIDING potentials + conjugate-gradient minimisation

**Stage 4b — Crosslink template** (user must configure):
- **Option A** (`fix bond/react`): template-based reactions with molecule pre/post files
- **Option B** (`fix bond/create`): simple distance-based bond formation

---

### `topon.simbox.workflow` — Full Workflow

The canonical high-level entry point.  Handles molecule building, packing, universal type-ID mapping, and writing in one call.

```python
from topon.simbox.workflow import run_workflow

files = run_workflow(
    output_dir="output/simbox_run",
    n_epoxy=600,
    n_amino=300,
    n_poss=0,
    density=0.85,
    seed=42,
)
```

Also exposes `UniversalTypeMapper` — a context manager that patches `topon.forcefield.dreiding` at write time to enforce stable DREIDING type IDs across all compositions, keeping pre-defined LAMMPS `fix bond/react` templates compatible.

**CLI equivalent:**

```bash
topon simbox --output output/simbox_run --n-epoxy 600 --n-amino 300 --n-poss 0 --seed 42
```

---

## Example Usage

```python
from topon.simbox.library import MoleculeLibrary
from topon.simbox.packer import BoxPacker
from topon.simbox.system import assemble
from topon.simbox.writer import write_lammps
from topon.simbox.inputs import write_inputs

lib = MoleculeLibrary()
epoxy = lib.epoxy_pdms(n_dms=2)
amino = lib.amino_pdms(n_dms=8)
poss  = lib.am0270_poss()

packer = BoxPacker(density=0.85, seed=42)
packed = packer.pack([(epoxy, 80), (amino, 40), (poss, 10)])

system = assemble(packed)
write_lammps(system, output_dir="output/02_Chemistry")
write_inputs(system, output_dir="output/04_Simulation", temperature=300.0)
```

Or use the CLI (recommended):

```bash
topon simbox --output output/simbox_run --n-epoxy 80 --n-amino 40 --n-poss 10 --seed 42
```

Or the Python workflow function:

```python
from topon.simbox.workflow import run_workflow
run_workflow("output/simbox_run", n_epoxy=80, n_amino=40, n_poss=10, seed=42)
```

---

## Output Structure

All files are written to a single flat output directory by `write_lammps` + `write_inputs`:

```
output_dir/
  system.data        ← LAMMPS data file (atom_style full, DREIDING)
  settings.in        ← pair_coeff, bond_coeff, angle_coeff, dihedral_coeff
  ff_coeffs.in       ← all *_coeff commands for use with fix/bond/react templates
  groups.txt         ← LAMMPS group definitions by reactive-group type
  1_minimize.in      ← soft push-off (pair/soft) → CG minimisation
  2_nvt.in           ← NVT thermalisation
  3_npt.in           ← NPT density equilibration
  4b_crosslink.in    ← crosslink template (fix bond/react or fix bond/create)
```

---

## Versioned Development History

The simbox sub-system has been developed through 15 iterations:

| Version | Key Change |
|---------|-----------|
| v2.0 | Initial simbox packing with Epoxy-PDMS + Amino-PDMS |
| v2.1–v2.9 | Refinements to packing algorithm, overlap detection, atom typing |
| v2.10 | `setup_crosslink.py` — manual crosslink setup helper |
| v2.11 | `gen_pre_template.py`, `check_template_match.py` — bond/react templates |
| v2.11 (poss_only) | `patch_data.py` — POSS-only system test |
| v2.12–v2.14 | Improved settings and atom type consistency |
| v2.15 | `verify_consistency.py` — data file validation |
| v3 | `generate_all.py` — full parameter sweep across POSS fractions (poss_0–poss_6) |
| v4 | `generate_simbox_crosslink.py` — canonical reference; fixed writer header bug (`max` not `len` for type counts); added `ff_coeffs.in` generation; UniversalTypeMapper enforces N_3=4, H_=5 ordering |
