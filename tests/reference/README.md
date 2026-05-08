# Regression Reference Outputs

These are the frozen reference outputs used to verify that the restructured
`topon` package produces exactly the same files as the original implementation.

**Do not regenerate these.** They were produced by the original code after months
of debugging and refinement. They are the ground truth.

---

## Reference Outputs by Use Case

### CG Polymer Network
- **Path**: `tests/output/v21_cg_combined/cg_combined/`
- **Workflow**: `tests/workflows/generate_cg_combined.py`
- **Key files for comparison**:
  - `02_Chemistry/system.data` — LAMMPS data file (atom types, bonds, beads)
  - `02_Chemistry/system.in.settings` — force field settings
  - `02_Chemistry/system_nodes.displace` — node displacement file
  - `02_Chemistry/system_beads.displace` — bead displacement file
  - `03_Conformation/system_conformed.data` — after coordinate placement
  - `04_Simulation/minimize_1_serial.in` — LAMMPS input script

### Atomistic PDMS Network (DREIDING)
- **Path**: `tests/output/v21_atomistic_combined/atomistic_combined/`
- **Workflow**: `tests/workflows/generate_atomistic_combined.py`
- **Key files for comparison**:
  - `02_Chemistry/system.data` — LAMMPS data file (atom types, DREIDING params)
  - `02_Chemistry/system.in.settings` — DREIDING force field settings
  - `02_Chemistry/system_nodes.displace`
  - `02_Chemistry/system_backbone.displace`
  - `02_Chemistry/system_pendant.displace`
  - `02_Chemistry/system_hydrogens.displace`
  - `03_Conformation/system_conformed.data`

### Atomistic with POSS (latest, v32)
- **Path**: `tests/output/v32_poss_sweep/poss_0/poss_frac_0.00/`
- **Workflow**: `tests/workflows/generate_poss_dataset.py` or `run_poss_sweep.py`
- **Key files for comparison**:
  - `02_Chemistry/system.data`
  - `03_Conformation/system_relaxed.data`
  - `04_Simulation/minimize_1_serial.in`

### SimBox Crosslink (POSS + Epoxy-PDMS + Amino-PDMS)
- **Path**: `tests/output/simbox_crosslink/v3/poss_0/`
- **Workflow**: `tests/workflows/generate_simbox_crosslink.py`
- **Key files for comparison**:
  - `system.data` — LAMMPS data file with universal type map
  - `groups.txt` — reactive group definitions
  - `1_minimize.in`, `2_nvt.in`, `3_npt.in` — LAMMPS protocols
  - `ff_coeffs.in` — force field coefficients
  - `settings.in` — simulation settings
  - `pre_react_epoxy_amine.mol`, `post_react_full.mol` — reaction templates

---

## What to Compare in Each Data File

When running regression tests, compare in this order of importance:

### 1. Header counts (must be exact)
```
N atoms
M bonds
K angles
J dihedrals
I impropers
X atom types
Y bond types
Z angle types
...
```

### 2. Coefficient sections (must be exact)
- `Pair Coeffs` — all types present, correct parameters
- `Bond Coeffs`
- `Angle Coeffs`
- `Dihedral Coeffs`
- `Improper Coeffs` (atomistic only)

### 3. Topology (must be exact)
- `Bonds` section — all bond pairs
- `Angles` section — all angle triples
- `Dihedrals` section

### 4. Atom types and charges (must be exact)
- Column 3 (type ID) in `Atoms` section
- Column 4 (charge) in `Atoms` section

### 5. Coordinates (tolerance: 1e-4 Å)
- Columns 5-7 (x, y, z) in `Atoms` section
- Exact match only required for pre-conformation `02_Chemistry/system.data`

### 6. LAMMPS input scripts (must be exact)
- `pair_style`, `pair_coeff` lines
- `bond_style`, `bond_coeff` lines
- `angle_style`, `angle_coeff` lines
- `dihedral_style`, `dihedral_coeff` lines
- Stage sequence (minimize → NVT → NPT)

---

## Comparison Logic (for test implementation)

```python
def compare_lammps_data(new_path, ref_path, coord_tol=1e-4):
    """
    Returns list of discrepancies. Empty list = pass.
    """
    # Parse both files into structured dicts
    # Compare header counts exactly
    # Compare all coeff sections exactly
    # Compare bonds/angles/dihedrals exactly
    # Compare atom type IDs and charges exactly
    # Compare coordinates within coord_tol
    ...
```

Test files go in `tests/regression/`:
- `test_cg_network.py`
- `test_atomistic_network.py`
- `test_simbox_crosslink.py`
- `test_lammps_writer.py`
