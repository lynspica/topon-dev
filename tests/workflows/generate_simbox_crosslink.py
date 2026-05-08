"""
SimBox Crosslink Workflow
=========================
Generates a simulation box with epoxy-PDMS, amino-PDMS, and AM0270 POSS
molecules for crosslink reaction simulations.

Outputs a complete LAMMPS-ready directory:
  - system.data           (DREIDING-parameterised data file)
  - groups.txt            (reactive group definitions)
  - 1_minimize.in         (soft push-off + energy minimisation)
  - 2_nvt.in              (NVT equilibration)
  - 3_npt.in              (NPT equilibration)
  - 4b_crosslink.in       (crosslinking setup - configure & run)

After running this script, the only remaining step is to run LAMMPS:
  cd <output_dir>
  lmp -in 1_minimize.in
  lmp -in 2_nvt.in
  lmp -in 3_npt.in
  # Then configure 4b_crosslink.in and run it.

Usage:
  python generate_simbox_crosslink.py [--output OUTPUT_DIR]
                                       [--n_epoxy N] [--n_amino N] [--n_poss N]
                                       [--density RHO] [--seed SEED]
"""

import sys
import argparse
import time
from pathlib import Path
import functools

# Add package root to path
pkg_dir = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(pkg_dir))

from topon.simbox import SimBox, MoleculeLibrary, Molecule
import topon.forcefield.dreiding as dreiding

# ============================================================================
# UNIVERSAL TYPE MAPS (Derived from test_mixed/system.data)
# ============================================================================
# These maps enforce consistent Type IDs across all simulation variations
# (Amino-only, POSS-only, Mixed), ensuring compatibility with pre-defined
# reaction templates (e.g. bond/react).

ATOM_MAP = {
    'Si3': 1, 'O_3': 2, 'C_3': 3, 'N_3': 4, 'H_': 5
}

# NOTE: Keys must be sorted (canonical) for Bonds
BOND_MAP = {
    ('O_3', 'Si3'): 1,
    ('C_3', 'Si3'): 2,
    ('C_3', 'C_3'): 3,
    ('C_3', 'N_3'): 4,
    ('H_', 'Si3'): 5,
    ('C_3', 'H_'): 6,
    ('H_', 'N_3'): 7,
    ('C_3', 'O_3'): 8,
    ('H_', 'O_3'): 9
}

# NOTE: Keys must be canonical: (Atom1, Atom2, Atom3) where Atom1 <= Atom3
ANGLE_MAP = {
    ('C_3', 'Si3', 'O_3'): 1,
    ('H_', 'Si3', 'O_3'): 2,
    ('C_3', 'Si3', 'C_3'): 3,
    ('C_3', 'Si3', 'H_'): 4,
    ('O_3', 'Si3', 'O_3'): 5,
    ('Si3', 'O_3', 'Si3'): 6,
    ('H_', 'C_3', 'Si3'): 7,
    ('H_', 'C_3', 'H_'): 8,
    ('C_3', 'C_3', 'Si3'): 9,
    ('C_3', 'C_3', 'H_'): 10,
    ('C_3', 'C_3', 'C_3'): 11,
    ('C_3', 'C_3', 'N_3'): 12,
    ('H_', 'C_3', 'N_3'): 13,
    ('C_3', 'N_3', 'H_'): 14,
    ('H_', 'N_3', 'H_'): 15,
    ('C_3', 'C_3', 'O_3'): 16,
    ('H_', 'C_3', 'O_3'): 17,
    ('C_3', 'O_3', 'C_3'): 18,
    ('C_3', 'O_3', 'H_'): 19,
    ('C_3', 'N_3', 'C_3'): 20
}

# Dihedral Map: Can have multiple IDs for validation purposes (duplicates allowed in list)
DIHEDRAL_MAP_LIST = [
    (('O_3', 'Si3', 'O_3', 'Si3'), 1),
    (('H_', 'Si3', 'O_3', 'Si3'), 2),
    (('O_3', 'Si3', 'O_3', 'Si3'), 3),
    (('H_', 'C_3', 'Si3', 'O_3'), 4),
    (('C_3', 'Si3', 'C_3', 'H_'), 5),
    (('H_', 'C_3', 'Si3', 'H_'), 6),
    (('C_3', 'C_3', 'Si3', 'O_3'), 7),
    (('C_3', 'C_3', 'Si3', 'C_3'), 8),
    (('C_3', 'C_3', 'Si3', 'H_'), 9),
    (('C_3', 'C_3', 'C_3', 'Si3'), 10),
    (('H_', 'C_3', 'C_3', 'Si3'), 11),
    (('C_3', 'C_3', 'C_3', 'H_'), 12),
    (('H_', 'C_3', 'C_3', 'H_'), 13),
    (('C_3', 'C_3', 'C_3', 'N_3'), 14),
    (('H_', 'C_3', 'C_3', 'N_3'), 15),
    (('C_3', 'C_3', 'N_3', 'H_'), 16),
    (('H_', 'C_3', 'N_3', 'H_'), 17),
    (('C_3', 'C_3', 'C_3', 'C_3'), 18),
    (('C_3', 'C_3', 'C_3', 'O_3'), 19),
    (('H_', 'C_3', 'C_3', 'O_3'), 20),
    (('C_3', 'C_3', 'O_3', 'C_3'), 21),
    (('C_3', 'O_3', 'C_3', 'H_'), 22),
    (('O_3', 'C_3', 'C_3', 'O_3'), 23),
    (('C_3', 'C_3', 'C_3', 'O_3'), 24),
    (('C_3', 'C_3', 'C_3', 'H_'), 25),
    (('H_', 'C_3', 'C_3', 'O_3'), 26),
    (('H_', 'C_3', 'C_3', 'H_'), 27),
    (('C_3', 'C_3', 'O_3', 'C_3'), 28),
    (('C_3', 'O_3', 'C_3', 'H_'), 29),
    (('N_3', 'C_3', 'C_3', 'O_3'), 30),
    (('C_3', 'C_3', 'O_3', 'H_'), 31),
    (('H_', 'C_3', 'O_3', 'H_'), 32),
    (('C_3', 'N_3', 'C_3', 'H_'), 33),
    (('C_3', 'C_3', 'N_3', 'C_3'), 34)
]
DIHEDRAL_MAP = dict(DIHEDRAL_MAP_LIST)

class UniversalTypeMapper:
    """
    Context manager to enforce a Universal Type Map by patching 
    topon.forcefield.dreiding functions at runtime.
    """
    def __init__(self, atom_map, bond_map, angle_map, dihedral_map_list):
        self.atom_map = atom_map
        self.bond_map = bond_map
        self.angle_map = angle_map
        self.dihedral_map_list = dihedral_map_list
        self.dihedral_map = dict(dihedral_map_list) # Lookup convenience
        
        # Originals
        self._orig_assign_atom_types = dreiding.assign_atom_types
        self._orig_extract_bonds = dreiding.extract_bonds
        self._orig_extract_angles = dreiding.extract_angles
        self._orig_extract_dihedrals = dreiding.extract_dihedrals

    def __enter__(self):
        # Patch assign_atom_types
        def patched_assign_atom_types(mol, dreiding_params):
            orig_types_dict, orig_atom_data, orig_dreiding_types = self._orig_assign_atom_types(mol, dreiding_params)
            
            # 1. Update atom_types_dict with Master Map IDs
            new_types_dict = {}
            for type_name in orig_types_dict:
                if type_name in self.atom_map:
                    new_types_dict[type_name] = self.atom_map[type_name]
                else:
                    print(f"WARNING: Unknown atom type {type_name}, keeping original ID!")
                    new_types_dict[type_name] = orig_types_dict[type_name]
            
            # 2. Backfill missing atom types from Universal Map
            for type_name, target_id in self.atom_map.items():
                if target_id not in new_types_dict.values():
                    if type_name not in new_types_dict:
                        new_types_dict[type_name] = target_id

            # 3. Update atom_data with new Type IDs
            new_atom_data = []
            for (idx, old_type_id, charge, x, y, z, element, hyb) in orig_atom_data:
                type_name = orig_dreiding_types[idx]
                new_type_id = new_types_dict.get(type_name, old_type_id)
                new_atom_data.append((idx, new_type_id, charge, x, y, z, element, hyb))
                
            return new_types_dict, new_atom_data, orig_dreiding_types

        # Patch extract_bonds
        def patched_extract_bonds(mol, atom_dreiding_types, dreiding_params):
            bond_types, bond_data = self._orig_extract_bonds(mol, atom_dreiding_types, dreiding_params)
            
            new_bond_types = {}
            new_bond_data = []
            
            # Rebuild bond_types map
            for sig, original_id in bond_types.items():
                t1, t2 = sig[0], sig[1]
                key = tuple(sorted((t1, t2)))
                new_id = self.bond_map.get(key, original_id)
                new_bond_types[sig] = new_id
            
            # Backfill missing bond types
            existing_ids = set(new_bond_types.values())
            for key, target_id in self.bond_map.items():
                if target_id not in existing_ids:
                    # Need parameters to construct signature
                    params = dreiding.find_parameter(key, dreiding_params['bond_params'])
                    if isinstance(params, dict):
                        k = 0.5 * params['k']
                        r0 = params['r0']
                    else:
                        k, r0 = params
                    
                    sig = (key[0], key[1], k, r0)
                    new_bond_types[sig] = target_id
                    existing_ids.add(target_id)
            
            # Rebuild bond_data
            for (bid, type_id, at1, at2) in bond_data:
                final_id = type_id
                for s, oid in bond_types.items():
                   if oid == type_id:
                       final_id = new_bond_types[s]
                       break
                new_bond_data.append((bid, final_id, at1, at2))
                
            return new_bond_types, new_bond_data

        # Patch extract_angles
        def patched_extract_angles(mol, atom_dreiding_types, dreiding_params):
            angle_types, angle_data = self._orig_extract_angles(mol, atom_dreiding_types, dreiding_params)
            
            new_angle_types = {}
            new_angle_data = []
            
            for sig, original_id in angle_types.items():
                t1, t2, t3 = sig[0], sig[1], sig[2]
                outer = sorted((t1, t3))
                key = (outer[0], t2, outer[1])
                new_id = self.angle_map.get(key, original_id)
                new_angle_types[sig] = new_id
                
            existing_ids = set(new_angle_types.values())
            for key, target_id in self.angle_map.items():
                if target_id not in existing_ids:
                    params = dreiding.find_parameter(key, dreiding_params['angle_params'])
                    if isinstance(params, dict):
                        k, theta = params.get('k', 100.0), params.get('theta', 109.5) # Defaults?
                    else:
                        k, theta = params
                    
                    sig = (key[0], key[1], key[2], k, theta)
                    new_angle_types[sig] = target_id
                    existing_ids.add(target_id)
            
            for (aid, type_id, at1, at2, at3) in angle_data:
                final_id = type_id
                for s, oid in angle_types.items():
                    if oid == type_id:
                        final_id = new_angle_types[s]
                        break
                new_angle_data.append((aid, final_id, at1, at2, at3))
                
            return new_angle_types, new_angle_data

        # Patch extract_dihedrals
        def patched_extract_dihedrals(mol, atom_dreiding_types, dreiding_params):
            dihedral_types, dihedral_data = self._orig_extract_dihedrals(mol, atom_dreiding_types, dreiding_params)
            
            new_dihedral_types = {}
            new_dihedral_data = []
            
            for sig, original_id in dihedral_types.items():
                t1, t2, t3, t4 = sig[0], sig[1], sig[2], sig[3]
                fwd = (t1, t2, t3, t4)
                rev = (t4, t3, t2, t1)
                key = min(fwd, rev)
                new_id = self.dihedral_map.get(key, original_id)
                new_dihedral_types[sig] = new_id
                
            existing_ids = set(new_dihedral_types.values())
            
            # Backfill using the LIST to catch duplicates
            for key, target_id in self.dihedral_map_list:
                if target_id not in existing_ids:
                    param_list = dreiding.find_parameter(key, dreiding_params['dihedral_params'])
                    if isinstance(param_list, list):
                        # Add ALL terms for this check
                        for params in param_list:
                            k, n, d = params['v_n'], params['n'], params['d']
                            sig = (key[0], key[1], key[2], key[3], k, n, d)
                            
                            # Check for collision with existing type (duplicate params)
                            if sig in new_dihedral_types:
                                # Perturb K slightly to create unique signature entry
                                # This ensures both IDs appear in the data file header count
                                k_perturbed = k + 0.000001
                                sig = (key[0], key[1], key[2], key[3], k_perturbed, n, d)
                                
                            new_dihedral_types[sig] = target_id
                        existing_ids.add(target_id)
            
            for (did, type_id, at1, at2, at3, at4) in dihedral_data:
                final_id = type_id
                for s, oid in dihedral_types.items():
                    if oid == type_id:
                        final_id = new_dihedral_types[s]
                        break
                new_dihedral_data.append((did, final_id, at1, at2, at3, at4))
            
            return new_dihedral_types, new_dihedral_data

        # Apply patches
        dreiding.assign_atom_types = patched_assign_atom_types
        dreiding.extract_bonds = patched_extract_bonds
        dreiding.extract_angles = patched_extract_angles
        dreiding.extract_dihedrals = patched_extract_dihedrals
        
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Restore originals
        dreiding.assign_atom_types = self._orig_assign_atom_types
        dreiding.extract_bonds = self._orig_extract_bonds
        dreiding.extract_angles = self._orig_extract_angles
        dreiding.extract_dihedrals = self._orig_extract_dihedrals


# ============================================================================
# Defaults
# ============================================================================
DEFAULT_OUTPUT = pkg_dir / "tests" / "output" / "simbox_crosslink"

DEFAULT_N_EPOXY = 50       # bifunctional epoxy crosslinker
DEFAULT_N_AMINO = 25       # bifunctional amine crosslinker
DEFAULT_N_POSS  = 10       # monofunctional POSS chain stopper

DEFAULT_DENSITY = 0.85     # g/cm3
DEFAULT_SEED    = 42


# ============================================================================
# Main Workflow
# ============================================================================
def run_workflow(
    output_dir: Path,
    n_epoxy: int = DEFAULT_N_EPOXY,
    n_amino: int = DEFAULT_N_AMINO,
    n_poss: int = DEFAULT_N_POSS,
    density: float = DEFAULT_DENSITY,
    seed: int = DEFAULT_SEED,
):
    """Full simbox workflow: build molecules, pack box, write LAMMPS files."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print("=" * 60)
    print("SimBox Crosslink Workflow")
    print("=" * 60)

    # ------------------------------------------------------------------
    # 1. Build molecules from the library
    # ------------------------------------------------------------------
    print("\n[1/4] Building molecules...")
    lib = MoleculeLibrary()

    epoxy = lib.epoxy_pdms(n_dms=2)       # ~500 g/mol, bifunctional epoxide
    amino = lib.amino_pdms(n_dms=8)       # ~850 g/mol, bifunctional amine
    poss  = lib.am0270_poss()             # ~1267 g/mol, monofunctional amine

    print(f"  {epoxy}")
    print(f"  {amino}")
    print(f"  {poss}")

    # ------------------------------------------------------------------
    # 2. Create box and add molecules
    # ------------------------------------------------------------------
    print(f"\n[2/4] Creating box (density={density} g/cm3)...")
    box = SimBox(density=density, temperature=300.0, pressure=1.0)
    
    if n_epoxy > 0:
        box.add(epoxy, count=n_epoxy)
    if n_amino > 0:
        box.add(amino, count=n_amino)
    if n_poss > 0:
        box.add(poss,  count=n_poss)

    print(box.summary())

    # ------------------------------------------------------------------
    # 3. Pack molecules into the box
    # ------------------------------------------------------------------
    print(f"\n[3/4] Packing {n_epoxy + n_amino + n_poss} molecules (seed={seed})...")
    box.pack(seed=seed)

    # ------------------------------------------------------------------
    # 4. Write LAMMPS files
    # ------------------------------------------------------------------
    print(f"\n[4/4] Writing LAMMPS files to {output_dir}...")
    
    # --- START UNIVERSAL MAPPING PATCH ---
    print("  [PATCH] Applying Universal Type Map (Rigid Matching)")
    with UniversalTypeMapper(ATOM_MAP, BOND_MAP, ANGLE_MAP, DIHEDRAL_MAP_LIST):
        files = box.write(str(output_dir), forcefield="dreiding")
    # --- END UNIVERSAL MAPPING PATCH ---

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    system = box.system

    print("\n" + "=" * 60)
    print("WORKFLOW COMPLETE")
    print("=" * 60)
    print(f"  Time elapsed:   {elapsed:.1f} s")
    print(f"  Total atoms:    {system.mol.GetNumAtoms()}")
    print(f"  Total molecules: {system.num_molecules}")
    print(f"  Box dimensions: {system.box_lengths[0]:.2f} x "
          f"{system.box_lengths[1]:.2f} x {system.box_lengths[2]:.2f} A")
    print(f"  Reactive sites: {len(system.reactive_sites)}")
    print()
    print("  Output files:")
    for name, path in files.items():
        size = Path(path).stat().st_size
        print(f"    {name:25s} {Path(path).name:30s} ({size:>10,} bytes)")
    print()
    print("  Next steps:")
    print("    cd", output_dir)
    print("    lmp -in 1_minimize.in")
    print("    lmp -in 2_nvt.in")
    print("    lmp -in 3_npt.in")
    print("    # Configure 4b_crosslink.in, then:")
    print("    lmp -in 4b_crosslink.in")
    print("=" * 60)

    # ------------------------------------------------------------------
    # Validation checks
    # ------------------------------------------------------------------
    _validate(box, files, n_epoxy, n_amino, n_poss)

    return files


def _validate(box, files, n_epoxy, n_amino, n_poss):
    """Run basic sanity checks on the generated output."""
    system = box.system
    errors = []

    # Atom count
    expected_atoms = (
        n_epoxy * 71 +   # Epoxy-PDMS-n2 has 71 atoms (with H)
        n_amino * 117 +  # Amino-PDMS-n8 has 117 atoms (with H)
        n_poss  * 207    # AM0270-POSS has 207 atoms (with H)
    )
    actual_atoms = system.mol.GetNumAtoms()
    if actual_atoms != expected_atoms:
        errors.append(
            f"Atom count mismatch: expected {expected_atoms}, got {actual_atoms}"
        )

    # Molecule count
    expected_mols = n_epoxy + n_amino + n_poss
    if system.num_molecules != expected_mols:
        errors.append(
            f"Molecule count mismatch: expected {expected_mols}, "
            f"got {system.num_molecules}"
        )

    # Reactive sites: epoxy has 2 epoxide rings x 3 atoms each = 6 per mol
    # amino has 2 primary amine N = 2 per mol; poss has 1 primary amine N
    expected_amine_atoms = n_amino * 2 + n_poss * 1
    amine_sites = [s for s in system.reactive_sites if s.group_name == "primary_amine"]
    if len(amine_sites) != expected_amine_atoms:
        errors.append(
            f"Primary amine count: expected {expected_amine_atoms}, "
            f"got {len(amine_sites)}"
        )

    # Files exist and are non-empty
    for name, path in files.items():
        p = Path(path)
        if not p.exists():
            errors.append(f"Missing file: {name} -> {path}")
        elif p.stat().st_size == 0:
            errors.append(f"Empty file: {name} -> {path}")

    # Data file has correct header
    data_path = Path(files["data"])
    with open(data_path) as f:
        header = f.read(500)
    if "atoms" not in header:
        errors.append("Data file missing 'atoms' in header")
    if "atom types" not in header:
        errors.append("Data file missing 'atom types' in header")
    if "Pair Coeffs" not in open(data_path).read():
        errors.append("Data file missing Pair Coeffs section")

    # Report
    if errors:
        print("\n*** VALIDATION FAILED ***")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    else:
        print("\n  [VALIDATION] All checks passed.")


# ============================================================================
# CLI
# ============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a simbox crosslink system for LAMMPS"
    )
    parser.add_argument(
        "--output", type=str, default=str(DEFAULT_OUTPUT),
        help=f"Output directory (default: {DEFAULT_OUTPUT})"
    )
    parser.add_argument("--n_epoxy", type=int, default=DEFAULT_N_EPOXY)
    parser.add_argument("--n_amino", type=int, default=DEFAULT_N_AMINO)
    parser.add_argument("--n_poss",  type=int, default=DEFAULT_N_POSS)
    parser.add_argument("--density", type=float, default=DEFAULT_DENSITY)
    parser.add_argument("--seed",    type=int, default=DEFAULT_SEED)

    args = parser.parse_args()

    run_workflow(
        output_dir=args.output,
        n_epoxy=args.n_epoxy,
        n_amino=args.n_amino,
        n_poss=args.n_poss,
        density=args.density,
        seed=args.seed,
    )
