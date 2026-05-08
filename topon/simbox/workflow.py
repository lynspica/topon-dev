"""
SimBox crosslink workflow — packaged entry point.

Contains the core ``run_workflow`` function and the ``UniversalTypeMapper``
context manager that enforces consistent DREIDING type IDs across all
simulation compositions (Amino-only, POSS-only, mixed).

This module is the canonical implementation.  The script at
``tests/workflows/generate_simbox_crosslink.py`` delegates to this.
"""

from __future__ import annotations

import time
from pathlib import Path

# dreiding is imported lazily inside UniversalTypeMapper.__enter__ and
# run_workflow to avoid a module-level rdkit import at CLI startup.
from topon.simbox import SimBox, MoleculeLibrary


# ---------------------------------------------------------------------------
# Universal type maps — keep IDs stable across all compositions so that
# pre-defined LAMMPS bond/react templates stay compatible.
# ---------------------------------------------------------------------------

ATOM_MAP = {
    'Si3': 1, 'O_3': 2, 'C_3': 3, 'N_3': 4, 'H_': 5
}

BOND_MAP = {
    ('O_3', 'Si3'): 1,
    ('C_3', 'Si3'): 2,
    ('C_3', 'C_3'): 3,
    ('C_3', 'N_3'): 4,
    ('H_', 'Si3'): 5,
    ('C_3', 'H_'): 6,
    ('H_', 'N_3'): 7,
    ('C_3', 'O_3'): 8,
    ('H_', 'O_3'): 9,
}

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
    ('C_3', 'N_3', 'C_3'): 20,
}

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
    (('C_3', 'C_3', 'N_3', 'C_3'), 34),
]


# ---------------------------------------------------------------------------
# UniversalTypeMapper
# ---------------------------------------------------------------------------

class UniversalTypeMapper:
    """
    Context manager that patches ``topon.forcefield.dreiding`` functions to
    enforce the universal DREIDING type-ID maps at write time.

    This guarantees that atom/bond/angle/dihedral IDs are identical across
    all compositions (Amino-only, POSS-only, mixed), so pre-defined LAMMPS
    bond/react templates remain compatible.
    """

    def __init__(self, atom_map=None, bond_map=None, angle_map=None, dihedral_map_list=None):
        self.atom_map = atom_map or ATOM_MAP
        self.bond_map = bond_map or BOND_MAP
        self.angle_map = angle_map or ANGLE_MAP
        self.dihedral_map_list = dihedral_map_list or DIHEDRAL_MAP_LIST
        self.dihedral_map = dict(self.dihedral_map_list)

        import topon.forcefield.dreiding as _dreiding
        self._dreiding = _dreiding
        self._orig_assign_atom_types = _dreiding.assign_atom_types
        self._orig_extract_bonds = _dreiding.extract_bonds
        self._orig_extract_angles = _dreiding.extract_angles
        self._orig_extract_dihedrals = _dreiding.extract_dihedrals

    def __enter__(self):
        atom_map = self.atom_map
        bond_map = self.bond_map
        angle_map = self.angle_map
        dihedral_map = self.dihedral_map
        dihedral_map_list = self.dihedral_map_list
        orig_assign = self._orig_assign_atom_types
        orig_bonds = self._orig_extract_bonds
        orig_angles = self._orig_extract_angles
        orig_dihedrals = self._orig_extract_dihedrals
        _dreiding = self._dreiding

        def patched_assign_atom_types(mol, dreiding_params):
            orig_types_dict, orig_atom_data, orig_dreiding_types = orig_assign(mol, dreiding_params)
            new_types_dict = {}
            for type_name in orig_types_dict:
                if type_name in atom_map:
                    new_types_dict[type_name] = atom_map[type_name]
                else:
                    print(f"WARNING: Unknown atom type {type_name}, keeping original ID!")
                    new_types_dict[type_name] = orig_types_dict[type_name]
            for type_name, target_id in atom_map.items():
                if target_id not in new_types_dict.values():
                    if type_name not in new_types_dict:
                        new_types_dict[type_name] = target_id
            new_atom_data = []
            for (idx, old_type_id, charge, x, y, z, element, hyb) in orig_atom_data:
                type_name = orig_dreiding_types[idx]
                new_type_id = new_types_dict.get(type_name, old_type_id)
                new_atom_data.append((idx, new_type_id, charge, x, y, z, element, hyb))
            return new_types_dict, new_atom_data, orig_dreiding_types

        def patched_extract_bonds(mol, atom_dreiding_types, dreiding_params):
            bond_types, bond_data = orig_bonds(mol, atom_dreiding_types, dreiding_params)
            new_bond_types = {}
            for sig, original_id in bond_types.items():
                t1, t2 = sig[0], sig[1]
                key = tuple(sorted((t1, t2)))
                new_bond_types[sig] = bond_map.get(key, original_id)
            existing_ids = set(new_bond_types.values())
            for key, target_id in bond_map.items():
                if target_id not in existing_ids:
                    params = _dreiding.find_parameter(key, dreiding_params['bond_params'])
                    if isinstance(params, dict):
                        k, r0 = 0.5 * params['k'], params['r0']
                    else:
                        k, r0 = params
                    new_bond_types[(key[0], key[1], k, r0)] = target_id
                    existing_ids.add(target_id)
            new_bond_data = []
            for (bid, type_id, at1, at2) in bond_data:
                final_id = type_id
                for s, oid in bond_types.items():
                    if oid == type_id:
                        final_id = new_bond_types[s]
                        break
                new_bond_data.append((bid, final_id, at1, at2))
            return new_bond_types, new_bond_data

        def patched_extract_angles(mol, atom_dreiding_types, dreiding_params):
            angle_types, angle_data = orig_angles(mol, atom_dreiding_types, dreiding_params)
            new_angle_types = {}
            for sig, original_id in angle_types.items():
                t1, t2, t3 = sig[0], sig[1], sig[2]
                outer = sorted((t1, t3))
                key = (outer[0], t2, outer[1])
                new_angle_types[sig] = angle_map.get(key, original_id)
            existing_ids = set(new_angle_types.values())
            for key, target_id in angle_map.items():
                if target_id not in existing_ids:
                    params = _dreiding.find_parameter(key, dreiding_params['angle_params'])
                    if isinstance(params, dict):
                        k = params.get('k', 100.0)
                        theta = params.get('theta0', params.get('theta', 109.5))
                    else:
                        k, theta = params
                    new_angle_types[(key[0], key[1], key[2], k, theta)] = target_id
                    existing_ids.add(target_id)
            new_angle_data = []
            for (aid, type_id, at1, at2, at3) in angle_data:
                final_id = type_id
                for s, oid in angle_types.items():
                    if oid == type_id:
                        final_id = new_angle_types[s]
                        break
                new_angle_data.append((aid, final_id, at1, at2, at3))
            return new_angle_types, new_angle_data

        def patched_extract_dihedrals(mol, atom_dreiding_types, dreiding_params):
            dihedral_types, dihedral_data = orig_dihedrals(mol, atom_dreiding_types, dreiding_params)
            new_dihedral_types = {}
            for sig, original_id in dihedral_types.items():
                t1, t2, t3, t4 = sig[0], sig[1], sig[2], sig[3]
                fwd, rev = (t1, t2, t3, t4), (t4, t3, t2, t1)
                key = min(fwd, rev)
                new_dihedral_types[sig] = dihedral_map.get(key, original_id)
            existing_ids = set(new_dihedral_types.values())
            for key, target_id in dihedral_map_list:
                if target_id not in existing_ids:
                    param_list = _dreiding.find_parameter(key, dreiding_params['dihedral_params'])
                    if isinstance(param_list, list):
                        for params in param_list:
                            k, n, d = params['v_n'], params['n'], params['d']
                            sig = (key[0], key[1], key[2], key[3], k, n, d)
                            if sig in new_dihedral_types:
                                sig = (key[0], key[1], key[2], key[3], k + 1e-6, n, d)
                            new_dihedral_types[sig] = target_id
                        existing_ids.add(target_id)
            new_dihedral_data = []
            for (did, type_id, at1, at2, at3, at4) in dihedral_data:
                final_id = type_id
                for s, oid in dihedral_types.items():
                    if oid == type_id:
                        final_id = new_dihedral_types[s]
                        break
                new_dihedral_data.append((did, final_id, at1, at2, at3, at4))
            return new_dihedral_types, new_dihedral_data

        _dreiding.assign_atom_types = patched_assign_atom_types
        _dreiding.extract_bonds = patched_extract_bonds
        _dreiding.extract_angles = patched_extract_angles
        _dreiding.extract_dihedrals = patched_extract_dihedrals
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._dreiding.assign_atom_types = self._orig_assign_atom_types
        self._dreiding.extract_bonds = self._orig_extract_bonds
        self._dreiding.extract_angles = self._orig_extract_angles
        self._dreiding.extract_dihedrals = self._orig_extract_dihedrals


# ---------------------------------------------------------------------------
# run_workflow
# ---------------------------------------------------------------------------

def run_workflow(
    output_dir: str | Path,
    n_epoxy: int = 50,
    n_amino: int = 25,
    n_poss: int = 10,
    density: float = 0.85,
    seed: int = 42,
    verbose: bool = True,
) -> dict:
    """
    Full SimBox crosslink workflow: build molecules, pack box, write LAMMPS files.

    Parameters
    ----------
    output_dir : path
        Directory where all output files will be written.
    n_epoxy : int
        Number of Epoxy-PDMS molecules.
    n_amino : int
        Number of Amino-PDMS molecules.
    n_poss : int
        Number of AM0270-POSS molecules.
    density : float
        Target packing density in g/cm³.
    seed : int
        Random seed for reproducible packing.
    verbose : bool
        Print progress and summary.

    Returns
    -------
    dict
        Mapping of file labels to absolute path strings (keys: ``data``,
        ``settings``, ``groups``, ``ff_coeffs``, ``minimize``, ``nvt``,
        ``npt``, ``crosslink``).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    if verbose:
        print("=" * 60)
        print("SimBox Crosslink Workflow")
        print("=" * 60)

    lib = MoleculeLibrary()
    epoxy = lib.epoxy_pdms(n_dms=2)
    amino = lib.amino_pdms(n_dms=8)
    poss = lib.am0270_poss()

    if verbose:
        print(f"\n[1/4] Molecules: {epoxy}  |  {amino}  |  {poss}")

    box = SimBox(density=density, temperature=300.0, pressure=1.0)
    if n_epoxy > 0:
        box.add(epoxy, count=n_epoxy)
    if n_amino > 0:
        box.add(amino, count=n_amino)
    if n_poss > 0:
        box.add(poss, count=n_poss)

    if verbose:
        print(f"[2/4] {box.summary()}")
        print(f"[3/4] Packing {n_epoxy + n_amino + n_poss} molecules (seed={seed})...")

    box.pack(seed=seed)

    if verbose:
        print(f"[4/4] Writing LAMMPS files to {output_dir}...")

    with UniversalTypeMapper():
        files = box.write(str(output_dir), forcefield="dreiding")

    if verbose:
        elapsed = time.time() - t0
        system = box.system
        print(f"\n{'=' * 60}")
        print("WORKFLOW COMPLETE")
        print(f"  Time:      {elapsed:.1f} s")
        print(f"  Atoms:     {system.mol.GetNumAtoms()}")
        print(f"  Molecules: {system.num_molecules}")
        bl = system.box_lengths
        print(f"  Box:       {bl[0]:.2f} x {bl[1]:.2f} x {bl[2]:.2f} Å")
        print(f"  Next: cd {output_dir} && lmp -in 1_minimize.in")
        print("=" * 60)

    return files
