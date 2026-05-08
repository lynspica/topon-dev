"""
topon.workflows.reactive_crosslink
====================================
Canonical workflow: reactive crosslink simulation box (Epoxy-PDMS / Amino-PDMS / POSS).

Four-stage pipeline
-------------------
1. Molecules  — build molecular definitions from MoleculeLibrary
2. Box        — create SimBox, add molecule counts, compute box size from density
3. Packing    — random placement of all molecules (packmol-style, seed-controlled)
4. Output     — write LAMMPS data + input scripts with Universal Type Map enforced

Universal Type Map
------------------
Enforces fixed atom/bond/angle/dihedral type IDs across all simulation variations
(amino-only, POSS-only, mixed) so that pre-defined bond/react templates remain valid.
IDs match the v3 reference (poss_0 composition: 600 epoxy, 300 amino, 0 POSS):
  Si3=1, O_3=2, C_3=3, N_3=4 (or 5 depending on H_ ordering), H_=5 (or 4)

Usage (Python)
--------------
    from topon.workflows.reactive_crosslink import run

    run(
        output_dir="output/crosslink_run",
        n_epoxy=50,
        n_amino=25,
        n_poss=0,
        density=0.85,
        seed=42,
    )

Usage (CLI)
-----------
    python -m topon.workflows.reactive_crosslink \\
        --output output/crosslink_run \\
        --n_epoxy 50 --n_amino 25 --n_poss 0 \\
        --density 0.85 --seed 42
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import topon.forcefield.dreiding as dreiding
from topon.simbox import SimBox, MoleculeLibrary


# ---------------------------------------------------------------------------
# Universal Type Maps
# Derived from v3 poss_0 reference (600 epoxy, 300 amino, 0 POSS, seed=42).
# These enforce consistent type IDs so bond/react templates stay valid.
# ---------------------------------------------------------------------------

ATOM_MAP = {
    "Si3": 1, "O_3": 2, "C_3": 3, "N_3": 4, "H_": 5
}

BOND_MAP = {
    ("O_3", "Si3"): 1,
    ("C_3", "Si3"): 2,
    ("C_3", "C_3"): 3,
    ("C_3", "N_3"): 4,
    ("H_",  "Si3"): 5,
    ("C_3", "H_"):  6,
    ("H_",  "N_3"): 7,
    ("C_3", "O_3"): 8,
    ("H_",  "O_3"): 9,
}

ANGLE_MAP = {
    ("C_3", "Si3", "O_3"): 1,
    ("H_",  "Si3", "O_3"): 2,
    ("C_3", "Si3", "C_3"): 3,
    ("C_3", "Si3", "H_"):  4,
    ("O_3", "Si3", "O_3"): 5,
    ("Si3", "O_3", "Si3"): 6,
    ("H_",  "C_3", "Si3"): 7,
    ("H_",  "C_3", "H_"):  8,
    ("C_3", "C_3", "Si3"): 9,
    ("C_3", "C_3", "H_"):  10,
    ("C_3", "C_3", "C_3"): 11,
    ("C_3", "C_3", "N_3"): 12,
    ("H_",  "C_3", "N_3"): 13,
    ("C_3", "N_3", "H_"):  14,
    ("H_",  "N_3", "H_"):  15,
    ("C_3", "C_3", "O_3"): 16,
    ("H_",  "C_3", "O_3"): 17,
    ("C_3", "O_3", "C_3"): 18,
    ("C_3", "O_3", "H_"):  19,
    ("C_3", "N_3", "C_3"): 20,
}

DIHEDRAL_MAP_LIST = [
    (("O_3", "Si3", "O_3", "Si3"), 1),
    (("H_",  "Si3", "O_3", "Si3"), 2),
    (("O_3", "Si3", "O_3", "Si3"), 3),
    (("H_",  "C_3", "Si3", "O_3"), 4),
    (("C_3", "Si3", "C_3", "H_"),  5),
    (("H_",  "C_3", "Si3", "H_"),  6),
    (("C_3", "C_3", "Si3", "O_3"), 7),
    (("C_3", "C_3", "Si3", "C_3"), 8),
    (("C_3", "C_3", "Si3", "H_"),  9),
    (("C_3", "C_3", "C_3", "Si3"), 10),
    (("H_",  "C_3", "C_3", "Si3"), 11),
    (("C_3", "C_3", "C_3", "H_"),  12),
    (("H_",  "C_3", "C_3", "H_"),  13),
    (("C_3", "C_3", "C_3", "N_3"), 14),
    (("H_",  "C_3", "C_3", "N_3"), 15),
    (("C_3", "C_3", "N_3", "H_"),  16),
    (("H_",  "C_3", "N_3", "H_"),  17),
    (("C_3", "C_3", "C_3", "C_3"), 18),
    (("C_3", "C_3", "C_3", "O_3"), 19),
    (("H_",  "C_3", "C_3", "O_3"), 20),
    (("C_3", "C_3", "O_3", "C_3"), 21),
    (("C_3", "O_3", "C_3", "H_"),  22),
    (("O_3", "C_3", "C_3", "O_3"), 23),
    (("C_3", "C_3", "C_3", "O_3"), 24),
    (("C_3", "C_3", "C_3", "H_"),  25),
    (("H_",  "C_3", "C_3", "O_3"), 26),
    (("H_",  "C_3", "C_3", "H_"),  27),
    (("C_3", "C_3", "O_3", "C_3"), 28),
    (("C_3", "O_3", "C_3", "H_"),  29),
    (("N_3", "C_3", "C_3", "O_3"), 30),
    (("C_3", "C_3", "O_3", "H_"),  31),
    (("H_",  "C_3", "O_3", "H_"),  32),
    (("C_3", "N_3", "C_3", "H_"),  33),
    (("C_3", "C_3", "N_3", "C_3"), 34),
]


# ---------------------------------------------------------------------------
# Universal Type Mapper (context manager)
# ---------------------------------------------------------------------------

class UniversalTypeMapper:
    """
    Context manager: patches topon.forcefield.dreiding at runtime to enforce
    fixed type IDs from the Universal Type Map. Restores originals on exit.
    """

    def __init__(self, atom_map, bond_map, angle_map, dihedral_map_list):
        self.atom_map = atom_map
        self.bond_map = bond_map
        self.angle_map = angle_map
        self.dihedral_map_list = dihedral_map_list
        self.dihedral_map = dict(dihedral_map_list)

        self._orig_assign_atom_types = dreiding.assign_atom_types
        self._orig_extract_bonds = dreiding.extract_bonds
        self._orig_extract_angles = dreiding.extract_angles
        self._orig_extract_dihedrals = dreiding.extract_dihedrals

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

        def patched_assign_atom_types(mol, dreiding_params):
            orig_types_dict, orig_atom_data, orig_dreiding_types = orig_assign(mol, dreiding_params)
            new_types_dict = {}
            for type_name in orig_types_dict:
                new_types_dict[type_name] = atom_map.get(type_name, orig_types_dict[type_name])
            for type_name, target_id in atom_map.items():
                if target_id not in new_types_dict.values() and type_name not in new_types_dict:
                    new_types_dict[type_name] = target_id
            new_atom_data = []
            for (idx, old_type_id, charge, x, y, z, element, hyb) in orig_atom_data:
                type_name = orig_dreiding_types[idx]
                new_atom_data.append((idx, new_types_dict.get(type_name, old_type_id),
                                      charge, x, y, z, element, hyb))
            return new_types_dict, new_atom_data, orig_dreiding_types

        def patched_extract_bonds(mol, atom_dreiding_types, dreiding_params):
            bond_types, bond_data = orig_bonds(mol, atom_dreiding_types, dreiding_params)
            new_bond_types = {}
            for sig, original_id in bond_types.items():
                key = tuple(sorted((sig[0], sig[1])))
                new_bond_types[sig] = bond_map.get(key, original_id)
            existing_ids = set(new_bond_types.values())
            for key, target_id in bond_map.items():
                if target_id not in existing_ids:
                    params = dreiding.find_parameter(key, dreiding_params["bond_params"])
                    k, r0 = (params["k"] * 0.5, params["r0"]) if isinstance(params, dict) else params
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
                outer = sorted((sig[0], sig[2]))
                key = (outer[0], sig[1], outer[1])
                new_angle_types[sig] = angle_map.get(key, original_id)
            existing_ids = set(new_angle_types.values())
            for key, target_id in angle_map.items():
                if target_id not in existing_ids:
                    params = dreiding.find_parameter(key, dreiding_params["angle_params"])
                    k, theta = (params.get("k", 100.0), params.get("theta", 109.5)) if isinstance(params, dict) else params
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
                fwd = (sig[0], sig[1], sig[2], sig[3])
                key = min(fwd, (fwd[3], fwd[2], fwd[1], fwd[0]))
                new_dihedral_types[sig] = dihedral_map.get(key, original_id)
            existing_ids = set(new_dihedral_types.values())
            for key, target_id in dihedral_map_list:
                if target_id not in existing_ids:
                    param_list = dreiding.find_parameter(key, dreiding_params["dihedral_params"])
                    if isinstance(param_list, list):
                        for params in param_list:
                            k, n, d = params["v_n"], params["n"], params["d"]
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

            # The simbox writer uses len(dict) for the "N dihedral types" header.
            # Multi-term DREIDING dihedrals can produce multiple sigs per type ID,
            # making len(dict) > max(values) and causing LAMMPS to expect coeffs
            # for types beyond max(values). Deduplicate: one sig per type ID.
            dedup = {}
            for sig, tid in sorted(new_dihedral_types.items(), key=lambda x: x[1]):
                if tid not in dedup.values():
                    dedup[sig] = tid
            return dedup, new_dihedral_data

        dreiding.assign_atom_types = patched_assign_atom_types
        dreiding.extract_bonds = patched_extract_bonds
        dreiding.extract_angles = patched_extract_angles
        dreiding.extract_dihedrals = patched_extract_dihedrals
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        dreiding.assign_atom_types = self._orig_assign_atom_types
        dreiding.extract_bonds = self._orig_extract_bonds
        dreiding.extract_angles = self._orig_extract_angles
        dreiding.extract_dihedrals = self._orig_extract_dihedrals


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run(
    output_dir: str | Path,
    n_epoxy: int = 50,
    n_amino: int = 25,
    n_poss: int = 0,
    density: float = 0.85,
    seed: int = 42,
) -> Path:
    """
    Run the full reactive crosslink workflow.

    Parameters
    ----------
    output_dir : output directory for LAMMPS files
    n_epoxy    : number of Epoxy-PDMS molecules (bifunctional epoxide)
    n_amino    : number of Amino-PDMS molecules (bifunctional amine)
    n_poss     : number of AM0270-POSS molecules (monofunctional amine)
    density    : target density in g/cm3
    seed       : random seed for packing

    Returns
    -------
    Path to output directory.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print("=" * 60)
    print("[Stage 1] Building molecules...")

    lib = MoleculeLibrary()
    epoxy = lib.epoxy_pdms(n_dms=2)    # bifunctional epoxide, ~71 atoms w/ H
    amino = lib.amino_pdms(n_dms=8)    # bifunctional amine,   ~117 atoms w/ H
    poss  = lib.am0270_poss()          # monofunctional amine, ~207 atoms w/ H

    print(f"  {epoxy}")
    print(f"  {amino}")
    print(f"  {poss}")

    print(f"[Stage 2] Creating box (density={density} g/cm3)...")
    box = SimBox(density=density, temperature=300.0, pressure=1.0)
    if n_epoxy > 0:
        box.add(epoxy, count=n_epoxy)
    if n_amino > 0:
        box.add(amino, count=n_amino)
    if n_poss > 0:
        box.add(poss,  count=n_poss)
    print(box.summary())

    print(f"[Stage 3] Packing {n_epoxy + n_amino + n_poss} molecules (seed={seed})...")
    box.pack(seed=seed)

    print(f"[Stage 4] Writing LAMMPS files to {output_dir}...")
    with UniversalTypeMapper(ATOM_MAP, BOND_MAP, ANGLE_MAP, DIHEDRAL_MAP_LIST):
        files = box.write(str(output_dir), forcefield="dreiding")

    elapsed = time.time() - t0
    system = box.system
    print(f"Reactive crosslink workflow complete -> {output_dir} ({elapsed:.1f}s)")
    print(f"  Atoms: {system.mol.GetNumAtoms()}, Molecules: {system.num_molecules}, "
          f"Reactive sites: {len(system.reactive_sites)}")

    return output_dir


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Reactive crosslink simulation box (Epoxy-PDMS / Amino-PDMS / POSS)"
    )
    p.add_argument("--output",   required=True, help="Output directory")
    p.add_argument("--n_epoxy",  type=int,   default=50)
    p.add_argument("--n_amino",  type=int,   default=25)
    p.add_argument("--n_poss",   type=int,   default=0)
    p.add_argument("--density",  type=float, default=0.85)
    p.add_argument("--seed",     type=int,   default=42)
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    run(
        output_dir=args.output,
        n_epoxy=args.n_epoxy,
        n_amino=args.n_amino,
        n_poss=args.n_poss,
        density=args.density,
        seed=args.seed,
    )
