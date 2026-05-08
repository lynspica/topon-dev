"""
LAMMPS data file writer for simbox.

Produces a complete LAMMPS data file (``atom_style full``) with DREIDING
force-field parameterisation.  Reuses the atom-typing and topology-extraction
functions from ``topon.forcefield.dreiding`` and adds per-molecule IDs, proper
box dimensions, and reactive-group definitions.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import numpy as np

from topon.simbox.system import AssembledSystem


def write_lammps(
    system: AssembledSystem,
    output_dir: str,
    data_filename: str = "system.data",
    groups_filename: str = "groups.txt",
) -> dict[str, str]:
    """Write LAMMPS data file and group definitions.

    Parameters
    ----------
    system : AssembledSystem
        Merged system from :func:`simbox.system.assemble`.
    output_dir : str
        Directory for output files (created if needed).
    data_filename : str
        LAMMPS data file name.
    groups_filename : str
        Reactive-group definitions file name.

    Returns
    -------
    dict[str, str]
        Mapping from logical name to written file path.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    data_path = out / data_filename
    groups_path = out / groups_filename
    settings_path = out / "settings.in"
    ff_coeffs_path = out / "ff_coeffs.in"

    _write_data_file(system, data_path)
    _write_groups_file(system, groups_path)
    _write_settings_file(system, settings_path)
    _write_ff_coeffs_file(system, ff_coeffs_path)

    written = {
        "data": str(data_path),
        "groups": str(groups_path),
        "settings": str(settings_path),
        "ff_coeffs": str(ff_coeffs_path),
    }

    print(f"[Writer] Data file     -> {data_path}")
    print(f"[Writer] Groups file   -> {groups_path}")
    print(f"[Writer] Settings file -> {settings_path}")
    print(f"[Writer] FF coeffs     -> {ff_coeffs_path}")
    return written


# ======================================================================
# Internal: LAMMPS data file
# ======================================================================

def _write_data_file(system: AssembledSystem, path: Path) -> None:
    """Write a complete LAMMPS data file with DREIDING parameterisation."""
    from topon.forcefield.dreiding import (
        assign_atom_types,
        extract_angles,
        extract_bonds,
        extract_dihedrals,
        extract_impropers,
        find_parameter,
        parse_dreiding_parameter_file,
    )

    # Locate the DREIDING parameter file shipped with topon
    dreiding_param_path = _find_dreiding_params()
    dreiding_params = parse_dreiding_parameter_file(dreiding_param_path)

    mol = system.mol

    # -- Atom typing --
    atom_types_dict, atom_data_raw, atom_dreiding_types = assign_atom_types(
        mol, dreiding_params
    )

    # -- Topology extraction --
    bond_types, bond_data = extract_bonds(mol, atom_dreiding_types, dreiding_params)
    angle_types, angle_data = extract_angles(mol, atom_dreiding_types, dreiding_params)
    dihedral_types, dihedral_data = extract_dihedrals(
        mol, atom_dreiding_types, dreiding_params
    )
    improper_types, improper_data = extract_impropers(
        mol, atom_dreiding_types, dreiding_params
    )

    # -- Coordinates from conformer --
    conf = mol.GetConformer()

    # -- Box --
    lx, ly, lz = system.box_lengths
    xlo, ylo, zlo = 0.0, 0.0, 0.0
    xhi, yhi, zhi = lx, ly, lz

    num_atoms = mol.GetNumAtoms()
    num_bonds = len(bond_data)
    num_angles = len(angle_data)
    num_dihedrals = len(dihedral_data)
    num_impropers = len(improper_data)

    with open(path, "w") as f:
        # ---- Header ----
        f.write("LAMMPS data file - simbox (topon) DREIDING parameterization\n\n")
        f.write(f"{num_atoms} atoms\n")
        f.write(f"{num_bonds} bonds\n")
        f.write(f"{num_angles} angles\n")
        f.write(f"{num_dihedrals} dihedrals\n")
        f.write(f"{num_impropers} impropers\n\n")

        n_atom_types = max(atom_types_dict.values()) if atom_types_dict else 0
        n_bond_types = max(bond_types.values()) if bond_types else 0
        n_angle_types = max(angle_types.values()) if angle_types else 0
        n_dihedral_types = max(dihedral_types.values()) if dihedral_types else 0
        n_improper_types = max(improper_types.values()) if improper_types else 0
        f.write(f"{n_atom_types} atom types\n")
        f.write(f"{n_bond_types} bond types\n")
        f.write(f"{n_angle_types} angle types\n")
        f.write(f"{n_dihedral_types} dihedral types\n")
        f.write(f"{n_improper_types} improper types\n\n")

        f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
        f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
        f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n")
        f.write("0.0 0.0 0.0 xy xz yz\n\n")

        # ---- Masses ----
        f.write("Masses\n\n")
        for type_name, type_id in sorted(atom_types_dict.items(), key=lambda x: x[1]):
            mass = _get_mass(type_name, dreiding_params)
            f.write(f"{type_id} {mass:.4f}  # {type_name}\n")
        f.write("\n")

        # ---- Pair Coeffs ----
        f.write("Pair Coeffs\n\n")
        for type_name, type_id in sorted(atom_types_dict.items(), key=lambda x: x[1]):
            eps, sig = _get_pair_params(type_name, dreiding_params)
            f.write(f"{type_id} {eps} {sig}  # {type_name}\n")
        f.write("\n")

        # ---- Bond Coeffs ----
        if bond_types:
            f.write("Bond Coeffs\n\n")
            for (t1, t2, k, r0), tid in sorted(bond_types.items(), key=lambda x: x[1]):
                f.write(f"{tid} {k} {r0}  # {t1}-{t2}\n")
            f.write("\n")

        # ---- Angle Coeffs ----
        if angle_types:
            f.write("Angle Coeffs\n\n")
            for (t1, t2, t3, k, th), tid in sorted(
                angle_types.items(), key=lambda x: x[1]
            ):
                f.write(f"{tid} {k} {th}  # {t1}-{t2}-{t3}\n")
            f.write("\n")

        # ---- Dihedral Coeffs ----
        if dihedral_types:
            f.write("Dihedral Coeffs\n\n")
            for key, tid in sorted(dihedral_types.items(), key=lambda x: x[1]):
                t1, t2, t3, t4, k, n, d = key
                f.write(
                    f"{tid} {k:.6f} {int(d)} {int(n)}"
                    f"  # {t1}-{t2}-{t3}-{t4}\n"
                )
            f.write("\n")

        # ---- Improper Coeffs ----
        if improper_types:
            f.write("Improper Coeffs\n\n")
            for key, tid in sorted(improper_types.items(), key=lambda x: x[1]):
                t1, t2, t3, t4, itype, k, chi0 = key
                # cvff: E = K[1 + d*cos(n*chi)]; DREIDING: d=-1, n=1 → E = K(1-cos chi)
                f.write(f"{tid} {k:.6f} -1 1  # {t1} {t2} {t3} {t4}\n")
            f.write("\n")

        # ---- Atoms (atom_style full) ----
        f.write("Atoms # full\n\n")
        for raw_entry in atom_data_raw:
            idx, type_id, charge, _x, _y, _z, element, hyb = raw_entry
            # idx is 1-based
            mol_id = system.molecule_ids[idx - 1]

            pos = conf.GetAtomPosition(idx - 1)
            x, y, z = pos.x, pos.y, pos.z

            f.write(
                f"{idx} {mol_id} {type_id} {charge:.6f} "
                f"{x:.6f} {y:.6f} {z:.6f}  # {element} ({hyb})\n"
            )
        f.write("\n")

        # ---- Bonds ----
        if bond_data:
            f.write("Bonds\n\n")
            for bond_id, type_id, a1, a2 in bond_data:
                f.write(f"{bond_id} {type_id} {a1} {a2}\n")
            f.write("\n")

        # ---- Angles ----
        if angle_data:
            f.write("Angles\n\n")
            for angle_id, type_id, a1, a2, a3 in angle_data:
                f.write(f"{angle_id} {type_id} {a1} {a2} {a3}\n")
            f.write("\n")

        # ---- Dihedrals ----
        if dihedral_data:
            f.write("Dihedrals\n\n")
            for did, type_id, a1, a2, a3, a4 in dihedral_data:
                f.write(f"{did} {type_id} {a1} {a2} {a3} {a4}\n")
            f.write("\n")

        # ---- Impropers ----
        if improper_data:
            f.write("Impropers\n\n")
            for iid, type_id, a1, a2, a3, a4 in improper_data:
                f.write(f"{iid} {type_id} {a1} {a2} {a3} {a4}\n")

    print(
        f"[Writer] Wrote {num_atoms} atoms, {num_bonds} bonds, "
        f"{num_angles} angles, {num_dihedrals} dihedrals, "
        f"{num_impropers} impropers"
    )


# ======================================================================
# Internal: groups file
# ======================================================================

def _write_groups_file(system: AssembledSystem, path: Path) -> None:
    """Write LAMMPS group definitions for reactive atoms.

    Produces commands like::

        group epoxide    id 15 18 234 237 ...
        group primary_amine id 45 312 ...
    """
    from collections import defaultdict

    groups: dict[str, list[int]] = defaultdict(list)
    for entry in system.reactive_sites:
        # Convert to 1-based LAMMPS atom ID
        groups[entry.group_name].append(entry.global_atom_idx + 1)

    with open(path, "w") as f:
        f.write("# Reactive-group definitions generated by topon.simbox\n")
        f.write("# Include in your LAMMPS input: include groups.txt\n\n")

        for group_name, atom_ids in sorted(groups.items()):
            ids_str = " ".join(str(i) for i in sorted(atom_ids))
            f.write(f"group {group_name} id {ids_str}\n")

        # Also write per-species molecule groups
        f.write("\n# Per-species molecule groups\n")
        species_mol_ids: dict[str, list[int]] = defaultdict(list)
        for mol_id_0, name in enumerate(system.species_names):
            species_mol_ids[name].append(mol_id_0 + 1)

        for sp_name, mol_ids in sorted(species_mol_ids.items()):
            safe_name = sp_name.replace(" ", "_").replace("-", "_").lower()
            ids_str = " ".join(str(i) for i in sorted(mol_ids))
            f.write(f"group {safe_name} molecule {ids_str}\n")

    print(
        f"[Writer] Groups: "
        + ", ".join(f"{k}({len(v)} atoms)" for k, v in groups.items())
    )


# ======================================================================
# Internal: settings file (pair_coeff for re-application after soft)
# ======================================================================

def _write_settings_file(system: AssembledSystem, path: Path) -> None:
    """Write LAMMPS-style pair_coeff commands to a separate include file.

    This file is needed to re-apply DREIDING pair coefficients after
    switching from ``pair_style soft`` back to ``pair_style lj/cut``.
    """
    from topon.forcefield.dreiding import (
        assign_atom_types,
        parse_dreiding_parameter_file,
    )

    dreiding_param_path = _find_dreiding_params()
    dreiding_params = parse_dreiding_parameter_file(dreiding_param_path)

    atom_types_dict, _, _ = assign_atom_types(system.mol, dreiding_params)

    with open(path, "w") as f:
        f.write("# DREIDING pair coefficients - generated by topon.simbox\n")
        f.write("# Include after (re-)declaring pair_style lj/cut 12.0\n\n")

        for type_name, type_id in sorted(atom_types_dict.items(), key=lambda x: x[1]):
            eps, sig = _get_pair_params(type_name, dreiding_params)
            f.write(f"pair_coeff {type_id} {type_id} {eps} {sig}  # {type_name}\n")

        f.write("\npair_modify mix arithmetic\n")


# ======================================================================
# Internal: ff_coeffs file (all coeff commands for LAMMPS include)
# ======================================================================

def _write_ff_coeffs_file(system: AssembledSystem, path: Path) -> None:
    """Write all force-field coefficients as LAMMPS *_coeff commands.

    Produces a self-contained ``ff_coeffs.in`` that can be included in any
    LAMMPS input script after read_data to set all pair/bond/angle/dihedral
    coefficients.  Equivalent to the Coeffs sections in system.data but in
    command form, which is required for fix/bond/react reaction templates.
    """
    from topon.forcefield.dreiding import (
        assign_atom_types,
        extract_angles,
        extract_bonds,
        extract_dihedrals,
        parse_dreiding_parameter_file,
    )

    dreiding_param_path = _find_dreiding_params()
    dreiding_params = parse_dreiding_parameter_file(dreiding_param_path)

    mol = system.mol
    atom_types_dict, _, _ = assign_atom_types(mol, dreiding_params)
    atom_dreiding_types = assign_atom_types(mol, dreiding_params)[2]
    bond_types, _ = extract_bonds(mol, atom_dreiding_types, dreiding_params)
    angle_types, _ = extract_angles(mol, atom_dreiding_types, dreiding_params)
    dihedral_types, _ = extract_dihedrals(mol, atom_dreiding_types, dreiding_params)

    with open(path, "w") as f:
        f.write("# Force field coefficients (DREIDING)\n")
        f.write("# Auto-generated by topon.simbox\n")
        f.write("# Include after read_data + pair_style/bond_style declarations\n\n")

        f.write("# --- Pair Coeffs (LJ) ---\n")
        for type_name, type_id in sorted(atom_types_dict.items(), key=lambda x: x[1]):
            eps, sig = _get_pair_params(type_name, dreiding_params)
            f.write(f"pair_coeff {type_id} {type_id} {eps} {sig}  # {type_name}\n")
        f.write("pair_modify mix arithmetic\n\n")

        if bond_types:
            f.write("# --- Bond Coeffs ---\n")
            for (t1, t2, k, r0), tid in sorted(bond_types.items(), key=lambda x: x[1]):
                f.write(f"bond_coeff {tid} {k} {r0}  # {t1}-{t2}\n")
            f.write("\n")

        if angle_types:
            f.write("# --- Angle Coeffs ---\n")
            for (t1, t2, t3, k, th), tid in sorted(
                angle_types.items(), key=lambda x: x[1]
            ):
                f.write(f"angle_coeff {tid} {k} {th}  # {t1}-{t2}-{t3}\n")
            f.write("\n")

        if dihedral_types:
            f.write("# --- Dihedral Coeffs ---\n")
            for key, tid in sorted(dihedral_types.items(), key=lambda x: x[1]):
                t1, t2, t3, t4, kd, n, d = key
                f.write(
                    f"dihedral_coeff {tid} {kd:.6f} {int(d)} {int(n)}"
                    f"  # {t1}-{t2}-{t3}-{t4}\n"
                )
            f.write("\n")


# ======================================================================
# Helpers
# ======================================================================

def _find_dreiding_params() -> str:
    """Locate ``DreidingX6parameters.txt`` shipped with topon."""
    # It lives next to topon/forcefield/dreiding.py
    ff_dir = Path(__file__).resolve().parent.parent / "forcefield"
    param_file = ff_dir / "DreidingX6parameters.txt"
    if param_file.exists():
        return str(param_file)
    # Fallback: look in utils (older location)
    utils_dir = Path(__file__).resolve().parent.parent / "utils"
    param_file2 = utils_dir / "DreidingX6parameters.txt"
    if param_file2.exists():
        return str(param_file2)
    raise FileNotFoundError(
        "Cannot find DreidingX6parameters.txt in topon/forcefield/ or topon/utils/"
    )


def _get_mass(type_name: str, dreiding_params: dict) -> float:
    """Get atomic mass for a DREIDING type."""
    if type_name in dreiding_params["atom_types"]:
        return dreiding_params["atom_types"][type_name]["mass"]
    # Fallback by element symbol
    element = "".join(c for c in type_name if c.isalpha() and c.isupper())
    return {
        "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
        "F": 18.998, "S": 32.065, "Si": 28.086, "P": 30.974,
    }.get(element, 1.0)


def _get_pair_params(type_name: str, dreiding_params: dict) -> tuple[float, float]:
    """Return (epsilon, sigma) for a DREIDING atom type."""
    if type_name in dreiding_params["vdw_params"]:
        p = dreiding_params["vdw_params"][type_name]
        if isinstance(p, dict):
            return (p.get("epsilon", 0.001), p.get("radius", 3.5))
    return (0.001, 3.5)
