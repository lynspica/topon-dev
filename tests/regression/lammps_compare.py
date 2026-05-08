"""
tests/regression/lammps_compare.py
===================================
Utility for comparing two LAMMPS data files section by section.

Used by all regression tests to verify that refactored code produces
output identical to the frozen reference files in tests/reference/.

Usage
-----
    diffs = compare_lammps_data(new_path, ref_path)
    assert diffs == [], "\\n".join(diffs)
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import NamedTuple


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

class LammpsData(NamedTuple):
    header: dict          # {'atoms': 1000, 'bonds': 990, ...}
    masses: dict          # {type_id: (mass, comment)}
    pair_coeffs: dict     # {type_id: [values...]}
    bond_coeffs: dict
    angle_coeffs: dict
    dihedral_coeffs: dict
    improper_coeffs: dict
    atoms: dict           # {atom_id: (mol, type, charge, x, y, z)}
    bonds: dict           # {bond_id: (type, a1, a2)}
    angles: dict          # {angle_id: (type, a1, a2, a3)}
    dihedrals: dict       # {dihedral_id: (type, a1, a2, a3, a4)}
    impropers: dict       # {improper_id: (type, a1, a2, a3, a4)}


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_lammps_data(path: str | Path) -> LammpsData:
    """Parse a LAMMPS data file into structured sections."""
    path = Path(path)
    text = path.read_text()
    lines = text.splitlines()

    header = {}
    masses = {}
    pair_coeffs = {}
    bond_coeffs = {}
    angle_coeffs = {}
    dihedral_coeffs = {}
    improper_coeffs = {}
    atoms = {}
    bonds = {}
    angles = {}
    dihedrals = {}
    impropers = {}

    # ---- Header ----
    header_keys = [
        "atoms", "bonds", "angles", "dihedrals", "impropers",
        "atom types", "bond types", "angle types",
        "dihedral types", "improper types",
    ]
    for line in lines:
        stripped = line.strip()
        for key in header_keys:
            m = re.match(rf"^(\d+)\s+{re.escape(key)}\s*$", stripped)
            if m:
                header[key] = int(m.group(1))

    # ---- Section reader ----
    current_section = None
    for line in lines:
        stripped = line.strip()

        # Blank lines don't change section
        if not stripped:
            continue

        # Section headers (must come before data lines)
        sec = stripped.lower()
        if sec == "masses":
            current_section = "masses"; continue
        elif sec in ("pair coeffs", "pair coeffs # lj/cut", "pair coeffs # lj/cut/coul/long"):
            current_section = "pair_coeffs"; continue
        elif sec == "bond coeffs":
            current_section = "bond_coeffs"; continue
        elif sec == "angle coeffs":
            current_section = "angle_coeffs"; continue
        elif sec == "dihedral coeffs":
            current_section = "dihedral_coeffs"; continue
        elif sec == "improper coeffs":
            current_section = "improper_coeffs"; continue
        elif sec.startswith("atoms"):
            current_section = "atoms"; continue
        elif sec == "bonds":
            current_section = "bonds"; continue
        elif sec == "angles":
            current_section = "angles"; continue
        elif sec == "dihedrals":
            current_section = "dihedrals"; continue
        elif sec == "impropers":
            current_section = "impropers"; continue
        elif sec == "velocities":
            current_section = "velocities"; continue

        # Skip comment-only lines
        if stripped.startswith("#"):
            continue

        # Strip inline comment
        data_part = stripped.split("#")[0].strip()
        if not data_part:
            continue

        parts = data_part.split()

        if current_section == "masses":
            # id mass [# comment]
            if len(parts) >= 2:
                tid = int(parts[0])
                mass = float(parts[1])
                comment = stripped.split("#")[1].strip() if "#" in stripped else ""
                masses[tid] = (mass, comment)

        elif current_section == "pair_coeffs":
            # id eps sig [cutoff]
            if len(parts) >= 2:
                tid = int(parts[0])
                pair_coeffs[tid] = [float(x) for x in parts[1:]]

        elif current_section == "bond_coeffs":
            # id k r0 [...]
            if len(parts) >= 2:
                tid = int(parts[0])
                bond_coeffs[tid] = [float(x) for x in parts[1:]]

        elif current_section == "angle_coeffs":
            # id k theta0 [...]
            if len(parts) >= 2:
                tid = int(parts[0])
                angle_coeffs[tid] = [float(x) for x in parts[1:]]

        elif current_section == "dihedral_coeffs":
            # id k n d [...]
            if len(parts) >= 2:
                tid = int(parts[0])
                dihedral_coeffs[tid] = [float(x) for x in parts[1:]]

        elif current_section == "improper_coeffs":
            if len(parts) >= 2:
                tid = int(parts[0])
                improper_coeffs[tid] = [float(x) for x in parts[1:]]

        elif current_section == "atoms":
            # atom_style full: id mol type charge x y z
            if len(parts) >= 7:
                aid = int(parts[0])
                mol = int(parts[1])
                atype = int(parts[2])
                charge = float(parts[3])
                x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                atoms[aid] = (mol, atype, charge, x, y, z)

        elif current_section == "bonds":
            if len(parts) >= 4:
                bid = int(parts[0])
                bonds[bid] = (int(parts[1]), int(parts[2]), int(parts[3]))

        elif current_section == "angles":
            if len(parts) >= 5:
                aid = int(parts[0])
                angles[aid] = (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))

        elif current_section == "dihedrals":
            if len(parts) >= 6:
                did = int(parts[0])
                dihedrals[did] = (int(parts[1]), int(parts[2]), int(parts[3]),
                                   int(parts[4]), int(parts[5]))

        elif current_section == "impropers":
            if len(parts) >= 6:
                iid = int(parts[0])
                impropers[iid] = (int(parts[1]), int(parts[2]), int(parts[3]),
                                   int(parts[4]), int(parts[5]))

    return LammpsData(
        header=header,
        masses=masses,
        pair_coeffs=pair_coeffs,
        bond_coeffs=bond_coeffs,
        angle_coeffs=angle_coeffs,
        dihedral_coeffs=dihedral_coeffs,
        improper_coeffs=improper_coeffs,
        atoms=atoms,
        bonds=bonds,
        angles=angles,
        dihedrals=dihedrals,
        impropers=impropers,
    )


# ---------------------------------------------------------------------------
# Comparator
# ---------------------------------------------------------------------------

def compare_lammps_data(
    new_path: str | Path,
    ref_path: str | Path,
    coord_tol: float = 1e-4,
    charge_tol: float = 1e-8,
    coeff_tol: float = 1e-6,
) -> list[str]:
    """
    Compare two LAMMPS data files. Returns a list of discrepancy strings.
    Empty list means the files are equivalent within tolerances.

    Comparison order (most critical first):
    1. Header counts   — exact
    2. Atom types      — exact (via masses section)
    3. Coeff sections  — within coeff_tol
    4. Bond/angle/dihedral connectivity — exact
    5. Atom type IDs and charges — exact / charge_tol
    6. Coordinates     — within coord_tol
    """
    new = parse_lammps_data(new_path)
    ref = parse_lammps_data(ref_path)
    diffs = []

    # 1. Header counts
    for key in ref.header:
        if new.header.get(key) != ref.header[key]:
            diffs.append(
                f"Header '{key}': new={new.header.get(key)}, ref={ref.header[key]}"
            )

    # 2. Masses / atom type count
    if set(new.masses.keys()) != set(ref.masses.keys()):
        diffs.append(
            f"Atom type IDs differ: new={sorted(new.masses)}, ref={sorted(ref.masses)}"
        )
    else:
        for tid in ref.masses:
            nm, _ = new.masses[tid]
            rm, _ = ref.masses[tid]
            if abs(nm - rm) > coeff_tol:
                diffs.append(f"Mass type {tid}: new={nm}, ref={rm}")

    # 3. Coefficient sections
    def _compare_coeffs(new_c, ref_c, section):
        if set(new_c.keys()) != set(ref_c.keys()):
            diffs.append(
                f"{section} type IDs differ: new={sorted(new_c)}, ref={sorted(ref_c)}"
            )
            return
        for tid in ref_c:
            nv, rv = new_c[tid], ref_c[tid]
            if len(nv) != len(rv):
                diffs.append(f"{section} type {tid} param count: new={len(nv)}, ref={len(rv)}")
                continue
            for i, (n, r) in enumerate(zip(nv, rv)):
                if abs(n - r) > coeff_tol:
                    diffs.append(f"{section} type {tid} param[{i}]: new={n}, ref={r}")

    _compare_coeffs(new.pair_coeffs, ref.pair_coeffs, "Pair Coeffs")
    _compare_coeffs(new.bond_coeffs, ref.bond_coeffs, "Bond Coeffs")
    _compare_coeffs(new.angle_coeffs, ref.angle_coeffs, "Angle Coeffs")
    _compare_coeffs(new.dihedral_coeffs, ref.dihedral_coeffs, "Dihedral Coeffs")
    _compare_coeffs(new.improper_coeffs, ref.improper_coeffs, "Improper Coeffs")

    # 4. Connectivity
    if new.bonds != ref.bonds:
        diffs.append(f"Bonds section differs ({len(new.bonds)} vs {len(ref.bonds)} entries)")
    if new.angles != ref.angles:
        diffs.append(f"Angles section differs ({len(new.angles)} vs {len(ref.angles)} entries)")
    if new.dihedrals != ref.dihedrals:
        diffs.append(f"Dihedrals section differs")
    if new.impropers != ref.impropers:
        diffs.append(f"Impropers section differs")

    # 5. Atom types and charges
    if set(new.atoms.keys()) != set(ref.atoms.keys()):
        diffs.append(f"Atom ID set differs: {len(new.atoms)} vs {len(ref.atoms)}")
    else:
        type_mismatches = []
        charge_mismatches = []
        for aid in ref.atoms:
            n_mol, n_type, n_charge, *_ = new.atoms[aid]
            r_mol, r_type, r_charge, *_ = ref.atoms[aid]
            if n_type != r_type:
                type_mismatches.append(aid)
            if abs(n_charge - r_charge) > charge_tol:
                charge_mismatches.append(aid)
        if type_mismatches:
            diffs.append(f"Atom type mismatch for {len(type_mismatches)} atoms (first: {type_mismatches[0]})")
        if charge_mismatches:
            diffs.append(f"Charge mismatch for {len(charge_mismatches)} atoms (first: {charge_mismatches[0]})")

    # 6. Coordinates (only if everything else matches — coords can drift after minimization)
    if not diffs:
        coord_mismatches = []
        for aid in ref.atoms:
            _, _, _, nx_, ny_, nz_ = new.atoms[aid]
            _, _, _, rx_, ry_, rz_ = ref.atoms[aid]
            if any(abs(a - b) > coord_tol for a, b in [(nx_, rx_), (ny_, ry_), (nz_, rz_)]):
                coord_mismatches.append(aid)
        if coord_mismatches:
            diffs.append(
                f"Coordinate mismatch for {len(coord_mismatches)} atoms "
                f"(tol={coord_tol} Å, first: atom {coord_mismatches[0]})"
            )

    return diffs


# ---------------------------------------------------------------------------
# Convenience assertion
# ---------------------------------------------------------------------------

def assert_lammps_identical(
    new_path: str | Path,
    ref_path: str | Path,
    coord_tol: float = 1e-4,
    label: str = "",
) -> None:
    """Raise AssertionError with a diff summary if the files differ."""
    diffs = compare_lammps_data(new_path, ref_path, coord_tol=coord_tol)
    if diffs:
        header = f"LAMMPS data mismatch{f' ({label})' if label else ''}:\n"
        raise AssertionError(header + "\n".join(f"  - {d}" for d in diffs))
