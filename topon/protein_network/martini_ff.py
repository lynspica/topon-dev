"""Loader and query API for the MARTINI 3 force field.

Reads `[ defaults ]`, `[ atomtypes ]`, `[ nonbond_params ]`, and per-molecule
blocks from one or more vendored ITP files. Designed for the pruned protein +
water ITPs under `topon/protein_network/data/`, but works on any GROMACS ITP.

Unit-conversion helpers convert GROMACS native units (nm, kJ/mol) to LAMMPS
``units real`` (Angstrom, kcal/mol). Conversion derivations:

* length     1 nm = 10 Angstrom                                  -> NM_TO_ANG
* energy     1 kJ/mol = 0.2390057361376673 kcal/mol              -> KJ_TO_KCAL
* LJ pair    eps unchanged unit-converted; sigma * 10
* harmonic bond
  GROMACS:  U = 0.5 * K_g * (r - r0)^2 ,  [K_g] = kJ/(mol*nm^2)
  LAMMPS:   U = K_l * (r - r0)^2 ,        [K_l] = kcal/(mol*A^2)
  K_l = 0.5 * K_g * KJ_TO_KCAL / 100  (since 1 nm^2 = 100 A^2)
       = K_g * 0.001195029 ...
* harmonic angle (funct 1)
  same prefactor handling; K_g in kJ/(mol*rad^2), K_l in kcal/(mol*rad^2)
  K_l = 0.5 * K_g * KJ_TO_KCAL
* cosine-squared angle (funct 2 in GROMACS, used by MARTINI sidechains)
  GROMACS:  U = 0.5 * K_g * (cos t - cos t0)^2  ,  [K_g] = kJ/mol
  LAMMPS angle_style cosine/squared:
            U = K_l * (cos t - cos t0)^2 ,         [K_l] = kcal/mol
  K_l = 0.5 * K_g * KJ_TO_KCAL
* restricted-bending angle (funct 10 in GROMACS, used by MARTINI 3 backbone)
  GROMACS:  U = 0.5 * K_g * (cos t - cos t0)^2 / sin^2 t
  No exact LAMMPS equivalent; angle_style cosine/squared is used as the first-
  cut surrogate (omits the 1/sin^2 t factor). Same K conversion as funct 2.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Iterator

NM_TO_ANG: float = 10.0
KJ_TO_KCAL: float = 0.2390057361376673


def gmx_bond_k_to_lammps(k_gmx: float) -> float:
    """GROMACS K [kJ/(mol*nm^2)] -> LAMMPS K [kcal/(mol*A^2)] (harmonic bond)."""
    return 0.5 * k_gmx * KJ_TO_KCAL / 100.0


def gmx_angle_k_kjmol_to_lammps(k_gmx: float) -> float:
    """GROMACS K [kJ/mol] -> LAMMPS K [kcal/mol] for cosine-squared angle styles.

    Also used as the first-cut conversion for funct 10 (restricted bending).
    """
    return 0.5 * k_gmx * KJ_TO_KCAL


def gmx_dihedral_k_to_lammps(k_gmx: float) -> float:
    """GROMACS dihedral K [kJ/mol] -> LAMMPS K [kcal/mol] (proper dihedrals)."""
    return k_gmx * KJ_TO_KCAL


def gmx_lj_to_lammps(sigma_nm: float, epsilon_kj: float) -> tuple[float, float]:
    """LJ (sigma_nm, epsilon_kJ_per_mol) -> (sigma_A, epsilon_kcal_per_mol)."""
    return sigma_nm * NM_TO_ANG, epsilon_kj * KJ_TO_KCAL


@dataclass(frozen=True)
class Atomtype:
    name: str
    mass: float            # amu
    charge_default: float  # e (overridden per-atom in molecule blocks)
    particle_type: str     # usually "A" for MARTINI atomtypes
    sigma_nm: float
    epsilon_kj: float


@dataclass
class MoleculeType:
    name: str
    nrexcl: int
    atoms: list[tuple[int, str, int, str, str, int, float]]  # (id, type, resnr, resname, atom_name, cgnr, charge)


def _strip_comment(line: str) -> str:
    return re.split(r"\s*;", line, maxsplit=1)[0]


def _iter_sections(text: str) -> Iterator[tuple[str, list[str]]]:
    """Yield (section_name, rows) tuples. Rows are non-empty stripped lines with
    comments removed, preserving order. Honours #ifdef FLEXIBLE / #ifndef FLEXIBLE
    by skipping the FLEXIBLE branch."""
    section = ""
    rows: list[str] = []
    skip = False
    for raw in text.splitlines():
        s = raw.strip()
        if s.startswith("#ifdef") and "FLEXIBLE" in s:
            skip = True
            continue
        if s.startswith("#ifndef") and "FLEXIBLE" in s:
            skip = False
            continue
        if s.startswith("#endif"):
            skip = False
            continue
        if skip:
            continue
        m = re.match(r"\[\s*(\w+)\s*\]", s)
        if m:
            if section:
                yield section, rows
            section = m.group(1)
            rows = []
            continue
        body = _strip_comment(raw).strip()
        if body and not body.startswith(";"):
            rows.append(body)
    if section:
        yield section, rows


class MartiniLibrary:
    """Container for atomtypes, nonbond pair table, and moleculetypes."""

    def __init__(self) -> None:
        self.defaults: dict[str, str] = {}
        self.atomtypes: dict[str, Atomtype] = {}
        self.nonbond_params: dict[tuple[str, str], tuple[float, float]] = {}
        self.moleculetypes: dict[str, MoleculeType] = {}

    @classmethod
    def from_files(cls, *paths: Path) -> "MartiniLibrary":
        lib = cls()
        for p in paths:
            lib.load_itp(Path(p))
        return lib

    @classmethod
    def from_package_data(cls) -> "MartiniLibrary":
        """Load the vendored protein + water ITPs shipped with topon."""
        with resources.as_file(resources.files("topon.protein_network.data")) as data_dir:
            return cls.from_files(
                data_dir / "martini_v3_protein.itp",
                data_dir / "martini_v3_water.itp",
            )

    def load_itp(self, path: Path) -> None:
        text = path.read_text(encoding="utf-8", errors="replace")
        current_mol: MoleculeType | None = None
        for name, rows in _iter_sections(text):
            if name == "defaults":
                for row in rows:
                    toks = row.split()
                    if len(toks) >= 2:
                        self.defaults["nbfunc"] = toks[0]
                        self.defaults["comb_rule"] = toks[1]
                    if len(toks) >= 4:
                        self.defaults["fudgeLJ"] = toks[2]
                        self.defaults["fudgeQQ"] = toks[3]
            elif name == "atomtypes":
                for row in rows:
                    toks = row.split()
                    if len(toks) < 6:
                        continue
                    try:
                        at = Atomtype(
                            name=toks[0],
                            mass=float(toks[1]),
                            charge_default=float(toks[2]),
                            particle_type=toks[3],
                            sigma_nm=float(toks[4]),
                            epsilon_kj=float(toks[5]),
                        )
                    except ValueError:
                        continue
                    self.atomtypes[at.name] = at
            elif name == "nonbond_params":
                for row in rows:
                    toks = row.split()
                    if len(toks) < 5:
                        continue
                    try:
                        a, b = toks[0], toks[1]
                        sigma, eps = float(toks[3]), float(toks[4])
                    except ValueError:
                        continue
                    self.nonbond_params[(a, b)] = (sigma, eps)
                    self.nonbond_params[(b, a)] = (sigma, eps)
            elif name == "moleculetype":
                for row in rows:
                    toks = row.split()
                    if len(toks) >= 2:
                        try:
                            current_mol = MoleculeType(name=toks[0], nrexcl=int(toks[1]), atoms=[])
                            self.moleculetypes[current_mol.name] = current_mol
                        except ValueError:
                            current_mol = None
                        break
            elif name == "atoms" and current_mol is not None:
                for row in rows:
                    toks = row.split()
                    if len(toks) < 7:
                        continue
                    try:
                        current_mol.atoms.append((
                            int(toks[0]), toks[1], int(toks[2]), toks[3], toks[4],
                            int(toks[5]), float(toks[6]),
                        ))
                    except ValueError:
                        continue

    def get_lj_pair(self, t1: str, t2: str) -> tuple[float, float]:
        """Return (sigma_nm, epsilon_kJ) for the given bead-type pair.

        Falls back to geometric mixing (sigma = sqrt(s1*s2), eps = sqrt(e1*e2))
        when the pair is not in the explicit MARTINI table.
        """
        if (t1, t2) in self.nonbond_params:
            return self.nonbond_params[(t1, t2)]
        if (t2, t1) in self.nonbond_params:
            return self.nonbond_params[(t2, t1)]
        a1 = self.atomtypes.get(t1)
        a2 = self.atomtypes.get(t2)
        if a1 is None or a2 is None:
            raise KeyError(f"unknown bead type(s): {t1!r}, {t2!r}")
        sigma = math.sqrt(a1.sigma_nm * a2.sigma_nm)
        eps = math.sqrt(a1.epsilon_kj * a2.epsilon_kj)
        return sigma, eps

    def get_lj_pair_lammps(self, t1: str, t2: str) -> tuple[float, float]:
        """Same as get_lj_pair but in LAMMPS units (Angstrom, kcal/mol)."""
        sigma, eps = self.get_lj_pair(t1, t2)
        return gmx_lj_to_lammps(sigma, eps)

    def get_mass(self, bead_type: str) -> float:
        return self.atomtypes[bead_type].mass

    def iter_atomtypes(self) -> Iterator[Atomtype]:
        return iter(self.atomtypes.values())

    def iter_unique_pairs(self, types: list[str]) -> Iterator[tuple[str, str, float, float]]:
        """Yield (t1, t2, sigma_nm, eps_kJ) for every i <= j pair drawn from `types`,
        deduplicating symmetric storage in nonbond_params. Order matches `types`.
        """
        seen: set[tuple[str, str]] = set()
        for i, ti in enumerate(types):
            for tj in types[i:]:
                key = (ti, tj) if ti <= tj else (tj, ti)
                if key in seen:
                    continue
                seen.add(key)
                sigma, eps = self.get_lj_pair(ti, tj)
                yield key[0], key[1], sigma, eps

    def water_bead_name(self) -> str:
        """Return the canonical regular-water moleculetype name (e.g. 'W')."""
        for name in ("W", "WL"):
            if name in self.moleculetypes:
                return name
        raise KeyError("no W water moleculetype found in loaded ITPs")
