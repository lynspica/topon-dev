"""In-memory dataclasses for a MARTINI protein-network system.

The builder, water packer, and crosslink resolver share these dataclasses; the
LAMMPS writer consumes a `System` and emits the data + settings + input files.

Units throughout are GROMACS-native (nm, kJ/mol). The writer applies LAMMPS
unit conversions.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Bead:
    atom_id: int                                # 1-based global index
    bead_type: str                              # MARTINI bead type, e.g. "P2"
    molecule_id: int                            # 1-based molecule (chain or water)
    residue_idx: int                            # 1-based residue number; 0 for water
    residue_name: str                           # 3-letter code; "W" for water
    atom_name: str                              # "BB", "SC1", "W", ...
    charge: float                               # e
    mass: float                                 # amu
    position: tuple[float, float, float]        # angstroms, may lie outside the box;
                                                 # LAMMPS auto-wraps with `boundary p p p`
                                                 # on `read_data` and assigns image flags.


@dataclass
class Bond:
    a: int                                      # atom_id
    b: int
    funct: int                                  # GROMACS function type (1 = harmonic)
    length_nm: float
    k_kj: float                                 # kJ/(mol*nm^2) for funct 1
    is_crosslink: bool = False                  # dityrosine SC4-SC4 (post-build)


@dataclass
class Constraint:
    a: int
    b: int
    length_nm: float                            # rigid distance; becomes `fix shake`


@dataclass
class Angle:
    a: int
    b: int                                      # apex
    c: int
    funct: int                                  # 1 harmonic / 2 cosine-sq / 10 restricted
    angle_deg: float
    k_kj: float


@dataclass
class Dihedral:
    a: int
    b: int
    c: int
    d: int
    funct: int                                  # 9 multi-term proper / 2 improper-harm
    angle_deg: float
    k_kj: float
    mult: int = 1                               # only relevant for funct=9
    is_improper: bool = False


@dataclass
class Exclusion:
    """One row from the GROMACS `[ exclusions ]` section: pivot atom + the
    additional neighbours it must not interact with non-bondedly."""
    atoms: tuple[int, ...]


@dataclass
class System:
    beads: list[Bead] = field(default_factory=list)
    bonds: list[Bond] = field(default_factory=list)
    constraints: list[Constraint] = field(default_factory=list)
    angles: list[Angle] = field(default_factory=list)
    dihedrals: list[Dihedral] = field(default_factory=list)
    exclusions: list[Exclusion] = field(default_factory=list)
    box_dims_ang: tuple[float, float, float] = (0.0, 0.0, 0.0)

    def n_atoms(self) -> int:
        return len(self.beads)

    def total_charge(self) -> float:
        return sum(b.charge for b in self.beads)

    def bead_types_in_use(self) -> list[str]:
        seen: dict[str, None] = {}
        for b in self.beads:
            seen.setdefault(b.bead_type, None)
        return list(seen)
