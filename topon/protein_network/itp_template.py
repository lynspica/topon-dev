"""Per-chain ITP template loader.

Loads a polyply-generated MARTINI 3 protein ITP and exposes it as a
`ChainTemplate` -- the full bonded / nonbonded topology of one chain at the
exact resolution polyply produced. The builder replicates this template per
chain with atom-ID offsets, so the user's system gets the full polyply FF
(every IDP-specific BBB angle, every BB-BB-BB-SC improper, every multi-term
dihedral) without any "canonical extraction" pattern-matching that would lose
context-dependent terms.

This is the recommended path for any sequence the user has a polyply ITP for.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path


@dataclass(frozen=True)
class ITPAtom:
    id: int
    bead_type: str
    resnr: int
    resname: str
    atom_name: str
    cgnr: int
    charge: float


@dataclass(frozen=True)
class ITPBond:
    i: int
    j: int
    funct: int
    length_nm: float
    k_kj: float


@dataclass(frozen=True)
class ITPConstraint:
    i: int
    j: int
    funct: int
    length_nm: float


@dataclass(frozen=True)
class ITPAngle:
    i: int
    j: int
    k: int
    funct: int
    angle_deg: float
    k_kj: float


@dataclass(frozen=True)
class ITPDihedral:
    i: int
    j: int
    k: int
    l: int
    funct: int
    angle_deg: float
    k_kj: float
    mult: int = 1
    is_improper: bool = False


@dataclass(frozen=True)
class ITPExclusion:
    atoms: tuple[int, ...]


@dataclass
class ChainTemplate:
    """Full ITP topology for a single chain. Atom IDs are 1-based as in the ITP."""
    name: str
    nrexcl: int
    atoms: list[ITPAtom] = field(default_factory=list)
    bonds: list[ITPBond] = field(default_factory=list)
    constraints: list[ITPConstraint] = field(default_factory=list)
    angles: list[ITPAngle] = field(default_factory=list)
    dihedrals_proper: list[ITPDihedral] = field(default_factory=list)
    dihedrals_improper: list[ITPDihedral] = field(default_factory=list)
    exclusions: list[ITPExclusion] = field(default_factory=list)

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def n_residues(self) -> int:
        return max(a.resnr for a in self.atoms) if self.atoms else 0

    def block_sequence(self) -> tuple[str, int]:
        """Detect the smallest repeating residue-name block. Returns (block, n_repeats)."""
        from .residues import THREE_TO_ONE
        seq3, last = [], -1
        for a in self.atoms:
            if a.resnr != last:
                seq3.append(a.resname); last = a.resnr
        seq1 = "".join(THREE_TO_ONE.get(r, "?") for r in seq3)
        n = len(seq1)
        for blen in range(1, n // 2 + 1):
            if n % blen != 0: continue
            block = seq1[:blen]
            if all(seq1[i:i + blen] == block for i in range(0, n, blen)):
                return block, n // blen
        return seq1, 1


def _strip_comment(line: str) -> str:
    return re.split(r"\s*;", line, maxsplit=1)[0]


def load_chain_template(path: str | Path) -> ChainTemplate:
    """Parse a polyply protein ITP into a `ChainTemplate`.

    Honours `#ifdef FLEXIBLE` (skipped) and `#ifndef FLEXIBLE` (kept) so that
    the returned template uses constraints (the default) instead of stiff
    surrogate bonds.
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")

    name, nrexcl = "", 1
    atoms: list[ITPAtom] = []
    bonds: list[ITPBond] = []
    constraints: list[ITPConstraint] = []
    angles: list[ITPAngle] = []
    propers: list[ITPDihedral] = []
    impropers: list[ITPDihedral] = []
    exclusions: list[ITPExclusion] = []

    section: str | None = None
    seen_dihedral_blocks = 0  # 1st [dihedrals] = propers (funct 9), 2nd = impropers (funct 2)
    skip = False

    for raw in text.splitlines():
        s = raw.strip()
        if s.startswith("#ifdef") and "FLEXIBLE" in s: skip = True; continue
        if s.startswith("#ifndef") and "FLEXIBLE" in s: skip = False; continue
        if s.startswith("#endif"): skip = False; continue
        if skip: continue
        m = re.match(r"\[\s*(\w+)\s*\]", s)
        if m:
            new_section = m.group(1)
            if new_section == "dihedrals":
                seen_dihedral_blocks += 1
            section = new_section
            continue
        if not s or s.startswith(";"): continue
        body = _strip_comment(raw).strip()
        if not body: continue
        toks = body.split()

        try:
            if section == "moleculetype":
                if len(toks) >= 2:
                    name = toks[0]; nrexcl = int(toks[1])
            elif section == "atoms":
                if len(toks) >= 7:
                    atoms.append(ITPAtom(
                        id=int(toks[0]), bead_type=toks[1], resnr=int(toks[2]),
                        resname=toks[3], atom_name=toks[4], cgnr=int(toks[5]),
                        charge=float(toks[6]),
                    ))
            elif section == "bonds":
                if len(toks) >= 5:
                    bonds.append(ITPBond(
                        i=int(toks[0]), j=int(toks[1]),
                        funct=int(toks[2]), length_nm=float(toks[3]), k_kj=float(toks[4]),
                    ))
            elif section == "constraints":
                if len(toks) >= 4:
                    constraints.append(ITPConstraint(
                        i=int(toks[0]), j=int(toks[1]),
                        funct=int(toks[2]), length_nm=float(toks[3]),
                    ))
            elif section == "angles":
                if len(toks) >= 6:
                    angles.append(ITPAngle(
                        i=int(toks[0]), j=int(toks[1]), k=int(toks[2]),
                        funct=int(toks[3]), angle_deg=float(toks[4]), k_kj=float(toks[5]),
                    ))
            elif section == "dihedrals":
                if len(toks) >= 7:
                    is_imp = (seen_dihedral_blocks >= 2)
                    mult = int(toks[7]) if (not is_imp and len(toks) >= 8) else 1
                    d = ITPDihedral(
                        i=int(toks[0]), j=int(toks[1]), k=int(toks[2]), l=int(toks[3]),
                        funct=int(toks[4]), angle_deg=float(toks[5]),
                        k_kj=float(toks[6]), mult=mult, is_improper=is_imp,
                    )
                    (impropers if is_imp else propers).append(d)
            elif section == "exclusions":
                ids = tuple(int(t) for t in toks)
                if len(ids) >= 2:
                    exclusions.append(ITPExclusion(atoms=ids))
        except (ValueError, IndexError):
            continue

    return ChainTemplate(
        name=name, nrexcl=nrexcl,
        atoms=atoms, bonds=bonds, constraints=constraints, angles=angles,
        dihedrals_proper=propers, dihedrals_improper=impropers, exclusions=exclusions,
    )


# Convenience: vendored ITPs keyed by canonical resilin block sequence.
_VENDORED: dict[str, str] = {
    "GGRPSDSYGAPGGGN": "nat_pro.itp",   # 18 PRO per chain (resilin consensus)
    "GPRPSDSYGAPGPGN": "high_pro.itp",  # 36 PRO per chain (high-PRO variant)
    "GGRGSDSYGAGGGGN": "no_pro.itp",    # 0 PRO (no-PRO variant)
}


def load_vendored_template(block_seq: str) -> ChainTemplate:
    """Load a vendored polyply ITP by canonical block sequence."""
    if block_seq not in _VENDORED:
        raise KeyError(
            f"No vendored ITP for block_seq={block_seq!r}. "
            f"Known: {sorted(_VENDORED)}. "
            f"Generate one via `polyply gen_params -seq <FASTA> -lib martini3` and "
            f"pass its path explicitly via load_chain_template(path)."
        )
    with resources.as_file(resources.files("topon.protein_network.data")) as data_dir:
        return load_chain_template(data_dir / _VENDORED[block_seq])
