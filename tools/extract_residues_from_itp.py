"""Extract MARTINI 3 residue-bead-and-bonded table from a polyply-generated protein ITP.

Reads a polyply-generated protein .itp (e.g. tests/_martini_extracted/Martini_Ahmet/
itp_files/nat_pro.itp) and emits topon/protein_network/residues.py: a generated
Python module mapping residue name -> bead pattern + intra-residue bonded terms,
plus the inter-residue backbone parameters and the canonical resilin block.

Optionally prunes the 16 MB master MARTINI 3 ITP down to just the bead types
referenced by the protein, plus their pairwise [ nonbond_params ] rows, and
copies the water solvent ITP into the package data folder.

Usage:
    python tools/extract_residues_from_itp.py \\
        --itp tests/_martini_extracted/Martini_Ahmet/itp_files/nat_pro.itp \\
        --master tests/_martini_extracted/Martini_Ahmet/ff/martini_v3.0.0.itp \\
        --water tests/_martini_extracted/Martini_Ahmet/ff/martini_v3.0.0_solvents_v1.itp \\
        --out-residues topon/protein_network/residues.py \\
        --out-protein-itp topon/protein_network/data/martini_v3_protein.itp \\
        --out-water-itp topon/protein_network/data/martini_v3_water.itp
"""
from __future__ import annotations

import argparse
import re
import shutil
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


@dataclass
class Atom:
    idx: int          # 1-based ITP atom id
    bead_type: str    # MARTINI bead type (Q5, P2, SP1, ...)
    resnr: int        # 1-based residue number
    resname: str      # 3-letter residue name
    atom_name: str    # BB, SC1, SC2, ...
    charge: float


@dataclass
class Section:
    name: str
    rows: list[list[str]] = field(default_factory=list)  # tokenized rows, comments stripped
    raw_lines: list[str] = field(default_factory=list)   # raw text for round-trip


def parse_itp(path: Path) -> dict[str, list[Section]]:
    """Parse an ITP file into ordered sections. Sections may repeat (e.g. dihedrals)."""
    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    sections: dict[str, list[Section]] = defaultdict(list)
    current: Section | None = None
    skip_block = False  # honour #ifdef FLEXIBLE / #ifndef FLEXIBLE

    for line in text:
        s = line.strip()
        # preprocessor: we want the constraints (default) and skip the FLEXIBLE bonds
        if s.startswith("#ifdef") and "FLEXIBLE" in s:
            skip_block = True
            continue
        if s.startswith("#ifndef") and "FLEXIBLE" in s:
            skip_block = False
            continue
        if s.startswith("#endif"):
            skip_block = False
            continue
        if skip_block:
            continue
        if not s or s.startswith(";"):
            continue
        m = re.match(r"\[\s*(\w+)\s*\]", s)
        if m:
            current = Section(name=m.group(1))
            sections[current.name].append(current)
            continue
        if current is None:
            continue
        # strip trailing comments
        body = re.split(r"\s+;", line, maxsplit=1)[0]
        toks = body.split()
        if toks:
            current.rows.append(toks)
            current.raw_lines.append(line.rstrip())
    return dict(sections)


def parse_atoms(section: Section) -> list[Atom]:
    out: list[Atom] = []
    for row in section.rows:
        if len(row) < 7:
            continue
        idx = int(row[0])
        bead_type = row[1]
        resnr = int(row[2])
        resname = row[3]
        atom_name = row[4]
        # row[5] is cgnr, row[6] is charge
        try:
            charge = float(row[6])
        except ValueError:
            charge = 0.0
        out.append(Atom(idx, bead_type, resnr, resname, atom_name, charge))
    return out


def group_by_residue(atoms: list[Atom]) -> dict[int, list[Atom]]:
    out: dict[int, list[Atom]] = defaultdict(list)
    for a in atoms:
        out[a.resnr].append(a)
    return out


def extract_canonical_residues(
    atoms: list[Atom],
) -> tuple[dict[str, list[tuple[str, str, float]]], list[tuple[str, str, float]] | None, list[tuple[str, str, float]] | None, dict[str, list[int]]]:
    """For each residue name, pick the first mid-chain occurrence as canonical.

    Returns (canonical, n_term_first, c_term_last, occurrences_by_name) where:
      canonical[resname] = [(atom_name, bead_type, charge), ...] for mid-chain.
      n_term_first: bead pattern of residue 1 (or None if matches canonical).
      c_term_last: bead pattern of last residue (or None if matches canonical).
      occurrences_by_name[resname] = list of resnr where this residue appears.
    """
    by_res = group_by_residue(atoms)
    sorted_resnrs = sorted(by_res.keys())
    first_resnr, last_resnr = sorted_resnrs[0], sorted_resnrs[-1]

    occ: dict[str, list[int]] = defaultdict(list)
    for r in sorted_resnrs:
        name = by_res[r][0].resname
        occ[name].append(r)

    canonical: dict[str, list[tuple[str, str, float]]] = {}
    for name, rnrs in occ.items():
        # pick first NON-terminal occurrence as canonical
        candidate = next((r for r in rnrs if r != first_resnr and r != last_resnr), rnrs[0])
        beads = [(a.atom_name, a.bead_type, a.charge) for a in by_res[candidate]]
        canonical[name] = beads

    first_residue = by_res[first_resnr]
    last_residue = by_res[last_resnr]
    n_term = [(a.atom_name, a.bead_type, a.charge) for a in first_residue]
    c_term = [(a.atom_name, a.bead_type, a.charge) for a in last_residue]

    n_term_diff = n_term if n_term != canonical.get(first_residue[0].resname) else None
    c_term_diff = c_term if c_term != canonical.get(last_residue[0].resname) else None
    return canonical, n_term_diff, c_term_diff, dict(occ)


def is_intra_residue(atom_idxs: list[int], by_idx: dict[int, Atom]) -> bool:
    rs = {by_idx[i].resnr for i in atom_idxs}
    return len(rs) == 1


def extract_intra_residue_bonded(
    section_rows: list[list[str]],
    by_idx: dict[int, Atom],
    n_atom_idxs: int,
    n_param_floats: int,
) -> dict[str, list[tuple]]:
    """Group rows by residue name; only keep intra-residue rows.

    Each row is (idx_a, idx_b, ...) followed by funct + n_param_floats numeric params.
    Returns {resname: [(atom_name_a, atom_name_b, ..., funct, *params), ...]} where
    only the FIRST occurrence's pattern per residue is kept (after de-dup by atom-name
    pattern + params).
    """
    out: dict[str, list[tuple]] = defaultdict(list)
    seen: dict[str, set[tuple]] = defaultdict(set)
    for row in section_rows:
        if len(row) < n_atom_idxs + 1:
            continue
        try:
            idxs = [int(row[i]) for i in range(n_atom_idxs)]
        except ValueError:
            continue
        if not all(i in by_idx for i in idxs):
            continue
        if not is_intra_residue(idxs, by_idx):
            continue
        try:
            funct = int(row[n_atom_idxs])
        except ValueError:
            continue
        params: list[float] = []
        for i in range(n_atom_idxs + 1, n_atom_idxs + 1 + n_param_floats):
            if i >= len(row):
                break
            try:
                params.append(float(row[i]))
            except ValueError:
                params.append(0.0)
        atoms = [by_idx[i] for i in idxs]
        resname = atoms[0].resname
        atom_names = tuple(a.atom_name for a in atoms)
        key = atom_names + (funct,) + tuple(params)
        if key in seen[resname]:
            continue
        seen[resname].add(key)
        out[resname].append(atom_names + (funct,) + tuple(params))
    return dict(out)


def extract_intra_residue_exclusions(
    section_rows: list[list[str]],
    by_idx: dict[int, Atom],
) -> dict[str, list[tuple[str, ...]]]:
    """Each row: pivot atom + atoms it is excluded from. Group by residue name,
    only intra-residue."""
    out: dict[str, list[tuple[str, ...]]] = defaultdict(list)
    seen: dict[str, set[tuple]] = defaultdict(set)
    for row in section_rows:
        try:
            idxs = [int(t) for t in row]
        except ValueError:
            continue
        if not all(i in by_idx for i in idxs):
            continue
        if not is_intra_residue(idxs, by_idx):
            continue
        atoms = [by_idx[i] for i in idxs]
        resname = atoms[0].resname
        names = tuple(a.atom_name for a in atoms)
        if names in seen[resname]:
            continue
        seen[resname].add(names)
        out[resname].append(names)
    return dict(out)


def extract_backbone_bonds(
    bond_rows: list[list[str]],
    by_idx: dict[int, Atom],
) -> dict[tuple[str, str, str], list[tuple[int, float, float]]]:
    """Backbone BB-BB bonds, grouped by (prev_resname, this_resname, next_resname).

    Returns map -> list of (funct, length_nm, k_kJ_per_mol_per_nm2). Most chains use
    one canonical pattern; PRO neighbours create variants.
    """
    out: dict[tuple[str, str, str], list[tuple[int, float, float]]] = defaultdict(list)
    seen: set[tuple] = set()
    # build a map of resnr -> ordered atoms
    by_res = group_by_residue(list(by_idx.values()))
    sorted_resnrs = sorted(by_res.keys())

    for row in bond_rows:
        if len(row) < 5:
            continue
        try:
            i, j = int(row[0]), int(row[1])
            funct = int(row[2])
            length = float(row[3])
            k = float(row[4])
        except ValueError:
            continue
        if i not in by_idx or j not in by_idx:
            continue
        a, b = by_idx[i], by_idx[j]
        if a.resnr == b.resnr:
            continue  # intra-residue handled elsewhere
        if a.atom_name != "BB" or b.atom_name != "BB":
            continue  # only canonical BB-BB
        # neighbour context
        rs = sorted([a.resnr, b.resnr])
        prev_r = rs[0] - 1
        next_r = rs[1] + 1
        prev_name = by_res[prev_r][0].resname if prev_r in by_res else "<NONE>"
        this_name_a = a.resname if a.resnr == rs[0] else b.resname
        this_name_b = a.resname if a.resnr == rs[1] else b.resname
        next_name = by_res[next_r][0].resname if next_r in by_res else "<NONE>"
        key = (this_name_a, this_name_b)
        rec = (funct, round(length, 4), round(k, 1))
        sk = key + rec
        if sk in seen:
            continue
        seen.add(sk)
        out[key + ("",)].append(rec)
    return dict(out)


def extract_bbb_angle(angle_rows: list[list[str]], by_idx: dict[int, Atom]) -> tuple[int, float, float] | None:
    """Find the canonical inter-residue BB-BB-BB angle parameters."""
    for row in angle_rows:
        if len(row) < 6:
            continue
        try:
            i, j, k = int(row[0]), int(row[1]), int(row[2])
            funct = int(row[3])
            angle = float(row[4])
            kc = float(row[5])
        except ValueError:
            continue
        if not all(x in by_idx for x in (i, j, k)):
            continue
        a, b, c = by_idx[i], by_idx[j], by_idx[k]
        if a.atom_name == b.atom_name == c.atom_name == "BB":
            return (funct, round(angle, 3), round(kc, 3))
    return None


def extract_bbbb_dihedrals(
    dihedral_rows: list[list[str]], by_idx: dict[int, Atom]
) -> dict[str, list[tuple]]:
    """Backbone BBBB dihedrals grouped by 4-residue sequence pattern (e.g. BBBB, GGGX).

    Reads comments to bucket terms; without comments, falls back to 'unlabelled'.
    Returns {label: [(funct, angle, k, mult), ...]}.
    """
    out: dict[str, list[tuple]] = defaultdict(list)
    seen: dict[str, set[tuple]] = defaultdict(set)
    return out  # populated downstream from raw_lines (need comment context)


def extract_bbbb_dihedrals_with_comments(
    section: Section, by_idx: dict[int, Atom]
) -> dict[str, list[tuple[int, float, float, int]]]:
    """Use the raw lines to capture the trailing comment label (e.g. 'BBBB-v1')."""
    out: dict[str, list[tuple[int, float, float, int]]] = defaultdict(list)
    seen: dict[str, set[tuple]] = defaultdict(set)
    for raw in section.raw_lines:
        body, sep, comment = raw.partition(";")
        toks = body.split()
        if len(toks) < 8:
            continue
        try:
            i, j, k, l = int(toks[0]), int(toks[1]), int(toks[2]), int(toks[3])
            funct = int(toks[4])
            angle = float(toks[5])
            kc = float(toks[6])
            mult = int(toks[7])
        except ValueError:
            continue
        if not all(x in by_idx for x in (i, j, k, l)):
            continue
        a, b, c, d = by_idx[i], by_idx[j], by_idx[k], by_idx[l]
        if not (a.atom_name == b.atom_name == c.atom_name == d.atom_name == "BB"):
            continue
        # bucket by trailing comment label up to the dash version
        label = comment.strip().split("-")[0] or "BBBB"
        rec = (funct, round(angle, 3), round(kc, 3), mult)
        if rec in seen[label]:
            continue
        seen[label].add(rec)
        out[label].append(rec)
    return dict(out)


def extract_resilin_block(atoms: list[Atom]) -> tuple[str, int]:
    """Detect the smallest repeating residue-name block in the chain.

    Returns (one_letter_block, n_repeats).
    """
    seq3 = []
    last_resnr = -1
    for a in atoms:
        if a.resnr != last_resnr:
            seq3.append(a.resname)
            last_resnr = a.resnr
    seq1 = "".join(THREE_TO_ONE.get(r, "?") for r in seq3)
    n = len(seq1)
    for blen in range(1, n // 2 + 1):
        if n % blen != 0:
            continue
        block = seq1[:blen]
        if all(seq1[i:i + blen] == block for i in range(0, n, blen)):
            return block, n // blen
    return seq1, 1


def collect_bead_types_used(atoms: list[Atom]) -> set[str]:
    return {a.bead_type for a in atoms}


def prune_master_itp(master_path: Path, types_used: set[str], out_path: Path) -> None:
    """Read master ITP and keep only [ defaults ], [ atomtypes ] for types_used,
    and [ nonbond_params ] for type pairs both in types_used. Writes a slim file."""
    text = master_path.read_text(encoding="utf-8", errors="replace").splitlines()
    section: str | None = None
    out_lines: list[str] = [
        f"; Pruned from {master_path.name} for topon/protein_network.",
        f"; Retained atomtypes ({len(types_used)}): " + ", ".join(sorted(types_used)),
        "",
    ]
    for line in text:
        s = line.strip()
        m = re.match(r"\[\s*(\w+)\s*\]", s)
        if m:
            section = m.group(1)
            out_lines.append(line)
            continue
        if not s or s.startswith(";") or s.startswith("#"):
            out_lines.append(line)
            continue
        if section == "atomtypes":
            toks = s.split()
            if toks and toks[0] in types_used:
                out_lines.append(line)
        elif section == "nonbond_params":
            toks = s.split()
            if len(toks) >= 2 and toks[0] in types_used and toks[1] in types_used:
                out_lines.append(line)
        elif section == "defaults":
            out_lines.append(line)
        else:
            # carry through other sections verbatim (rare in master ITP)
            out_lines.append(line)
    out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def emit_residues_module(
    canonical: dict[str, list[tuple[str, str, float]]],
    n_term: list[tuple[str, str, float]] | None,
    c_term: list[tuple[str, str, float]] | None,
    intra_bonds: dict[str, list[tuple]],
    intra_constraints: dict[str, list[tuple]],
    intra_angles: dict[str, list[tuple]],
    intra_dihedrals_propers: dict[str, list[tuple]],
    intra_dihedrals_impropers: dict[str, list[tuple]],
    intra_exclusions: dict[str, list[tuple[str, ...]]],
    bbb_angle: tuple[int, float, float] | None,
    bbbb_dihedrals_by_label: dict[str, list[tuple]],
    bb_bond_variants: list[tuple[str, str, int, float, float]],
    resilin_block: str,
    n_repeats: int,
    occurrences: dict[str, list[int]],
    source_itp: Path,
    out_path: Path,
) -> None:
    types_used = sorted({bt for beads in canonical.values() for (_n, bt, _c) in beads}
                        | ({bt for (_n, bt, _c) in (n_term or [])})
                        | ({bt for (_n, bt, _c) in (c_term or [])}))

    def fmt_beads(b):
        return "[" + ", ".join(f"({n!r}, {t!r}, {c!r})" for n, t, c in b) + "]"

    def fmt_records(records):
        return "[" + ", ".join(repr(r) for r in records) + "]"

    lines = [
        '"""MARTINI 3 protein residue table.',
        "",
        f"Auto-generated from {source_itp.name} by tools/extract_residues_from_itp.py.",
        "Edit the extractor and re-run; do not hand-edit this file.",
        "",
        f"Detected repeating block: {resilin_block!r} x {n_repeats}.",
        f"Bead types referenced ({len(types_used)}): " + ", ".join(types_used) + ".",
        '"""',
        "from __future__ import annotations",
        "",
        f"SOURCE_ITP = {source_itp.name!r}",
        f"REFERENCE_BLOCK = {resilin_block!r}",
        f"REFERENCE_N_REPEATS = {n_repeats}",
        "",
        "THREE_TO_ONE = {",
    ]
    for k, v in sorted(THREE_TO_ONE.items()):
        lines.append(f"    {k!r}: {v!r},")
    lines += ["}", "ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}", ""]

    # Canonical residues
    lines.append("# Canonical mid-chain bead pattern per residue: [(atom_name, bead_type, charge), ...]")
    lines.append("RESIDUES: dict[str, dict] = {")
    for resname in sorted(canonical):
        lines.append(f"    {resname!r}: {{")
        lines.append(f'        "beads": {fmt_beads(canonical[resname])},')
        lines.append(f'        "intra_bonds": {fmt_records(intra_bonds.get(resname, []))},')
        lines.append(f'        "intra_constraints": {fmt_records(intra_constraints.get(resname, []))},')
        lines.append(f'        "intra_angles": {fmt_records(intra_angles.get(resname, []))},')
        lines.append(f'        "intra_dihedrals_proper": {fmt_records(intra_dihedrals_propers.get(resname, []))},')
        lines.append(f'        "intra_dihedrals_improper": {fmt_records(intra_dihedrals_impropers.get(resname, []))},')
        lines.append(f'        "intra_exclusions": {fmt_records(intra_exclusions.get(resname, []))},')
        lines.append(f'        "occurrences_in_reference": {occurrences.get(resname, [])},')
        lines.append("    },")
    lines.append("}")
    lines.append("")

    # Terminal patches
    lines.append("TERMINAL_PATCHES = {")
    if n_term:
        lines.append(f'    "N_term": {fmt_beads(n_term)},')
    if c_term:
        lines.append(f'    "C_term": {fmt_beads(c_term)},')
    lines.append("}")
    lines.append("")

    # Backbone params
    lines.append("# Backbone BB-BB bonds extracted as (resname_a, resname_b, funct, length_nm, k_kJ_mol_nm2).")
    lines.append("BACKBONE_BB_BONDS = [")
    for v in bb_bond_variants:
        lines.append(f"    {v!r},")
    lines.append("]")
    lines.append("")

    if bbb_angle:
        lines.append(f"BACKBONE_BBB_ANGLE = {bbb_angle!r}  # (funct, angle_deg, k_kJ_mol_rad2)")
    else:
        lines.append("BACKBONE_BBB_ANGLE = None")
    lines.append("")

    lines.append("# Backbone BBBB dihedrals bucketed by polyply-comment label (e.g. 'BBBB', 'GGGX').")
    lines.append("BACKBONE_BBBB_DIHEDRALS: dict[str, list[tuple]] = {")
    for label in sorted(bbbb_dihedrals_by_label):
        terms = bbbb_dihedrals_by_label[label]
        lines.append(f"    {label!r}: [")
        for t in terms:
            lines.append(f"        {t!r},")
        lines.append("    ],")
    lines.append("}")
    lines.append("")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--itp", type=Path, required=True, help="Input polyply protein ITP")
    ap.add_argument("--master", type=Path, default=None, help="Master MARTINI 3 ITP (for pruning)")
    ap.add_argument("--water", type=Path, default=None, help="MARTINI solvents ITP (copied verbatim)")
    ap.add_argument("--out-residues", type=Path, required=True, help="Output residues.py")
    ap.add_argument("--out-protein-itp", type=Path, default=None, help="Output pruned protein ITP")
    ap.add_argument("--out-water-itp", type=Path, default=None, help="Output water ITP")
    args = ap.parse_args()

    print(f"[extract] reading {args.itp}")
    sections = parse_itp(args.itp)
    atom_section = sections["atoms"][0]
    atoms = parse_atoms(atom_section)
    by_idx = {a.idx: a for a in atoms}
    print(f"[extract] parsed {len(atoms)} atoms across {len({a.resnr for a in atoms})} residues")

    canonical, n_term, c_term, occurrences = extract_canonical_residues(atoms)
    print(f"[extract] residues found: {sorted(canonical.keys())}")
    if n_term:
        print(f"[extract] N-terminal patch detected: {n_term}")
    if c_term:
        print(f"[extract] C-terminal patch detected: {c_term}")

    bond_rows = sections["bonds"][0].rows if "bonds" in sections else []
    constraint_rows = sections["constraints"][0].rows if "constraints" in sections else []
    angle_rows = sections["angles"][0].rows if "angles" in sections else []
    dihedral_sections = sections.get("dihedrals", [])
    excl_rows = sections["exclusions"][0].rows if "exclusions" in sections else []

    intra_bonds = extract_intra_residue_bonded(bond_rows, by_idx, n_atom_idxs=2, n_param_floats=2)
    intra_constraints = extract_intra_residue_bonded(constraint_rows, by_idx, n_atom_idxs=2, n_param_floats=1)
    intra_angles = extract_intra_residue_bonded(angle_rows, by_idx, n_atom_idxs=3, n_param_floats=2)

    propers, impropers = {}, {}
    if dihedral_sections:
        # First [ dihedrals ] section is propers (funct 9), second is impropers (funct 2)
        propers = extract_intra_residue_bonded(dihedral_sections[0].rows, by_idx, n_atom_idxs=4, n_param_floats=3)
        if len(dihedral_sections) > 1:
            impropers = extract_intra_residue_bonded(dihedral_sections[1].rows, by_idx, n_atom_idxs=4, n_param_floats=2)

    intra_exclusions = extract_intra_residue_exclusions(excl_rows, by_idx)

    # Inter-residue backbone params
    bb_bond_variants_raw: dict[tuple[str, str], set[tuple[int, float, float]]] = defaultdict(set)
    by_res = group_by_residue(atoms)
    for row in bond_rows:
        if len(row) < 5:
            continue
        try:
            i, j = int(row[0]), int(row[1])
            funct = int(row[2])
            length = float(row[3])
            k = float(row[4])
        except ValueError:
            continue
        if i not in by_idx or j not in by_idx:
            continue
        a, b = by_idx[i], by_idx[j]
        if a.resnr == b.resnr:
            continue
        if a.atom_name != "BB" or b.atom_name != "BB":
            continue
        bb_bond_variants_raw[(a.resname, b.resname)].add((funct, round(length, 4), round(k, 1)))
    bb_bond_variants = []
    for (ra, rb), recs in sorted(bb_bond_variants_raw.items()):
        for rec in sorted(recs):
            bb_bond_variants.append((ra, rb) + rec)
    print(f"[extract] BB-BB bond pair patterns: {len(bb_bond_variants)} (residue-pair, funct, length, k)")

    bbb_angle = extract_bbb_angle(angle_rows, by_idx)
    if bbb_angle:
        print(f"[extract] backbone BBB angle: {bbb_angle}")

    bbbb_dihedrals = {}
    if dihedral_sections:
        bbbb_dihedrals = extract_bbbb_dihedrals_with_comments(dihedral_sections[0], by_idx)
        print(f"[extract] backbone BBBB dihedral labels: {sorted(bbbb_dihedrals.keys())}")

    block, n_repeats = extract_resilin_block(atoms)
    print(f"[extract] detected repeat block: {block!r} x {n_repeats}")

    args.out_residues.parent.mkdir(parents=True, exist_ok=True)
    emit_residues_module(
        canonical=canonical,
        n_term=n_term,
        c_term=c_term,
        intra_bonds=intra_bonds,
        intra_constraints=intra_constraints,
        intra_angles=intra_angles,
        intra_dihedrals_propers=propers,
        intra_dihedrals_impropers=impropers,
        intra_exclusions=intra_exclusions,
        bbb_angle=bbb_angle,
        bbbb_dihedrals_by_label=bbbb_dihedrals,
        bb_bond_variants=bb_bond_variants,
        resilin_block=block,
        n_repeats=n_repeats,
        occurrences=occurrences,
        source_itp=args.itp,
        out_path=args.out_residues,
    )
    print(f"[extract] wrote {args.out_residues}")

    if args.master and args.out_protein_itp:
        types_used = collect_bead_types_used(atoms)
        # also include W/SW/TW for water + Q5/SQ5n etc already included
        args.out_protein_itp.parent.mkdir(parents=True, exist_ok=True)
        # also include water (W/SW/TW) and ion (TQ5 used by NA/CL) types
        prune_master_itp(args.master, types_used | {"W", "SW", "TW", "TQ5"}, args.out_protein_itp)
        before = args.master.stat().st_size
        after = args.out_protein_itp.stat().st_size
        print(f"[extract] pruned master ITP: {before:,} B -> {after:,} B "
              f"({100 * after / before:.2f}%)")

    if args.water and args.out_water_itp:
        args.out_water_itp.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(args.water, args.out_water_itp)
        print(f"[extract] copied water ITP -> {args.out_water_itp}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
