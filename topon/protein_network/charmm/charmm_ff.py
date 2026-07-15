"""
topro.protein.charmm_ff — CHARMM RTF/PRM force-field parser.

Parses topology (RTF) and parameter (PRM) files in CHARMM36 format.
Hard-codes TIP3P/SOD/CLA masses that may be absent in some RTF distributions.
"""

import itertools


class CHARMMForceField:
    """Parse CHARMM36 RTF and PRM files."""

    def __init__(self, prm_path, rtf_path):
        self.masses = {}         # atom_type → mass
        self.bonds_prm = {}      # (t1,t2) sorted → (Kb, b0)
        self.angles_prm = {}     # (t1,t2,t3)    → (Ktheta, theta0, Kub, S0)
        self.dihedrals_prm = {}  # (t1,t2,t3,t4) → [(Kchi, n, delta), …]
        self.impropers_prm = {}  # (t1,t2,t3,t4) → (Kpsi, psi0)
        self.vdw_prm = {}        # atom_type     → (epsilon, Rmin/2)
        self.residues = {}       # res_name       → data dict
        self.patches = {}        # patch_name     → data dict

        self._parse_prm(prm_path)
        self._parse_rtf(rtf_path)

        # Guarantee TIP3P and ion entries even if absent from the RTF
        self.masses.setdefault("OT", 15.9994)
        self.masses.setdefault("HT", 1.0080)
        self.masses.setdefault("SOD", 22.9898)
        self.masses.setdefault("CLA", 35.4500)

    # ── PRM parser ────────────────────────────────────────────────────────────

    def _parse_prm(self, path):
        section = None
        with open(path, "r") as f:
            for raw in f:
                line = raw.split("!")[0].strip()
                if not line:
                    continue
                parts = line.split()
                first = parts[0] if parts else ""

                if first in ("BONDS", "ANGLES", "DIHEDRALS", "IMPROPER",
                             "NONBONDED", "CMAP", "HBOND", "NBFIX", "END"):
                    section = first
                    continue

                if first == "MASS" and len(parts) >= 4:
                    self.masses[parts[2]] = float(parts[3])
                    continue

                if section == "BONDS" and len(parts) >= 4:
                    try:
                        key = tuple(sorted([parts[0], parts[1]]))
                        self.bonds_prm[key] = (float(parts[2]), float(parts[3]))
                    except (ValueError, IndexError):
                        pass

                elif section == "ANGLES" and len(parts) >= 5:
                    try:
                        t1, t2, t3 = parts[0], parts[1], parts[2]
                        ktheta = float(parts[3])
                        theta0 = float(parts[4])
                        kub = float(parts[5]) if len(parts) > 5 else 0.0
                        s0 = float(parts[6]) if len(parts) > 6 else 0.0
                        self.angles_prm[(t1, t2, t3)] = (ktheta, theta0, kub, s0)
                    except (ValueError, IndexError):
                        pass

                elif section == "DIHEDRALS" and len(parts) >= 7:
                    try:
                        key = (parts[0], parts[1], parts[2], parts[3])
                        entry = (float(parts[4]), int(parts[5]), float(parts[6]))
                        self.dihedrals_prm.setdefault(key, []).append(entry)
                    except (ValueError, IndexError):
                        pass

                elif section == "IMPROPER" and len(parts) >= 7:
                    try:
                        key = (parts[0], parts[1], parts[2], parts[3])
                        self.impropers_prm[key] = (float(parts[4]), float(parts[6]))
                    except (ValueError, IndexError):
                        pass

                elif section == "NONBONDED" and len(parts) >= 4:
                    try:
                        self.vdw_prm[parts[0]] = (float(parts[2]), float(parts[3]))
                    except (ValueError, IndexError):
                        pass

    # ── RTF parser ────────────────────────────────────────────────────────────

    def _parse_rtf(self, path):
        current_data = None
        with open(path, "r") as f:
            for raw in f:
                line = raw.split("!")[0].strip()
                if not line:
                    continue
                parts = line.split()

                if parts[0] in ("RESI", "PRES"):
                    block_type = parts[0]
                    name = parts[1]
                    current_data = {
                        "atoms": {},
                        "bonds": [],
                        "impropers": [],
                        "deletes": [],
                        "ics": [],
                        "charge": float(parts[2]) if len(parts) > 2 else 0.0,
                    }
                    if block_type == "RESI":
                        self.residues[name] = current_data
                    else:
                        self.patches[name] = current_data

                elif parts[0] == "ATOM" and current_data and len(parts) >= 4:
                    current_data["atoms"][parts[1]] = (parts[2], float(parts[3]))

                elif parts[0] in ("BOND", "DOUBLE") and current_data:
                    for i in range(1, len(parts) - 1, 2):
                        current_data["bonds"].append((parts[i], parts[i + 1]))

                elif parts[0] == "IMPR" and current_data:
                    for i in range(1, len(parts) - 3, 4):
                        if i + 3 < len(parts):
                            current_data["impropers"].append(
                                (parts[i], parts[i + 1], parts[i + 2], parts[i + 3])
                            )

                elif parts[0] == "DELETE" and len(parts) >= 3 and current_data:
                    if parts[1] == "ATOM":
                        current_data["deletes"].append(parts[2])

                elif parts[0] == "IC" and current_data and len(parts) >= 10:
                    # Internal-coordinate entry:
                    #   IC a1 a2 a3 a4  R(1-2) A(1-2-3) D(1-2-3-4) A(2-3-4) R(3-4)
                    # A leading '*' on a3 marks an improper IC (a3 central); the
                    # geometry to place a4 from a1,a2,a3 is the same either way
                    # (dihedral a1-a2-a3-a4, angle a2-a3-a4, bond a3-a4), so we
                    # keep the flag only for reference.
                    a3 = parts[3]
                    improper = a3.startswith("*")
                    a3 = a3[1:] if improper else a3
                    try:
                        current_data["ics"].append({
                            "atoms": (parts[1], parts[2], a3, parts[4]),
                            "r12": float(parts[5]), "a123": float(parts[6]),
                            "d1234": float(parts[7]), "a234": float(parts[8]),
                            "r34": float(parts[9]), "improper": improper,
                        })
                    except ValueError:
                        pass

    # ── Look-up helpers ───────────────────────────────────────────────────────

    def lookup_bond(self, t1, t2):
        return self.bonds_prm.get(tuple(sorted([t1, t2])))

    def lookup_angle(self, t1, t2, t3):
        return (self.angles_prm.get((t1, t2, t3)) or
                self.angles_prm.get((t3, t2, t1)))

    def lookup_dihedral(self, t1, t2, t3, t4):
        """Look up a proper dihedral, honouring CHARMM wildcard terms.

        CHARMM .prm dihedrals are matched in priority order:
          1. fully specified term (all four atom types), forward or reverse;
          2. wildcard term ``X t2 t3 X`` (wildcards only ever sit on the two
             OUTER atoms in CHARMM36), forward or reverse on the inner pair.

        The previous implementation tried only step (1), so every dihedral
        whose CHARMM parameter is defined as a wildcard term -- which is the
        majority of proline-ring (CP1/CP2/CP3) and many sidechain torsions --
        fell through to the writer's K=0 ``(DEFAULT)`` fallback, silently
        zeroing those torsion barriers. Mirrors the wildcard fallback that
        ``lookup_improper`` already performs.
        """
        # 1. Fully specified term (forward or reverse).
        hit = (self.dihedrals_prm.get((t1, t2, t3, t4)) or
               self.dihedrals_prm.get((t4, t3, t2, t1)))
        if hit:
            return hit
        # 2. Wildcard on the outer atoms (inner pair forward or reverse).
        return (self.dihedrals_prm.get(("X", t2, t3, "X")) or
                self.dihedrals_prm.get(("X", t3, t2, "X")))

    def lookup_improper(self, t1, t2, t3, t4):
        """CHARMM impropers, matched bidirectionally with X wildcards.

        CHARMM improper parameters may list the central atom either first
        (e.g. ``C NC2 NC2 NC2``) or last (e.g. ``NC2 X X C`` for the
        guanidinium, ``OC X X CC`` for the carboxylate), and use X
        wildcards on the two non-central positions. The lookup therefore
        tries the improper in BOTH directions (t1->t4 and t4->t1), fixing
        the leading atom as central and permuting the other three, against
        exact parameter keys first and wildcard keys second.

        The previous implementation fixed only ``t1`` as central, so
        impropers whose CHARMM parameter places the central atom last
        (ARG guanidinium, ASP/GLU/C-term carboxylate) fell through to the
        writer's K=20 ``(DEFAULT)`` fallback, applying a wrong planarity
        constant (should be 45 for guanidinium, 96 for carboxylate).
        """
        directions = [(t1, t2, t3, t4), (t4, t3, t2, t1)]
        # 1. Exact: leading atom central, permute the other three.
        for c0, a, b, c in directions:
            for perm in itertools.permutations([a, b, c]):
                key = (c0, *perm)
                if key in self.impropers_prm:
                    return self.impropers_prm[key]
        # 2. Wildcard fallback (X stands in for any non-central atom).
        for c0, a, b, c in directions:
            for perm in itertools.permutations([a, b, c]):
                for key in self.impropers_prm:
                    if key[0] in (c0, "X") and all(
                        k in (p, "X") for k, p in zip(key[1:], perm)
                    ):
                        return self.impropers_prm[key]
        return None
