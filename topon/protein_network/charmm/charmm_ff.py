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

    # ── Look-up helpers ───────────────────────────────────────────────────────

    def lookup_bond(self, t1, t2):
        return self.bonds_prm.get(tuple(sorted([t1, t2])))

    def lookup_angle(self, t1, t2, t3):
        return (self.angles_prm.get((t1, t2, t3)) or
                self.angles_prm.get((t3, t2, t1)))

    def lookup_dihedral(self, t1, t2, t3, t4):
        return (self.dihedrals_prm.get((t1, t2, t3, t4)) or
                self.dihedrals_prm.get((t4, t3, t2, t1)))

    def lookup_improper(self, t1, t2, t3, t4):
        """CHARMM impropers: central atom first.  Try all orderings of (t2,t3,t4)."""
        for perm in itertools.permutations([t2, t3, t4]):
            key = (t1, *perm)
            if key in self.impropers_prm:
                return self.impropers_prm[key]
        # Wildcard fallback
        for perm in itertools.permutations([t2, t3, t4]):
            for key in self.impropers_prm:
                if key[0] in (t1, "X") and all(
                    k in (p, "X") for k, p in zip(key[1:], perm)
                ):
                    return self.impropers_prm[key]
        return None
