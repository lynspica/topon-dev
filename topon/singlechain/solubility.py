"""
Hansen Solubility Parameter estimation for polymers and solvents.

Uses the Hoftyzer–Van Krevelen group-contribution method to estimate the
three HSP components (δ_d, δ_p, δ_h) from SMILES.  Also provides chain-level
estimation accounting for repeat-unit composition and end-group effects.

Reference
---------
D.W. van Krevelen & K. te Nijenhuis, *Properties of Polymers*, 4th ed.,
Elsevier, 2009, Ch. 7.

Usage::

    from topon.singlechain.solubility import estimate_hsp, compute_ra

    hsp_toluene = estimate_hsp("Cc1ccccc1")
    hsp_pdms    = estimate_hsp("[Si](C)(C)O")
    ra = compute_ra(hsp_pdms, hsp_toluene)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


# ---------------------------------------------------------------------------
# Group-contribution tables (Hoftyzer–Van Krevelen)
#
# Units:  F_di  — (J cm³)^½ mol⁻¹
#         F_pi² — J cm³ mol⁻²   (squared polar attraction)
#         E_hi  — J mol⁻¹
#         V_i   — cm³ mol⁻¹
#
# SMARTS pattern → (F_di, F_pi², E_hi, V_i)
# ---------------------------------------------------------------------------
_GROUP_CONTRIBUTIONS: list[tuple[str, str, tuple[float, float, float, float]]] = [
    # (SMARTS, label, (F_di, F_pi_sq, E_hi, V_i))
    #
    # --- Silicon ---
    ("[Si](-[#6])(-[#6])(-O)-[#6]",  "Si(CH3)2-O (siloxane repeat)",  (450.0, 0.0, 0.0, 55.0)),
    ("[Si](-[#6])(-O)",              "Si-O linkage",                   (270.0, 0.0, 0.0, 30.0)),

    # --- Fluorine ---
    ("[CX4](F)(F)(F)",               "-CF3",                           (560.0, 112000.0, 0.0, 57.5)),
    ("[CX4](F)(F)",                  "-CF2- (di-fluoro)",              (420.0, 75000.0, 0.0, 40.0)),
    ("F",                            "-F (generic)",                   (160.0, 25000.0, 0.0, 18.0)),

    # --- Nitrogen ---
    ("C#N",                          "-C≡N (nitrile)",                 (725.0, 1350000.0, 3100.0, 24.0)),

    # --- Oxygen functional groups ---
    ("[CX3](=O)[OX2]",              "-COO- (ester)",                  (668.0, 490000.0, 7000.0, 33.5)),
    ("[CX3](=O)[#6]",              "-C(=O)- (ketone)",               (563.0, 750000.0, 2000.0, 22.3)),
    ("[OX2H]",                      "-OH (hydroxyl)",                 (210.0, 250000.0, 20000.0, 10.0)),
    ("[OX2]([#6])[#6]",            "-O- (ether)",                    (235.0, 100000.0, 3000.0, 6.5)),
    ("[OX2H2]",                     "H2O",                            (210.0, 250000.0, 20000.0, 18.0)),

    # --- Chlorine ---
    ("Cl",                           "-Cl",                            (420.0, 490000.0, 400.0, 24.0)),

    # --- Aromatic ---
    ("c1ccccc1",                     "phenyl ring (C6H4=)",            (1503.0, 37000.0, 0.0, 71.4)),
    ("c1cccc1",                      "5-ring aromatic",                (1350.0, 37000.0, 0.0, 60.0)),

    # --- Hydrocarbon backbone groups ---
    ("[CH3]",                        "-CH3",                           (420.0, 0.0, 0.0, 33.5)),
    ("[CH2]",                        "-CH2-",                          (270.0, 0.0, 0.0, 16.1)),
    ("[CH1]",                        ">CH-",                           (80.0, 0.0, 0.0, -1.0)),
    ("[CX4;H0]",                     ">C< (quaternary)",               (-70.0, 0.0, 0.0, -19.2)),
    ("[CX3]=[CX3]",                  "-CH=CH- (olefin)",               (444.0, 0.0, 0.0, 28.0)),
]


@dataclass
class HSP:
    """Hansen Solubility Parameters in MPa^½."""
    delta_d: float
    delta_p: float
    delta_h: float

    @property
    def delta_total(self) -> float:
        """Total (Hildebrand) solubility parameter."""
        return (self.delta_d**2 + self.delta_p**2 + self.delta_h**2) ** 0.5

    def as_dict(self) -> dict[str, float]:
        return {
            "delta_d": round(self.delta_d, 2),
            "delta_p": round(self.delta_p, 2),
            "delta_h": round(self.delta_h, 2),
            "delta_total": round(self.delta_total, 2),
        }

    def __repr__(self) -> str:
        return (f"HSP(δ_d={self.delta_d:.1f}, δ_p={self.delta_p:.1f}, "
                f"δ_h={self.delta_h:.1f}, δ_t={self.delta_total:.1f} MPa^½)")


# ---------------------------------------------------------------------------
# Core estimation
# ---------------------------------------------------------------------------

def estimate_hsp(smiles: str, verbose: bool = False) -> HSP:
    """Estimate Hansen Solubility Parameters from a SMILES string.

    Uses the Hoftyzer–Van Krevelen group-contribution method with SMARTS-based
    fragment identification via RDKit.

    Parameters
    ----------
    smiles : str
        SMILES of the molecule or polymer repeat unit.
    verbose : bool
        Print matched groups.

    Returns
    -------
    HSP
        Estimated solubility parameters in MPa^½.
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    # Track which atoms have been assigned to a group to avoid double-counting
    assigned_atoms: set[int] = set()

    sum_Fdi = 0.0
    sum_Fpi_sq = 0.0
    sum_Ehi = 0.0
    sum_Vi = 0.0

    for smarts_str, label, (Fdi, Fpi_sq, Ehi, Vi) in _GROUP_CONTRIBUTIONS:
        pat = Chem.MolFromSmarts(smarts_str)
        if pat is None:
            continue
        matches = mol.GetSubstructMatches(pat)
        for match in matches:
            # Only count this match if at least one heavy atom hasn't been
            # assigned yet (simple greedy assignment).
            heavy_in_match = {idx for idx in match
                              if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1}
            if heavy_in_match and heavy_in_match <= assigned_atoms:
                continue  # All heavy atoms already assigned
            assigned_atoms.update(heavy_in_match)

            sum_Fdi += Fdi
            sum_Fpi_sq += Fpi_sq
            sum_Ehi += Ehi
            sum_Vi += Vi

            if verbose:
                print(f"  Matched: {label:30s}  Fdi={Fdi:7.0f}  Vi={Vi:6.1f}")

    # Fallback: use molar volume from molecular weight if group V is too small
    # (happens when not all atoms matched)
    mw = Descriptors.ExactMolWt(mol)
    if sum_Vi <= 0:
        # Rough estimate: V ≈ MW / 1.0 for liquids, MW / 1.1 for polymers
        sum_Vi = mw / 1.0
        if verbose:
            print(f"  Using fallback molar volume: {sum_Vi:.1f} cm³/mol")

    # Hoftyzer–Van Krevelen equations
    delta_d = sum_Fdi / sum_Vi
    delta_p = (sum_Fpi_sq ** 0.5) / sum_Vi
    delta_h = (max(0.0, sum_Ehi) / sum_Vi) ** 0.5

    return HSP(delta_d=delta_d, delta_p=delta_p, delta_h=delta_h)


def estimate_chain_hsp(
    repeat_smiles: str,
    dp: int,
    end_smiles: Optional[str] = None,
    copolymer_composition: Optional[list[dict]] = None,
    verbose: bool = False,
) -> HSP:
    """Estimate HSP for a full polymer chain, accounting for copolymer
    composition and end-group dilution.

    For homopolymers, this is simply ``estimate_hsp(repeat_smiles)`` (end-
    group effects are negligible at high DP).  For copolymers, contributions
    are weighted by mole fraction.

    Parameters
    ----------
    repeat_smiles : str
        SMILES for the primary repeat unit.
    dp : int
        Degree of polymerization.
    end_smiles : str, optional
        SMILES for the end group (if different from the repeat unit).
    copolymer_composition : list[dict], optional
        Same format as ``run_workflow``:
        ``[{"monomer": "M0", "smiles": "...", "fraction": 0.6}, ...]``
    verbose : bool
        Print diagnostics.

    Returns
    -------
    HSP
    """
    if copolymer_composition and len(copolymer_composition) > 1:
        # Weighted average of monomer contributions
        total_d = total_p = total_h = 0.0
        total_frac = sum(e["fraction"] for e in copolymer_composition)

        for entry in copolymer_composition:
            mon_smiles = entry.get("smiles", repeat_smiles)
            frac = entry["fraction"] / total_frac
            hsp_i = estimate_hsp(mon_smiles, verbose=verbose)
            total_d += frac * hsp_i.delta_d
            total_p += frac * hsp_i.delta_p
            total_h += frac * hsp_i.delta_h

        return HSP(delta_d=total_d, delta_p=total_p, delta_h=total_h)

    # Homopolymer (or single-component copolymer)
    hsp_repeat = estimate_hsp(repeat_smiles, verbose=verbose)

    # End-group correction (only meaningful at low DP)
    if end_smiles and dp < 50:
        hsp_end = estimate_hsp(end_smiles, verbose=verbose)
        # Weight: (dp-2) repeat units + 2 end groups
        w_rep = max(dp - 2, 0) / dp
        w_end = min(2, dp) / dp
        return HSP(
            delta_d=w_rep * hsp_repeat.delta_d + w_end * hsp_end.delta_d,
            delta_p=w_rep * hsp_repeat.delta_p + w_end * hsp_end.delta_p,
            delta_h=w_rep * hsp_repeat.delta_h + w_end * hsp_end.delta_h,
        )

    return hsp_repeat


# ---------------------------------------------------------------------------
# Hansen distance
# ---------------------------------------------------------------------------

def compute_ra(hsp1: HSP, hsp2: HSP) -> float:
    """Compute the Hansen distance R_a between two sets of HSP.

    .. math::

        R_a = \\sqrt{4(\\delta_{d1}-\\delta_{d2})^2
                     + (\\delta_{p1}-\\delta_{p2})^2
                     + (\\delta_{h1}-\\delta_{h2})^2}

    Smaller R_a → more compatible (higher swelling).
    """
    return (
        4.0 * (hsp1.delta_d - hsp2.delta_d) ** 2
        + (hsp1.delta_p - hsp2.delta_p) ** 2
        + (hsp1.delta_h - hsp2.delta_h) ** 2
    ) ** 0.5


# ---------------------------------------------------------------------------
# Summary printer
# ---------------------------------------------------------------------------

def print_hsp_matrix(
    polymers: dict[str, HSP],
    solvents: dict[str, HSP],
) -> None:
    """Print a formatted R_a matrix for polymers × solvents."""
    solv_names = list(solvents.keys())
    header = f"{'Polymer':>12s}" + "".join(f"{s:>14s}" for s in solv_names)
    print(header)
    print("-" * len(header))
    for pname, phsp in polymers.items():
        row = f"{pname:>12s}"
        for sname in solv_names:
            ra = compute_ra(phsp, solvents[sname])
            row += f"{ra:14.2f}"
        print(row)
