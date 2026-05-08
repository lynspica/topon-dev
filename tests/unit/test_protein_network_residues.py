"""Sanity tests for the auto-extracted MARTINI residue table.

Pure structural assertions — no FF parsing, no chain building. Confirms the
extractor produced a well-formed table and key invariants from nat_pro.itp
(the resilin reference) survive the round-trip.
"""
from __future__ import annotations

import pytest

from topon.protein_network import residues as R


def test_eight_resilin_residues_present():
    expected = {"ALA", "ARG", "ASN", "ASP", "GLY", "PRO", "SER", "TYR"}
    assert set(R.RESIDUES) == expected


def test_repeat_block_is_resilin_consensus():
    assert R.REFERENCE_BLOCK == "GGRPSDSYGAPGGGN"
    assert R.REFERENCE_N_REPEATS == 18
    assert len(R.REFERENCE_BLOCK) * R.REFERENCE_N_REPEATS == 270


def test_three_letter_to_one_letter_round_trip():
    for three, one in R.THREE_TO_ONE.items():
        assert R.ONE_TO_THREE[one] == three


@pytest.mark.parametrize(
    "resname,expected_atoms",
    [
        ("GLY", ["BB"]),
        ("ALA", ["BB", "SC1"]),
        ("ARG", ["BB", "SC1", "SC2"]),
        ("PRO", ["BB", "SC1"]),
        ("SER", ["BB", "SC1"]),
        ("ASP", ["BB", "SC1"]),
        ("TYR", ["BB", "SC1", "SC2", "SC3", "SC4"]),
        ("ASN", ["BB", "SC1"]),
    ],
)
def test_per_residue_bead_atom_names(resname, expected_atoms):
    actual = [a for (a, _bt, _q) in R.RESIDUES[resname]["beads"]]
    assert actual == expected_atoms


def test_backbone_bb_bond_default_and_pro_variants():
    # Default GLY-GLY backbone bond: funct 1, length 0.360 nm, k 8000 kJ/mol/nm^2
    gly_gly = next(b for b in R.BACKBONE_BB_BONDS if b[0] == "GLY" and b[1] == "GLY")
    assert gly_gly[2:] == (1, 0.360, 8000.0)
    # Bond out of PRO (PRO -> GLY or PRO -> SER) is shorter (0.305) and stiffer (10000)
    pro_outs = [b for b in R.BACKBONE_BB_BONDS if b[0] == "PRO"]
    assert pro_outs, "expected PRO BB outgoing bonds in the reference"
    for b in pro_outs:
        assert b[2:] == (1, 0.305, 10000.0)


def test_backbone_bbb_angle_idp_restricted_bending():
    # MARTINI 3 IDP backbone: GROMACS angle funct 10, theta0=137 deg, k=25
    assert R.BACKBONE_BBB_ANGLE == (10, 137.0, 25.0)


def test_idp_backbone_dihedral_label_set():
    # polyply emits 5 BBBB labels for an IDP chain (default + Gly-rich variants)
    expected = {"BBBB", "GGGX", "GGXG", "GXGG", "XGGG"}
    assert set(R.BACKBONE_BBBB_DIHEDRALS) == expected


def test_tyr_has_ring_constraints_and_improper():
    tyr = R.RESIDUES["TYR"]
    # 4 SC beads + ring planarity improper
    assert len(tyr["beads"]) == 5
    assert any(rec[0] == "SC4" and rec[1] == "TN6" for rec in tyr["beads"])
    # Ring constraints: 3 sides at 0.300 nm + 2 at 0.285 nm
    lengths = sorted(c[3] for c in tyr["intra_constraints"])
    assert lengths == [0.285, 0.285, 0.3, 0.3, 0.3]
    # Planarity improper
    assert len(tyr["intra_dihedrals_improper"]) == 1
    impr = tyr["intra_dihedrals_improper"][0]
    assert impr[5:] == (180.0, 50.0)


def test_terminal_patches():
    # N-term first residue (GLY in the resilin reference): BB becomes Q5 +1
    assert R.TERMINAL_PATCHES["N_term"] == [("BB", "Q5", 1.0)]
    # C-term last residue (ASN): BB becomes Q5 -1, SC1 stays SP5
    assert R.TERMINAL_PATCHES["C_term"] == [("BB", "Q5", -1.0), ("SC1", "SP5", 0.0)]


def test_charge_neutrality_for_neutral_residues():
    # GLY/ALA/ASN/PRO/SER/TYR have zero net charge per canonical pattern
    for n in ("GLY", "ALA", "ASN", "PRO", "SER", "TYR"):
        net = sum(q for (_a, _t, q) in R.RESIDUES[n]["beads"])
        assert net == 0.0, f"{n} should be neutral"
    # ARG carries +1 (SC2 = SQ3p), ASP carries -1 (SC1 = SQ5n)
    assert sum(q for (_a, _t, q) in R.RESIDUES["ARG"]["beads"]) == 1.0
    assert sum(q for (_a, _t, q) in R.RESIDUES["ASP"]["beads"]) == -1.0
