"""MartiniLibrary parser + unit conversion tests."""
from __future__ import annotations

import math

import pytest

from topon.protein_network.martini_ff import (
    KJ_TO_KCAL,
    NM_TO_ANG,
    MartiniLibrary,
    gmx_angle_k_kjmol_to_lammps,
    gmx_bond_k_to_lammps,
    gmx_lj_to_lammps,
)


@pytest.fixture(scope="module")
def lib() -> MartiniLibrary:
    return MartiniLibrary.from_package_data()


def test_unit_conversion_constants():
    assert NM_TO_ANG == 10.0
    # 1 kcal = 4.184 kJ exactly; check round-trip to the precision we publish
    assert math.isclose(KJ_TO_KCAL * 4.184, 1.0, rel_tol=0, abs_tol=1e-9)


def test_bond_k_conversion_against_known_reference():
    # MARTINI 3 backbone BB-BB bond: K_GMX = 8000 kJ/(mol*nm^2)
    # Expected K_LAMMPS = 0.5 * 8000 * 0.2390057361376673 / 100 = 9.5602 kcal/(mol*A^2)
    k_lammps = gmx_bond_k_to_lammps(8000.0)
    assert math.isclose(k_lammps, 9.560229445506693, rel_tol=1e-9)


def test_angle_k_conversion():
    # MARTINI 3 backbone restricted-bending: K_GMX = 25 kJ/mol
    # Expected K_LAMMPS = 0.5 * 25 * 0.2390057... = 2.9876
    k_lammps = gmx_angle_k_kjmol_to_lammps(25.0)
    assert math.isclose(k_lammps, 2.9875717017208417, rel_tol=1e-9)


def test_lj_unit_conversion():
    sig_a, eps_kcal = gmx_lj_to_lammps(0.470, 4.65)  # W-W
    assert math.isclose(sig_a, 4.70, rel_tol=1e-12)
    assert math.isclose(eps_kcal, 4.65 * KJ_TO_KCAL, rel_tol=1e-12)


def test_pruned_protein_itp_loaded(lib):
    # All resilin bead types should be present
    regular = {"P2", "Q5"}                          # bare names: regular size, 72 amu
    small = {"SC3", "SP1", "SP2", "SP2a", "SP5", "SQ3p", "SQ5n", "TP1"}  # 'S' prefix: small, 54 amu
    tiny = {"TC3", "TC4", "TC5", "TN6"}             # 'T' prefix: tiny, 36 amu
    expected = regular | small | tiny
    assert expected.issubset(set(lib.atomtypes.keys()))
    for t in regular:
        assert lib.get_mass(t) == pytest.approx(72.0), f"{t} expected regular-size mass"
    # MARTINI 3 size hierarchy: regular > small > tiny.
    for t in small | tiny:
        assert lib.get_mass(t) < 72.0, f"{t} expected small or tiny size"


def test_water_loaded(lib):
    assert "W" in lib.moleculetypes
    w_mol = lib.moleculetypes["W"]
    assert len(w_mol.atoms) == 1
    bead = w_mol.atoms[0]
    assert bead[1] == "W"  # bead type column


def test_get_lj_pair_explicit_table_hit(lib):
    # P2-P2 should be in the explicit nonbond_params table
    sigma_nm, eps_kj = lib.get_lj_pair("P2", "P2")
    # MARTINI 3 regular-regular sigma = 0.470 nm; P2-P2 eps from the master table is 4.06.
    assert sigma_nm == pytest.approx(0.470, rel=1e-6)
    assert 3.5 < eps_kj < 4.5, f"P2-P2 eps {eps_kj} outside expected MARTINI 3 range"


def test_get_lj_pair_symmetric(lib):
    a = lib.get_lj_pair("P2", "TN6")
    b = lib.get_lj_pair("TN6", "P2")
    assert a == b


def test_get_lj_pair_lammps_units(lib):
    sigma_nm, eps_kj = lib.get_lj_pair("W", "W")
    sigma_a, eps_kcal = lib.get_lj_pair_lammps("W", "W")
    assert sigma_a == pytest.approx(sigma_nm * 10.0)
    assert eps_kcal == pytest.approx(eps_kj * KJ_TO_KCAL)


def test_iter_unique_pairs_dedups(lib):
    types = ["P2", "Q5", "SC3"]
    pairs = list(lib.iter_unique_pairs(types))
    # 3 types -> 3 self-pairs + 3 cross-pairs = 6
    assert len(pairs) == 6
    # No duplicates in either direction
    seen_keys = {tuple(sorted([p[0], p[1]])) for p in pairs}
    assert len(seen_keys) == 6


def test_water_bead_name(lib):
    assert lib.water_bead_name() == "W"
