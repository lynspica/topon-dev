"""Tests for the LAMMPS writer."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from topon.protein_network import bfm, builder, lammps_writer, sequence, water
from topon.protein_network.martini_ff import MartiniLibrary
from topon.protein_network.system import System


@pytest.fixture(scope="module")
def lib() -> MartiniLibrary:
    return MartiniLibrary.from_package_data()


@pytest.fixture()
def small_system_no_water(lib):
    topo = bfm.generate_topology(n_chains=2, n_repeats=3, segs_per_block=2,
                                 equil_steps=0, seed=11, verbose=False)
    snap = topo["snapshots"][0]
    seq3 = sequence.build_full_sequence("GGRPSDSYGAPGGGN", 3)
    return builder.build_protein_system(snap, seq3, lib, block_seq="GGRPSDSYGAPGGGN", seed=42)


def _read_section(path: Path, section: str) -> list[str]:
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    out: list[str] = []
    in_section = False
    for line in lines:
        if line.strip() == section:
            in_section = True
            continue
        if in_section:
            stripped = line.strip()
            if not stripped:
                if out:
                    break
                continue
            if stripped.endswith("types") or stripped.endswith("xhi") or "atom" in stripped:
                # we hit the next section header / count line
                break
            if stripped[0].isalpha() and ":" not in stripped:
                # next named section
                break
            out.append(line)
    return out


def test_writer_produces_data_settings_groups_and_three_stages(lib, small_system_no_water, tmp_path: Path):
    paths = lammps_writer.write_lammps(small_system_no_water, lib, tmp_path)
    assert set(paths) == {"data", "settings", "groups", "stage1", "stage2", "stage3"}
    for kind, p in paths.items():
        assert p.exists() and p.stat().st_size > 0, f"{kind} file empty: {p}"
    # The three relaxation stages live in a relaxation/ subdirectory
    assert (tmp_path / "relaxation").is_dir()
    for stage in ("stage1", "stage2", "stage3"):
        assert paths[stage].parent.name == "relaxation"


def test_data_header_counts_match_system(lib, small_system_no_water, tmp_path: Path):
    paths = lammps_writer.write_lammps(small_system_no_water, lib, tmp_path)
    text = paths["data"].read_text(encoding="utf-8")
    n_atoms = small_system_no_water.n_atoms()
    n_constraints = len(small_system_no_water.constraints)
    n_real_bonds = len(small_system_no_water.bonds)
    n_bonds = n_real_bonds + n_constraints
    n_angles = len(small_system_no_water.angles)
    n_propers = sum(1 for d in small_system_no_water.dihedrals if not d.is_improper)
    n_impropers = sum(1 for d in small_system_no_water.dihedrals if d.is_improper)
    assert f"{n_atoms} atoms" in text
    assert f"{n_bonds} bonds" in text
    assert f"{n_angles} angles" in text
    assert f"{n_propers} dihedrals" in text
    assert f"{n_impropers} impropers" in text


def test_atoms_section_uses_full_style_wrap_only(lib, small_system_no_water, tmp_path: Path):
    """Wrap-only convention (matches core topon writer): 7 columns, no image flags.

    LAMMPS handles bonds across periodic boundaries via min-image through the
    neighbor/ghost system. Writing image flags was the source of phantom bonds
    at BFM-merged crosslink endpoints whose chain-walks accumulated different
    image counts.
    """
    paths = lammps_writer.write_lammps(small_system_no_water, lib, tmp_path)
    text = paths["data"].read_text(encoding="utf-8")
    assert "Atoms  # full" in text
    # atom_style full WITHOUT image flags = 7 columns: id mol type q x y z
    after = text.split("Atoms  # full", 1)[1]
    rows = [l for l in after.splitlines() if l.strip() and not l.strip().startswith("#")]
    first = rows[0]
    body = first.split("#", 1)[0]
    cols = body.split()
    assert len(cols) == 7, f"expected 7 cols (full, wrap-only no image flags), got {len(cols)}"
    int(cols[0]); int(cols[1]); int(cols[2])              # id, mol, type
    float(cols[3]); float(cols[4]); float(cols[5]); float(cols[6])  # q, x, y, z


def test_pair_coeff_includes_self_and_cross_pairs(lib, small_system_no_water, tmp_path: Path):
    paths = lammps_writer.write_lammps(small_system_no_water, lib, tmp_path)
    text = paths["settings"].read_text(encoding="utf-8")
    n_types = len(small_system_no_water.bead_types_in_use())
    pair_lines = [l for l in text.splitlines() if l.startswith("pair_coeff")]
    expected = n_types * (n_types + 1) // 2
    assert len(pair_lines) == expected


def test_bond_types_include_constraint_and_real(lib, small_system_no_water, tmp_path: Path):
    paths = lammps_writer.write_lammps(small_system_no_water, lib, tmp_path)
    text = paths["settings"].read_text(encoding="utf-8")
    bond_coeff_lines = [l for l in text.splitlines() if l.startswith("bond_coeff")]
    # Must include at least one constraint line (TYR ring) and one regular bond
    assert any("constraint" in l for l in bond_coeff_lines)
    assert any("funct 1" in l for l in bond_coeff_lines)


def test_stage_scripts_have_expected_physics(lib, small_system_no_water, tmp_path: Path):
    paths = lammps_writer.write_lammps(small_system_no_water, lib, tmp_path)
    s1 = paths["stage1"].read_text(encoding="utf-8")
    s2 = paths["stage2"].read_text(encoding="utf-8")
    s3 = paths["stage3"].read_text(encoding="utf-8")

    # Stage 1 = soft-push overlap removal, ramped CG min + brief nve/limit.
    for kw in ("pair_style      soft", "fix             soft_push all adapt",
               "min_style       cg", "fix             nve_limit all nve/limit",
               "write_data      system_after_soft.data",
               f"read_data       ../protein_network.data"):
        assert kw in s1, f"stage1 missing: {kw!r}"

    # Stage 2 = epsilon ramp 0.001 -> 1.0 with fix adapt under nve/limit.
    for kw in ("read_data       system_after_soft.data",
               "variable        scale equal ramp(0.001,1.0)",
               "fix             1 all adapt 1 pair lj/cut/coul/cut epsilon * * v_scale",
               "fix             fxnve all nve/limit 0.1",
               "write_data      system_ramped.data"):
        assert kw in s2, f"stage2 missing: {kw!r}"

    # Stage 3 = tight CG min + brief NVT/NPT at 310 K -> system_equilibrated.data
    for kw in ("read_data       system_ramped.data",
               "minimize        1.0e-6 1.0e-8",
               "fix             1 all nvt temp ${T} ${T}",
               "fix             1 all npt temp ${T} ${T}",
               "write_data      ../system_equilibrated.data"):
        assert kw in s3, f"stage3 missing: {kw!r}"

    # All three stages use the MARTINI 3 + RF-approx style headers
    for s, name in ((s1, "stage1"), (s2, "stage2"), (s3, "stage3")):
        assert "atom_style      full" in s, f"{name} missing atom_style full"
        assert "angle_style     cosine/squared" in s, f"{name} missing cosine/squared angles"
        assert "dielectric      15.0" in s, f"{name} missing dielectric 15.0 (RF approx)"


def test_groups_file_has_protein_and_water_sections(lib, tmp_path: Path):
    topo = bfm.generate_topology(n_chains=2, n_repeats=3, segs_per_block=2,
                                 equil_steps=0, seed=11, verbose=False)
    seq3 = sequence.build_full_sequence("GGRPSDSYGAPGGGN", 3)
    sys_ = builder.build_protein_system(topo["snapshots"][0], seq3, lib,
                                        block_seq="GGRPSDSYGAPGGGN", seed=42)
    water.pack_water(sys_, lib, density_w_per_nm3=4.0, seed=99)
    paths = lammps_writer.write_lammps(sys_, lib, tmp_path)
    text = paths["groups"].read_text(encoding="utf-8")
    assert "group protein molecule" in text
    assert "group water molecule" in text
