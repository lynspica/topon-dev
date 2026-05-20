"""Tests for the LAMMPS writer."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from topon.protein_network import bfm, builder, lammps_writer, sequence, water
from topon.protein_network.martini_ff import MartiniLibrary
from topon.protein_network.system import Bead, Bond, Constraint, System


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


def test_atoms_section_emits_10_columns_with_image_flags(
        lib, small_system_no_water, tmp_path: Path):
    """Image-flag convention (2026-05; replaces prior wrap-only writer).

    Each Atoms row is 10 columns:
        atom_id  mol_id  type_id  charge  x  y  z  ix  iy  iz

    Image flags are computed by a priority-weighted MST over the bond
    graph (see ``lammps_writer._kruskal_image_flags_and_drop``): non-
    crosslink bonds (priority 0) are processed before crosslinks
    (priority 1) so every tree-edge bond is minimum-image. Winding-
    cycle back-edges that cannot be made MIC are dropped; with the
    priority key, a hard assertion guarantees that real funct-int
    non-crosslink bonds (backbone, sidechain) can never drop — only
    crosslinks, plus rare TYR-ring constraints whose sidechain
    straddles a box boundary.

    Note: this writer's earlier 7-column "wrap-only" convention was
    adopted after v33-v38 per-chain image walks produced phantom 240
    A bonds at crosslinks. The current MST-over-global-bond-graph
    approach is structurally different from the v33-v38 anti-pattern:
    tree edges are MIC by construction regardless of which chain they
    belong to.
    """
    paths = lammps_writer.write_lammps(small_system_no_water, lib, tmp_path)
    text = paths["data"].read_text(encoding="utf-8")
    assert "Atoms  # full" in text
    after = text.split("Atoms  # full", 1)[1]
    rows = [l for l in after.splitlines() if l.strip() and not l.strip().startswith("#")]
    first = rows[0]
    body = first.split("#", 1)[0]
    cols = body.split()
    assert len(cols) == 10, (
        f"expected 10 cols (full + image flags), got {len(cols)}"
    )
    int(cols[0]); int(cols[1]); int(cols[2])                 # id, mol, type
    float(cols[3]); float(cols[4]); float(cols[5]); float(cols[6])  # q, x, y, z
    int(cols[7]); int(cols[8]); int(cols[9])                 # ix, iy, iz


def test_bonds_are_minimum_image_under_emitted_image_flags(
        lib, small_system_no_water, tmp_path: Path):
    """After the MST image-flag pass, every emitted bond should be
    minimum-image when its endpoints' wrapped positions are unwrapped
    via the file's image flags. This is the structural invariant of the
    fix: it guarantees LAMMPS' parallel ghost-shell construction has
    consistent bookkeeping for bonded interactions across MPI ranks.
    """
    import re

    paths = lammps_writer.write_lammps(small_system_no_water, lib, tmp_path)
    text = paths["data"].read_text(encoding="utf-8")

    # Box bounds
    Lx = Ly = Lz = None
    for line in text.splitlines():
        m = re.match(r"\s*([\-\d.eE]+)\s+([\-\d.eE]+)\s+xlo xhi", line)
        if m: Lx = float(m.group(2)) - float(m.group(1))
        m = re.match(r"\s*([\-\d.eE]+)\s+([\-\d.eE]+)\s+ylo yhi", line)
        if m: Ly = float(m.group(2)) - float(m.group(1))
        m = re.match(r"\s*([\-\d.eE]+)\s+([\-\d.eE]+)\s+zlo zhi", line)
        if m: Lz = float(m.group(2)) - float(m.group(1))
    assert None not in (Lx, Ly, Lz)

    atoms = {}     # id -> (x, y, z, ix, iy, iz)
    bonds = []     # (a, b)
    section = None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("Atoms"):
            section = "atoms"; continue
        if s.startswith("Bonds"):
            section = "bonds"; continue
        if s.startswith(("Velocities", "Angles", "Dihedrals", "Impropers",
                         "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
                         "Dihedral Coeffs", "Improper Coeffs", "Masses")):
            section = None
        if not s or s.startswith("#"):
            continue
        toks = s.split("#", 1)[0].split()
        if section == "atoms" and len(toks) >= 10:
            atoms[int(toks[0])] = (
                float(toks[4]), float(toks[5]), float(toks[6]),
                int(toks[7]), int(toks[8]), int(toks[9]),
            )
        elif section == "bonds" and len(toks) >= 4:
            bonds.append((int(toks[2]), int(toks[3])))

    # Every bond's unwrapped length, computed via the file's image flags,
    # must equal its wrapped minimum-image length within numerical noise.
    for (a, b) in bonds:
        xa, ya, za, ixa, iya, iza = atoms[a]
        xb, yb, zb, ixb, iyb, izb = atoms[b]
        dx_unwrap = (xa + ixa*Lx) - (xb + ixb*Lx)
        dy_unwrap = (ya + iya*Ly) - (yb + iyb*Ly)
        dz_unwrap = (za + iza*Lz) - (zb + izb*Lz)
        # Wrapped MIC reference
        dx_mic = xa - xb;  dx_mic -= Lx * round(dx_mic / Lx)
        dy_mic = ya - yb;  dy_mic -= Ly * round(dy_mic / Ly)
        dz_mic = za - zb;  dz_mic -= Lz * round(dz_mic / Lz)
        # The two displacements should agree (each component within 1 Å —
        # MARTINI bonds + 0.05 Å perturbation can't push the difference
        # past 0 since both are computed deterministically from the same
        # atom positions).
        assert abs(dx_unwrap - dx_mic) < 1e-6, (
            f"bond {a}-{b}: unwrapped Δx = {dx_unwrap}, MIC Δx = {dx_mic}"
        )
        assert abs(dy_unwrap - dy_mic) < 1e-6
        assert abs(dz_unwrap - dz_mic) < 1e-6


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


# --------------------------------------------------------------------------
# Priority-MST corner cases (2026-05-19): real-bond assertion & constraint
# drop. See lammps_writer._kruskal_image_flags_and_drop docstring and
# docs/JOURNAL.md 2026-05-19 (later) for the rationale.
# --------------------------------------------------------------------------

def _bead(aid: int, x: float, y: float, z: float, bead_type: str = "P2") -> Bead:
    """Minimal Bead for image-flag tests; mass/charge/molecule/residue
    fields are placeholders — only atom_id, bead_type, and position
    are exercised by the writer's MST."""
    return Bead(
        atom_id=aid, bead_type=bead_type, molecule_id=1, residue_idx=1,
        residue_name="GLY", atom_name="BB", charge=0.0, mass=72.0,
        position=(x, y, z),
    )


def test_corrupted_real_bond_winding_cycle_triggers_assertion(
        lib, tmp_path: Path):
    """Three real (funct=1, non-crosslink, non-constraint) bonds forming
    a triangle that straddles a box face create a winding cycle. With
    the priority-MST sort, the longest of those bonds (a real bond)
    would normally be the back-edge — the writer's hard assertion must
    fire instead of silently corrupting the topology.

    This guards the "BB-BB drop" regression: prior to 2026-05-19 (later)
    the writer used a length-only sort and silently dropped 16 real
    backbone bonds on the resilin_martini_highpro v2 system. The
    priority-MST + restored assertion together guarantee that scenario
    can never recur without a loud failure.
    """
    sys_ = System(box_dims_ang=(10.0, 10.0, 10.0))
    # Three atoms straddling x=0 / x=10 boundary. Bond 1-2 is short via
    # wrap (~0.2 A); bonds 2-3 and 3-1 are ~5 A each within the box.
    # Triangle traversal accumulates one image-x crossing → winding=1.
    sys_.beads.extend([
        _bead(1, 0.1, 5.0, 5.0),
        _bead(2, 9.9, 5.0, 5.0),
        _bead(3, 5.0, 5.0, 5.0),
    ])
    sys_.bonds.extend([
        Bond(a=1, b=2, funct=1, length_nm=0.36, k_kj=8000.0),
        Bond(a=2, b=3, funct=1, length_nm=0.36, k_kj=8000.0),
        Bond(a=3, b=1, funct=1, length_nm=0.36, k_kj=8000.0),
    ])

    with pytest.raises(AssertionError, match="real .*non-crosslink"):
        lammps_writer.write_lammps(
            sys_, lib, tmp_path, coord_perturbation_ang=0.0,
        )


def test_winding_constraint_triangle_drops_silently(lib, tmp_path: Path):
    """A constraint-only triangle (e.g. an isolated TYR-ring fragment
    straddling a box face) closes a winding cycle that no image-flag
    assignment can resolve. With the priority-MST, the constraint that
    is the longest in that cycle is silently dropped (priority 0 ==
    same as real bonds, but `funct=='constraint'` routes the drop to
    ``n_constraint_dropped`` instead of the assertion path).

    The writer must succeed and the drop must be reported under the
    "constraints" category, not the assertion.
    """
    sys_ = System(box_dims_ang=(10.0, 10.0, 10.0))
    sys_.beads.extend([
        _bead(1, 0.1, 5.0, 5.0, bead_type="TC4"),
        _bead(2, 9.9, 5.0, 5.0, bead_type="TC5"),
        _bead(3, 5.0, 5.0, 5.0, bead_type="TC5"),
    ])
    sys_.constraints.extend([
        Constraint(a=1, b=2, length_nm=0.30),
        Constraint(a=2, b=3, length_nm=0.30),
        Constraint(a=3, b=1, length_nm=0.30),
    ])

    paths = lammps_writer.write_lammps(
        sys_, lib, tmp_path, coord_perturbation_ang=0.0,
    )
    # 3 constraints in, 1 winding-cycle back-edge dropped → 2 emitted.
    header = paths["data"].read_text(encoding="utf-8").splitlines()[:10]
    n_bonds_line = next(line for line in header if line.strip().endswith("bonds"))
    n_bonds = int(n_bonds_line.split()[0])
    assert n_bonds == 2, (
        f"expected 2 bonds emitted (3 constraints - 1 winding drop), got {n_bonds}"
    )


def test_priority_mst_keeps_real_bond_when_crosslink_in_cycle(
        lib, tmp_path: Path):
    """A real bond + a crosslink forming a 2-edge cycle (parallel bonds
    between the same atom pair) should keep the real bond as the tree
    edge regardless of length. The priority key promotes the real bond
    to priority 0; the crosslink at priority 1 becomes the back-edge.

    Verifies the structural invariant the priority key was added for.
    """
    sys_ = System(box_dims_ang=(10.0, 10.0, 10.0))
    sys_.beads.extend([
        _bead(1, 2.0, 5.0, 5.0),
        _bead(2, 4.0, 5.0, 5.0),
    ])
    # Real backbone bond, long actual length (2 A).
    sys_.bonds.append(Bond(a=1, b=2, funct=1, length_nm=0.36, k_kj=8000.0))
    # Crosslink at the same pair — BFM-merged TYR scenario where MIC
    # distance is the same as the real bond (no merging cheat here for
    # the test, just a parallel bond). Even though both are length 2 A,
    # the priority key forces the crosslink to lose.
    sys_.bonds.append(Bond(
        a=1, b=2, funct=1, length_nm=0.27, k_kj=8000.0, is_crosslink=True,
    ))

    paths = lammps_writer.write_lammps(
        sys_, lib, tmp_path, coord_perturbation_ang=0.0,
    )
    # The 2-edge cycle (atoms 1-2 connected twice) has no winding (both
    # bonds have the same MIC delta), so neither edge is dropped — both
    # are emitted. The assertion path (real-bond drop) must not fire.
    n_bonds_line = next(
        line for line in paths["data"].read_text(encoding="utf-8").splitlines()
        if line.strip().endswith("bonds")
    )
    assert int(n_bonds_line.split()[0]) == 2
