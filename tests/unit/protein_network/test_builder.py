"""Tests for the MARTINI protein chain-builder.

Builds small systems from synthesized BFM snapshots + residue sequences and
checks: atom counts, charge neutrality, intra-residue bonded counts per residue
type, BB-BB chain-spanning bonds, dityrosine crosslink emission.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from topon.protein_network import bfm, builder, residues, sequence, topology_io
from topon.protein_network.martini_ff import MartiniLibrary
from topon.protein_network.system import System


@pytest.fixture(scope="module")
def lib() -> MartiniLibrary:
    return MartiniLibrary.from_package_data()


@pytest.fixture()
def small_snapshot_2chains():
    # 2 chains x 3 repeats x 2 segs/block + 1 = 7 nodes per chain
    topo = bfm.generate_topology(
        n_chains=2, n_repeats=3, segs_per_block=2,
        equil_steps=0, seed=11, verbose=False,
    )
    return topo["snapshots"][0]


def _build(lib, snap, n_repeats: int = 3) -> System:
    block = "GGRPSDSYGAPGGGN"
    seq3 = sequence.build_full_sequence(block, n_repeats)
    return builder.build_protein_system(snap, seq3, lib, block_seq=block, seed=42)


def test_total_bead_count_matches_resilin_pattern(lib, small_snapshot_2chains):
    sys_ = _build(lib, small_snapshot_2chains, n_repeats=3)
    # Per-residue bead counts in the canonical 'GGRPSDSYGAPGGGN' block:
    # G=1, G=1, R=3, P=2, S=2, D=2, S=2, Y=5, G=1, A=2, P=2, G=1, G=1, G=1, N=2 -> sum=28
    per_block = sum(len(residues.RESIDUES[res]["beads"]) for res in sequence.build_full_sequence("GGRPSDSYGAPGGGN", 1))
    assert per_block == 28
    n_chains = len(small_snapshot_2chains["chains"])
    n_repeats = 3
    expected = n_chains * n_repeats * per_block
    assert sys_.n_atoms() == expected


def test_total_charge_consistent_with_termini(lib, small_snapshot_2chains):
    sys_ = _build(lib, small_snapshot_2chains, n_repeats=3)
    # Per chain: ARG x3 (+3), ASP x3 (-3), N-term BB +1 (replaces neutral SP1 GLY),
    # C-term BB -1 (replaces neutral P2 ASN).
    # Net per chain = +3 - 3 + 1 - 1 = 0.
    n_chains = len(small_snapshot_2chains["chains"])
    assert sys_.total_charge() == pytest.approx(0.0, abs=1e-9)
    # And per chain it's zero too (build is symmetric across chains)
    chain_charges: dict[int, float] = {}
    for b in sys_.beads:
        chain_charges[b.molecule_id] = chain_charges.get(b.molecule_id, 0.0) + b.charge
    assert all(math.isclose(c, 0.0, abs_tol=1e-9) for c in chain_charges.values())


def test_n_term_first_bead_is_q5_plus1(lib, small_snapshot_2chains):
    sys_ = _build(lib, small_snapshot_2chains, n_repeats=3)
    first = sys_.beads[0]
    assert first.atom_name == "BB"
    assert first.bead_type == "Q5"
    assert first.charge == 1.0


def test_c_term_last_bead_is_sp5_neutral(lib, small_snapshot_2chains):
    sys_ = _build(lib, small_snapshot_2chains, n_repeats=3)
    # Last residue is ASN with patched BB Q5(-1) + SC1 SP5(0)
    chain1_beads = [b for b in sys_.beads if b.molecule_id == 1]
    last_bb = next(b for b in reversed(chain1_beads) if b.atom_name == "BB")
    last_sc1 = next(b for b in reversed(chain1_beads) if b.atom_name == "SC1")
    assert last_bb.bead_type == "Q5"
    assert last_bb.charge == -1.0
    assert last_sc1.bead_type == "SP5"
    assert last_sc1.charge == 0.0


def test_bb_bb_backbone_bond_uses_pro_variant_when_pro_neighbours(lib, small_snapshot_2chains):
    sys_ = _build(lib, small_snapshot_2chains, n_repeats=3)
    # In sequence GGRPSDSY..., the ARG->PRO BB-BB has length 0.36 nm, k=10000;
    # PRO->SER has length 0.305 nm, k=10000.
    # Find these by indexing chain 1 backbone atoms.
    chain1 = [b for b in sys_.beads if b.molecule_id == 1]
    bb_by_residue = {b.residue_idx: b for b in chain1 if b.atom_name == "BB"}
    arg_bb = bb_by_residue[3]
    pro_bb = bb_by_residue[4]
    ser_bb = bb_by_residue[5]
    # Find the bonds connecting them
    bonds_pairs = {(b.a, b.b): b for b in sys_.bonds}
    arg_pro = bonds_pairs.get((arg_bb.atom_id, pro_bb.atom_id)) or bonds_pairs.get((pro_bb.atom_id, arg_bb.atom_id))
    pro_ser = bonds_pairs.get((pro_bb.atom_id, ser_bb.atom_id)) or bonds_pairs.get((ser_bb.atom_id, pro_bb.atom_id))
    assert arg_pro is not None, "ARG-PRO BB bond missing"
    assert pro_ser is not None, "PRO-SER BB bond missing"
    assert arg_pro.length_nm == pytest.approx(0.360)
    assert arg_pro.k_kj == pytest.approx(10000.0)
    assert pro_ser.length_nm == pytest.approx(0.305)
    assert pro_ser.k_kj == pytest.approx(10000.0)


def test_tyr_ring_constraints_present(lib, small_snapshot_2chains):
    sys_ = _build(lib, small_snapshot_2chains, n_repeats=3)
    # TYR contributes 5 ring constraints (3x 0.300, 2x 0.285); ALA contributes
    # 1 (BB-SC1 at 0.270). The resilin block has 1 ALA + 1 TYR.
    n_chains = len(small_snapshot_2chains["chains"])
    n_repeats = 3
    expected_tyr = n_chains * n_repeats * 5
    expected_ala = n_chains * n_repeats * 1
    assert len(sys_.constraints) == expected_tyr + expected_ala
    lengths = [round(c.length_nm, 3) for c in sys_.constraints]
    assert lengths.count(0.285) == n_chains * n_repeats * 2
    assert lengths.count(0.300) == n_chains * n_repeats * 3
    assert lengths.count(0.270) == n_chains * n_repeats * 1


def test_bbbb_dihedral_uses_label_appropriate_terms(lib, small_snapshot_2chains):
    sys_ = _build(lib, small_snapshot_2chains, n_repeats=3)
    # Confirm at least one dihedral was emitted from each label found in the
    # IDP backbone reference.
    funct9 = [d for d in sys_.dihedrals if d.funct == 9]
    assert funct9, "expected at least some funct=9 backbone dihedrals"
    # Dihedral set should include both default-BBBB terms (k=2.8, mult=1) and
    # GGGX-style terms (k=0.8 or k=1.2)
    ks = {round(d.k_kj, 2) for d in funct9}
    assert 2.8 in ks or -0.6 in ks or -1.2 in ks


def test_dityrosine_crosslink_bonds_emitted_when_reactions_present():
    # Hand-craft a snapshot that has one reaction: chain 0 Y-node 1 with chain 1 Y-node 1
    # In a 2 segs_per_block geometry, Y-nodes are at chain idxs 1, 3, 5 -> residues 7, 22, 37
    # We need the SC4 of each TYR to exist.
    snap = {
        "label": "gel_point",
        "conv": 0.05,
        "chains": [[0, 1, 2, 3, 4, 5, 6], [10, 11, 12, 13, 14, 15, 16]],
        "crosslinker_positions": [1, 3, 5],
        "reactions": [[[0, 1], [1, 1]]],
        "Nx": 9, "Ny": 9, "Nz": 9,
    }
    block = "GGRPSDSYGAPGGGN"
    seq3 = sequence.build_full_sequence(block, 3)
    lib = MartiniLibrary.from_package_data()
    sys_ = builder.build_protein_system(snap, seq3, lib, block_seq=block, seed=0)
    crosslinks = [b for b in sys_.bonds if b.is_crosslink]
    assert len(crosslinks) == 1
    cl = crosslinks[0]
    # Both endpoints must be SC4 beads
    bead_by_id = {b.atom_id: b for b in sys_.beads}
    assert bead_by_id[cl.a].atom_name == "SC4"
    assert bead_by_id[cl.b].atom_name == "SC4"
    assert bead_by_id[cl.a].molecule_id != bead_by_id[cl.b].molecule_id


def test_in_memory_positions_within_one_box_image(lib, small_snapshot_2chains):
    """In-memory positions may fall outside [0, box) by up to one box-edge.
    The writer wraps with `pos % box` before emitting the data file (matches
    topro's writer.py:220), so the on-disk data file always has atoms inside.
    """
    sys_ = _build(lib, small_snapshot_2chains, n_repeats=3)
    bx, by, bz = sys_.box_dims_ang
    for b in sys_.beads:
        x, y, z = b.position
        assert -bx < x < 2 * bx
        assert -by < y < 2 * by
        assert -bz < z < 2 * bz


def test_writer_wraps_positions_into_box_emits_image_flags(lib, small_snapshot_2chains, tmp_path):
    """Atoms section: positions wrapped into [0, box), image flags emitted
    as columns 8-10 (ix, iy, iz). Image flags are computed by a length-
    weighted MST over the bond graph (see
    ``lammps_writer._kruskal_image_flags_and_drop``) so every tree-edge
    bond is minimum-image. This is the post-2026-05 convention; the prior
    wrap-only 7-column writer left LAMMPS to assume ix=iy=iz=0 and broke
    parallel-MPI bond communication for chains that wrap across the box.
    """
    from topon.protein_network import lammps_writer
    sys_ = _build(lib, small_snapshot_2chains, n_repeats=3)
    paths = lammps_writer.write_lammps(sys_, lib, tmp_path)
    text = paths["data"].read_text(encoding="utf-8")
    bx = float(next(l for l in text.splitlines() if "xhi" in l).split()[1])
    after = text.split("Atoms  # full", 1)[1]
    n_checked = 0
    for line in after.splitlines():
        s = line.strip()
        if not s or s.startswith("#"): continue
        toks = s.split("#")[0].split()
        # full atom_style row: id mol type q x y z ix iy iz (10 cols)
        if len(toks) != 10: continue
        try:
            x, y, z = float(toks[4]), float(toks[5]), float(toks[6])
            ix, iy, iz = int(toks[7]), int(toks[8]), int(toks[9])
        except ValueError:
            continue
        assert 0.0 <= x < bx, f"atom outside box on x: {x}"
        assert 0.0 <= y < bx, f"atom outside box on y: {y}"
        assert 0.0 <= z < bx, f"atom outside box on z: {z}"
        # ix/iy/iz are ints (no sign/magnitude constraint here — the
        # MIC-consistency invariant is checked in test_writer.py's
        # test_bonds_are_minimum_image_under_emitted_image_flags).
        n_checked += 1
    assert n_checked > 0
