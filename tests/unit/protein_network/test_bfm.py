"""Tests for the BFM 6-neighbour cubic-lattice topology fork.

Covers determinism (single-seed reproducibility), gel-point detection, snapshot
schema, JSON round-trip via topology_io, and the 6-neighbour adjacency rule
(which the snapshot-keyed flat-index encoding depends on).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from topon.protein_network import bfm, sequence, topology_io


def test_lattice_size_rounds_to_odd_min_nine():
    assert bfm.compute_lattice_size(1, 1) == 9
    L = bfm.compute_lattice_size(16, 25, target_packing=0.45)
    assert L >= 9
    assert L % 2 == 1


def test_flat_index_round_trip_consistency():
    Nx = Ny = Nz = 9
    for x in (0, 4, 8):
        for y in (0, 4, 8):
            for z in (0, 4, 8):
                flat = bfm.xyz_to_flat(x, y, z, Nx, Ny)
                assert bfm.flat_to_xyz(flat, Nx, Ny) == (x, y, z)


def test_six_neighbours_periodic():
    Nx = Ny = Nz = 5
    nbs = bfm.get_neighbors(0, 0, 0, Nx, Ny, Nz)
    assert len(nbs) == 6
    # corner site wraps periodically
    assert bfm.xyz_to_flat(4, 0, 0, Nx, Ny) in nbs  # x-1 wrap
    assert bfm.xyz_to_flat(0, 4, 0, Nx, Ny) in nbs  # y-1 wrap
    assert bfm.xyz_to_flat(0, 0, 4, Nx, Ny) in nbs  # z-1 wrap


def test_unionfind_basic():
    uf = bfm.UnionFind(5)
    assert uf.n_components() == 5
    uf.union(0, 1)
    uf.union(1, 2)
    assert uf.n_components() == 3
    assert uf.find(0) == uf.find(2)


def test_topology_deterministic_with_seed():
    """Two runs with the same seed must produce bitwise-identical snapshots."""
    a = bfm.generate_topology(
        n_chains=4, n_repeats=4, segs_per_block=2,
        equil_steps=200, seed=12345, verbose=False,
    )
    b = bfm.generate_topology(
        n_chains=4, n_repeats=4, segs_per_block=2,
        equil_steps=200, seed=12345, verbose=False,
    )
    assert a["config"] == b["config"]
    assert a["snapshots"] == b["snapshots"]


def test_topology_different_seeds_differ():
    a = bfm.generate_topology(
        n_chains=4, n_repeats=4, segs_per_block=2, equil_steps=0,
        seed=1, verbose=False,
    )
    b = bfm.generate_topology(
        n_chains=4, n_repeats=4, segs_per_block=2, equil_steps=0,
        seed=2, verbose=False,
    )
    assert a["snapshots"][0]["chains"] != b["snapshots"][0]["chains"]


def test_snapshot_keys_match_topro_format():
    topo = bfm.generate_topology(
        n_chains=4, n_repeats=4, segs_per_block=2,
        equil_steps=100, seed=7, verbose=False,
    )
    assert topo["snapshots"], "expected at least one snapshot"
    for snap in topo["snapshots"]:
        assert set(snap) == {"label", "conv", "chains", "crosslinker_positions",
                             "reactions", "Nx", "Ny", "Nz"}
        assert isinstance(snap["label"], str)
        assert isinstance(snap["conv"], float)
        assert isinstance(snap["chains"], list)
        # reactions is list-of-list-of-list-of-int (NOT tuples)
        for rxn in snap["reactions"]:
            assert isinstance(rxn, list) and len(rxn) == 2
            for half in rxn:
                assert isinstance(half, list) and len(half) == 2
                assert all(isinstance(v, int) for v in half)


def test_gel_point_label_present_in_dense_system():
    # Larger system + more equilibration usually gels
    topo = bfm.generate_topology(
        n_chains=8, n_repeats=8, segs_per_block=2,
        target_packing=0.45, equil_steps=2000, seed=42, verbose=False,
    )
    labels = [s["label"] for s in topo["snapshots"]]
    # Either it gelled (gel_point present) or it didn't (no_gel present); both
    # are valid outcomes for a small test system. Confirm the LABEL set is sane.
    assert any(label in {"gel_point", "no_gel"} for label in labels)
    assert all(label.startswith(("gel_point", "post_gel_", "no_gel")) for label in labels)


def test_topology_json_round_trip(tmp_path: Path):
    topo = bfm.generate_topology(
        n_chains=4, n_repeats=4, segs_per_block=2,
        equil_steps=100, seed=7, verbose=False,
    )
    path = tmp_path / "topo.json"
    topology_io.save_topology(topo, str(path))
    loaded = topology_io.load_topology(str(path))
    # config and snapshots survive the JSON round-trip
    assert loaded["config"] == topo["config"]
    assert loaded["snapshots"] == topo["snapshots"]
    # File envelope has the version and timestamp
    raw = json.loads(path.read_text())
    assert raw["version"] == topology_io.TOPOLOGY_VERSION
    assert "created" in raw


def test_get_snapshot_by_label_and_index():
    topo = bfm.generate_topology(
        n_chains=4, n_repeats=4, segs_per_block=2,
        equil_steps=100, seed=7, verbose=False,
    )
    first_label = topo["snapshots"][0]["label"]
    a = topology_io.get_snapshot(topo, first_label)
    b = topology_io.get_snapshot(topo, 0)
    assert a == b


# Sequence module tests =====================================================

def test_build_full_sequence_resilin():
    seq = sequence.build_full_sequence("GGRPSDSYGAPGGGN", 2)
    # Should be 30 residues
    assert len(seq) == 30
    assert seq[0] == "GLY" and seq[1] == "GLY"
    assert seq[7] == "TYR"  # position 7 in first repeat
    assert seq[22] == "TYR"  # position 7 in second repeat (15 + 7)


def test_build_full_sequence_rejects_unknown_letter():
    with pytest.raises(ValueError):
        sequence.build_full_sequence("ABZ", 1)


def test_get_node_residue_mapping_anchors_at_ends():
    n_repeats, segs = 3, 2
    n_nodes = n_repeats * segs + 1  # 7 nodes
    block_size = 15
    n_residues = n_repeats * block_size
    mapping = sequence.get_node_residue_mapping(
        n_repeats=n_repeats, segs_per_block=segs, block_size=block_size, y_in_block=7,
    )
    assert mapping[0] == 0
    assert mapping[n_nodes - 1] == n_residues - 1


def test_get_node_residue_mapping_y_anchors_at_tyr():
    """Y-node residues must land on the crosslinker letter index."""
    block = "GGRPSDSYGAPGGGN"  # Y at index 7
    n_repeats = 3
    mapping = sequence.get_node_residue_mapping(
        n_repeats=n_repeats, segs_per_block=2, block_seq=block,
    )
    # Y-nodes at chain indices 1, 3, 5 -> residues 7, 22, 37
    assert mapping[1] == 7
    assert mapping[3] == 22
    assert mapping[5] == 37


def test_get_tyr_node_indices_consistent_with_y_positions():
    n_repeats, segs = 6, 2
    via_seq = sequence.get_tyr_node_indices(n_repeats, segs)
    via_bfm = set(bfm.get_y_positions(n_repeats, segs))
    assert via_seq == via_bfm


# --------------------------------------------------------------------------
# winding-safe crosslinking (2026-05-19): a candidate crosslink is rejected
# if it would close a cycle in the chain-plus-crosslink bond graph whose
# accumulated lattice image-delta is non-zero around the periodic box.
# Goal: the projected geometry has zero winding-cycle drops at write time.
# --------------------------------------------------------------------------

def test_winding_safe_produces_zero_writer_drops():
    """End-to-end: a `winding_safe` BFM topology projected via the
    builder and written via the LAMMPS writer must produce zero
    winding-cycle drops. The new method's whole point is to filter the
    candidate set at reaction time so the writer never has to drop.

    Uses a system small enough to run fast (~3s) but large enough that
    a few candidates trigger the winding rejection path (verified via
    the ``n_rejected_winding`` field in the snapshot).
    """
    from topon.protein_network import builder, lammps_writer, sequence
    from topon.protein_network.martini_ff import MartiniLibrary
    import tempfile

    topo = bfm.generate_topology(
        n_chains=20, n_repeats=12, segs_per_block=2,
        target_packing=0.45, equil_steps=50_000,
        crosslink_method="winding_safe", seed=99, verbose=False,
    )
    snap = topo["snapshots"][0]  # gel_point if reached, else no_gel
    assert snap["n_rejected_winding"] >= 1, (
        f"expected at least one winding rejection for seed=99 (got "
        f"{snap['n_rejected_winding']}); this canary test relies on the "
        f"realisation actually exercising the rejection path."
    )

    lib = MartiniLibrary.from_package_data()
    seq3 = sequence.build_full_sequence("GGRPSDSYGAPGGGN", 12)
    sys_ = builder.build_protein_system(
        snap, seq3, lib, block_seq="GGRPSDSYGAPGGGN", seed=99,
    )

    # The writer prints a drop report to stdout if it dropped anything,
    # AND the kept_bonds count would be less than the total. Easier
    # check: parse the header.
    with tempfile.TemporaryDirectory() as tmp:
        paths = lammps_writer.write_lammps(sys_, lib, tmp)
        data_text = paths["data"].read_text(encoding="utf-8")

    # Total bonds in the System (input to writer)
    n_input = len(sys_.bonds) + len(sys_.constraints)
    # Emitted bonds count (writer header)
    n_bonds_line = next(
        l for l in data_text.splitlines() if l.strip().endswith("bonds")
    )
    n_emitted = int(n_bonds_line.split()[0])
    assert n_emitted == n_input, (
        f"winding_safe should produce zero drops, but writer emitted "
        f"{n_emitted} of {n_input} bonds ({n_input - n_emitted} dropped). "
        f"Snapshot conversion was {snap['conv']:.4f} with "
        f"{snap['n_rejected_winding']} winding-rejections."
    )


def test_winding_safe_matches_adjacent_on_no_winding_seed():
    """For a small / sparse system where no candidates wind, the
    `winding_safe` method must produce the same accepted crosslink set
    as the default `adjacent` method (modulo nothing — same shuffle
    order, same rng seed). This guards against unintended divergence
    in the easy case.
    """
    common_args = dict(
        n_chains=6, n_repeats=4, segs_per_block=2,
        target_packing=0.45, equil_steps=2000, seed=3, verbose=False,
    )
    topo_adj = bfm.generate_topology(**common_args, crosslink_method="adjacent")
    topo_ws = bfm.generate_topology(**common_args, crosslink_method="winding_safe")

    # No rejections for this small/sparse setup
    last_ws = topo_ws["snapshots"][-1]
    if last_ws.get("n_rejected_winding", 0) == 0:
        rxns_adj = topo_adj["snapshots"][-1]["reactions"]
        rxns_ws = topo_ws["snapshots"][-1]["reactions"]
        assert rxns_adj == rxns_ws, (
            "with zero winding-rejections, winding_safe should produce "
            "an identical reaction list to the adjacent method"
        )


def test_crosslink_method_none_emits_single_uncrosslinked_snapshot():
    """`crosslink_method="none"` must skip crosslinking entirely and emit one
    conv=0 snapshot labelled `uncrosslinked` with an empty reaction list."""
    topo = bfm.generate_topology(
        n_chains=4, n_repeats=6, segs_per_block=3, y_offset_in_block=1,
        equil_steps=20_000, crosslink_method="none", seed=501, verbose=False,
    )
    snaps = topo["snapshots"]
    assert len(snaps) == 1
    snap = snaps[0]
    assert snap["label"] == "uncrosslinked"
    assert snap["conv"] == 0.0
    assert snap["reactions"] == []
    assert topo["config"]["crosslink_method"] == "none"
    # schema parity with the crosslinked methods
    for key in ("chains", "crosslinker_positions", "Nx", "Ny", "Nz"):
        assert key in snap


def test_crosslink_method_none_leaves_every_y_node_on_its_own_site():
    """The CHARMM builder infers a crosslink from two Y nodes sharing a lattice
    site (the BFM merges a reacted pair onto one site). An uncrosslinked
    snapshot must therefore have NO duplicated site at all -- otherwise the
    builder would silently stitch a CE2-CE2 bond, or leave two tyrosines
    superimposed at r~0.
    """
    topo = bfm.generate_topology(
        n_chains=4, n_repeats=6, segs_per_block=3, y_offset_in_block=1,
        equil_steps=20_000, crosslink_method="none", seed=7, verbose=False,
    )
    snap = topo["snapshots"][0]
    y_pos = set(snap["crosslinker_positions"])

    all_sites, y_sites = [], []
    for chain in snap["chains"]:
        for ni, flat in enumerate(chain):
            all_sites.append(flat)
            if ni in y_pos:
                y_sites.append(flat)

    assert len(all_sites) == len(set(all_sites)), "lattice site double-occupied"
    assert len(y_sites) == len(set(y_sites)), "two Y nodes share a site"
    assert len(y_sites) == 4 * 6, "expected one Y per repeat per chain"
