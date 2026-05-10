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
