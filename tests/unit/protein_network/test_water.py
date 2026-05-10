"""Tests for the voxel-grid W water packer."""
from __future__ import annotations

import math

import numpy as np
import pytest

from topon.protein_network import bfm, builder, sequence, water
from topon.protein_network.martini_ff import MartiniLibrary
from topon.protein_network.system import Bead, System


@pytest.fixture(scope="module")
def lib() -> MartiniLibrary:
    return MartiniLibrary.from_package_data()


def test_pack_into_empty_box(lib):
    sys_ = System(box_dims_ang=(20.0, 20.0, 20.0))
    n = water.pack_water(sys_, lib, density_w_per_nm3=10.0, seed=42)
    # 8 nm^3 box at 10 W/nm^3 -> ~80 beads
    assert 60 <= n <= 100, f"expected ~80 W beads, got {n}"
    assert all(b.bead_type == "W" for b in sys_.beads)
    assert all(b.atom_name == "W" for b in sys_.beads)


def test_each_water_gets_unique_molecule_id(lib):
    sys_ = System(box_dims_ang=(15.0, 15.0, 15.0))
    n = water.pack_water(sys_, lib, density_w_per_nm3=10.0, seed=1)
    mol_ids = {b.molecule_id for b in sys_.beads}
    assert len(mol_ids) == n, "each water bead should be its own molecule"


def test_exclusion_zone_keeps_water_away(lib):
    # One protein bead at the box centre; pack water with a 5 A exclusion.
    sys_ = System(box_dims_ang=(20.0, 20.0, 20.0))
    sys_.beads.append(Bead(atom_id=1, bead_type="P2", molecule_id=1, residue_idx=1,
                           residue_name="X", atom_name="BB", charge=0.0,
                           mass=72.0, position=(10.0, 10.0, 10.0)))
    water.pack_water(sys_, lib, density_w_per_nm3=10.0, exclusion_radius_ang=5.0, seed=42)
    box = np.array([20.0, 20.0, 20.0])
    centre = np.array([10.0, 10.0, 10.0])
    for b in sys_.beads:
        if b.bead_type != "W":
            continue
        d = np.array(b.position) - centre
        d -= box * np.round(d / box)
        assert np.linalg.norm(d) >= 5.0 - 1e-6


def test_atom_ids_continuous_after_existing_protein(lib):
    sys_ = System(box_dims_ang=(15.0, 15.0, 15.0))
    sys_.beads.append(Bead(atom_id=1, bead_type="P2", molecule_id=1, residue_idx=1,
                           residue_name="X", atom_name="BB", charge=0.0,
                           mass=72.0, position=(7.5, 7.5, 7.5)))
    water.pack_water(sys_, lib, density_w_per_nm3=5.0, seed=0)
    atom_ids = [b.atom_id for b in sys_.beads]
    assert atom_ids == list(range(1, len(sys_.beads) + 1))


def test_water_pack_after_protein_build_is_deterministic(lib):
    topo = bfm.generate_topology(n_chains=2, n_repeats=3, segs_per_block=2,
                                 equil_steps=0, seed=11, verbose=False)
    snap = topo["snapshots"][0]
    seq3 = sequence.build_full_sequence("GGRPSDSYGAPGGGN", 3)

    sys_a = builder.build_protein_system(snap, seq3, lib, block_seq="GGRPSDSYGAPGGGN", seed=42)
    sys_b = builder.build_protein_system(snap, seq3, lib, block_seq="GGRPSDSYGAPGGGN", seed=42)
    n_a = water.pack_water(sys_a, lib, density_w_per_nm3=8.0, seed=99)
    n_b = water.pack_water(sys_b, lib, density_w_per_nm3=8.0, seed=99)
    assert n_a == n_b
    assert [b.position for b in sys_a.beads if b.bead_type == "W"] == \
           [b.position for b in sys_b.beads if b.bead_type == "W"]


def test_max_beads_caps_packing(lib):
    sys_ = System(box_dims_ang=(20.0, 20.0, 20.0))
    n = water.pack_water(sys_, lib, density_w_per_nm3=10.0, seed=42, max_beads=10)
    assert n == 10
