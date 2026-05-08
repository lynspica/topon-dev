"""
Unit tests for the simbox sub-system.

Tests cover:
  - Molecule creation (from_smiles), reactive-site detection
  - BoxPacker: atom counts, box sizing, determinism
  - assemble(): molecule_ids length, reactive-site global indices
"""

import pytest
import numpy as np


# ---------------------------------------------------------------------------
# Helpers — tiny fast molecules (no library conformers needed)
# ---------------------------------------------------------------------------

def _oxirane():
    """Oxirane (ethylene oxide): one epoxide, 7 atoms total with H."""
    from topon.simbox.molecule import Molecule
    return Molecule.from_smiles("oxirane", "C1CO1")


def _methylamine():
    """Methylamine: one primary amine, 7 atoms total with H."""
    from topon.simbox.molecule import Molecule
    return Molecule.from_smiles("methylamine", "CN")


def _propane():
    """Propane: inert, no reactive sites, 11 atoms with H."""
    from topon.simbox.molecule import Molecule
    return Molecule.from_smiles("propane", "CCC")


# ---------------------------------------------------------------------------
# Molecule tests
# ---------------------------------------------------------------------------

class TestMolecule:
    def test_from_smiles_atom_count(self):
        mol = _oxirane()
        # C1CO1 with explicit H: 2C + 1O + 4H = 7
        assert mol.num_atoms == 7

    def test_from_smiles_mw_positive(self):
        mol = _oxirane()
        assert mol.mw > 0

    def test_from_smiles_name(self):
        mol = _oxirane()
        assert mol.name == "oxirane"

    def test_from_smiles_has_3d_conformer(self):
        mol = _oxirane()
        coords = mol.get_coordinates()
        assert coords.shape == (mol.num_atoms, 3)
        # Not all zero
        assert not np.allclose(coords, 0)

    def test_get_centroid_shape(self):
        mol = _oxirane()
        c = mol.get_centroid()
        assert c.shape == (3,)

    def test_repr_contains_name(self):
        mol = _oxirane()
        assert "oxirane" in repr(mol)

    def test_invalid_smiles_raises(self):
        from topon.simbox.molecule import Molecule
        with pytest.raises(ValueError):
            Molecule.from_smiles("bad", "not_a_smiles!!!")


class TestReactiveSiteDetection:
    def test_epoxide_detected(self):
        mol = _oxirane()
        assert "epoxide" in mol.reactive_sites
        # Epoxide SMARTS [C]1[O][C]1 matches: 2C + 1O
        assert len(mol.reactive_sites["epoxide"]) == 3

    def test_no_epoxide_in_methylamine(self):
        mol = _methylamine()
        assert "epoxide" not in mol.reactive_sites

    def test_primary_amine_detected(self):
        mol = _methylamine()
        assert "primary_amine" in mol.reactive_sites
        assert len(mol.reactive_sites["primary_amine"]) >= 1

    def test_no_amine_in_propane(self):
        mol = _propane()
        assert "primary_amine" not in mol.reactive_sites
        assert "secondary_amine" not in mol.reactive_sites

    def test_inert_molecule_no_reactive_sites(self):
        mol = _propane()
        assert mol.reactive_sites == {}

    def test_custom_smarts_override(self):
        from topon.simbox.molecule import Molecule
        # Detect any carbon as "carbon_site"
        mol = Molecule.from_smiles(
            "propane_custom", "CCC",
            reactive_smarts={"carbon_site": "[#6]"}
        )
        assert "carbon_site" in mol.reactive_sites
        assert len(mol.reactive_sites["carbon_site"]) == 3  # 3 carbons


# ---------------------------------------------------------------------------
# BoxPacker tests
# ---------------------------------------------------------------------------

class TestBoxPacker:
    @pytest.fixture(scope="class")
    def packed_small(self):
        """Pack 3 oxirane + 2 methylamine with fixed seed."""
        from topon.simbox.packer import BoxPacker
        packer = BoxPacker(density=0.85, seed=0, max_attempts=500)
        mol_ox = _oxirane()
        mol_ma = _methylamine()
        return packer.pack([(mol_ox, 3), (mol_ma, 2)])

    def test_total_molecules(self, packed_small):
        assert packed_small.total_molecules == 5

    def test_total_atoms(self, packed_small):
        mol_ox = _oxirane()
        mol_ma = _methylamine()
        expected = 3 * mol_ox.num_atoms + 2 * mol_ma.num_atoms
        assert packed_small.total_atoms == expected

    def test_box_lengths_positive(self, packed_small):
        assert all(l > 0 for l in packed_small.box_lengths)

    def test_box_lengths_shape(self, packed_small):
        assert packed_small.box_lengths.shape == (3,)

    def test_placement_coordinates_shape(self, packed_small):
        mol_ox = _oxirane()
        for pm in packed_small.placements:
            assert pm.coordinates.shape[1] == 3

    def test_deterministic_with_same_seed(self):
        """Same seed → identical box lengths and coordinates."""
        from topon.simbox.packer import BoxPacker
        mol = _propane()
        packed_a = BoxPacker(density=0.85, seed=7).pack([(mol, 3)])
        packed_b = BoxPacker(density=0.85, seed=7).pack([(mol, 3)])
        np.testing.assert_allclose(packed_a.box_lengths, packed_b.box_lengths)
        for pa, pb in zip(packed_a.placements, packed_b.placements):
            np.testing.assert_allclose(pa.coordinates, pb.coordinates)

    def test_different_seeds_differ(self):
        """Different seeds → different placements."""
        from topon.simbox.packer import BoxPacker
        mol = _propane()
        packed_a = BoxPacker(density=0.85, seed=1).pack([(mol, 3)])
        packed_b = BoxPacker(density=0.85, seed=2).pack([(mol, 3)])
        coords_a = packed_a.placements[0].coordinates
        coords_b = packed_b.placements[0].coordinates
        assert not np.allclose(coords_a, coords_b)

    def test_box_size_scales_with_density(self):
        """Lower density → larger box."""
        from topon.simbox.packer import BoxPacker
        mol = _propane()
        p_hi = BoxPacker(density=1.5, seed=0).pack([(mol, 5)])
        p_lo = BoxPacker(density=0.5, seed=0).pack([(mol, 5)])
        assert p_lo.box_lengths[0] > p_hi.box_lengths[0]

    def test_min_dist_respected(self):
        """No pair of atoms from different placements closer than min_dist."""
        from topon.simbox.packer import BoxPacker
        min_dist = 2.0
        mol = _propane()
        packed = BoxPacker(density=0.5, seed=99, min_dist=min_dist).pack([(mol, 5)])
        # Check inter-molecular distances only (intra-molecular bonds are shorter)
        box = packed.box_lengths
        placements = packed.placements
        for i in range(len(placements)):
            for j in range(i + 1, len(placements)):
                ci = placements[i].coordinates
                cj = placements[j].coordinates
                for ai in ci:
                    for aj in cj:
                        delta = ai - aj
                        delta -= box * np.round(delta / box)
                        dist = np.linalg.norm(delta)
                        assert dist >= min_dist - 1e-6, (
                            f"Placement ({i},{j}) inter-molecular distance "
                            f"{dist:.3f} < min_dist {min_dist}"
                        )


# ---------------------------------------------------------------------------
# Assembler tests
# ---------------------------------------------------------------------------

class TestAssemble:
    @pytest.fixture(scope="class")
    def system(self):
        from topon.simbox.packer import BoxPacker
        from topon.simbox.system import assemble
        mol_ox = _oxirane()
        mol_ma = _methylamine()
        packed = BoxPacker(density=0.5, seed=1).pack([(mol_ox, 2), (mol_ma, 1)])
        return assemble(packed)

    def test_total_atom_count(self, system):
        mol_ox = _oxirane()
        mol_ma = _methylamine()
        expected = 2 * mol_ox.num_atoms + 1 * mol_ma.num_atoms
        assert system.mol.GetNumAtoms() == expected

    def test_molecule_ids_length(self, system):
        assert len(system.molecule_ids) == system.mol.GetNumAtoms()

    def test_molecule_ids_range(self, system):
        ids = system.molecule_ids
        assert min(ids) == 1
        assert max(ids) == system.num_molecules

    def test_num_molecules(self, system):
        assert system.num_molecules == 3  # 2 oxirane + 1 methylamine

    def test_species_names_count(self, system):
        assert len(system.species_names) == system.num_molecules

    def test_reactive_sites_global_indices_in_range(self, system):
        n_atoms = system.mol.GetNumAtoms()
        for site in system.reactive_sites:
            assert 0 <= site.global_atom_idx < n_atoms

    def test_reactive_sites_have_epoxide_and_amine(self, system):
        groups = {s.group_name for s in system.reactive_sites}
        assert "epoxide" in groups
        assert "primary_amine" in groups

    def test_box_lengths_preserved(self, system):
        assert system.box_lengths.shape == (3,)
        assert all(l > 0 for l in system.box_lengths)

    def test_molecule_id_contiguous_per_molecule(self, system):
        """Each molecule's atoms must all have the same molecule ID."""
        mol_ox = _oxirane()
        n = mol_ox.num_atoms
        # First molecule: atoms 0..n-1 should all be mol_id=1
        ids = system.molecule_ids
        assert all(ids[i] == 1 for i in range(n))
        # Second molecule: atoms n..2n-1 should all be mol_id=2
        assert all(ids[i] == 2 for i in range(n, 2 * n))
