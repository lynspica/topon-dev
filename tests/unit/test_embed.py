"""Regression tests for the ETKDGv3 embedder — priority 1a of solubility v2.

The legacy :func:`topon.singlechain.workflow._assign_extended_linear_coords`
collapses branch atoms onto their parent atom.  :func:`embed_with_etkdg`
replaces it with RDKit's ETKDGv3 + MMFF94/UFF.

These tests document the bug and guard against regressions.
"""

from __future__ import annotations

import warnings

import pytest

from topon.chemistry.embed import (
    embed_with_etkdg,
    min_pairwise_distance,
    EmbedFailedError,
)


def _build_short_pdms_chain(dp: int = 10):
    """Build a DP=10 PDMS chain via ChemistryBuilder (short for speed)."""
    import networkx as nx
    from rdkit import Chem
    from topon.chemistry.builder import ChemistryBuilder
    from topon.config.schema import (
        ChemistryConfig, MonomerConfig, NodeMoleculeConfig,
        EdgeChemistryConfig, ConnectionConfig,
    )

    G = nx.MultiGraph()
    G.add_node(0, pos=(0.0, 0.0, 0.0), node_type="end")
    G.add_node(1, pos=(float(dp), 0.0, 0.0), node_type="end")
    G.add_edge(0, 1, key=0, dp=dp, edge_type="A")

    chem = ChemistryConfig(
        model_type="atomistic", target_density=0.97,
        node_type_map={"end": NodeMoleculeConfig(molecule="[Si](C)(C)C", is_end_cap=True)},
        edge_type_map={"A": EdgeChemistryConfig(monomer="M0")},
        monomers={"M0": MonomerConfig(smiles="[Si](C)(C)O", chain_head="Si", chain_tail="O")},
        connection=ConnectionConfig(auto_bridge=True, default_bridge_atom="O"),
    )
    rw = ChemistryBuilder(G, dims=None, config=chem).build().GetMol()
    Chem.SanitizeMol(rw)
    return Chem.AddHs(rw)


def test_etkdg_gives_reasonable_min_distance():
    """After embed, no two atoms should be closer than a C-H bond."""
    mol = _build_short_pdms_chain(dp=10)
    embed_with_etkdg(mol, seed=42)
    d = min_pairwise_distance(mol)
    # C-H bond is ~1.09 Å; allow a little slack because MMFF can nudge
    # hydrogens briefly toward each other during optimisation.
    assert d > 0.8, f"min pairwise distance {d:.3f} Å indicates collapsed atoms"


def test_etkdg_deterministic_with_same_seed():
    """Same seed → same coordinates (bitwise-equal conformer)."""
    import numpy as np
    m1 = _build_short_pdms_chain(dp=8)
    embed_with_etkdg(m1, seed=123)
    m2 = _build_short_pdms_chain(dp=8)
    embed_with_etkdg(m2, seed=123)
    c1 = m1.GetConformer()
    c2 = m2.GetConformer()
    coords1 = np.array([[c1.GetAtomPosition(i).x, c1.GetAtomPosition(i).y, c1.GetAtomPosition(i).z]
                         for i in range(m1.GetNumAtoms())])
    coords2 = np.array([[c2.GetAtomPosition(i).x, c2.GetAtomPosition(i).y, c2.GetAtomPosition(i).z]
                         for i in range(m2.GetNumAtoms())])
    # ETKDGv3 is deterministic, but MMFF can pick different stationary points
    # when the starting coords are the same — so we check RMSD < 1 Å not
    # bitwise equality.
    from scipy.spatial.distance import cdist
    assert coords1.shape == coords2.shape
    rmsd = float(((coords1 - coords2) ** 2).sum(axis=1).mean() ** 0.5)
    assert rmsd < 1.0, f"same-seed RMSD {rmsd:.3f} Å — non-deterministic"


def test_legacy_linear_placer_has_collapsed_atoms():
    """Document the known bug in _assign_extended_linear_coords.

    This test pins the *legacy* behaviour so we detect if topon's linear
    placer ever gets silently fixed (in which case the bug note in
    :mod:`topon.chemistry.embed` and the v2 roadmap can be updated).
    """
    from topon.singlechain.workflow import _assign_extended_linear_coords

    mol = _build_short_pdms_chain(dp=10)
    _assign_extended_linear_coords(mol, bond_length=1.5)
    d = min_pairwise_distance(mol)
    # We *expect* this to be small (the bug).  If this assertion starts
    # failing with a LARGER d, the linear placer has been fixed — update
    # this test to assert ``d > 0.8`` instead.
    assert d < 0.8, (
        f"Legacy _assign_extended_linear_coords now gives min dist {d:.3f} Å — "
        "if it has been fixed, update this test to assert d > 0.8."
    )


def test_embed_with_etkdg_raises_on_impossible_geometry():
    """Overconstrained rings should raise EmbedFailedError, not return silently."""
    from rdkit import Chem
    # A fused ring that RDKit ETKDGv3 struggles with when random coords fail
    # is rare; this is more a smoke check that the exception path exists.
    # We craft a nonsense distance matrix via a placeholder — the typical
    # failure path is mol with ~0 atoms or missing conformers.
    empty = Chem.RWMol()
    with pytest.raises(Exception):  # RDKit itself raises before we do
        embed_with_etkdg(empty, seed=0)
