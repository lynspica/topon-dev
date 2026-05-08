"""Regression tests for the ``node_type`` attribute fallback (v2 priority 1b).

Legacy callers (e.g. ``generate_matrix.py`` in DOW/solvent_effects/)
annotate graph nodes with ``type="end"`` instead of ``node_type="end"``.
Before v2, topon silently ignored the legacy attribute and fell through
to its default ``"A"`` node type, which maps to the Si-atom fallback —
leaking Si atoms into hydrocarbon polymers.

After v2, topon reads ``type`` as a fallback and emits a
DeprecationWarning.
"""

from __future__ import annotations

import warnings

import networkx as nx
import pytest
from rdkit import Chem

from topon.chemistry.builder import ChemistryBuilder
from topon.config.schema import (
    ChemistryConfig, MonomerConfig, NodeMoleculeConfig,
    EdgeChemistryConfig, ConnectionConfig,
)


def _build_butyl_with_legacy_type_attr():
    """Build a DP=4 polyisobutylene chain using ``type="end"`` (legacy)."""
    G = nx.MultiGraph()
    # NOTE: *only* 'type', no 'node_type' — this is the legacy pattern
    G.add_node(0, pos=(0.0, 0.0, 0.0), type="end")
    G.add_node(1, pos=(4.0, 0.0, 0.0), type="end")
    G.add_edge(0, 1, key=0, dp=4, type="A")

    chem = ChemistryConfig(
        model_type="atomistic", target_density=0.91,
        node_type_map={"end": NodeMoleculeConfig(molecule="C", is_end_cap=True)},
        edge_type_map={"A": EdgeChemistryConfig(monomer="M0")},
        monomers={"M0": MonomerConfig(
            smiles="CC(C)(C)CC", chain_head="C", chain_tail="C")},
        connection=ConnectionConfig(auto_bridge=True, default_bridge_atom="C"),
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rw = ChemistryBuilder(G, dims=None, config=chem).build().GetMol()
    return rw, caught


def test_legacy_type_attribute_is_honoured():
    """topon must read ``type`` and place the methyl end-cap, NOT default Si."""
    rw, caught = _build_butyl_with_legacy_type_attr()
    Chem.SanitizeMol(rw)
    # Polyisobutylene contains only C and H — no Si.  If the legacy
    # attribute were ignored, Si end-caps would leak in.
    n_si = sum(1 for a in rw.GetAtoms() if a.GetAtomicNum() == 14)
    assert n_si == 0, f"Si leaked into hydrocarbon chain: {n_si} Si atoms"


def test_legacy_type_attribute_raises_deprecation_warning():
    """A DeprecationWarning should be emitted when only ``type`` is set."""
    _, caught = _build_butyl_with_legacy_type_attr()
    msgs = [str(w.message) for w in caught if issubclass(w.category, DeprecationWarning)]
    assert any("node_type" in m for m in msgs), (
        f"expected DeprecationWarning about 'node_type', got: {msgs}"
    )


def test_new_code_using_node_type_is_silent():
    """No warning when ``node_type`` is set (new code path)."""
    G = nx.MultiGraph()
    G.add_node(0, pos=(0.0, 0.0, 0.0), node_type="end")
    G.add_node(1, pos=(4.0, 0.0, 0.0), node_type="end")
    G.add_edge(0, 1, key=0, dp=4, type="A")

    chem = ChemistryConfig(
        model_type="atomistic", target_density=0.91,
        node_type_map={"end": NodeMoleculeConfig(molecule="C", is_end_cap=True)},
        edge_type_map={"A": EdgeChemistryConfig(monomer="M0")},
        monomers={"M0": MonomerConfig(
            smiles="CC(C)(C)CC", chain_head="C", chain_tail="C")},
        connection=ConnectionConfig(auto_bridge=True, default_bridge_atom="C"),
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ChemistryBuilder(G, dims=None, config=chem).build()
    msgs = [str(w.message) for w in caught if issubclass(w.category, DeprecationWarning)]
    legacy_warnings = [m for m in msgs if "node_type" in m]
    assert not legacy_warnings, (
        f"unexpected deprecation warning for new-style graph: {legacy_warnings}"
    )
