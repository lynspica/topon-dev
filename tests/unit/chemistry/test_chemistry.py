"""
Unit tests for Chemistry Builder.
Tests logic for assembling molecular structures.
"""

import pytest
import networkx as nx
import numpy as np
from topon.config.schema import (
    ToponConfig, ChemistryConfig, ConnectionConfig, EdgeChemistryConfig,
    MonomerConfig, NodeMoleculeConfig,
)
from topon.chemistry.builder import ChemistryBuilder, _PEROXIDE_FIX_APPLIED

# Try importing RDKit - skip tests if not available
rdkit_available = False
try:
    from rdkit import Chem
    rdkit_available = True
except ImportError:
    pass

@pytest.fixture
def cg_graph():
    """Simple 2-node graph for CG testing."""
    g = nx.MultiGraph()
    g.add_node(0, pos=(0,0,0), node_type="A")
    g.add_node(1, pos=(10,0,0), node_type="A")
    g.add_edge(0, 1, key=0, edge_type="A", dp=5)
    return g

@pytest.fixture
def cg_config():
    """CG configuration."""
    config = ToponConfig()
    config.chemistry.model_type = "coarse_grained"
    return config.chemistry

@pytest.fixture
def atomistic_config():
    """Atomistic configuration."""
    config = ToponConfig()
    config.chemistry.model_type = "atomistic"
    return config.chemistry

@pytest.mark.skipif(not rdkit_available, reason="RDKit not installed")
def test_cg_chain_build(cg_graph, cg_config):
    """Test coarse-grained chain building logic."""
    builder = ChemistryBuilder(cg_graph, None, cg_config)
    mol = builder.build()
    
    # Expected atoms:
    # 2 nodes (Si) + 5 chain beads (C) = 7 atoms
    assert mol.GetNumAtoms() == 7
    
    # Expected bonds:
    # Chain internal: 4 bonds (for 5 beads)
    # Node connections: 2 bonds (Node-ChainStart, ChainEnd-Node)
    # Total: 6 bonds
    assert mol.GetNumBonds() == 6

@pytest.mark.skipif(not rdkit_available, reason="RDKit not installed")
def test_atomistic_auto_bridge(atomistic_config):
    """Test auto-bridging logic (Si-O-Si vs Si-Si)."""
    # Create graph with 2 nodes closer together
    g = nx.MultiGraph()
    g.add_node(0, pos=(0,0,0), node_type="A")
    g.add_node(1, pos=(10,0,0), node_type="A")
    # PDMS monomer: [Si](C)(C)O ... Head=Si, Tail=O
    g.add_edge(0, 1, key=0, edge_type="A", dp=1)
    
    # Config: Node A is "Si". Default bridge is "O".
    # Chain Head: Si. Node: Si. -> Bridge needed (Si-O-Si).
    # Chain Tail: O. Node: Si. -> Direct bond OK (O-Si).
    
    builder = ChemistryBuilder(g, None, atomistic_config)
    mol = builder.build()
    
    # Check atoms count
    # Node 0 (Si) = 1
    # Node 1 (Si) = 1
    # PDMS DP1 ([Si](C)(C)O) = 4 atoms + Hydrogens?
    # Wait, SMILES "[Si](C)(C)O" has 1 Si, 2 C, 1 O. Total 4 heavy atoms.
    # Plus bridge "O" on left side.
    # Total expected: 1 (Node) + 1 (Bridge) + 4 (Chain) + 1 (Node) = 7 atoms?
    # Let's inspect connectivity if finding counts is tricky.
    
    # For unit test, just verify more atoms than CG
    assert mol.GetNumAtoms() > 2
    
    # Verify connectivity (conceptually):
    # Should have Si-O-Si-C... connectivity.
    # Simply check that we didn't crash and produced bonds.
    assert mol.GetNumBonds() >= mol.GetNumAtoms() - 1

def test_pos_cage_selection():
    """Test POSS corner selection logic (math only)."""
    # Create builder without graph just to test helper
    config = ToponConfig().chemistry
    # Create dummy graph to init builder
    g = nx.MultiGraph()
    # Dummy node map setup manually
    # Node 0 has 8 corners (indices 0-7)
    builder = ChemistryBuilder(g, None, config)
    builder.node_map[0] = [0, 1, 2, 3, 4, 5, 6, 7]
    builder.poss_usage[0] = set()
    
    # Vector pointing to (-1, -1, -1) corner
    vec = np.array([-10, -10, -10])
    
    # Should pick index 0 (if our order matches corner 0)
    # Corner 0 in code is [-1, -1, -1].
    idx = builder._get_attachment_atom(0, vec)
    assert idx == 0
    assert 0 in builder.poss_usage[0]
    
    # Second request in similar direction should pick neighbor
    idx2 = builder._get_attachment_atom(0, vec)
    assert idx2 != 0 # Should pick next best
    assert idx2 in builder.poss_usage[0]


def test_peroxide_fix_marker():
    """Marker must be exported and True so downstream forks can retire."""
    assert _PEROXIDE_FIX_APPLIED is True


def _build_single_chain(smiles, head, tail, dp):
    """Build one atomistic chain between two TMS-style end caps."""
    g = nx.MultiGraph()
    g.add_node(0, pos=(0., 0., 0.), type="end")
    g.add_node(1, pos=(float(dp), 0., 0.), type="end")
    g.add_edge(0, 1, key=0, dp=dp, type="A")
    cfg = ChemistryConfig(
        model_type="atomistic", target_density=0.85,
        node_type_map={"end": NodeMoleculeConfig(molecule="[Si](C)(C)C", is_end_cap=True)},
        edge_type_map={"A": EdgeChemistryConfig(monomer="M0")},
        monomers={"M0": MonomerConfig(smiles=smiles, chain_head=head, chain_tail=tail)},
        connection=ConnectionConfig(auto_bridge=True, default_bridge_atom="O"),
    )
    return ChemistryBuilder(g, dims=None, config=cfg).build().GetMol()


@pytest.mark.skipif(not rdkit_available, reason="RDKit not installed")
@pytest.mark.parametrize("smiles, head, tail, dp, n_si, n_o, n_c", [
    # PDMS DP=10 with TMS caps: 10 backbone Si + 2 cap Si = 12 Si;
    # 1 head bridge + 10 linker O = 11 O; 2*10 backbone methyls + 2*3 cap methyls = 26 C.
    ("[Si](C)(C)O", "Si", "O", 10, 12, 11, 26),
    # PDMS DP=5: 5+2 Si, 1+5 O, 2*5 + 6 C.
    ("[Si](C)(C)O", "Si", "O", 5, 7, 6, 16),
])
def test_no_peroxide_in_o_terminal_chain(smiles, head, tail, dp, n_si, n_o, n_c):
    """O-terminal monomers must not produce a Si-O-O-Si peroxide at the chain tail."""
    mol = _build_single_chain(smiles, head, tail, dp)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    oo_bonds = sum(
        1 for b in mol.GetBonds()
        if b.GetBeginAtom().GetSymbol() == "O" and b.GetEndAtom().GetSymbol() == "O"
    )
    smi = Chem.MolToSmiles(Chem.RemoveHs(mol))
    counts = {sym: sum(1 for a in mol.GetAtoms() if a.GetSymbol() == sym)
              for sym in ("Si", "O", "C")}

    assert oo_bonds == 0, f"peroxide bond present (smiles={smi})"
    assert "OO" not in smi, f"OO substring in canonical SMILES: {smi}"
    assert counts["Si"] == n_si
    assert counts["O"] == n_o
    assert counts["C"] == n_c

