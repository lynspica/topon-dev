"""
Generate Atomistic (PDMS) Polymer Network
=========================================

This script demonstrates how to generate a complete set of LAMMPS input files
for an Atomistic (PDMS) polymer network using the `topon` package.

Workflow:
1.  **Chemistry**: Creates a fully atomistic PDMS structure from a graph topology.
    -   Adds side chains (Methyls) and hydrogens.
    -   Generates `system.data` (placeholders) and 4 displacement files.
2.  **Conformation**: Applies displacements to map topology to 3D space.
    -   Resolves hard overlaps using iterative relaxation (~50 steps).
3.  **Simulation**: Generates LAMMPS input scripts.
    -   `minimize_1_serial.in`: Soft potential minimization.
    -   `minimize_2_parallel.in`: Soft-to-Real Ramp.
    -   `minimize_3_parallel.in`: Equilibration.

Usage:
    python generate_atomistic.py
"""

import sys
from pathlib import Path

# Ensure 'topon' package is reachable if run from examples/ folder
pkg_dir = Path(__file__).parent.parent
sys.path.insert(0, str(pkg_dir))

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import shutil

from topon.topology.loader import load_graph
from topon.writers import DreidingWriter, LammpsInputGenerator
from topon.conformation import ConformationManager
from topon.utils import write_lammps_displacement_file, generate_approximate_side_chain_coords

# === CONFIGURATION ===
DP = 5                 # Degree of Polymerization (Monomers per chain)
DENSITY = 0.90          # Target Density (g/cm^3)
OUTPUT_DIR = Path("output_atomistic")

def main():
    print("--- Genering Atomistic PDMS Network ---")
    
    # 1. Load Topology (Sample 5x5x5 Graph)
    graph_path = pkg_dir / "tests/sample_graphs/network_N5x5x5_trial3.nodes"
    edges_path = pkg_dir / "tests/sample_graphs/network_N5x5x5_trial3.edges"
    
    print(f"Loading graph: {graph_path.name}")
    G, dims = load_graph(nodes_path=str(graph_path), edges_path=str(edges_path))
    
    # Update DP on edges
    for u, v, data in G.edges(data=True):
        data['dp'] = DP
    
    # 2. Setup Output Directory
    study_name = "atomistic_polynet"
    root = OUTPUT_DIR / study_name
    chem_dir = root / "02_Chemistry"
    chem_dir.mkdir(parents=True, exist_ok=True)
    
    # 3. STAGE 1: Chemistry Generation
    print("\nStage 1: Chemistry Generation (Molecule Building)")
    
    mol = Chem.RWMol()
    node_map = {} # node_id -> Si idx
    edge_backbone = {} # edge_idx -> [Si idx]
    
    # Add Nodes (Crosslinkers)
    for node in sorted(G.nodes()):
        idx = mol.AddAtom(Chem.Atom('Si'))
        node_map[node] = idx
        
    # Add PDMS Chains
    edges = list(G.edges(data=True))
    for i, (u, v, data) in enumerate(edges):
        dp_val = data.get('dp', DP)
        backbone_si = []
        
        # Link u -> O -> Si ... -> O -> v
        o_start = mol.AddAtom(Chem.Atom('O'))
        mol.AddBond(node_map[u], o_start, Chem.BondType.SINGLE)
        prev = o_start
        
        for _ in range(dp_val):
            si = mol.AddAtom(Chem.Atom('Si'))
            backbone_si.append(si)
            mol.AddBond(prev, si, Chem.BondType.SINGLE)
            
            # Side Methyls
            for _ in range(2):
                c = mol.AddAtom(Chem.Atom('C'))
                mol.AddBond(si, c, Chem.BondType.SINGLE)
            
            o = mol.AddAtom(Chem.Atom('O'))
            mol.AddBond(si, o, Chem.BondType.SINGLE)
            prev = o
            
        mol.AddBond(prev, node_map[v], Chem.BondType.SINGLE)
        edge_backbone[i] = backbone_si
        
    # Cap Node Valences
    for node, idx in node_map.items():
        atom = mol.GetAtomWithIdx(idx)
        caps = 4 - atom.GetDegree()
        for _ in range(max(0, caps)):
            c = mol.AddAtom(Chem.Atom('C'))
            mol.AddBond(idx, c, Chem.BondType.SINGLE)
            
    # Add Hydrogens & Charges
    print("  Sanitizing and adding Hydrogens...")
    final_mol = mol.GetMol()
    Chem.SanitizeMol(final_mol)
    mol_h = Chem.AddHs(final_mol)
    AllChem.ComputeGasteigerCharges(mol_h)
    
    n_total = mol_h.GetNumAtoms()
    print(f"  Total Atoms: {n_total}")
    
    # Calculate Scaling
    mol_weight = sum(a.GetMass() for a in mol_h.GetAtoms())
    vol = (mol_weight / DENSITY) * 1.66054
    scale = (vol / np.prod(dims)) ** (1/3.0)
    sx, sy, sz = scale, scale, scale
    print(f"  Scaling factor: {scale:.4f}")
    
    # Write Data & Settings
    data_file = chem_dir / "system.data"
    writer = DreidingWriter(mol_h, str(data_file), use_charges=True)
    writer.write() # Writes system.data AND system.in.settings

    # Displacements
    print("  Generating displacement files...")
    # Nodes
    node_coords = {idx: G.nodes[node].get('pos', (0,0,0)) for node, idx in node_map.items()}
    write_lammps_displacement_file(node_coords, sx, sy, sz, str(chem_dir / "system_nodes.displace"), "nodes")
    
    # Backbone
    backbone_coords = {}
    for i, atoms in edge_backbone.items():
        u, v, _ = edges[i]
        pos_u = np.array(G.nodes[u].get('pos', (0,0,0)))
        pos_v = np.array(G.nodes[v].get('pos', (0,0,0)))
        vec = pos_v - pos_u
        mic = vec - dims * np.round(vec/dims)
        for j, a_idx in enumerate(atoms):
            frac = (j+1)/(len(atoms)+1)
            backbone_coords[a_idx] = tuple(pos_u + frac*mic)
    write_lammps_displacement_file(backbone_coords, sx, sy, sz, str(chem_dir / "system_backbone.displace"), "backbone")
    
    # Side Chains
    known = {**node_coords, **backbone_coords}
    side_coords = generate_approximate_side_chain_coords(mol_h, known)
    
    h_coords = {k:v for k,v in side_coords.items() if mol_h.GetAtomWithIdx(k).GetSymbol() == 'H'}
    p_coords = {k:v for k,v in side_coords.items() if mol_h.GetAtomWithIdx(k).GetSymbol() != 'H'}
    
    write_lammps_displacement_file(p_coords, sx, sy, sz, str(chem_dir / "system_pendant.displace"), "pendant")
    write_lammps_displacement_file(h_coords, sx, sy, sz, str(chem_dir / "system_hydrogens.displace"), "hydrogens")
    
    # Groups
    with open(chem_dir / "system.groups", 'w') as f:
        f.write("# Groups\n")
        ids = [idx+1 for idx in sorted(node_map.values())]
        f.write(f"group nodes id {' '.join(str(x) for x in ids)}\n")

    # 4. STAGE 2: Conformation
    print("\nStage 2: Conformation Mapping")
    cm = ConformationManager(str(OUTPUT_DIR), study_name)
    conformed_file, atom_roles = cm.apply_displacements("system.data")
    
    # Relax (longer due to many atoms)
    relaxed_file = cm.resolve_overlaps(conformed_file, atom_roles, cutoff=0.85, max_iters=50)
    
    # 5. STAGE 3: LAMMPS Scripts
    print("\nStage 3: LAMMPS Input Generation")
    gen = LammpsInputGenerator(str(OUTPUT_DIR), study_name)
    
    gen.write_serial_soft_minimization(
        input_data="system_relaxed.data",
        groups_file="system.groups",
        settings_file="system.in.settings",
        model_type="atomistic"
    )
    
    gen.write_parallel_production(
        settings_file="system.in.settings",
        model_type="atomistic"
    )
    
    print(f"\nDone! Output available in: {root.resolve()}")
    print("To run:")
    print(f"  cd {root / '04_Simulation'}")
    print("  lmp_serial -in minimize_1_serial.in")
    print("  lmp_mpi -np 4 -in minimize_2_parallel.in")

if __name__ == "__main__":
    main()
