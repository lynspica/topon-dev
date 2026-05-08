"""
Generate Coarse-Grained (CG) Polymer Network
============================================

This script demonstrates how to generate a complete set of LAMMPS input files
for a Coarse-Grained (Bead-Spring) polymer network using the `topon` package.

Workflow:
1.  **Chemistry**: Creates a bead-spring representation from a graph topology.
    -   Generates `system.data` (placeholders) and displacement files.
2.  **Conformation**: Applies displacements to map topology to 3D space.
    -   Applies random noise to break symmetry.
    -   Resolves severe overlaps using a grid-based check.
3.  **Simulation**: Generates LAMMPS input scripts.
    -   `minimize_1_serial.in`: Soft potential minimization (creates `1.restart`).
    -   `minimize_2_cg.in`: Harmonic ramp (robust minimization).
    -   `minimize_3_cg.in`: Equilibration using FENE bonds.

Usage:
    python generate_cg.py
"""

import sys
from pathlib import Path

# Ensure 'topon' package is reachable if run from examples/ folder
pkg_dir = Path(__file__).parent.parent
sys.path.insert(0, str(pkg_dir))

from rdkit import Chem
import numpy as np

from topon.topology.loader import load_graph
from topon.writers import CGWriter, LammpsInputGenerator
from topon.conformation import ConformationManager
from topon.conformation import ConformationManager
from topon.utils import write_lammps_displacement_file
from topon.simulation import SimulationRunner

import json
import argparse

# === CONFIGURATION ===
CONFIG_FILE = pkg_dir / "examples/config_cg.json"

def load_config(path):
    with open(path, 'r') as f:
        return json.load(f)

def main(args):
    print("--- Genering CG Polymer Network ---")
    
    # Load Config
    config = load_config(CONFIG_FILE)
    
    # Extract Chemistry Params
    DP = config['chemistry']['degree_of_polymerization']
    DENSITY = config['chemistry']['bead_density']
    INCLUDE_ANGLES = config['simulation'].get('include_angles', True)
    
    OUTPUT_DIR = Path("output_cg")
    
    # 1. Load Topology (Sample 5x5x5 Graph)
    graph_path = pkg_dir / "tests/sample_graphs/network_N5x5x5_trial3.nodes"
    edges_path = pkg_dir / "tests/sample_graphs/network_N5x5x5_trial3.edges"
    
    print(f"Loading graph: {graph_path.name}")
    G, dims = load_graph(nodes_path=str(graph_path), edges_path=str(edges_path))
    
    # Update DP on edges
    for u, v, data in G.edges(data=True):
        data['dp'] = DP
    
    # 2. Setup Output Directory
    study_name = "cg_polynet"
    root = OUTPUT_DIR / study_name
    chem_dir = root / "02_Chemistry"
    chem_dir.mkdir(parents=True, exist_ok=True)
    
    # 3. STAGE 1: Chemistry (Placeholders)
    print("\nStage 1: Chemistry Generation")
    
    mol = Chem.RWMol()
    node_map = {}
    edge_atom_map = {}
    
    # Add Junctions (Nodes)
    for node in sorted(G.nodes()):
        atom = Chem.Atom('Si')
        atom.SetProp("bead_type", "J")
        idx = mol.AddAtom(Chem.Atom('Si'))
        mol.GetAtomWithIdx(idx).SetProp("bead_type", "J")
        node_map[node] = idx
        
    # Add Chains
    edges = list(G.edges(data=True))
    for i, (u, v, data) in enumerate(edges):
        dp_val = data.get('dp', DP)
        chain_atoms = []
        prev = node_map[u]
        
        for _ in range(dp_val):
            atom = Chem.Atom('Si')
            atom = Chem.Atom('Si')
            idx = mol.AddAtom(atom)
            mol.GetAtomWithIdx(idx).SetProp("bead_type", "A")
            chain_atoms.append(idx)
            mol.AddBond(prev, idx, Chem.BondType.SINGLE)
            prev = idx
        
        mol.AddBond(prev, node_map[v], Chem.BondType.SINGLE)
        edge_atom_map[i] = chain_atoms
        
    n_beads = mol.GetNumAtoms()
    print(f"  Total Beads: {n_beads}")
    
    # Calculate box scale
    vol = n_beads / DENSITY
    scale = (vol / np.prod(dims)) ** (1/3.0)
    sx, sy, sz = scale, scale, scale
    print(f"  Scaling factor: {scale:.4f}")
    
    # Write Data File
    data_file = chem_dir / "system.data"
    # Note: angles enabled by default in writer, can disable if needed
    writer = CGWriter(mol, str(data_file), include_angles=INCLUDE_ANGLES)
    writer.write()
    
    # Dummy Settings (Required by generic minimization script)
    with open(chem_dir / "system.in.settings", "w") as f:
        f.write("# Dummy settings for CG (Coeffs in data file)\n")

    # Generate Displacement Files
    # Nodes
    node_coords = {idx: G.nodes[node].get('pos', (0,0,0)) for node, idx in node_map.items()}
    write_lammps_displacement_file(node_coords, sx, sy, sz, str(chem_dir / "system_nodes.displace"), "nodes")
    
    # Chains (Interpolated)
    chain_coords = {}
    for i, atoms in edge_atom_map.items():
        u, v, _ = edges[i]
        pos_u = np.array(G.nodes[u].get('pos', (0,0,0)))
        pos_v = np.array(G.nodes[v].get('pos', (0,0,0)))
        raw_vec = pos_v - pos_u
        mic_vec = raw_vec - dims * np.round(raw_vec / dims)
        
        for j, atom_idx in enumerate(atoms):
            frac = (j + 1) / (len(atoms) + 1)
            chain_coords[atom_idx] = tuple(pos_u + frac * mic_vec)
            
    write_lammps_displacement_file(chain_coords, sx, sy, sz, str(chem_dir / "system_beads.displace"), "beads")
    
    # Groups File
    with open(chem_dir / "system.groups", 'w') as f:
        f.write("# Groups\n")
        node_ids = [idx+1 for idx in sorted(node_map.values())]
        f.write(f"group nodes id {' '.join(str(x) for x in node_ids)}\n")
        f.write("group beads subtract all nodes\n")
        
    # 4. STAGE 2: Conformation
    print("\nStage 2: Conformation Mapping")
    cm = ConformationManager(str(OUTPUT_DIR), study_name)
    conformed_file, atom_roles = cm.apply_displacements("system.data")
    
    # Apply Noise (Symmetry Breaking)
    noisy_file = cm.apply_noise(conformed_file, magnitude=1e-4) # LJ units
    
    # Check Overlaps (Redundant check)
    cutoff = config['conformation']['overlap_cutoff']
    max_iters = config['conformation']['overlap_max_iters']
    relaxed_file = cm.resolve_overlaps(noisy_file, atom_roles, cutoff=cutoff, max_iters=max_iters)
    
    # 5. STAGE 3: LAMMPS Scripts
    print("\nStage 3: LAMMPS Input Generation")
    gen = LammpsInputGenerator(str(OUTPUT_DIR), study_name, config=config['simulation'])
    
    # Script 1: Serial Soft Min (Creates 1.restart)
    gen.write_serial_soft_minimization(
        input_data="system_relaxed.data",
        groups_file="system.groups",
        settings_file="system.in.settings",
        model_type="cg"
    )
    
    # Script 2 & 3: Parallel Min/Equil (Reads 1.restart)
    gen.write_parallel_production(
        settings_file="system.in.settings",
        model_type="cg"
    )
    
    print("  lmp_mpi -np 4 -in minimize_2_cg.in")
    
    # Run if requested
    if args.run or config.get('execution', {}).get('auto_run', False):
        exec_conf = config.get('execution', {'executable': 'lmp', 'n_procs': 1})
        runner = SimulationRunner(
            sim_dir=root/"04_Simulation",
            executable=exec_conf.get('executable', 'lmp'),
            n_procs=exec_conf.get('n_procs', 1),
            use_mpi=False # Force serial for now as lmp_mpi is missing
        )
        
        scripts = ["minimize_1_serial.in", "minimize_2_parallel.in", "minimize_3_parallel.in"]
        runner.run_sequence(scripts)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="Execute LAMMPS simulations after generation")
    args = parser.parse_args()
    
    main(args)
