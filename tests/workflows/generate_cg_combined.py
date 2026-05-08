"""
CG Workflow: Combined Entanglements + Grafts
============================================
Generates 5x5x5 CG Network with:
1. Entanglements (Kinked Backbones).
2. Grafts (Attached to Kinked Backbones).
Uses v20 Dynamic Scaling for Grafts.
"""

import sys
import argparse
from pathlib import Path
import json
import numpy as np
import random
import networkx as nx

# Add package root to path
pkg_dir = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(pkg_dir))

from rdkit import Chem

from topon.topology.loader import load_graph
from topon.writers import CGWriter, LammpsInputGenerator
from topon.conformation import ConformationManager
from topon.utils import write_lammps_displacement_file
from topon.simulation import SimulationRunner
from topon.assignment.attributor import EntanglementsConfig
from topon.assignment.entanglements import select_entanglements
from topon.utils.network_helpers import calculate_entangled_kink

# Load Configs
CONFIG_FILE = pkg_dir / "examples/config_cg_combined.json"
EXPERIMENTAL_FILE = pkg_dir / "examples/experimental_test.json"

def load_config(path):
    with open(path, 'r') as f:
        return json.load(f)

def load_experimental(path):
    with open(path, 'r') as f:
        exp = json.load(f)
    return exp

def run_workflow(output_dir):
    print(f"--- Running CG Combined Workflow -> {output_dir} ---")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Topology
    graph_path = pkg_dir / "tests/sample_graphs/network_N5x5x5_trial3.nodes"
    edges_path = pkg_dir / "tests/sample_graphs/network_N5x5x5_trial3.edges"
    
    print(f"Loading graph: {graph_path.name}")
    G, dims = load_graph(nodes_path=str(graph_path), edges_path=str(edges_path))
    
    # Load Config
    config = load_config(CONFIG_FILE)
    DP = config['chemistry']['degree_of_polymerization']
    DENSITY = config['chemistry']['bead_density']
    ASSIGNMENT = config['assignment']
    
    # Entanglements
    ENT_CONF_DICT = ASSIGNMENT.get('entanglements', {})
    if ENT_CONF_DICT.get('enabled'):
        print(f"Applying entanglements...")
        ent_config = EntanglementsConfig(**ENT_CONF_DICT)
        select_entanglements(G, ent_config, dims)
    
    # Grafts
    GRAFTS = ASSIGNMENT.get('grafts', {})
    grafts_enabled = GRAFTS.get('enabled', False)
    per_edge_grafts = GRAFTS.get('per_edge_type', {})
    
    # Experimental
    experimental = load_experimental(EXPERIMENTAL_FILE)
    extension_factor = experimental.get('cg', {}).get('graft_extension_factor', 0.5)
    print(f"Grafts: {grafts_enabled}, Factor: {extension_factor}")
    
    # 2. Chemistry
    study_name = "cg_combined"
    root = output_dir / study_name
    chem_dir = root / "02_Chemistry"
    chem_dir.mkdir(parents=True, exist_ok=True)
    
    mol = Chem.RWMol()
    node_map = {}
    
    # Nodes
    for node in sorted(G.nodes()):
        idx = mol.AddAtom(Chem.Atom('Si'))
        mol.GetAtomWithIdx(idx).SetProp("bead_type", "J")
        node_map[node] = idx
        
    edges = list(G.edges(data=True)) 
    edge_atom_map = {}
    graft_atom_map = {}
    
    total_grafts = 0
    
    for i, (u, v, data) in enumerate(edges):
        dp_val = data.get('dp', DP)
        edge_type = data.get('type', 'A')
        
        # Graft Config
        graft_conf = per_edge_grafts.get(edge_type)
        should_graft = grafts_enabled and graft_conf is not None
        if should_graft:
            density = graft_conf.get('graft_density', 0.0)
            side_dp = graft_conf.get('side_chain_dp', 5)
        
        chain_atoms = []
        prev = node_map[u]
        edge_grafts = []
        
        for k in range(dp_val):
            # Backbone
            atom = Chem.Atom('Si')
            idx = mol.AddAtom(atom)
            mol.GetAtomWithIdx(idx).SetProp("bead_type", "A") 
            chain_atoms.append(idx)
            mol.AddBond(prev, idx, Chem.BondType.SINGLE)
            
            # Graft
            if should_graft and random.random() < density:
                g_prev = idx
                graft_chain = []
                for _ in range(side_dp):
                    g_atom = Chem.Atom('Si')
                    g_idx = mol.AddAtom(g_atom)
                    mol.GetAtomWithIdx(g_idx).SetProp("bead_type", "G") 
                    graft_chain.append(g_idx)
                    mol.AddBond(g_prev, g_idx, Chem.BondType.SINGLE)
                    g_prev = g_idx
                edge_grafts.append((k, graft_chain)) 
                total_grafts += 1
                
            prev = idx
            
        mol.AddBond(prev, node_map[v], Chem.BondType.SINGLE)
        edge_atom_map[i] = chain_atoms
        graft_atom_map[i] = edge_grafts
        
    print(f"Total Grafts Added: {total_grafts}")
    
    # Count Entanglements in Loop
    entangled_count = 0
    for i, (_, _, data) in enumerate(edges):
        if data.get('entangled_with') is not None:
            entangled_count += 1
    print(f"Total Entangled Edges Processed: {entangled_count}")
    
    n_beads = mol.GetNumAtoms()
    vol = n_beads / DENSITY
    scale = (vol / np.prod(dims)) ** (1/3.0)
    sx, sy, sz = scale, scale, scale
    
    writer = CGWriter(mol, str(chem_dir / "system.data"), include_angles=True)
    writer.write()
    
    with open(chem_dir / "system.in.settings", "w") as f:
        f.write("# Dummy settings\n")
        
    # Displacements
    node_coords = {idx: G.nodes[node].get('pos', (0,0,0)) for node, idx in node_map.items()}
    write_lammps_displacement_file(node_coords, sx, sy, sz, str(chem_dir / "system_nodes.displace"), "nodes")
    
    chain_coords = {}
    graft_coords = {}
    
    for i, atoms in edge_atom_map.items():
        u, v, data = edges[i]
        pos_u = np.array(G.nodes[u].get('pos', (0,0,0)))
        pos_v = np.array(G.nodes[v].get('pos', (0,0,0)))
        vec = pos_v - pos_u
        mic = vec - dims * np.round(vec/dims)
        edge_len = np.linalg.norm(mic)
        
        # Determine Perp (Used for orientation of Kink AND Grafts)
        unit_vec = mic / (edge_len + 1e-9)
        rand_vec = np.random.randn(3)
        perp = np.cross(unit_vec, rand_vec)
        if np.linalg.norm(perp) < 1e-6: perp = np.cross(unit_vec, np.array([1,0,0]))
        perp_unit = perp / np.linalg.norm(perp)
        
        # Calculate Backbone Coords (Linear or Kinked)
        # Note: calculate_entangled_kink returns UNWRAPPED coordinates starting from 0,0,0 relative to start
        # We need to shift them to pos_u.
        
        entangled_partner_key = data.get('entangled_with')
        backbone_xyz = []
        
        if entangled_partner_key is not None:
             # Find partner properties for orientation
             # entangled_partner_key is (u, v, k) or (u, v)
             p_u, p_v = entangled_partner_key[0], entangled_partner_key[1]
             # p_data lookup not needed for geometry
             p_pos_u = np.array(G.nodes[p_u]['pos'])
             p_pos_v = np.array(G.nodes[p_v]['pos'])
             p_vec = p_pos_v - p_pos_u
             p_mic = p_vec - dims * np.round(p_vec/dims)
             
             # Midpoint logic (v15.2 fix)
             my_mid = pos_u + 0.5 * mic
             p_mid = p_pos_u + 0.5 * p_mic
             # Wrap p_mid to be close to my_mid
             delta = p_mid - my_mid
             delta -= dims * np.round(delta/dims)
             p_mid_wrapped = my_mid + delta
             
             orient_vec = p_mid_wrapped - my_mid
             if np.linalg.norm(orient_vec) < 0.01: orient_vec = perp_unit # Fallback
             
             # Calculate Kink
             # num_atoms = len(atoms) + 2 (start/end nodes) ?
             # calculate_entangled_kink takes num_atoms for full chain including nodes?
             # Let's check signature. usually num_atoms = DP+2 or internal beads?
             # Implementation in network_helpers uses parameter `num_atoms`.
             # If we pass `len(atoms)`, it generates that many points?
             # Let's generate `len(atoms)` points between start and end.
             # Actually `calculate_entangled_kink` returns points for the beads.
             
             kink_coords_dict = calculate_entangled_kink(
                 start_pos=np.zeros(3),
                 end_pos=mic,
                 num_atoms=len(atoms) + 2,    # Generate points for StartNode + Beads + EndNode
                 orientation_vec=orient_vec,
                 z_phase=1.0 # v15.2 fix
             )
             # Convert dict to list
             full_kink_path = [kink_coords_dict[k] for k in sorted(kink_coords_dict.keys())]
             
             # Use only internal points for beads
             backbone_xyz = [pos_u + np.array(pt) for pt in full_kink_path[1:-1]]
        else:
            # Linear
            for j in range(len(atoms)):
                frac = (j + 1) / (len(atoms) + 1)
                backbone_xyz.append(pos_u + frac * mic)
                
        # Assign Backbone Coords
        for j, atom_idx in enumerate(atoms):
            chain_coords[atom_idx] = tuple(backbone_xyz[j])
            
        # Assign Grafts (Anchored to Backbone)
        grafts = graft_atom_map.get(i, [])
        for k, g_atoms in grafts:
            # Anchor is backbone_xyz[k]
            # Wait, k is linear index in the edge loop (0 to DP-1).
            # So it corresponds to atoms[k].
            
            anchor_pos = backbone_xyz[k]
            
            # v20 Dynamic Logic
            graft_dp = len(g_atoms)
            backbone_dp = len(atoms)
            eff_factor = min(extension_factor, graft_dp / backbone_dp)
            target_len = edge_len * eff_factor # Scale by Edge Length (Mesh Size)
            
            # Graft Vector
            # We use `perp_unit` (Global Perpendicular).
            # Should we adjust it if we are entangled?
            # If entangled, `perp_unit` might point INTO the kink plane?
            # Ideally `perp` should be orthogonal to the LOCAL tangent.
            # But global `perp` is likely okay for topology generation.
            
            graft_vec = perp_unit * target_len
            
            for m, g_idx in enumerate(g_atoms):
                g_frac = (m + 1) / len(g_atoms)
                g_pos = anchor_pos + g_frac * graft_vec
                graft_coords[g_idx] = tuple(g_pos)
                
    write_lammps_displacement_file(chain_coords, sx, sy, sz, str(chem_dir / "system_beads.displace"), "beads")
    write_lammps_displacement_file(graft_coords, sx, sy, sz, str(chem_dir / "system_grafts.displace"), "grafts")
    
    with open(chem_dir / "system.groups", 'w') as f:
        f.write("# Groups\n")
        ids = [idx+1 for idx in sorted(node_map.values())]
        f.write(f"group nodes id {' '.join(str(x) for x in ids)}\n")
        f.write("group beads subtract all nodes\n")
        
    cm = ConformationManager(str(output_dir), study_name)
    conformed, roles = cm.apply_displacements("system.data")
    noisy = cm.apply_noise(conformed, magnitude=1e-4) 
    
    # 3. Conformation
    # We have entanglements + grafts. Overlaps are likely.
    cm.resolve_overlaps(noisy, roles, cutoff=config['conformation']['overlap_cutoff'], max_iters=20)
    
    gen = LammpsInputGenerator(str(output_dir), study_name, config=config['simulation'], experimental=experimental)
    gen.write_serial_soft_minimization(settings_file="system.in.settings", model_type="cg")
    gen.write_parallel_production(settings_file="system.in.settings", model_type="cg")
    
    print(f"CG Combined Workflow Completed. Output: {root}")
    
    if config.get('execution', {}).get('auto_run', False):
        runner = SimulationRunner(
            sim_dir=root/"04_Simulation",
            executable='lmp',
            n_procs=1,
            use_mpi=False
        )
        scripts = ["minimize_1_serial.in", "minimize_2_parallel.in", "minimize_3_parallel.in"]
        runner.run_sequence(scripts)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, help="Output directory")
    args = parser.parse_args()
    
    run_workflow(args.output)
