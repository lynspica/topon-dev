"""
Atomistic Workflow: Combined Entanglements + Grafts
===================================================
Generates 5x5x5 Atomistic Network with:
1. Entanglements (Kinked Backbones).
2. Grafts (Attached to Kinked Backbones).
Uses v20 Dynamic Scaling and v21.1 Geometry Fixes.
"""

import sys
import argparse
from pathlib import Path
import json
import numpy as np
import random

# Add package root to path
pkg_dir = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(pkg_dir))

from rdkit import Chem
from rdkit.Chem import AllChem

from topon.topology.loader import load_graph
from topon.writers import DreidingWriter, LammpsInputGenerator
from topon.conformation import ConformationManager
from topon.utils import write_lammps_displacement_file, generate_approximate_side_chain_coords
from topon.simulation import SimulationRunner
from topon.assignment.attributor import EntanglementsConfig
from topon.assignment.entanglements import select_entanglements
from topon.utils.network_helpers import calculate_entangled_kink

# Load Configs
CONFIG_FILE = pkg_dir / "examples/config_atomistic_combined.json"
EXPERIMENTAL_FILE = pkg_dir / "examples/experimental_test.json" 

def load_config(path):
    with open(path, 'r') as f:
        return json.load(f)

def load_experimental(path):
    with open(path, 'r') as f:
        exp = json.load(f)
    return exp

def run_workflow(output_dir):
    print(f"--- Running Atomistic Combined Workflow -> {output_dir} ---")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load Config
    config = load_config(CONFIG_FILE)
    CHEM = config['chemistry']
    DP = CHEM['degree_of_polymerization']
    DENSITY = CHEM.get('target_density', 0.9)
    
    ASSIGNMENT = config['assignment']
    
    # 1. Topology
    graph_path = pkg_dir / "tests/sample_graphs/network_N5x5x5_trial3.nodes"
    edges_path = pkg_dir / "tests/sample_graphs/network_N5x5x5_trial3.edges"
    
    print(f"Loading graph: {graph_path.name}")
    G, dims = load_graph(nodes_path=str(graph_path), edges_path=str(edges_path))
    
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
    
    # Experimental Config
    experimental = load_experimental(EXPERIMENTAL_FILE)
    extension_factor = experimental.get('atomistic', {}).get('graft_extension_factor', 0.5)
    print(f"Config: DP={DP}, Density={DENSITY}, Grafts={grafts_enabled}, Factor={extension_factor}")
    
    study_name = "atomistic_combined"
    root = output_dir / study_name
    chem_dir = root / "02_Chemistry"
    chem_dir.mkdir(parents=True, exist_ok=True)
    
    mol = Chem.RWMol()
    node_map = {}
    edge_backbone = {} 
    graft_map = {} 
    
    # Nodes (Si)
    for node in sorted(G.nodes()):
        idx = mol.AddAtom(Chem.Atom('Si'))
        node_map[node] = idx
        
    edges = list(G.edges(data=True))
    total_grafts = 0
    
    for i, (u, v, data) in enumerate(edges):
        dp_val = data.get('dp', DP)
        edge_type = data.get('type', 'A')
        
        # Check graft config
        graft_conf = per_edge_grafts.get(edge_type)
        should_graft = grafts_enabled and graft_conf is not None
        if should_graft:
            density = graft_conf.get('graft_density', 0.0)
            side_dp = graft_conf.get('side_chain_dp', 5)
            # cap_atom = None  -> backbone-matched cap (extra methyl on
            # terminal Si for PDMS, giving trimethylsilyl -Si(CH3)3 instead
            # of an H-cap). cap_atom = 'H' preserves the legacy RDKit
            # implicit-H fill on the terminal Si.
            cap_atom = graft_conf.get('cap_atom', None)
        
        backbone_si = []
        edge_id_grafts = []
        
        # Add O start
        o_start = mol.AddAtom(Chem.Atom('O'))
        mol.AddBond(node_map[u], o_start, Chem.BondType.SINGLE)
        prev = o_start
        
        for k in range(dp_val):
            si = mol.AddAtom(Chem.Atom('Si'))
            backbone_si.append(si)
            mol.AddBond(prev, si, Chem.BondType.SINGLE)
            
            # Graft Logic
            is_grafted = should_graft and random.random() < density
            
            if is_grafted:
                c1 = mol.AddAtom(Chem.Atom('C'))
                mol.AddBond(si, c1, Chem.BondType.SINGLE)
                
                g_o = mol.AddAtom(Chem.Atom('O'))
                mol.AddBond(si, g_o, Chem.BondType.SINGLE)
                
                g_prev = g_o
                g_atoms = [g_prev] # O
                
                for _ in range(side_dp):
                    g_si = mol.AddAtom(Chem.Atom('Si'))
                    g_atoms.append(g_si) # Si
                    mol.AddBond(g_prev, g_si, Chem.BondType.SINGLE)

                    g_c1 = mol.AddAtom(Chem.Atom('C'))
                    g_c2 = mol.AddAtom(Chem.Atom('C'))
                    mol.AddBond(g_si, g_c1, Chem.BondType.SINGLE)
                    mol.AddBond(g_si, g_c2, Chem.BondType.SINGLE)

                    if _ < side_dp - 1:
                        g_next_o = mol.AddAtom(Chem.Atom('O'))
                        g_atoms.append(g_next_o) # O
                        mol.AddBond(g_si, g_next_o, Chem.BondType.SINGLE)
                        g_prev = g_next_o
                    else:
                        # Terminal Si of the side chain. Default behavior
                        # (cap_atom=None) is to add one extra methyl so
                        # the cap is trimethylsilyl -Si(CH3)3, matching
                        # the backbone PDMS chemistry. cap_atom='H' falls
                        # through (no extra atom) -> RDKit fills the
                        # missing fourth bond with an implicit H.
                        if cap_atom is None:
                            g_c3 = mol.AddAtom(Chem.Atom('C'))
                            mol.AddBond(g_si, g_c3, Chem.BondType.SINGLE)
                        
                edge_id_grafts.append(((k+1)/(dp_val+1), g_atoms))
                total_grafts += 1
                
            else:
                c1 = mol.AddAtom(Chem.Atom('C'))
                c2 = mol.AddAtom(Chem.Atom('C'))
                mol.AddBond(si, c1, Chem.BondType.SINGLE)
                mol.AddBond(si, c2, Chem.BondType.SINGLE)
            
            o = mol.AddAtom(Chem.Atom('O'))
            mol.AddBond(si, o, Chem.BondType.SINGLE)
            prev = o
            
        mol.AddBond(prev, node_map[v], Chem.BondType.SINGLE)
        edge_backbone[i] = backbone_si
        if edge_id_grafts:
            graft_map[i] = edge_id_grafts
        
    print(f"Total Grafts Added: {total_grafts}")

    # Cap Nodes
    for node, idx in node_map.items():
        atom = mol.GetAtomWithIdx(idx)
        caps = 4 - atom.GetDegree()
        for _ in range(max(0, caps)):
            c = mol.AddAtom(Chem.Atom('C'))
            mol.AddBond(idx, c, Chem.BondType.SINGLE)
            
    Chem.SanitizeMol(mol)
    mol_h = Chem.AddHs(mol)
    AllChem.ComputeGasteigerCharges(mol_h)
    
    # Scale
    mass = sum(a.GetMass() for a in mol_h.GetAtoms())
    vol = (mass / DENSITY) * 1.66054
    scale = (vol / np.prod(dims)) ** (1/3.0)
    sx, sy, sz = scale, scale, scale
    
    writer = DreidingWriter(mol_h, str(chem_dir / "system.data"), use_charges=True)
    writer.write()

    # Displacements
    node_coords = {idx: G.nodes[node].get('pos', (0,0,0)) for node, idx in node_map.items()}
    write_lammps_displacement_file(node_coords, sx, sy, sz, str(chem_dir / "system_nodes.displace"), "nodes")
    
    backbone_coords = {}
    graft_coords = {}
    
    for i, atoms in edge_backbone.items():
        u, v, data = edges[i]
        pos_u = np.array(G.nodes[u].get('pos', (0,0,0)))
        pos_v = np.array(G.nodes[v].get('pos', (0,0,0)))
        vec = pos_v - pos_u
        mic = vec - dims * np.round(vec/dims)
        edge_len = np.linalg.norm(mic)
        
        # Determine Perp
        unit_vec = mic / (edge_len + 1e-9)
        rand_vec = np.random.randn(3)
        perp = np.cross(unit_vec, rand_vec)
        if np.linalg.norm(perp) < 1e-6: perp = np.cross(unit_vec, np.array([1,0,0]))
        perp_unit = perp / np.linalg.norm(perp)
        
        # Calculate Backbone Coords (v21.1 Fix: N+2 logic)
        entangled_partner_key = data.get('entangled_with')
        backbone_xyz = []
        
        if entangled_partner_key is not None:
             p_u, p_v = entangled_partner_key[0], entangled_partner_key[1]
             p_pos_u = np.array(G.nodes[p_u]['pos'])
             p_pos_v = np.array(G.nodes[p_v]['pos'])
             p_vec = p_pos_v - p_pos_u
             p_mic = p_vec - dims * np.round(p_vec/dims)
             
             my_mid = pos_u + 0.5 * mic
             p_mid = p_pos_u + 0.5 * p_mic
             delta = p_mid - my_mid
             delta -= dims * np.round(delta/dims)
             p_mid_wrapped = my_mid + delta
             orient_vec = p_mid_wrapped - my_mid
             if np.linalg.norm(orient_vec) < 0.01: orient_vec = perp_unit
             
             # Calculate Kink with N+2 fix
             kink_coords_dict = calculate_entangled_kink(
                 start_pos=np.zeros(3),
                 end_pos=mic,
                 num_atoms=len(atoms) + 2, # N+2
                 orientation_vec=orient_vec,
                 z_phase=1.0
             )
             full_kink_path = [kink_coords_dict[k] for k in sorted(kink_coords_dict.keys())]
             
             # Slice internal [1:-1]
             backbone_xyz = [pos_u + np.array(pt) for pt in full_kink_path[1:-1]]
        else:
             # Linear
             for j in range(len(atoms)):
                 frac = (j + 1) / (len(atoms) + 1)
                 backbone_xyz.append(pos_u + frac * mic)
                 
        # Assign Backbone Coords (Si)
        for j, a_idx in enumerate(atoms):
            backbone_coords[a_idx] = tuple(backbone_xyz[j])
            
        # Graft Coords (Revised)
        if i in graft_map:
            for frac, g_atoms in graft_map[i]:
                # Anchor is backbone_xyz[k].
                # We need to find k corresponding to frac?
                # frac was (k+1)/(dp+1).
                # k = index in loop.
                # In graft_map construction: edge_id_grafts.append(((k+1)..., g_atoms))
                # BUT graft_map keys (values) are list of (frac, g_atoms).
                # I lost the integer 'k' in this structure in atomistic_grafts script?
                # Check line 140: ((k+1)/(dp_val+1), g_atoms).
                # I can recover k from frac?
                # k = integer roughly frac * (dp+1) - 1.
                # Or better, store k in graft_map tuple?
                # Or just assume linear spacing and match closest?
                # Wait, internal beads are exactly at (k+1)/(dp+1).
                # The 'backbone_xyz' has dp_val elements.
                # k goes 0..dp-1.
                # So backbone_xyz[k] is the bead.
                
                # I'll calculate k from frac.
                k_float = frac * (len(atoms) + 1) - 1
                k = int(round(k_float))
                
                anchor_pos = backbone_xyz[k]
                
                # v20 Dynamic Logic
                graft_dp = len(g_atoms) / 2.0
                backbone_dp = len(atoms)
                
                eff_factor = min(extension_factor, graft_dp / backbone_dp)
                target_len = edge_len * eff_factor
                graft_vec = perp_unit * target_len
                
                for m, g_idx in enumerate(g_atoms):
                    g_frac = (m + 1) / len(g_atoms)
                    pos = anchor_pos + g_frac * graft_vec
                    graft_coords[g_idx] = tuple(pos)
                    
    write_lammps_displacement_file(backbone_coords, sx, sy, sz, str(chem_dir / "system_backbone.displace"), "backbone")
    write_lammps_displacement_file(graft_coords, sx, sy, sz, str(chem_dir / "system_grafts.displace"), "grafts")
    
    # Side Groups
    known = {**node_coords, **backbone_coords, **graft_coords}
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

    cm = ConformationManager(str(output_dir), study_name)
    conformed, roles = cm.apply_displacements("system.data")
    cm.resolve_overlaps(conformed, roles, cutoff=config['conformation']['overlap_cutoff'], max_iters=20)
    
    gen = LammpsInputGenerator(str(output_dir), study_name, config=config['simulation'], experimental=experimental)
    gen.write_serial_soft_minimization(settings_file="system.in.settings", model_type="atomistic")
    gen.write_parallel_production(settings_file="system.in.settings", model_type="atomistic")
    
    print(f"Atomistic Combined Workflow Completed. Output: {root}")
    
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
