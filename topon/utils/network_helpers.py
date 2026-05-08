import os
import re
import pickle
import random
import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist
from . import network_config

# =======================================================================
# --- ENTANGLEMENT & GEOMETRY FUNCTIONS ---
# =======================================================================

def calculate_entangled_kink(start_pos, end_pos, num_atoms, params=None, orientation_vec=None, z_phase=1.0, num_entanglements=1.0):
    """
    Generates UNWRAPPED coordinates for a chain using the "Gaussian Overshoot" logic.
    Supports repeating the entanglement motif N times along the chain.
    """
    if params is None: params = {}
    overshoot = params.get('overshoot', 0.2) 
    z_amp = params.get('z_amp', 0.5)
    sigma = params.get('sigma', 0.15)
    
    # Allow num_entanglements from params if not explicitly passed (or default)
    if num_entanglements == 1.0 and 'num_entanglements' in params:
        num_entanglements = params['num_entanglements']

    p1 = np.array(start_pos)
    p2 = np.array(end_pos)
    vec_bond = p2 - p1
    length = np.linalg.norm(vec_bond)
    
    if length < 1e-6: return {k: tuple(p1) for k in range(num_atoms)}

    # --- 1. Define Local Basis ---
    u_x = vec_bond / length
    
    if orientation_vec is not None:
        # orientation_vec should already be the MIC vector from Builder
        proj = np.dot(orientation_vec, u_x)
        v_perp = orientation_vec - proj * u_x
        if np.linalg.norm(v_perp) > 1e-6:
            u_y = v_perp / np.linalg.norm(v_perp)
        else:
            u_y = _get_arbitrary_perp(u_x)
    else:
        u_y = _get_arbitrary_perp(u_x)
    
    u_z = np.cross(u_x, u_y)
    coords = {}
    
    # Determine integer number of segments
    N = max(1, int(round(num_entanglements)))
    
    # --- 2. Generate Points ---
    for k in range(num_atoms):
        t = k / (num_atoms - 1) if num_atoms > 1 else 0.5
        
        # Determine which segment (0 to N-1) this atom belongs to
        scaled_t = t * N
        if k == num_atoms - 1:
            # Handle end of chain explicitly
            segment_idx = N - 1
            t_local = 1.0
        else:
            segment_idx = int(scaled_t)
            # t_local is the fractional progress within the segment [0, 1)
            t_local = scaled_t - segment_idx
            
        # Parametric Equations
        # X is linear along the ENTIRE chain
        val_x = length * t
        
        # Y and Z are calculated based on LOCAL segment progress (t_local)
        # Gaussian peaks at t_local = 0.5
        gaussian = np.exp(-((t_local - 0.5)**2) / (2 * sigma**2))
        
        # Y: Reach out towards partner
        val_y = length * (0.5 + overshoot) * gaussian
        
        # Z: Oscillate (Sine wave)
        # Original: val_z = length * z_amp * z_phase * np.sin(2 * np.pi * t) 
        # Here we use t_local to repeat the sine wave N times
        val_z = length * z_amp * z_phase * np.sin(2 * np.pi * t_local)
        
        # Global Pos (Unwrapped)
        global_pos = p1 + (val_x * u_x) + (val_y * u_y) + (val_z * u_z)
        coords[k] = tuple(global_pos)
        
    return coords

def _get_arbitrary_perp(u_x):
    v_temp = np.array([0, 0, 1])
    if np.abs(np.dot(u_x, v_temp)) > 0.9:
        v_temp = np.array([0, 1, 0])
    u_y = np.cross(u_x, v_temp)
    return u_y / np.linalg.norm(u_y)

def _get_mic_vector(v, dims):
    """Calculates shortest vector respecting PBC."""
    if dims is not None:
        return v - dims * np.round(v / dims)
    return v

def _get_mic_distance(p1, p2, dims):
    v = p1 - p2
    if dims is not None:
        v = v - dims * np.round(v / dims)
    return np.linalg.norm(v)

def find_crossing_candidates(G, node_positions, dims=None):
    """
    Identifies edge pairs for entanglement using '2nd Unique Distance' logic with MIC.
    """
    print("  -> Identifying crossing candidates (Nearest Disjoint Neighbor)...")
    if dims is not None:
        print(f"     (Using Periodic Boundaries: {dims})")

    edges = list(G.edges(data=True))
    valid_indices = []
    midpoints = []
    edge_nodes_map = {}
    
    # 1. Calculate Midpoints (using MIC from origin approximation or just geometric center)
    # Note: Midpoint of u,v in PBC is u + 0.5*MIC(v-u)
    for i, (u, v, data) in enumerate(edges):
        if G.degree(u) > 1 and G.degree(v) > 1:
            try:
                pos_u = np.array(node_positions[u])
                pos_v = np.array(node_positions[v])
                
                # Correct Midpoint Calculation for PBC
                vec_uv = _get_mic_vector(pos_v - pos_u, dims)
                mid = pos_u + 0.5 * vec_uv
                
                # Wrap midpoint to be inside box for consistent distance checking
                if dims is not None:
                    mid = mid - dims * np.floor(mid / dims)
                
                valid_indices.append(i)
                midpoints.append(mid)
                edge_nodes_map[i] = {u, v}
            except KeyError:
                continue

    if len(valid_indices) < 2: return []
    
    midpoints = np.array(midpoints)
    candidates = []
    processed_pairs = set()
    
    # 2. Distance Analysis
    for k in range(len(valid_indices)):
        idx_a = valid_indices[k]
        nodes_a = edge_nodes_map[idx_a]
        mid_a = midpoints[k]
        
        # Manual distance loop to apply MIC correctly between midpoints
        dists = []
        for m in range(len(valid_indices)):
            d = _get_mic_distance(mid_a, midpoints[m], dims)
            dists.append(d)
        dists = np.array(dists)
        
        # Find nearest disjoint
        sorted_indices = np.argsort(dists)
        
        for target_local in sorted_indices:
            if target_local == k: continue
            
            idx_b = valid_indices[target_local]
            nodes_b = edge_nodes_map[idx_b]
            
            if nodes_a.isdisjoint(nodes_b):
                # Found closest disjoint neighbor
                # We want to group all neighbors at roughly this distance (tolerance)
                min_dist = dists[target_local]
                
                # Get all indices that are disjoint AND within tolerance of min_dist
                # This handles symmetry where multiple edges are at the same distance
                potential_matches = []
                for p_idx in sorted_indices:
                    if p_idx == k: continue
                    if not nodes_a.isdisjoint(edge_nodes_map[valid_indices[p_idx]]): continue
                    
                    if abs(dists[p_idx] - min_dist) < 0.01:
                        potential_matches.append(p_idx)
                    elif dists[p_idx] > min_dist + 0.01:
                        break # Sorted, so we can stop
                
                if potential_matches:
                    final_local = random.choice(potential_matches)
                    idx_final = valid_indices[final_local]
                    
                    pair = tuple(sorted((idx_a, idx_final)))
                    if pair not in processed_pairs:
                        candidates.append(pair)
                        processed_pairs.add(pair)
                break # Done for this edge
    
    print(f"  -> Found {len(candidates)} candidate pairs.")
    return candidates

# =======================================================================
# --- CHEMISTRY HELPERS (Unchanged) ---
# =======================================================================

def resolve_smiles(input_str):
    if input_str in network_config.PENDANT_GROUP_SMILES:
        return network_config.PENDANT_GROUP_SMILES[input_str]
    mol = Chem.MolFromSmiles(input_str)
    if mol: return input_str
    else: raise ValueError(f"Input '{input_str}' is neither a Library key nor valid SMILES.")

def graft_side_chain(backbone, side_chain):
    graft = f"({side_chain})"
    if "(C)" in backbone: return backbone.replace("(C)", graft, 1)
    elif "([H])" in backbone: return backbone.replace("([H])", graft, 1)
    elif "(H)" in backbone: return backbone.replace("(H)", graft, 1)
    elif "(F)" in backbone: return backbone.replace("(F)", graft, 1)
    elif "(Cl)" in backbone: return backbone.replace("(Cl)", graft, 1)
    else:
        match = re.search(r"[A-Z][a-z]?", backbone)
        if match: return backbone[:match.end()] + graft + backbone[match.end():]
        else: return backbone + graft

def generate_chain_string(dp, monomer_input, bottlebrush_config=None):
    base_smiles = resolve_smiles(monomer_input)
    if base_smiles.endswith('O'): unit, linker = base_smiles[:-1], "O"
    else: unit, linker = base_smiles, "" 

    native = unit
    grafted = unit
    use_graft = False
    if bottlebrush_config and bottlebrush_config.get('enabled'):
        use_graft = True
        density = bottlebrush_config.get('graft_density', 0.5)
        sc_mono = bottlebrush_config.get('side_chain_monomer')
        sc_dp = bottlebrush_config.get('side_chain_dp', 5)
        sc_str = generate_chain_string(sc_dp, sc_mono, None).replace(":1", "").replace(":2", "")
        grafted = graft_side_chain(native, sc_str)

    seq = []
    for _ in range(dp):
        seq.append(grafted if use_graft and random.random() < density else native)
            
    if linker == "O":
        # Create body: Si...O - Si...O
        body = "O".join(seq)
        
        # We need the Head (:1) to be on the first Silicon, and Tail (:2) on the last Oxygen.
        # We append the final Oxygen cap [O:2].
        full = body + "[O:2]"
        
        # Inject the Head mapping (:1) into the first atom (assuming it is [Si])
        # This prevents the double-oxygen issue (O-Si-O-O-Si)
        if full.startswith("[Si]"):
            full = full.replace("[Si]", "[Si:1]", 1)
        else:
            # Fallback: If monomer doesn't start with Si, prepend a generic mapping or handle differently
            # For PDMS/FPDMS this safe.
            pass
            
    else: 
        full = "".join(seq)
        
    return full

def create_chain_mol(dp, monomer_input, bottlebrush_config=None):
    if dp < 1: raise ValueError("DP must be >= 1")
    smiles = generate_chain_string(dp, monomer_input, bottlebrush_config)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError("RDKit could not parse generated chain.")
    return Chem.RemoveHs(mol)

# =======================================================================
# --- FILE WRITING (Unchanged) ---
# =======================================================================

def write_lammps_displacement_file(atom_coords, sx, sy, sz, output_file, comment):
    with open(output_file, "w") as f:
        f.write(f"# LAMMPS displacement file for {comment}\n\n")
        f.write("# Define scaling factors\n")
        f.write(f"variable scale_x equal {sx:.6f}\n")
        f.write(f"variable scale_y equal {sy:.6f}\n")
        f.write(f"variable scale_z equal {sz:.6f}\n\n")
        f.write(f"# Displace {comment}\n")
        for aid in sorted(atom_coords.keys()):
            x, y, z = atom_coords[aid]
            lid = aid + 1
            f.write(f"variable dx_{lid} equal v_scale_x*{x:.6f}\n")
            f.write(f"variable dy_{lid} equal v_scale_y*{y:.6f}\n")
            f.write(f"variable dz_{lid} equal v_scale_z*{z:.6f}\n")
            f.write(f"group current_atom id {lid}\n")
            f.write(f"displace_atoms current_atom move v_dx_{lid} v_dy_{lid} v_dz_{lid}\n")
            f.write("group current_atom delete\n\n")
    print(f"Successfully wrote {comment} displacement to '{os.path.basename(output_file)}'")

def write_group_definitions_to_file(mol, node_ids, scale, periodicity, sys_type, outfile, model_type='atomistic', atom_type_map=None):
    label = "nodes (Si crosslinkers)" if sys_type == 'crosslinked' else "nodes (terminal Si)"
    sx, sy, sz = scale
    if not periodicity[0]: sx += 1
    if not periodicity[1]: sy += 1
    if not periodicity[2]: sz += 1

    with open(outfile, "w") as f:
        f.write("# LAMMPS group definitions\n\n")
        if node_ids: f.write(f"group nodes id {' '.join(map(str, sorted(node_ids)))}\n\n")
        if model_type == 'atomistic' and atom_type_map:
            f.write("# Element Groups\n")
            for sym, types in sorted(atom_type_map.items()):
                f.write(f"group {sym.lower()}_atoms type {' '.join(map(str, sorted(types)))}\n")
        elif model_type == 'coarse_grained':
            f.write("group beads type 1\n")
    print(f"Group definitions written to '{os.path.basename(outfile)}'")

# =======================================================================
# --- ATOM PLACEMENT HELPERS ---
# =======================================================================

def generate_approximate_side_chain_coords(mol, known_coords, seed_offset=0):
    """
    Generates approximate coordinates for side chain atoms (non-backbone)
    by placing them near their already-placed neighbors with a small random offset.
    
    Args:
        mol: RDKit molecule (with Hydrogens).
        known_coords: Dict {atom_idx: (x, y, z)} of known positions (graph units).
        seed_offset: Integer offset for random seed to ensure reproducibility.
        
    Returns:
        Dict {atom_idx: (x, y, z)} of new coordinates for side chain atoms.
    """
    new_coords = {}
    
    # Identify all atoms that need placement
    all_indices = {a.GetIdx() for a in mol.GetAtoms()}
    known_indices = set(known_coords.keys())
    unknown_indices = all_indices - known_indices
    
    # We may need multiple passes if side chains are multiple atoms long (e.g. propyl)
    # But for PDMS (Methyl) and Hydrogens, 1 pass usually covers most, 
    # except H on Methyl which needs Methyl to be placed first.
    
    max_passes = 3
    for _ in range(max_passes):
        progress = False
        current_unknowns = list(unknown_indices) # Copy to iterate
        
        for idx in current_unknowns:
            atom = mol.GetAtomWithIdx(idx)
            
            # Find a neighbor that is known
            parent_pos = None
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in known_coords:
                    parent_pos = np.array(known_coords[nbr_idx])
                    break
                elif nbr_idx in new_coords:
                    parent_pos = np.array(new_coords[nbr_idx])
                    break
            
            if parent_pos is not None:
                # Place with small random offset
                # 0.05 lattice units is roughly 0.5-1.0 Angstroms depending on scale
                random.seed(idx + seed_offset)
                offset = (np.random.rand(3) - 0.5) * 0.05
                pos = tuple(parent_pos + offset)
                
                new_coords[idx] = pos
                unknown_indices.remove(idx)
                progress = True
        
        if not progress:
            break
            
    return new_coords


def generate_poss_coordinates(poss_structure, node_positions, dims, arm_length_fraction=0.2):
    """
    Generates 3D coordinates for POSS cage and arm atoms.
    
    Args:
        poss_structure: Dict from ChemistryBuilder.poss_structure
            {node_id: {
                'corner_si_ids': [8 atom indices],
                'cage_oxygen_ids': [12 atom indices],
                'propyl_arm': {'corner_idx': 0, 'atom_ids': [3 indices]},
                'isooctyl_arms': {corner_idx: [8 indices], ...}
            }}
        node_positions: Dict {node_id: (x, y, z)} in graph units.
        dims: Box dimensions array for edge length scaling.
        arm_length_fraction: Fraction of average edge length for arm terminal position.
        
    Returns:
        Dict {atom_idx: (x, y, z)} of coordinates for all POSS atoms.
    """
    coords = {}
    
    # Cube corner unit vectors (normalized space diagonals)
    # Corner i at position (sign_x, sign_y, sign_z) * 0.5 * cage_size from center
    # The diagonals point outward from the cage center
    corner_dirs = np.array([
        [-1, -1, -1],  # 0
        [-1, -1, +1],  # 1
        [-1, +1, -1],  # 2
        [-1, +1, +1],  # 3
        [+1, -1, -1],  # 4
        [+1, -1, +1],  # 5
        [+1, +1, -1],  # 6
        [+1, +1, +1],  # 7
    ], dtype=float)
    # Normalize to unit vectors
    corner_dirs = corner_dirs / np.linalg.norm(corner_dirs, axis=1, keepdims=True)
    
    # Cage edge connections for oxygen placement
    cage_edges = [
        (0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3), (2, 6),
        (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)
    ]
    
    # Estimate average edge length
    avg_edge = np.mean(dims) if dims is not None else 1.0
    
    # Cage parameters (in graph units)
    cage_half_size = 0.05  # Half the cage edge length
    arm_length = arm_length_fraction * avg_edge
    
    for node_id, structure in poss_structure.items():
        if node_id not in node_positions:
            continue
            
        center = np.array(node_positions[node_id])
        
        corner_si_ids = structure['corner_si_ids']
        cage_oxygen_ids = structure['cage_oxygen_ids']
        propyl_arm = structure['propyl_arm']
        isooctyl_arms = structure['isooctyl_arms']
        
        # --- 1. Corner Si positions (cube vertices) ---
        for i, si_idx in enumerate(corner_si_ids):
            pos = center + corner_dirs[i] * cage_half_size
            coords[si_idx] = tuple(pos)
        
        # --- 2. Oxygen positions (edge midpoints) ---
        for o_idx, (a, b) in zip(cage_oxygen_ids, cage_edges):
            pos_a = center + corner_dirs[a] * cage_half_size
            pos_b = center + corner_dirs[b] * cage_half_size
            coords[o_idx] = tuple((pos_a + pos_b) / 2)
        
        # --- 3. Propyl arm (corner 0, direction inward toward network) ---
        # This arm connects to the network, so it points INWARD (opposite to corner 0 diagonal)
        # But actually, the propyl linker extends FROM the Si, so the direction depends on connectivity.
        # For simplicity: place propyl atoms along the corner 0 diagonal, extending outward.
        # The terminal C (propyl_arm['atom_ids'][-1]) is the attachment point, which should be
        # close to the network chain, so we don't extend it too far.
        propyl_dir = -corner_dirs[0]  # Point toward cage center / network
        propyl_atoms = propyl_arm['atom_ids']
        n_prop = len(propyl_atoms)
        
        corner_0_pos = center + corner_dirs[0] * cage_half_size
        for j, atom_idx in enumerate(propyl_atoms):
            # Interpolate from corner to terminal at 0.1 * edge
            frac = (j + 1) / (n_prop + 1)
            pos = corner_0_pos + propyl_dir * (frac * 0.1 * avg_edge)
            coords[atom_idx] = tuple(pos)
        
        # --- 4. IsoOctyl arms (corners 1-7, extending outward along diagonals) ---
        for corner_idx, arm_atoms in isooctyl_arms.items():
            arm_dir = corner_dirs[corner_idx]  # Outward
            corner_pos = center + arm_dir * cage_half_size
            n_arm = len(arm_atoms)
            
            for j, atom_idx in enumerate(arm_atoms):
                # Interpolate from corner to terminal at arm_length
                frac = (j + 1) / (n_arm + 1)
                pos = corner_pos + arm_dir * (frac * arm_length)
                coords[atom_idx] = tuple(pos)
    
    return coords

