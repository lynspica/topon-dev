# make_kremer_grest.py (FINAL, FULLY AUTOMATED VERSION)
import itertools
from rdkit import Chem

def get_params_with_default(param_dict, key, default_key='default'):
    """Helper to find specific parameters or fall back to the default set."""
    if key in param_dict: return param_dict[key]
    if isinstance(key, tuple) and key[::-1] in param_dict: return param_dict[key[::-1]]
    if default_key in param_dict: return param_dict[default_key]
    return None

def create_lammps_data_file(mol, filename="polymer_network_cg.lammps"):
    """
    Generate a LAMMPS data file with a fully automated parameter system.
    New monomer types are automatically defined by cloning Type 'A'.
    """
    # =======================================================================
    # --- 1. SCRIPT TOGGLES & BASE PARAMETERS ---
    # =======================================================================
    ASSIGN_JUNCTION_TYPE = True
    PROCESS_ANGLES = True

    # You only need to define your base types and any special overrides.
    cg_params = {
        'atom_types': {
            'A': {'mass': 1.0, 'sigma': 1.0, 'epsilon': 1.0},
            'B': {'mass': 1.0, 'sigma': 1.0, 'epsilon': 1.0},
            'J': {'mass': 1.0, 'sigma': 1.0, 'epsilon': 1.0},
            'POSS_CORE': {'mass': 1.0, 'sigma': 1.0, 'epsilon': 1.0}, # <-- ADD THIS LINE
            #'J': {'mass': 1.0, 'sigma': 1.2, 'epsilon': 1.0},
        },
        'bond_types': {
            'default': {'K': 30.0, 'R0': 1.5, 'epsilon': 1.0, 'sigma': 1.0},
            ('A', 'J'): {'K': 30.0, 'R0': 1.5, 'epsilon': 1.0, 'sigma': 1.0},
            #('A', 'J'): {'K': 35.0, 'R0': 1.5, 'epsilon': 1.0, 'sigma': 1.1},
        },
        'angle_types': {
            'default': {'K': 1.5, 'Theta0': 180.0},
            ('A', 'J', 'A'): {'K': 1.5, 'Theta0': 180.0},
            #('A', 'J', 'A'): {'K': 5.0, 'Theta0': 120.0},
        },
        'cutoff': 2**(1/6)
    }

    # =======================================================================
    # --- 2. MAP RDKIT SYMBOLS & EXTRACT TOPOLOGY ---
    # =======================================================================
    #symbol_to_abstract_type = {
    #   'Ge': 'J', 'Si': 'A', 'C': 'B', 'O': 'C_type', 'N': 'D_type', 
    #   'S': 'E_type', 'P': 'F_type', 'F': 'G_type', 'Cl': 'H_type',
    #   
    #}

    #symbol_to_abstract_type = {
    #    'Ge': 'J', 'Si': 'A', 'C': 'B', 'O': 'C', 'N': 'D', # <-- CHANGE 'C_type' to 'C', etc.
    #    'S': 'E', 'P': 'F', 'F': 'G', 'Cl': 'H'
    #}
    
    #if not ASSIGN_JUNCTION_TYPE: symbol_to_abstract_type['Ge'] = 'A'

    # Extract all bead types present in the molecule
    print("Mapping RDKit atoms to abstract types using 'bead_type' property...") # Optional: Add print
    idx_to_type = {}
    all_discovered_types_from_prop = set() # Store types found via property

    for atom in mol.GetAtoms():
        try:
            # --- USE GetProp("bead_type") INSTEAD of GetSymbol() ---
            abstract_type = atom.GetProp("bead_type")
            idx_to_type[atom.GetIdx()] = abstract_type
            all_discovered_types_from_prop.add(abstract_type) # Track found types
        except KeyError:
            # Fallback if "bead_type" is missing (shouldn't happen with current build script)
            print(f"Warning: Atom {atom.GetIdx()} missing 'bead_type' property. Defaulting to 'A'.")
            idx_to_type[atom.GetIdx()] = 'A'
            all_discovered_types_from_prop.add('A')

    # --- Auto-populate missing atom types (using the newly discovered types) ---
    defined_atom_types = set(cg_params['atom_types'].keys())
    # Use the set derived from GetProp
    undefined_types = all_discovered_types_from_prop - defined_atom_types
    
    if undefined_types:
        print(f"INFO: Auto-defining parameters for new monomer types: {', '.join(undefined_types)}")
        if 'A' in cg_params['atom_types']:
            template_params = cg_params['atom_types']['A']
            for new_type in undefined_types:
                cg_params['atom_types'][new_type] = template_params.copy()
        else:
            print("Warning: Template type 'A' not found. Cannot auto-define new types.")
    # --- <<< END NEW SECTION >>>

    # --- Dynamically create integer mappings (unchanged) ---
    atom_type_map = { t: i + 1 for i, t in enumerate(sorted(list(cg_params['atom_types'].keys()))) }
    bond_type_map, angle_type_map = {}, {}

    # --- Verify the final map (Optional Debug Print) ---
    print(f"Final LAMMPS Atom Type Map: {atom_type_map}")

    # Extract bond and angle data (this logic is now robust to new types)
    # ... (bond and angle extraction loops are the same as the previous version) ...
    bond_data = []
    for bond in mol.GetBonds():
        type1, type2 = idx_to_type[bond.GetBeginAtomIdx()], idx_to_type[bond.GetEndAtomIdx()]
        bond_key = tuple(sorted((type1, type2)))
        params = get_params_with_default(cg_params['bond_types'], bond_key)
        if not params: continue

        if bond_key not in bond_type_map:
            bond_type_map[bond_key] = len(bond_type_map) + 1
        lammps_type_id = bond_type_map[bond_key]
        bond_data.append((lammps_type_id, bond.GetBeginAtomIdx()+1, bond.GetEndAtomIdx()+1))
    
    angle_data = []
    if PROCESS_ANGLES:
        for i in range(mol.GetNumAtoms()):
            neighbors = mol.GetAtomWithIdx(i).GetNeighbors()
            if len(neighbors) < 2: continue
            
            for j_atom, k_atom in itertools.combinations(neighbors, 2):
                idx_j, idx_i, idx_k = j_atom.GetIdx(), i, k_atom.GetIdx()
                type_j, type_i, type_k = idx_to_type[idx_j], idx_to_type[idx_i], idx_to_type[idx_k]
                angle_key = (type_j, type_i, type_k)
                params = get_params_with_default(cg_params['angle_types'], angle_key)
                if not params: continue

                if angle_key not in angle_type_map and angle_key[::-1] not in angle_type_map:
                    angle_type_map[angle_key] = len(angle_type_map) + 1
                
                stored_key = angle_key if angle_key in angle_type_map else angle_key[::-1]
                lammps_type_id = angle_type_map[stored_key]
                angle_data.append((lammps_type_id, idx_j + 1, idx_i + 1, idx_k + 1))
    
    # =======================================================================
    # --- 3. WRITE LAMMPS DATA FILE ---
    # =======================================================================
    # The file writing logic is unchanged and will work with the auto-populated types.
    # ... (code for writing the file is the same as the previous version) ...
    num_beads, num_bonds, num_angles = len(idx_to_type), len(bond_data), len(angle_data)
    num_atom_types, num_bond_types, num_angle_types = len(atom_type_map), len(bond_type_map), len(angle_type_map)
    try:
        with open(filename, 'w') as f:
            f.write("LAMMPS data file for Kremer-Grest Copolymer Network\n\n")
            f.write(f"{num_beads} atoms\n{num_bonds} bonds\n{num_angles} angles\n\n")
            f.write(f"{num_atom_types} atom types\n")
            f.write(f"{num_bond_types if num_bond_types > 0 else 1} bond types\n")
            f.write(f"{num_angle_types if num_angle_types > 0 else 1} angle types\n\n")
            f.write("-1000.0 1000.0 xlo xhi\n-1000.0 1000.0 ylo yhi\n-1000.0 1000.0 zlo zhi\n\n")
            
            f.write("Masses\n\n")
            for abstract_type, lammps_id in sorted(atom_type_map.items(), key=lambda item: item[1]):
                mass = cg_params['atom_types'][abstract_type]['mass']
                f.write(f"{lammps_id} {mass:.1f}\n")
            f.write("\n")
            
            f.write("Pair Coeffs # lj/cut\n\n")
            for abstract_type, lammps_id in sorted(atom_type_map.items(), key=lambda item: item[1]):
                p = cg_params['atom_types'][abstract_type]
                f.write(f"{lammps_id} {p['epsilon']:.2f} {p['sigma']:.2f}\n")
            f.write("\n")

            if num_bond_types > 0:
                f.write("Bond Coeffs # fene\n\n")
                for bond_key, lammps_id in sorted(bond_type_map.items(), key=lambda item: item[1]):
                    p = get_params_with_default(cg_params['bond_types'], bond_key)
                    f.write(f"{lammps_id} {p['K']:.2f} {p['R0']:.2f} {p['epsilon']:.2f} {p['sigma']:.2f}\n")
                f.write("\n")
            
            if num_angle_types > 0:
                f.write("Angle Coeffs # harmonic\n\n")
                for angle_key, lammps_id in sorted(angle_type_map.items(), key=lambda item: item[1]):
                    p = get_params_with_default(cg_params['angle_types'], angle_key)
                    f.write(f"{lammps_id} {p['K']:.2f} {p['Theta0']:.1f}\n")
                f.write("\n")
            
            f.write("Atoms # full\n\n")
            # We need to re-map bead_data here since it was built before auto-population
            for i in range(num_beads):
                bead_id = i + 1
                mol_id = 1
                type_id = atom_type_map[idx_to_type[i]]
                x, y, z = 0.0, 0.0, 0.0
                f.write(f"{bead_id} {mol_id} {type_id} 0.0 {x:.6f} {y:.6f} {z:.6f}\n")
            f.write("\n")
            
            if num_bonds > 0:
                f.write("Bonds\n\n")
                for i, (type_id, bead1, bead2) in enumerate(bond_data):
                    f.write(f"{i + 1} {type_id} {bead1} {bead2}\n")
                f.write("\n")
            
            if num_angles > 0:
                f.write("Angles\n\n")
                for i, (type_id, bead1, bead2, bead3) in enumerate(angle_data):
                    f.write(f"{i + 1} {type_id} {bead1} {bead2} {bead3}\n")
            
        print(f"Generated Kremer-Grest CG LAMMPS data file: {filename}")
        return {'success': True}

    except (IOError, KeyError) as e:
        print(f"Error during file generation: {e}")
        import traceback
        traceback.print_exc()
        return {'success': False}