import itertools
from pathlib import Path
from rdkit import Chem
import math

# Bundled parameter file shipped with the package.
_BUNDLED_PARAM_FILE = Path(__file__).parent / "DreidingX6parameters.txt"


def create_lammps_data_file(mol, filename="polymer_network.lammps", dreiding_param_file=None):
    """
    Generate a LAMMPS data file with automated DREIDING force field parameterization
    
    Args:
        mol: RDKit molecule object (chemical_space_withH)
        filename: Output LAMMPS data file name
        dreiding_param_file: Path to DreidingX6parameters.txt file.
            If None, uses the file bundled with the topon package.
    """
    if dreiding_param_file is None:
        dreiding_param_file = str(_BUNDLED_PARAM_FILE)

    # Parse the DREIDING parameter file
    dreiding_params = parse_dreiding_parameter_file(dreiding_param_file)

    # Count atoms and bonds
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()

    # Map RDKit atoms to DREIDING atom types
    atom_types_dict, atom_data, atom_dreiding_types = assign_atom_types(mol, dreiding_params)
    # Extract bonded interactions with parameters from DreidingX6parameters.txt
    bond_types, bond_data = extract_bonds(mol, atom_dreiding_types, dreiding_params)
    angle_types, angle_data = extract_angles(mol, atom_dreiding_types, dreiding_params)
    dihedral_types, dihedral_data = extract_dihedrals(mol, atom_dreiding_types, dreiding_params)
    improper_types, improper_data = extract_impropers(mol, atom_dreiding_types, dreiding_params)
    
    # Write LAMMPS data file
    with open(filename, 'w') as f:
        # Header
        f.write("LAMMPS data file for polymer network with DREIDING parameters\n\n")
        
        # Counts
        f.write(f"{num_atoms} atoms\n")
        f.write(f"{len(bond_data)} bonds\n")
        f.write(f"{len(angle_data)} angles\n")
        f.write(f"{len(dihedral_data)} dihedrals\n")
        f.write(f"{len(improper_data)} impropers\n\n")
        
        f.write(f"{len(atom_types_dict)} atom types\n")
        f.write(f"{len(bond_types)} bond types\n")
        f.write(f"{len(angle_types)} angle types\n")
        f.write(f"{len(dihedral_types)} dihedral types\n")
        f.write(f"{len(improper_types)} improper types\n\n")
        
        # Box dimensions
        calc_simulation_box(f, atom_data)
        
        # Write force field parameters
        write_masses_section(f, atom_types_dict, dreiding_params)
        write_pair_coeffs_section(f, atom_types_dict, dreiding_params)
        write_bond_coeffs_section(f, bond_types, dreiding_params)
        write_angle_coeffs_section(f, angle_types, dreiding_params)
        write_dihedral_coeffs_section(f, dihedral_types, dreiding_params)
        write_improper_coeffs_section(f, improper_types, dreiding_params)
        
        # Write atoms and topology
        write_atoms_section(f, atom_data)
        write_bonds_section(f, bond_data)
        write_angles_section(f, angle_data)
        write_dihedrals_section(f, dihedral_data)
        write_impropers_section(f, improper_data)
    
    print(f"Generated LAMMPS data file: {filename}")
    print(f"  Atoms: {num_atoms}")
    print(f"  Bonds: {len(bond_data)}")
    print(f"  Angles: {len(angle_data)}")
    print(f"  Dihedrals: {len(dihedral_data)}")
    print(f"  Impropers: {len(improper_data)}")
    
    return filename

def parse_dreiding_parameter_file(param_file):
    """
    Parse the DreidingX6parameters.txt file into structured dictionaries
    """
    params = {
        'atom_types': {},
        'vdw_params': {},
        'bond_params': {},
        'angle_params': {},
        'dihedral_params': {},
        'improper_params': {},
        'element_to_type': {}
    }
    
    current_section = None
    
    with open(param_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
            
            # Check for section headers (ending with 'END' marks the end of a section)
            if line == 'END':
                current_section = None
                continue
            
            # Check for section headers
            if line == 'ATOMTYPES':
                current_section = "atom_types"
                continue
            elif line == 'DIAGONAL_VDW':
                current_section = "vdw_params"
                continue
            elif line == 'BOND_STRETCH':
                current_section = "bond_params"
                continue
            elif line == 'ANGLE_BEND':
                current_section = "angle_params"
                continue
            elif line == 'TORSIONS' or line == 'TORSION_STRETCH':  # Look for torsion section
                current_section = "dihedral_params"
                continue
            elif line == 'INVERSIONS' or line == 'OUT_OF_PLANE':  # Look for improper section
                current_section = "improper_params"
                continue
            
            # Skip comment lines
            if line.startswith('#'):
                continue
                
            # Parse data based on current section
            if current_section == "atom_types":
                parts = line.split()
                if len(parts) >= 5:
                    atom_type = parts[0]
                    element = parts[1]
                    mass = float(parts[2])
                    charge = float(parts[3])
                    
                    params['atom_types'][atom_type] = {
                        'element': element,
                        'mass': mass,
                        'charge': charge,
                        'hybridization': parts[4:]
                    }
                    
                    # Extract element for mapping
                    if element not in params['element_to_type']:
                        params['element_to_type'][element] = []
                    params['element_to_type'][element].append(atom_type)
            
            elif current_section == "vdw_params":
                parts = line.split()
                if len(parts) >= 4:
                    atom_type = parts[0]
                    vdw_type = parts[1]  # e.g., EXPO_6
                    vdw_radius = float(parts[2])
                    epsilon = float(parts[3])
                    
                    params['vdw_params'][atom_type] = {
                        'type': vdw_type,
                        'radius': vdw_radius,
                        'epsilon': epsilon
                    }
                    # Add additional parameters if present
                    if len(parts) > 4:
                        params['vdw_params'][atom_type]['additional'] = parts[4:]
            
            elif current_section == "bond_params":
                parts = line.split()
                if len(parts) >= 5:
                    type1, type2 = parts[0], parts[1]
                    bond_type = parts[2]  # e.g., HARMONIC
                    k = float(parts[3])
                    r0 = float(parts[4])
                    
                    key = (type1, type2)
                    params['bond_params'][key] = {
                        'type': bond_type,
                        'k': k,
                        'r0': r0
                    }
                    # Add reverse order too
                    params['bond_params'][(type2, type1)] = {
                        'type': bond_type,
                        'k': k,
                        'r0': r0
                    }
            
            elif current_section == "angle_params":
                parts = line.split()
                if len(parts) >= 6:
                    type1, type2, type3 = parts[0], parts[1], parts[2]
                    angle_type = parts[3]  # e.g., COS_HARMON
                    k = float(parts[4])
                    theta0 = float(parts[5])
                    
                    key = (type1, type2, type3)
                    params['angle_params'][key] = {
                        'type': angle_type,
                        'k': k,
                        'theta0': theta0
                    }
                    # Add reverse order too
                    params['angle_params'][(type3, type2, type1)] = {
                        'type': angle_type,
                        'k': k,
                        'theta0': theta0
                    }
    
            elif current_section == "dihedral_params":
                # --- MODIFIED SECTION TO HANDLE MULTI-TERM DIHEDRALS ---
                parts = line.split()
                if len(parts) >= 7:
                    type1, type2, type3, type4 = parts[0], parts[1], parts[2], parts[3]
                    torsion_type = parts[4]
                    key = (type1, type2, type3, type4)
                    
                    # Initialize as a list to hold multiple terms
                    if key not in params['dihedral_params']:
                        params['dihedral_params'][key] = []
                    
                    param_values = [float(p) for p in parts[5:]]
                    # Iterate through parameter sets in chunks of 3 (Vn, n, d)
                    for i in range(0, len(param_values), 3):
                        try:
                            term_params = {
                                'type': torsion_type,
                                'v_n': param_values[i],
                                'n': int(param_values[i+1]),
                                'd': int(param_values[i+2])
                            }
                            params['dihedral_params'][key].append(term_params)
                        except (IndexError, ValueError):
                            # This handles cases where a line might be malformed
                            # or doesn't have a full set of 3 parameters at the end.
                            break
                    
                    # Apply the same list of parameter dictionaries to the reverse key
                    params['dihedral_params'][(type4, type3, type2, type1)] = params['dihedral_params'][key]
            
            elif current_section == "improper_params":
                # New parser code for impropers (inversions)
                parts = line.split()
                if len(parts) >= 6:
                    central_type = parts[0]  # Central atom type
                    type1, type2, type3 = parts[1], parts[2], parts[3]
                    improper_type = parts[4]  # e.g., HARMONIC or UMBRELLA
                    k = float(parts[5])  # Force constant
                    
                    # Some improper types might have additional parameters
                    additional = {}
                    if len(parts) > 6:
                        additional['chi0'] = float(parts[6])  # For equilibrium value
                    
                    key = (central_type, type1, type2, type3)
                    params['improper_params'][key] = {
                        'type': improper_type,
                        'k': k,
                        'additional': additional
                    }
    
    # Add some general DREIDING wildcards for missing parameters
    add_wildcard_parameters(params)
    
    return params

def add_wildcard_parameters(params):
    """Add wildcard parameters for cases not explicitly defined.

    All wildcards use the same dict format as the parser so that consumers
    can access parameters uniformly (e.g. ``p['k']``) regardless of whether
    the key came from the file or from a wildcard fallback.
    """
    # Generic/wildcard bond parameters (X-X)
    if ('X', 'X') not in params['bond_params']:
        params['bond_params'][('X', 'X')] = {'type': 'HARMONIC', 'k': 700.0, 'r0': 1.5}

    # Generic angle parameters (X-X-X)
    if ('X', 'X', 'X') not in params['angle_params']:
        params['angle_params'][('X', 'X', 'X')] = {'type': 'COS_HARMON', 'k': 100.0, 'theta0': 109.5}

    # Generic dihedral parameters (X-X-X-X) — list of term dicts (matches parser)
    if ('X', 'X', 'X', 'X') not in params['dihedral_params']:
        params['dihedral_params'][('X', 'X', 'X', 'X')] = [{'type': 'TORSION', 'v_n': 0.0, 'n': 0, 'd': 0}]

    # Generic improper parameters (X-X-X-X)
    if ('X', 'X', 'X', 'X') not in params['improper_params']:
        params['improper_params'][('X', 'X', 'X', 'X')] = {'type': 'HARMONIC', 'k': 40.0, 'additional': {'chi0': 0.0}}

def assign_atom_types(mol, dreiding_params):
    """Map RDKit atoms to DREIDING atom types based on element and hybridization"""
    atom_types_dict = {}
    atom_data = []
    atom_dreiding_types = {}
    
    # Hybridization mapping to DREIDING types
    hybridization_map = {
        'SP3': '3',
        'SP2': '2',
        'SP': '1',
        'S': '_'
    }
    
    for atom in mol.GetAtoms():
        idx = atom.GetIdx() + 1  # LAMMPS uses 1-based indexing
        element = atom.GetSymbol()
        if element == "H":
            element = element + "_"
        #if element == "F":
        #    element = element + "_"
        hyb = str(atom.GetHybridization())
        
        # Map to DREIDING type
        if hyb in hybridization_map:
            hyb_suffix = hybridization_map[hyb]
            possible_type = f"{element}{hyb_suffix}" if element != "H" else "H_"
        else:
            possible_type = f"{element}_UNSPECIFIED"
            if element == "H_":
                possible_type = "H_"
        
        # Check if this type exists in parameters
        dreiding_type = None
        if possible_type in dreiding_params['atom_types']:
            dreiding_type = possible_type
        else:
            # Create version with underscore before first digit
            underscore_type = None
            for i, char in enumerate(possible_type):
                if char.isdigit():
                    underscore_type = possible_type[:i] + '_' + possible_type[i:]
                    break
            
            # Check if underscore version exists
            if underscore_type and underscore_type in dreiding_params['atom_types']:
                dreiding_type = underscore_type
        
            elif element in dreiding_params['atom_types']:
                dreiding_type = element
        
            else:
                underscore_type = element + '_'
                dreiding_type = underscore_type

        # Get position from RDKit
        x, y, z = 0.0, 0.0, 0.0
        
        # Get charge if available
        charge = 0.0
        if atom.HasProp('_GasteigerCharge'):
            charge = atom.GetDoubleProp('_GasteigerCharge')
        
        # Add to atom types
        if dreiding_type not in atom_types_dict:
            atom_types_dict[dreiding_type] = len(atom_types_dict) + 1
        
        # Store atom data
        atom_data.append((idx, atom_types_dict[dreiding_type], charge, x, y, z, element, hyb))
        atom_dreiding_types[idx] = dreiding_type
    
    return atom_types_dict, atom_data, atom_dreiding_types

def find_parameter(key_tuple, param_dict, wildcard='X'):
    """Find the most specific parameter match, falling back to wildcards"""
    # Try exact match
    if key_tuple in param_dict:
        return param_dict[key_tuple]
    
    # Try with wildcards (progressively replacing specific types with X)
    n = len(key_tuple)
    
    # Create all possible partial wildcard combinations
    for num_wildcards in range(1, n+1):
        for positions in itertools.combinations(range(n), num_wildcards):
            wildcard_key = list(key_tuple)
            for pos in positions:
                wildcard_key[pos] = wildcard
            wildcard_key = tuple(wildcard_key)
            
            if wildcard_key in param_dict:
                return param_dict[wildcard_key]
    
    # If no match found, return defaults
    if n == 2:  # Bond
        return (700.0, 1.5)
    elif n == 3:  # Angle
        return (100.0, 109.5)

def extract_bonds(mol, atom_dreiding_types, dreiding_params):
    """Extract bonds with parameters from DreidingX6parameters.txt"""
    bond_types = {}
    bond_data = []
    
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx() + 1
        end_idx = bond.GetEndAtomIdx() + 1
        
        type1 = atom_dreiding_types[begin_idx]
        type2 = atom_dreiding_types[end_idx]
        
        # --- START FIX ---
        # Create a canonical key by sorting the atom types.
        # This ensures ('Si3', 'C_3') is treated the same as ('C_3', 'Si3').
        canonical_key = tuple(sorted((type1, type2)))
        # --- END FIX ---
        
        # Use the canonical key for parameter lookup
        params = find_parameter(canonical_key, dreiding_params['bond_params'])
        
        if isinstance(params, dict):
            k = 0.5 * params['k']
            r0 = params['r0']
        else:
            k, r0 = params
        
        # --- START FIX ---
        # Use the canonical key to define the bond type.
        type_key = (canonical_key[0], canonical_key[1], k, r0)
        # --- END FIX ---

        if type_key not in bond_types:
            bond_types[type_key] = len(bond_types) + 1
        
        bond_data.append((len(bond_data) + 1, bond_types[type_key], begin_idx, end_idx))
    
    return bond_types, bond_data

def extract_angles(mol, atom_dreiding_types, dreiding_params):
    """Extract angles with parameters from DreidingX6parameters.txt"""
    angle_types = {}
    angle_data = []
    
    for atom2 in mol.GetAtoms():
        atom2_idx = atom2.GetIdx() + 1
        neighbors = [n.GetIdx() + 1 for n in atom2.GetNeighbors()]
        
        if len(neighbors) < 2:
            continue
        
        for i, atom1_idx in enumerate(neighbors):
            for atom3_idx in neighbors[i+1:]:
                type1 = atom_dreiding_types[atom1_idx]
                type2 = atom_dreiding_types[atom2_idx]
                type3 = atom_dreiding_types[atom3_idx]
                
                # --- START FIX ---
                # Create a canonical key by sorting the outer atom types.
                # This ensures ('C_3', 'Si3', 'O_3') is the same as ('O_3', 'Si3', 'C_3').
                outer_types = sorted((type1, type3))
                canonical_key = (outer_types[0], type2, outer_types[1])
                # --- END FIX ---

                # Use the canonical key for parameter lookup
                params = find_parameter(canonical_key, dreiding_params['angle_params'])

                if isinstance(params, dict):
                    k_harmonic = params['k']
                    theta0 = params['theta0']
                else:
                    k_harmonic, theta0 = params
                
                k_harmonic = 0.5 * k_harmonic
                
                # --- START FIX ---
                # Use the canonical key to define the angle type.
                type_key = (canonical_key[0], canonical_key[1], canonical_key[2], k_harmonic, theta0)
                # --- END FIX ---

                if type_key not in angle_types:
                    angle_types[type_key] = len(angle_types) + 1
                
                angle_data.append((len(angle_data) + 1, angle_types[type_key], 
                                   atom1_idx, atom2_idx, atom3_idx))
    
    return angle_types, angle_data

def extract_dihedrals(mol, atom_dreiding_types, dreiding_params):
    """Extract dihedrals with parameters from DreidingX6parameters.txt"""
    dihedral_types = {}
    dihedral_data = []
    
    central_bond_count = {}
    dihedral_list = []
    
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        
        begin_idx = begin_atom.GetIdx() + 1
        end_idx = end_atom.GetIdx() + 1
        
        if begin_atom.GetDegree() == 1 or end_atom.GetDegree() == 1:
            continue
        
        begin_neighbors = [n.GetIdx() + 1 for n in begin_atom.GetNeighbors() if n.GetIdx() + 1 != end_idx]
        end_neighbors = [n.GetIdx() + 1 for n in end_atom.GetNeighbors() if n.GetIdx() + 1 != begin_idx]
        
        central_bond_key = tuple(sorted([begin_idx, end_idx]))
        
        count = 0
        for atom1_idx in begin_neighbors:
            for atom4_idx in end_neighbors:
                # Skip degenerate dihedrals (e.g. in 3-membered rings like epoxides)
                if atom1_idx == atom4_idx:
                    continue
                count += 1
                dihedral_list.append((atom1_idx, begin_idx, end_idx, atom4_idx, central_bond_key))
        central_bond_count[central_bond_key] = count
    
    for atom1_idx, begin_idx, end_idx, atom4_idx, central_bond_key in dihedral_list:
        type1 = atom_dreiding_types[atom1_idx]
        type2 = atom_dreiding_types[begin_idx]
        type3 = atom_dreiding_types[end_idx]
        type4 = atom_dreiding_types[atom4_idx]
        
        # --- START FIX ---
        # Create a canonical key by comparing the forward and reverse atom type sequences.
        forward_key = (type1, type2, type3, type4)
        reverse_key = (type4, type3, type2, type1)
        canonical_key = min(forward_key, reverse_key)
        # --- END FIX ---
        
        # Use the canonical key for parameter lookup
        param_list = find_parameter(canonical_key, dreiding_params['dihedral_params'])
        
        if isinstance(param_list, list):
            for params in param_list:
                v_n = params.get('v_n', 0.0)
                n_value = int(params.get('n', 0))
                d_value = int(params.get('d', 0))
                
                num_possibilities = central_bond_count[central_bond_key]
                k_value = (0.5 * v_n) / num_possibilities if num_possibilities > 0 else 0.5 * v_n
                
                lammps_params = (k_value, n_value, d_value)
                
                # --- START FIX ---
                # Use the canonical key to define the dihedral type.
                type_key = canonical_key + lammps_params
                # --- END FIX ---
                
                if type_key not in dihedral_types:
                    dihedral_types[type_key] = len(dihedral_types) + 1
                
                dihedral_data.append((len(dihedral_data) + 1, dihedral_types[type_key],
                                      atom1_idx, begin_idx, end_idx, atom4_idx))

    return dihedral_types, dihedral_data

def extract_impropers(mol, atom_dreiding_types, dreiding_params):
    """Extract impropers with parameters from DreidingX6parameters.txt"""
    improper_types = {}
    improper_data = []
    
    # Look for atoms that might need impropers (trigonal centers, etc.)
    for atom in mol.GetAtoms():
        if atom.GetDegree() < 3:
            continue
            
        central_idx = atom.GetIdx() + 1
        central_type = atom_dreiding_types[central_idx]
        
        # Get neighbors
        neighbors = [n.GetIdx() + 1 for n in atom.GetNeighbors()]
        
        # SP2 centers (trigonal) - one improper for planarity
        if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2 and len(neighbors) == 3:
            # For SP2 centers, the key is typically just the central atom type with wildcards
            key = (central_type, 'X', 'X', 'X')
            params = find_parameter(key, dreiding_params['improper_params'])
            
            # Extract parameters based on the DREIDING format shown in the image
            if isinstance(params, dict):
                improper_type = params.get('type', 'UMBRELLA')
                k = params.get('k', 40.0)
                chi0 = params.get('additional', {}).get('chi0', 0.0)
            else:
                # Assuming params is a tuple if not a dict
                k, chi0 = params
                improper_type = 'UMBRELLA'  # Default type
            
            # Add to improper types
            type_key = (central_type, 'X', 'X', 'X', improper_type, k, chi0)
            if type_key not in improper_types:
                improper_types[type_key] = len(improper_types) + 1
            
            # Store improper data
            improper_data.append((len(improper_data) + 1, improper_types[type_key],
                                central_idx, neighbors[0], neighbors[1], neighbors[2]))

    return improper_types, improper_data


def calc_simulation_box(f, atom_data):
    """Calculate and write simulation box dimensions"""
    f.write("-1000.0 1000.0 xlo xhi\n")
    f.write("-1000.0 1000.0 ylo yhi\n")
    f.write("-1000.0 1000.0 zlo zhi\n")
    f.write("0.0 0.0 0.0 xy xz yz\n\n")
    
def write_masses_section(f, atom_types, dreiding_params):
    """Write the Masses section to LAMMPS data file"""
    f.write("Masses\n\n")
    for type_name, type_id in sorted(atom_types.items(), key=lambda x: x[1]):
        # Get mass from parameters or default
        if type_name in dreiding_params['atom_types']:
            mass = dreiding_params['atom_types'][type_name]['mass']
        else:
            # Extract element symbol and use standard mass
            element = ''.join([c for c in type_name if c.isalpha() and c.isupper()])
            if element == 'H':
                mass = 1.008
            elif element == 'C':
                mass = 12.011
            elif element == 'O':
                mass = 15.999
            elif element == 'Si':
                mass = 28.086
            elif element == 'F':
                mass = 18.998  # Added fluorine mass
            else:
                mass = 0.0  # Default
        
        f.write(f"{type_id} {mass:.4f}  # {type_name}\n")
    f.write("\n")

def write_pair_coeffs_section(f, atom_types_dict, dreiding_params):
    """Write pair coefficients section for LAMMPS data file"""
    f.write("\nPair Coeffs\n\n")
    
    for type_name, type_idx in atom_types_dict.items():
        # Check if the atom type exists in parameters
        if type_name in dreiding_params['vdw_params']:
            atom_params = dreiding_params['vdw_params'][type_name]
            
            # Get epsilon with fallback to default
            if isinstance(atom_params, dict) and 'epsilon' in atom_params:
                epsilon = atom_params['epsilon']
            else:
                # Default value for epsilon
                epsilon = 0.001  # Typical small default
                print(f"Warning: epsilon not found for {type_name}, using default value")
            
            # Get sigma with fallback to default
            if isinstance(atom_params, dict) and 'radius' in atom_params:
                sigma = atom_params['radius']
            else:
                # Default value for sigma
                sigma = 3.5  # Typical carbon-like default
                print(f"Warning: sigma not found for {type_name}, using default value")
        else:
            # Atom type not found at all
            print(f"Warning: No parameters found for atom type {type_name}, using defaults")
            epsilon = 0.001
            sigma = 3.5
        
        # Write the pair coefficients
        f.write(f"{type_idx} {epsilon} {sigma}\n")


def write_bond_coeffs_section(f, bond_types, dreiding_params):
    """Write bond coefficients section for LAMMPS data file"""
    f.write("\nBond Coeffs\n\n")
    
    for (type1, type2, k, r0), type_idx in bond_types.items():
        # We already have k and r0 from extract_bonds function
        f.write(f"{type_idx} {k} {r0}\n")

def write_angle_coeffs_section(f, angle_types, dreiding_params):
    """Write angle coefficients section for LAMMPS data file"""
    f.write("\nAngle Coeffs\n\n")
    
    for (type1, type2, type3, k, theta0), type_idx in angle_types.items():
        # We already have k and theta0 from extract_angles function
        f.write(f"{type_idx} {k} {theta0}\n")

def write_dihedral_coeffs_section(f, dihedral_types, dreiding_params):
    """Write dihedral coefficients section for LAMMPS data file"""
    f.write("\nDihedral Coeffs\n\n")
    
    # Write header line specifying dihedral style
    #f.write("# Harmonic dihedral style: K [1 + d*cos(n*phi)]\n")
    
    for (type1, type2, type3, type4, k, n, d), type_idx in dihedral_types.items():
        # Format: coeff_id k_value n_value d_value
        di = int(d)
        f.write(f"{type_idx} {k:.6f} {di} {n} # {type1} {type2} {type3} {type4}\n")

def write_improper_coeffs_section(f, improper_types, dreiding_params):
    """Write improper coefficients section for LAMMPS data file"""
    f.write("\nImproper Coeffs\n\n")
    
    # Write header line specifying improper style
    #f.write("# cvff improper style: K(1 + d*cos(n*phi))\n")
    
    for (type1, type2, type3, type4, improper_type, k, chi0), type_idx in improper_types.items():
        # For UMBRELLA type in DREIDING, use cvff style in LAMMPS
        # Convert parameters if needed
        # LAMMPS cvff format: K d n  →  E = K[1 + d*cos(n*chi)]
        # DREIDING inversion: E = K(1 - cos(chi))  →  d=-1, n=1
        f.write(f"{type_idx} {k:.6f} -1 1  # {type1} {type2} {type3} {type4} ({improper_type})\n")

def write_atoms_section(f, atom_data):
    """Write the Atoms section to LAMMPS data file"""
    f.write("\n")
    f.write("Atoms # full\n\n")
    for idx, type_id, charge, x, y, z, element, hyb in atom_data:
        f.write(f"{idx} 1 {type_id} {charge:.6f} {x:.6f} {y:.6f} {z:.6f}  # {element} ({hyb})\n")
    f.write("\n")

def write_bonds_section(f, bond_data):
    """Write the Bonds section to LAMMPS data file"""
    if not bond_data:
        return
        
    f.write("Bonds\n\n")
    for bond_id, type_id, atom1, atom2 in bond_data:
        f.write(f"{bond_id} {type_id} {atom1} {atom2}\n")
    f.write("\n")

def write_angles_section(f, angle_data):
    """Write the Angles section to LAMMPS data file"""
    if not angle_data:
        return
        
    f.write("Angles\n\n")
    for angle_id, type_id, atom1, atom2, atom3 in angle_data:
        f.write(f"{angle_id} {type_id} {atom1} {atom2} {atom3}\n")
    f.write("\n")

def write_dihedrals_section(f, dihedral_data):
    """Write the Dihedrals section to LAMMPS data file"""
    if not dihedral_data:
        return
        
    f.write("Dihedrals\n\n")
    for dihedral_id, type_id, atom1, atom2, atom3, atom4 in dihedral_data:
        f.write(f"{dihedral_id} {type_id} {atom1} {atom2} {atom3} {atom4}\n")
    f.write("\n")

def write_impropers_section(f, improper_data):
    """Write the Impropers section to LAMMPS data file"""
    if not improper_data:
        return
        
    f.write("Impropers\n\n")
    for improper_id, type_id, atom1, atom2, atom3, atom4 in improper_data:
        f.write(f"{improper_id} {type_id} {atom1} {atom2} {atom3} {atom4}\n")
