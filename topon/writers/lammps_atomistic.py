import os
import itertools
from rdkit import Chem
from ..utils import network_config

class DreidingWriter:
    """
    A modular Writer that converts an RDKit molecule into LAMMPS Atomistic input files
    using the EXACT logic from the legacy 'make_dreiding.py' script.
    """

    def __init__(self, mol, output_file, use_charges=False):
        self.mol = mol
        self.output_file = output_file
        self.use_charges = use_charges
        
        # Paths
        self.output_dir = os.path.dirname(output_file)
        self.basename = os.path.splitext(os.path.basename(output_file))[0]
        self.settings_file = f"{os.path.splitext(output_file)[0]}.in.settings"
        
        # Locate parameter file
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.param_file = os.path.join(base_dir, "utils", network_config.DREIDING_PARAM_FILE)
        
        # Containers matching legacy script structure
        self.params = {
            'atom_types': {},
            'vdw_params': {},
            'bond_params': {},
            'angle_params': {},
            'dihedral_params': {},
            'improper_params': {},
            'element_to_type': {}
        }
        
        # Mappings
        self.atom_types_dict = {}      # {'C_3': 1}
        self.atom_dreiding_types = {}  # {atom_idx: 'C_3'}
        self.atom_data = []            # List of atom lines
        
        self.bond_types = {}
        self.bond_data = []
        
        self.angle_types = {}
        self.angle_data = []
        
        self.dihedral_types = {}
        self.dihedral_data = []
        
        self.improper_types = {}
        self.improper_data = []

    def write(self):
        """Main execution method."""
        if not os.path.exists(self.param_file):
            raise FileNotFoundError(f"Dreiding parameter file not found at: {self.param_file}")

        print(f"Parsing forcefield parameters from {os.path.basename(self.param_file)}...")
        self._parse_dreiding_params()
        
        print(f"Assigning atom types (Charges included: {self.use_charges})...")
        self._assign_atom_types()
        
        print("Extracting topology...")
        self._extract_bonds()
        self._extract_angles()
        self._extract_dihedrals()
        self._extract_impropers()
        
        print(f"Writing LAMMPS Data File: {self.output_file}")
        self._write_data_file()
        
        print(f"Writing LAMMPS Settings File: {self.settings_file}")
        self._write_settings_file()
        
        print("Write complete.")

    # =========================================================================
    # --- 1. PARAMETER PARSING ---
    # =========================================================================

    def _parse_dreiding_params(self):
        current_section = None
        with open(self.param_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line == 'END': 
                    current_section = None
                    continue
                
                # Headers
                if line == 'ATOMTYPES': current_section = "atom_types"; continue
                elif line == 'DIAGONAL_VDW': current_section = "vdw_params"; continue
                elif line == 'BOND_STRETCH': current_section = "bond_params"; continue
                elif line == 'ANGLE_BEND': current_section = "angle_params"; continue
                elif line in ['TORSIONS', 'TORSION_STRETCH']: current_section = "dihedral_params"; continue
                elif line in ['INVERSIONS', 'OUT_OF_PLANE']: current_section = "improper_params"; continue
                if line.startswith('#'): continue
                
                parts = line.split()
                
                if current_section == "atom_types" and len(parts) >= 5:
                    self.params['atom_types'][parts[0]] = {
                        'element': parts[1], 
                        'mass': float(parts[2]), 
                        'charge': float(parts[3])
                    }
                    if parts[1] not in self.params['element_to_type']:
                        self.params['element_to_type'][parts[1]] = []
                    self.params['element_to_type'][parts[1]].append(parts[0])

                elif current_section == "vdw_params" and len(parts) >= 4:
                    self.params['vdw_params'][parts[0]] = {
                        'type': parts[1], 'radius': float(parts[2]), 'epsilon': float(parts[3])
                    }

                elif current_section == "bond_params" and len(parts) >= 5:
                    # Store both forward and reverse
                    key = (parts[0], parts[1])
                    val = {'type': parts[2], 'k': float(parts[3]), 'r0': float(parts[4])}
                    self.params['bond_params'][key] = val
                    self.params['bond_params'][(parts[1], parts[0])] = val

                elif current_section == "angle_params" and len(parts) >= 6:
                    key = (parts[0], parts[1], parts[2])
                    val = {'type': parts[3], 'k': float(parts[4]), 'theta0': float(parts[5])}
                    self.params['angle_params'][key] = val
                    self.params['angle_params'][(parts[2], parts[1], parts[0])] = val # Reverse lookup

                elif current_section == "dihedral_params" and len(parts) >= 7:
                    key = (parts[0], parts[1], parts[2], parts[3])
                    if key not in self.params['dihedral_params']: self.params['dihedral_params'][key] = []
                    
                    # Parse multi-term dihedrals (Vn, n, d)
                    try:
                        # Dreiding file structure: I J K L TYPE [Vn n d]...
                        # params start at index 5
                        param_values = [float(p) for p in parts[5:]]
                        for i in range(0, len(param_values), 3):
                            self.params['dihedral_params'][key].append({
                                'type': parts[4], 
                                'v_n': param_values[i],
                                'n': int(param_values[i+1]), 
                                'd': int(param_values[i+2])
                            })
                    except (IndexError, ValueError): pass
                    
                    # Mirror
                    self.params['dihedral_params'][(parts[3], parts[2], parts[1], parts[0])] = self.params['dihedral_params'][key]

                elif current_section == "improper_params" and len(parts) >= 6:
                    key = (parts[0], parts[1], parts[2], parts[3])
                    chi0 = float(parts[6]) if len(parts) > 6 else 0.0
                    self.params['improper_params'][key] = {
                        'type': parts[4], 'k': float(parts[5]), 'chi0': chi0
                    }

        # Wildcard Defaults
        self.params['bond_params'].setdefault(('X', 'X'), {'k': 700.0, 'r0': 1.5})
        self.params['angle_params'].setdefault(('X', 'X', 'X'), {'k': 100.0, 'theta0': 109.5})
        self.params['dihedral_params'].setdefault(('X', 'X', 'X', 'X'), [{'v_n': 0.0, 'n': 0, 'd': 1}])
        self.params['improper_params'].setdefault(('X', 'X', 'X', 'X'), {'type': 'UMBRELLA', 'k': 40.0, 'chi0': 0.0})

    # =========================================================================
    # --- 2. ATOM TYPING ---
    # =========================================================================

    def _assign_atom_types(self):
        hybridization_map = {'SP3': '3', 'SP2': '2', 'SP': '1', 'S': '_'}
        
        for atom in self.mol.GetAtoms():
            idx = atom.GetIdx() + 1
            element = atom.GetSymbol()
            if element == "H": element = "H_"
            
            hyb = str(atom.GetHybridization())
            
            # 1. Construct possible type string
            possible_type = f"{element}_UNSPECIFIED"
            if hyb in hybridization_map:
                suffix = hybridization_map[hyb]
                possible_type = f"{element}{suffix}" if element != "H_" else "H_"

            # 2. Matching Logic
            dreiding_type = None
            
            # A. Direct Match
            if possible_type in self.params['atom_types']:
                dreiding_type = possible_type
            else:
                # B. Try Underscore insertion (C3 -> C_3)
                underscore_type = None
                for i, char in enumerate(possible_type):
                    if char.isdigit():
                        underscore_type = possible_type[:i] + '_' + possible_type[i:]
                        break
                
                if underscore_type and underscore_type in self.params['atom_types']:
                    dreiding_type = underscore_type
                # C. Try Element Symbol
                elif element in self.params['atom_types']:
                    dreiding_type = element
                # D. Try Element_
                else:
                    dreiding_type = f"{element}_"

            # 3. Store
            if dreiding_type not in self.atom_types_dict:
                self.atom_types_dict[dreiding_type] = len(self.atom_types_dict) + 1
            
            self.atom_dreiding_types[idx] = dreiding_type
            
            # Charge logic
            charge = 0.0
            if self.use_charges and atom.HasProp('_GasteigerCharge'):
                try: charge = float(atom.GetProp('_GasteigerCharge'))
                except: pass
                
            self.atom_data.append({
                'id': idx, 
                'type_idx': self.atom_types_dict[dreiding_type],  # Standardized Key
                'charge': charge, 
                'x': 0.0, 'y': 0.0, 'z': 0.0,
                'label': dreiding_type # Standardized Key
            })

    # =========================================================================
    # --- 3. PARAMETER EXTRACTION ---
    # =========================================================================

    def _find_param(self, key_tuple, param_dict):
        if key_tuple in param_dict: return param_dict[key_tuple]
        n = len(key_tuple)
        for num_wildcards in range(1, n+1):
            for positions in itertools.combinations(range(n), num_wildcards):
                wildcard_key = list(key_tuple)
                for pos in positions: wildcard_key[pos] = 'X'
                if tuple(wildcard_key) in param_dict:
                    return param_dict[tuple(wildcard_key)]
        return None

    def _extract_bonds(self):
        for bond in self.mol.GetBonds():
            idx1, idx2 = bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1
            t1, t2 = self.atom_dreiding_types[idx1], self.atom_dreiding_types[idx2]
            
            # Lookup Key (Original order available in dict)
            key = tuple(sorted((t1, t2)))
            p = self._find_param(key, self.params['bond_params'])
            if not p: p = self.params['bond_params'][('X', 'X')]
            
            k_val = 0.5 * p['k']
            
            # Canonical Key for LAMMPS Type Definition
            canon_t1, canon_t2 = sorted((t1, t2))
            type_sig = (canon_t1, canon_t2, k_val, p['r0'])
            
            if type_sig not in self.bond_types:
                self.bond_types[type_sig] = len(self.bond_types) + 1
            self.bond_data.append((len(self.bond_data)+1, self.bond_types[type_sig], idx1, idx2))

    def _extract_angles(self):
        for atom in self.mol.GetAtoms():
            if atom.GetDegree() < 2: continue
            center_idx = atom.GetIdx() + 1
            t_center = self.atom_dreiding_types[center_idx]
            
            neighbors = [n.GetIdx()+1 for n in atom.GetNeighbors()]
            for n1, n2 in itertools.combinations(neighbors, 2):
                t1, t2 = self.atom_dreiding_types[n1], self.atom_dreiding_types[n2]
                
                # Lookup Key
                outer_sorted = sorted((t1, t2))
                key = (outer_sorted[0], t_center, outer_sorted[1])
                
                p = self._find_param(key, self.params['angle_params'])
                if not p: p = self.params['angle_params'][('X', 'X', 'X')]
                
                k_val = 0.5 * p['k']
                
                # Canonical Type Key
                type_sig = (key[0], t_center, key[2], k_val, p['theta0'])
                
                if type_sig not in self.angle_types:
                    self.angle_types[type_sig] = len(self.angle_types) + 1
                self.angle_data.append((len(self.angle_data)+1, self.angle_types[type_sig], n1, center_idx, n2))

    def _extract_dihedrals(self):
        # Pre-calculate multiplicity for renormalization
        central_bond_counts = {}
        for bond in self.mol.GetBonds():
            b, e = bond.GetBeginAtom(), bond.GetEndAtom()
            if b.GetDegree() > 1 and e.GetDegree() > 1:
                count = (len(b.GetNeighbors()) - 1) * (len(e.GetNeighbors()) - 1)
                idx1, idx2 = sorted((b.GetIdx()+1, e.GetIdx()+1))
                central_bond_counts[(idx1, idx2)] = count

        for bond in self.mol.GetBonds():
            b, e = bond.GetBeginAtom(), bond.GetEndAtom()
            if b.GetDegree() < 2 or e.GetDegree() < 2: continue
            
            b_idx, e_idx = b.GetIdx()+1, e.GetIdx()+1
            neighs_b = [n.GetIdx()+1 for n in b.GetNeighbors() if n.GetIdx()+1 != e_idx]
            neighs_e = [n.GetIdx()+1 for n in e.GetNeighbors() if n.GetIdx()+1 != b_idx]
            
            bond_key = tuple(sorted((b_idx, e_idx)))
            multiplicity = central_bond_counts.get(bond_key, 1)
            if multiplicity < 1: multiplicity = 1
            
            for a_idx in neighs_b:
                for d_idx in neighs_e:
                    ta, tb, tc, td = [self.atom_dreiding_types[x] for x in (a_idx, b_idx, e_idx, d_idx)]
                    
                    # Lookup
                    key = (ta, tb, tc, td)
                    params_list = self._find_param(key, self.params['dihedral_params'])
                    if not params_list:
                         # Try reverse if _find_param didn't handle it (it handles 'X' but not explicit reverse if dict missing)
                         params_list = self._find_param((td, tc, tb, ta), self.params['dihedral_params'])
                    if not params_list:
                        params_list = self.params['dihedral_params'][('X', 'X', 'X', 'X')]
                    
                    for p in params_list:
                        v_n = p['v_n']
                        # Renormalize: K = (0.5 * V) / Multiplicity
                        k_val = (0.5 * v_n) / multiplicity
                        
                        # Canonical Type Key
                        fwd = (ta, tb, tc, td)
                        rev = (td, tc, tb, ta)
                        canon_t = min(fwd, rev)
                        
                        type_sig = (canon_t, k_val, p['n'], p['d'])
                        
                        if type_sig not in self.dihedral_types:
                            self.dihedral_types[type_sig] = len(self.dihedral_types) + 1
                            
                        self.dihedral_data.append((len(self.dihedral_data)+1, self.dihedral_types[type_sig], a_idx, b_idx, e_idx, d_idx))

    def _extract_impropers(self):
        for atom in self.mol.GetAtoms():
            if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2 and atom.GetDegree() == 3:
                c_idx = atom.GetIdx() + 1
                neighbors = [n.GetIdx()+1 for n in atom.GetNeighbors()]
                neighbors.sort()
                
                c_type = self.atom_dreiding_types[c_idx]
                key = (c_type, 'X', 'X', 'X')
                p = self._find_param(key, self.params['improper_params'])
                if not p: p = self.params['improper_params'][('X', 'X', 'X', 'X')]
                
                type_sig = (c_type, p['type'], p['k'], p['chi0'])
                if type_sig not in self.improper_types:
                    self.improper_types[type_sig] = len(self.improper_types) + 1
                    
                self.improper_data.append((len(self.improper_data)+1, self.improper_types[type_sig], c_idx, neighbors[0], neighbors[1], neighbors[2]))

    # =========================================================================
    # --- 4. FILE WRITING ---
    # =========================================================================

    def _write_data_file(self):
        with open(self.output_file, 'w') as f:
            f.write("LAMMPS data file (Dreiding)\n\n")
            f.write(f"{len(self.atom_data)} atoms\n")
            f.write(f"{len(self.bond_data)} bonds\n")
            f.write(f"{len(self.angle_data)} angles\n")
            f.write(f"{len(self.dihedral_data)} dihedrals\n")
            f.write(f"{len(self.improper_data)} impropers\n\n")
            
            f.write(f"{len(self.atom_types_dict)} atom types\n")
            f.write(f"{len(self.bond_types)} bond types\n")
            f.write(f"{len(self.angle_types)} angle types\n")
            f.write(f"{len(self.dihedral_types)} dihedral types\n")
            f.write(f"{len(self.improper_types)} improper types\n\n")
            
            f.write("-1000.0 1000.0 xlo xhi\n-1000.0 1000.0 ylo yhi\n-1000.0 1000.0 zlo zhi\n\n")
            
            f.write("Masses\n\n")
            for name, tid in sorted(self.atom_types_dict.items(), key=lambda x: x[1]):
                m = 1.0
                if name in self.params['atom_types']:
                    m = self.params['atom_types'][name]['mass']
                f.write(f"{tid} {m:.4f} # {name}\n")
            f.write("\n")
            
            f.write("Atoms # full\n\n")
            for a in self.atom_data:
                # Use standardized keys 'type_idx' and 'label'
                f.write(f"{a['id']} 1 {a['type_idx']} {a['charge']:.12f} {a['x']} {a['y']} {a['z']} # {a['label']}\n")
            f.write("\n")
            
            if self.bond_data:
                f.write("Bonds\n\n")
                for b in self.bond_data: f.write(f"{b[0]} {b[1]} {b[2]} {b[3]}\n")
                f.write("\n")
            
            if self.angle_data:
                f.write("Angles\n\n")
                for a in self.angle_data: f.write(f"{a[0]} {a[1]} {a[2]} {a[3]} {a[4]}\n")
                f.write("\n")
            
            if self.dihedral_data:
                f.write("Dihedrals\n\n")
                for d in self.dihedral_data: f.write(f"{d[0]} {d[1]} {d[2]} {d[3]} {d[4]} {d[5]}\n")
                f.write("\n")
                
            if self.improper_data:
                f.write("Impropers\n\n")
                for i in self.improper_data: f.write(f"{i[0]} {i[1]} {i[2]} {i[3]} {i[4]} {i[5]}\n")
                f.write("\n")

    def _write_settings_file(self):
        with open(self.settings_file, 'w') as f:
            f.write("# Forcefield Settings (Auto-Generated)\n\n")
            
            f.write("# Pair Coeffs\n")
            id_to_name = {v: k for k, v in self.atom_types_dict.items()}
            for tid in sorted(id_to_name.keys()):
                name = id_to_name[tid]
                if name in self.params['vdw_params']:
                    p = self.params['vdw_params'][name]
                    eps = p['epsilon']
                    sig = p['radius'] / 1.122462
                else:
                    eps, sig = 0.001, 3.5
                f.write(f"pair_coeff {tid} {tid} {eps:.4f} {sig:.4f} # {name}\n")
            f.write("\n")
            
            if self.bond_data:
                f.write("# Bond Coeffs\n")
                for type_sig, tid in sorted(self.bond_types.items(), key=lambda x: x[1]):
                    # type_sig = (canon_t1, canon_t2, k_val, r0)
                    t1, t2, k, r0 = type_sig
                    f.write(f"bond_coeff {tid} {k:.4f} {r0:.4f} # {t1}-{t2}\n")
                f.write("\n")
            
            if self.angle_data:
                f.write("# Angle Coeffs\n")
                for type_sig, tid in sorted(self.angle_types.items(), key=lambda x: x[1]):
                    # type_sig = (t1, t_center, t2, k_val, theta0)
                    t1, t2, t3, k, t0 = type_sig
                    f.write(f"angle_coeff {tid} {k:.4f} {t0:.4f} # {t1}-{t2}-{t3}\n")
                f.write("\n")
            
            if self.dihedral_data:
                f.write("# Dihedral Coeffs\n")
                for type_sig, tid in sorted(self.dihedral_types.items(), key=lambda x: x[1]):
                    # type_sig = (canon_t, k_val, n, d)
                    types, k, n, d = type_sig
                    label = f"{types[0]}-{types[1]}-{types[2]}-{types[3]}"
                    f.write(f"dihedral_coeff {tid} {k:.6f} {int(d)} {int(n)} # {label}\n")
                f.write("\n")
                
            if self.improper_data:
                f.write("# Improper Coeffs\n")
                for type_sig, tid in sorted(self.improper_types.items(), key=lambda x: x[1]):
                    c_type, _, k, _ = type_sig
                    f.write(f"improper_coeff {tid} {k:.4f} -1 0 # {c_type}\n")
                f.write("\n")