"""
LAMMPS Data File Writer for Topon.

Fixed version with:
- Correct atom_style full (not molecular)
- DREIDING force field integration for atomistic
- Proper box coordinate handling (placeholder first, then displace)
"""

from typing import Optional
from pathlib import Path
import numpy as np
import shutil


def write_lammps_data(
    molecule,
    coords: dict,
    box_dims: tuple,
    output_path: str,
    model_type: str = "coarse_grained",
    atom_type_map: Optional[dict] = None,
    dreiding_param_file: Optional[str] = None,
    comment: str = "Topon generated structure"
) -> None:
    """
    Write LAMMPS data file.
    
    Args:
        molecule: RDKit RWMol object.
        coords: Dict mapping atom_idx to (x, y, z) NORMALIZED coordinates [0,1).
        box_dims: Box dimensions (lx, ly, lz).
        output_path: Output file path.
        model_type: 'atomistic' or 'coarse_grained'.
        atom_type_map: Optional element -> type mapping.
        dreiding_param_file: Path to DREIDING parameter file (for atomistic).
        comment: File header comment.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if model_type == "atomistic" and dreiding_param_file:
        # Use DREIDING parameterization
        _write_atomistic_data(molecule, output_path, dreiding_param_file, comment)
    else:
        # Use simple CG or generic format
        _write_cg_data(molecule, coords, box_dims, output_path, comment)


def _write_cg_data(
    molecule,
    coords: dict,
    box_dims: tuple,
    output_path: Path,
    comment: str
) -> None:
    """Write coarse-grained LAMMPS data file."""
    
    num_atoms = molecule.GetNumAtoms()
    num_bonds = molecule.GetNumBonds()
    lx, ly, lz = box_dims
    
    # Get unique atom types (for CG, typically just 1 or few types)
    atom_types = {}
    type_counter = 1
    atom_type_list = []
    
    for i in range(num_atoms):
        atom = molecule.GetAtomWithIdx(i)
        bead_type = atom.GetProp("bead_type") if atom.HasProp("bead_type") else "B"
        
        if bead_type not in atom_types:
            atom_types[bead_type] = type_counter
            type_counter += 1
        
        atom_type_list.append(atom_types[bead_type])
    
    with open(output_path, "w") as f:
        f.write(f"# {comment}\n\n")
        f.write(f"{num_atoms} atoms\n")
        f.write(f"{num_bonds} bonds\n")
        f.write(f"0 angles\n")
        f.write(f"0 dihedrals\n")
        f.write(f"0 impropers\n\n")
        
        f.write(f"{len(atom_types)} atom types\n")
        f.write(f"1 bond types\n\n")
        
        # Box dimensions - use the actual scaled box
        f.write(f"0.0 {lx:.6f} xlo xhi\n")
        f.write(f"0.0 {ly:.6f} ylo yhi\n")
        f.write(f"0.0 {lz:.6f} zlo zhi\n\n")
        
        # Masses
        f.write("Masses\n\n")
        for sym, type_id in sorted(atom_types.items(), key=lambda x: x[1]):
            mass = 1.0  # CG bead mass
            f.write(f"{type_id} {mass:.4f} # {sym}\n")
        f.write("\n")
        
        # Pair Coeffs (for CG Kremer-Grest style)
        f.write("Pair Coeffs # lj/cut\n\n")
        for sym, type_id in sorted(atom_types.items(), key=lambda x: x[1]):
            f.write(f"{type_id} 1.0 1.0 # epsilon sigma\n")
        f.write("\n")
        
        # Bond Coeffs (FENE for Kremer-Grest)
        f.write("Bond Coeffs # fene\n\n")
        f.write("1 30.0 1.5 1.0 1.0 # K R0 epsilon sigma\n\n")
        
        # Atoms section - use atom_style full
        f.write("Atoms # full\n\n")
        for i in range(num_atoms):
            atom_type = atom_type_list[i]
            mol_id = 1
            charge = 0.0
            
            if i in coords:
                # Coords are normalized [0,1), scale to box
                nx, ny, nz = coords[i]
                x, y, z = nx * lx, ny * ly, nz * lz
            else:
                x, y, z = 0.0, 0.0, 0.0
            
            # atom_style full: atom-ID molecule-ID atom-type q x y z
            f.write(f"{i+1} {mol_id} {atom_type} {charge:.6f} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write("\n")
        
        # Bonds section
        f.write("Bonds\n\n")
        bond_id = 1
        for bond in molecule.GetBonds():
            atom1 = bond.GetBeginAtomIdx() + 1
            atom2 = bond.GetEndAtomIdx() + 1
            f.write(f"{bond_id} 1 {atom1} {atom2}\n")
            bond_id += 1
    
    print(f"  Wrote CG LAMMPS data file: {output_path}")


def _write_atomistic_data(
    molecule,
    output_path: Path,
    dreiding_param_file: str,
    comment: str
) -> None:
    """
    Write atomistic LAMMPS data file with DREIDING parameterization.
    
    This writes atoms at origin with placeholder box, 
    to be followed by displacement files for positioning.
    """
    # Import existing make_dreiding module if available
    try:
        from topon.forcefield import dreiding
        dreiding.create_lammps_data_file(molecule, str(output_path), dreiding_param_file)
        print(f"  Wrote atomistic LAMMPS data file: {output_path}")
    except ImportError:
        # Fallback: write simple format without full parameterization
        _write_simple_atomistic_data(molecule, output_path, comment)


def _write_simple_atomistic_data(
    molecule,
    output_path: Path,
    comment: str
) -> None:
    """Fallback simple atomistic writer without DREIDING."""
    
    num_atoms = molecule.GetNumAtoms()
    num_bonds = molecule.GetNumBonds()
    
    # Get unique atom types by element
    atom_types = {}
    type_counter = 1
    atom_type_list = []
    
    for i in range(num_atoms):
        atom = molecule.GetAtomWithIdx(i)
        symbol = atom.GetSymbol()
        
        if symbol not in atom_types:
            atom_types[symbol] = type_counter
            type_counter += 1
        
        atom_type_list.append(atom_types[symbol])
    
    with open(output_path, "w") as f:
        f.write(f"# {comment}\n")
        f.write("# WARNING: Simplified atomistic format - no DREIDING parameterization\n\n")
        
        f.write(f"{num_atoms} atoms\n")
        f.write(f"{num_bonds} bonds\n")
        f.write(f"0 angles\n")
        f.write(f"0 dihedrals\n")
        f.write(f"0 impropers\n\n")
        
        f.write(f"{len(atom_types)} atom types\n")
        f.write(f"1 bond types\n\n")
        
        # Placeholder box (large) - will be scaled after displacements
        f.write("-1000.0 1000.0 xlo xhi\n")
        f.write("-1000.0 1000.0 ylo yhi\n")
        f.write("-1000.0 1000.0 zlo zhi\n\n")
        
        # Masses
        f.write("Masses\n\n")
        masses = {"H": 1.008, "C": 12.011, "O": 15.999, "Si": 28.086, "F": 18.998, "N": 14.007}
        for sym, type_id in sorted(atom_types.items(), key=lambda x: x[1]):
            mass = masses.get(sym, 1.0)
            f.write(f"{type_id} {mass:.4f} # {sym}\n")
        f.write("\n")
        
        # Atoms section - atom_style full with atoms at origin
        f.write("Atoms # full\n\n")
        for i in range(num_atoms):
            atom_type = atom_type_list[i]
            mol_id = 1
            charge = 0.0
            x, y, z = 0.0, 0.0, 0.0  # At origin, will be displaced later
            
            f.write(f"{i+1} {mol_id} {atom_type} {charge:.6f} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write("\n")
        
        # Bonds section
        f.write("Bonds\n\n")
        bond_id = 1
        for bond in molecule.GetBonds():
            atom1 = bond.GetBeginAtomIdx() + 1
            atom2 = bond.GetEndAtomIdx() + 1
            f.write(f"{bond_id} 1 {atom1} {atom2}\n")
            bond_id += 1
    
    print(f"  Wrote simple atomistic LAMMPS data file: {output_path}")
    print(f"  NOTE: This file needs DREIDING parameters added manually")


def write_displacement_file(
    atom_coords: dict,
    scale: tuple,
    output_path: str,
    group_name: str = "atoms"
) -> None:
    """
    Write LAMMPS displacement file.
    
    Args:
        atom_coords: Dict mapping atom_idx (0-based) to (x, y, z) NORMALIZED coords.
        scale: Scaling factors (sx, sy, sz) for final box size.
        output_path: Output file path.
        group_name: Name for displacement group.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    sx, sy, sz = scale
    
    with open(output_path, "w") as f:
        f.write(f"# LAMMPS displacement file for {group_name}\n\n")
        f.write("# Define scaling factors\n")
        f.write(f"variable scale_x equal {sx:.6f}\n")
        f.write(f"variable scale_y equal {sy:.6f}\n")
        f.write(f"variable scale_z equal {sz:.6f}\n\n")
        
        f.write(f"# Displace {group_name}\n")
        for atom_idx in sorted(atom_coords.keys()):
            x, y, z = atom_coords[atom_idx]
            lid = atom_idx + 1  # LAMMPS 1-indexed
            
            f.write(f"variable dx_{lid} equal v_scale_x*{x:.6f}\n")
            f.write(f"variable dy_{lid} equal v_scale_y*{y:.6f}\n")
            f.write(f"variable dz_{lid} equal v_scale_z*{z:.6f}\n")
            f.write(f"group current_atom id {lid}\n")
            f.write(f"displace_atoms current_atom move v_dx_{lid} v_dy_{lid} v_dz_{lid}\n")
            f.write("group current_atom delete\n\n")
    
    print(f"  Wrote displacement file: {output_path}")


def write_group_definitions(
    molecule,
    node_ids: list,
    output_path: str,
    model_type: str = "coarse_grained",
    scale_factors: tuple = (1.0, 1.0, 1.0),
    periodicity: tuple = (True, True, True),
    atom_type_map: Optional[dict] = None
) -> None:
    """
    Write LAMMPS group definitions and box scaling commands.
    
    Args:
        molecule: RDKit molecule.
        node_ids: List of atom indices that are nodes (0-based).
        output_path: Output file path.
        model_type: 'atomistic' or 'coarse_grained'.
        scale_factors: Box scaling (sx, sy, sz).
        periodicity: Periodic boundaries (px, py, pz).
        atom_type_map: Element -> list of type IDs mapping.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    scale_x, scale_y, scale_z = scale_factors
    
    # Adjust scale for non-periodic dimensions
    if not periodicity[0]: scale_x += 1
    if not periodicity[1]: scale_y += 1
    if not periodicity[2]: scale_z += 1
    
    with open(output_path, "w") as f:
        f.write("# LAMMPS group definitions and box settings\n\n")
        
        # Node group
        f.write("# Group definition for nodes (crosslinkers)\n")
        if node_ids:
            node_str = " ".join(str(i + 1) for i in sorted(node_ids))
            f.write(f"group nodes id {node_str}\n\n")
        
        # Element groups for atomistic
        if model_type == "atomistic":
            f.write("# Group definitions for ATOMISTIC model based on element types\n")
            if atom_type_map:
                for symbol, type_list in sorted(atom_type_map.items()):
                    types_str = ' '.join(map(str, sorted(type_list)))
                    f.write(f"group {symbol.lower()}_atoms type {types_str}\n")
            else:
                # Fallback: group by element from molecule
                by_element = {}
                for i in range(molecule.GetNumAtoms()):
                    sym = molecule.GetAtomWithIdx(i).GetSymbol()
                    if sym not in by_element:
                        by_element[sym] = []
                    by_element[sym].append(i + 1)
                
                for sym in sorted(by_element.keys()):
                    ids = " ".join(str(i) for i in by_element[sym])
                    f.write(f"group {sym.lower()}_atoms id {ids}\n")
        else:
            f.write("# Group definitions for COARSE-GRAINED model\n")
            f.write("group beads type 1\n")
        
        # Box scaling commands
        f.write("\n# Change boundary for scaling, then set to periodic\n")
        f.write("change_box all boundary s s s\n")
        f.write(f"change_box all x scale {scale_x} y scale {scale_y} z scale {scale_z}\n")
        f.write("change_box all boundary p p p\n")
    
    print(f"  Wrote group definitions: {output_path}")


def get_mass(symbol: str) -> float:
    """Get atomic mass for element symbol."""
    masses = {
        "H": 1.008,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "F": 18.998,
        "Si": 28.086,
        "S": 32.065,
        "Cl": 35.453,
        "B": 1.0,  # CG bead
        "J": 1.0,  # Junction bead
    }
    return masses.get(symbol, 1.0)
