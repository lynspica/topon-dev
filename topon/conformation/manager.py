import os
import re
import math
import numpy as np
from collections import defaultdict

class ConformationManager:
    def __init__(self, output_dir, study_name):
        self.root_dir = os.path.join(output_dir, study_name)
        self.chem_dir = os.path.join(self.root_dir, "02_Chemistry")
        self.conf_dir = os.path.join(self.root_dir, "03_Conformation")
        
        if not os.path.exists(self.conf_dir):
            os.makedirs(self.conf_dir)

        self.xyz_indices = [4, 5, 6] 
            
    def _detect_xyz_indices(self, sample_line):
        parts = sample_line.split()
        n = len(parts)
        if n >= 7: return [4, 5, 6] 
        elif n == 6: return [3, 4, 5]
        elif n == 5: return [2, 3, 4]
        else: return [4, 5, 6]

    def apply_displacements(self, base_filename="system.data"):
        """
        Applies displacements, calculates the IDEAL lattice box, 
        and wraps all coordinates into it.
        """
        base_path = os.path.join(self.chem_dir, base_filename)
        print(f"Applying displacements to {base_filename}...")
        
        with open(base_path, 'r') as f: lines = f.readlines()
            
        atom_lines, header_lines, footer_lines = [], [], []
        section = 'header'
        for line in lines:
            s = line.strip()
            if section == 'header':
                if s.startswith('Atoms'): section = 'atoms_start'; header_lines.append(line)
                else: header_lines.append(line)
            elif section == 'atoms_start':
                if s == '': header_lines.append(line)
                else: section = 'atoms_body'; atom_lines.append(line)
            elif section == 'atoms_body':
                if s == '' or s.startswith(('Bonds', 'Velocities', 'Masses')):
                    section = 'footer'; footer_lines.append(line)
                else: atom_lines.append(line)
            elif section == 'footer': footer_lines.append(line)

        if atom_lines:
            for line in atom_lines:
                if line.strip() and not line.strip().startswith('#'):
                    self.xyz_indices = self._detect_xyz_indices(line)
                    break

        atoms_data = {}
        for line in atom_lines:
            parts = line.split()
            if parts: atoms_data[int(parts[0])] = parts

        import glob
        displace_files = glob.glob(os.path.join(self.chem_dir, "*.displace"))
        updates = {} 
        atom_roles = {} 
        
        # --- BOX DIMENSION DISCOVERY ---
        # We need to find the max lattice coordinate to determine the box size.
        # We look specifically at '_nodes.displace' for this.
        lattice_max = {'x': 0.0, 'y': 0.0, 'z': 0.0}
        scales = {'x': 1.0, 'y': 1.0, 'z': 1.0}
        
        for d_file in displace_files:
            fname = os.path.basename(d_file)
            is_node_file = "_nodes" in fname
            
            if "_nodes" in fname: priority = 0      
            elif "_backbone" in fname: priority = 1 
            elif "_beads" in fname: priority = 1    
            else: priority = 2                      
            
            with open(d_file, 'r') as f:
                for line in f:
                    if line.startswith('variable scale_'):
                        parts = line.split()
                        scales[parts[1].split('_')[1]] = float(parts[3])
                    elif line.startswith('variable d'):
                        match = re.match(r'variable d([xyz])_(\d+) equal v_scale_[xyz]\*(-?\d+\.?\d*)', line)
                        if match:
                            dim, atom_id, val_str = match.groups()
                            val = float(val_str)
                            atom_id = int(atom_id)
                            
                            if atom_id not in updates: updates[atom_id] = {}
                            updates[atom_id][dim] = val * scales[dim]
                            
                            curr_p = atom_roles.get(atom_id, 99)
                            atom_roles[atom_id] = min(curr_p, priority)
                            
                            # Track max lattice coordinate from Nodes file
                            if is_node_file:
                                lattice_max[dim] = max(lattice_max[dim], val)

        # --- CALCULATE IDEAL BOX ---
        # Box Size = (MaxNodeCoord + 1.0) * Scale
        # This recovers the 'dims' from the Builder
        box_x = (lattice_max['x'] + 1.0) * scales['x']
        box_y = (lattice_max['y'] + 1.0) * scales['y']
        box_z = (lattice_max['z'] + 1.0) * scales['z']
        
        print(f"  - Detected Lattice Extents (Graph Units): [{lattice_max['x']:.1f}, {lattice_max['y']:.1f}, {lattice_max['z']:.1f}]")
        print(f"  - Set Simulation Box: [{0.0:.2f}, {box_x:.2f}] x [{0.0:.2f}, {box_y:.2f}] x [{0.0:.2f}, {box_z:.2f}]")

        ix, iy, iz = self.xyz_indices
        count = 0
        
        # Apply Updates with WRAPPING
        for atom_id, coords in updates.items():
            if atom_id in atoms_data:
                if 'x' in coords:
                    val = coords['x'] % box_x # WRAP X
                    atoms_data[atom_id][ix] = f"{val:.6f}"
                if 'y' in coords: 
                    val = coords['y'] % box_y # WRAP Y
                    atoms_data[atom_id][iy] = f"{val:.6f}"
                if 'z' in coords: 
                    val = coords['z'] % box_z # WRAP Z
                    atoms_data[atom_id][iz] = f"{val:.6f}"
                count += 1
                
        print(f"  - Updated and wrapped {count} atoms.")

        # Rewrite Header with Ideal Box
        new_header = []
        for line in header_lines:
            if 'xlo xhi' in line: new_header.append(f"{0.0:.6f} {box_x:.6f} xlo xhi\n")
            elif 'ylo yhi' in line: new_header.append(f"{0.0:.6f} {box_y:.6f} ylo yhi\n")
            elif 'zlo zhi' in line: new_header.append(f"{0.0:.6f} {box_z:.6f} zlo zhi\n")
            else: new_header.append(line)

        output_path = os.path.join(self.conf_dir, "system_conformed.data")
        with open(output_path, 'w') as f:
            f.writelines(new_header)
            for atom_id in sorted(atoms_data.keys()):
                f.write(" ".join(atoms_data[atom_id]) + "\n")
            f.writelines(footer_lines)
            
        return output_path, atom_roles

    def resolve_overlaps(self, input_file, atom_roles, cutoff=0.85, max_iters=50):
        """
        Iterative Smart Relax.
        Reads box from input file and enforces wrapping during relaxation.
        """
        print(f"Checking for hard overlaps (cutoff < {cutoff} Å)...")
        
        with open(input_file, 'r') as f: lines = f.readlines()
        
        header, atoms, footer = [], [], []
        stage = 0 
        
        # Parse and Extract Box Dims from Header
        box_dims = {'x': 100.0, 'y': 100.0, 'z': 100.0} # Defaults
        
        for line in lines:
            if stage == 0:
                header.append(line)
                if 'xlo xhi' in line:
                    parts = line.split()
                    box_dims['x'] = float(parts[1]) - float(parts[0])
                elif 'ylo yhi' in line:
                    parts = line.split()
                    box_dims['y'] = float(parts[1]) - float(parts[0])
                elif 'zlo zhi' in line:
                    parts = line.split()
                    box_dims['z'] = float(parts[1]) - float(parts[0])
                    
                if line.strip().startswith('Atoms'): stage = 1
            elif stage == 1:
                if line.strip() == '': header.append(line)
                else: atoms.append(line); stage = 2
            elif stage == 2:
                if line.strip() == '' or line.strip()[0].isalpha(): footer.append(line); stage = 3
                else: atoms.append(line)
            elif stage == 3: footer.append(line)

        if atoms: self.xyz_indices = self._detect_xyz_indices(atoms[0])
        ix, iy, iz = self.xyz_indices

        coords = {}
        atom_lines_map = {}
        for i, line in enumerate(atoms):
            parts = line.split()
            atom_id = int(parts[0])
            pos = np.array([float(parts[ix]), float(parts[iy]), float(parts[iz])])
            coords[atom_id] = pos
            atom_lines_map[atom_id] = (i, parts) 

        # Relaxation Loop
        for iteration in range(max_iters):
            grid = defaultdict(list)
            for aid, pos in coords.items():
                cell = tuple((pos // cutoff).astype(int))
                grid[cell].append(aid)

            moved_count = 0
            for cell, cell_atoms in grid.items():
                cx, cy, cz = cell
                neighbors = []
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        for dz in [-1, 0, 1]:
                            neighbors.extend(grid.get((cx+dx, cy+dy, cz+dz), []))
                
                for id_a in cell_atoms:
                    pos_a = coords[id_a]
                    role_a = atom_roles.get(id_a, 2) 
                    
                    for id_b in neighbors:
                        if id_a >= id_b: continue
                        pos_b = coords[id_b]
                        
                        # Minimum Image Distance Check
                        diff = pos_a - pos_b
                        # Apply MIC for distance check (assuming orthogonal)
                        diff[0] -= box_dims['x'] * round(diff[0] / box_dims['x'])
                        diff[1] -= box_dims['y'] * round(diff[1] / box_dims['y'])
                        diff[2] -= box_dims['z'] * round(diff[2] / box_dims['z'])
                        
                        dist_sq = np.sum(diff**2)
                        
                        if dist_sq < (cutoff * cutoff):
                            dist = math.sqrt(dist_sq)
                            role_b = atom_roles.get(id_b, 2)
                            
                            mover_id = None
                            static_pos = None
                            # Determine mover...
                            if role_a == role_b: mover_id = id_b; static_pos = pos_a
                            elif role_a > role_b: mover_id = id_a; static_pos = pos_b
                            else: mover_id = id_b; static_pos = pos_a
                                
                            if mover_id is not None:
                                current_pos = coords[mover_id]
                                # Push Vector (normalized diff)
                                if dist < 1e-6: vec = np.random.rand(3) - 0.5
                                else: vec = diff / dist 
                                
                                # If we are moving 'a', vec points a->b (bad), we want a away from b.
                                # diff was a - b. So vec points b -> a. Correct.
                                # If we are moving 'b', we want b away from a. -vec.
                                
                                direction = 1.0 if mover_id == id_a else -1.0
                                push_dist = cutoff - dist + 0.05
                                new_pos = current_pos + (vec * direction * push_dist)
                                
                                # WRAP BACK INTO BOX
                                new_pos[0] %= box_dims['x']
                                new_pos[1] %= box_dims['y']
                                new_pos[2] %= box_dims['z']
                                
                                coords[mover_id] = new_pos
                                moved_count += 1

            if moved_count == 0:
                print(f"  - Converged in {iteration} iterations. System is clean.")
                break
            else:
                print(f"  - Iteration {iteration+1}: Resolved {moved_count} overlaps.")

        if moved_count > 0:
            print(f"  - Warning: Max iterations ({max_iters}) reached with overlaps remaining.")

        for aid, pos in coords.items():
            idx, parts = atom_lines_map[aid]
            parts[ix] = f"{pos[0]:.6f}"
            parts[iy] = f"{pos[1]:.6f}"
            parts[iz] = f"{pos[2]:.6f}"
            atoms[idx] = " ".join(parts) + "\n"

        # Write file (Header preserved from input, which had the correct box)
        final_path = os.path.join(self.conf_dir, "system_relaxed.data")
        with open(final_path, 'w') as f:
            f.writelines(header + atoms + footer)
            
        return final_path

    def apply_noise(self, input_file, magnitude=0.0001, output_name="system_relaxed.data"):
        """
        Applies a negligible random perturbation to all atoms to break symmetry.
        Used when overlap resolution is skipped or minimal.
        """
        print(f"Applying random noise (magnitude +/- {magnitude})...")
        
        with open(input_file, 'r') as f: lines = f.readlines()
        
        header, atoms, footer = [], [], []
        stage = 0 
        
        for line in lines:
            if stage == 0:
                header.append(line)
                if line.strip().startswith('Atoms'): stage = 1
            elif stage == 1:
                if line.strip() == '': header.append(line)
                else: atoms.append(line); stage = 2
            elif stage == 2:
                if line.strip() == '' or line.strip()[0].isalpha(): footer.append(line); stage = 3
                else: atoms.append(line)
            elif stage == 3: footer.append(line)

        if atoms: self.xyz_indices = self._detect_xyz_indices(atoms[0])
        ix, iy, iz = self.xyz_indices

        new_atoms = []
        for line in atoms:
            parts = line.split()
            x = float(parts[ix]) + (np.random.rand() - 0.5) * 2 * magnitude
            y = float(parts[iy]) + (np.random.rand() - 0.5) * 2 * magnitude
            z = float(parts[iz]) + (np.random.rand() - 0.5) * 2 * magnitude
            
            parts[ix] = f"{x:.6f}"
            parts[iy] = f"{y:.6f}"
            parts[iz] = f"{z:.6f}"
            new_atoms.append(" ".join(parts) + "\n")

        final_path = os.path.join(self.conf_dir, output_name)
        with open(final_path, 'w') as f:
            f.writelines(header + new_atoms + footer)
            
        print(f"  - Perturbed {len(new_atoms)} atoms.")
        return final_path