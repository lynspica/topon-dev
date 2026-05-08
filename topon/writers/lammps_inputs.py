import os
import json

class LammpsInputGenerator:
    def __init__(self, output_dir, study_name, config=None, experimental=None):
        self.root_dir = os.path.join(output_dir, study_name)
        self.config = config or {}
        self.experimental = experimental or {}
        self.sim_dir = os.path.join(self.root_dir, "04_Simulation")
        self.chem_dir = os.path.join(self.root_dir, "02_Chemistry")
        self.conf_dir = os.path.join(self.root_dir, "03_Conformation")
        
        # Determine if in test mode
        self.test_mode = self.experimental.get('test_mode', False)
        
        if not os.path.exists(self.sim_dir):
            os.makedirs(self.sim_dir)
    
    def _get_cg_param(self, *keys, default=None):
        """Get CG parameter from experimental config."""
        d = self.experimental.get('cg', {})
        for k in keys:
            if isinstance(d, dict):
                d = d.get(k, default if k == keys[-1] else {})
            else:
                return default
        return d
    
    def _get_run_steps(self, model_type='cg'):
        """Get run steps from experimental config."""
        return self.experimental.get(model_type, {}).get('dynamics', {}).get('run_steps', 10000)

    def write_serial_soft_minimization(self, input_data="system_relaxed.data", groups_file="system.groups", settings_file="system.in.settings", model_type="atomistic"):
        """
        Stage 1: Serial Soft Minimization.
        Freezes nodes, resolves hard overlaps using soft potential.
        """
        script_path = os.path.join(self.sim_dir, "minimize_1_serial.in")
        
        data_path = os.path.relpath(os.path.join(self.conf_dir, input_data), self.sim_dir).replace("\\", "/")
        groups_path = os.path.relpath(os.path.join(self.chem_dir, groups_file), self.sim_dir).replace("\\", "/")
        settings_path = os.path.relpath(os.path.join(self.chem_dir, settings_file), self.sim_dir).replace("\\", "/")
        
        with open(script_path, 'w') as f:
            f.write(f"# LAMMPS Stage 1: Serial Soft Minimization ({model_type.upper()})\n\n")
            
            if model_type == 'cg':
                f.write("units           lj\n")
                f.write("atom_style      full\n")
                f.write("boundary        p p p\n")
                f.write("comm_modify     mode single cutoff 5.0\n")
                # User requested Harmonic start
                f.write("bond_style      harmonic\n")
                if self.config.get('include_angles', True):
                    f.write("angle_style     harmonic\n")
                # Pair style: repulsive (2^(1/6)) or attractive (2.5)
                pair_style = self.config.get('pair_style', 'attractive')
                pair_cutoff = 1.122462 if pair_style == 'repulsive' else 2.5
                f.write(f"pair_style      lj/cut {pair_cutoff}\n")
                f.write("special_bonds   lj 0.0 1.0 1.0\n\n")
            else:
                f.write("units           real\n")
                f.write("atom_style      full\n")
                f.write("boundary        p p p\n")
                f.write("bond_style      harmonic\n")
                f.write("angle_style     harmonic\n")
                f.write("dihedral_style  harmonic\n")
                f.write("improper_style  cvff\n")
                f.write("special_bonds   dreiding\n")
                f.write("pair_style      lj/cut/coul/long 12.0\n")
                f.write("kspace_style    pppm 1.0e-4\n\n")
            
            f.write(f"read_data       {data_path}\n")
            f.write(f"include         {settings_path}\n")
            f.write(f"include         {groups_path}\n\n")
            
            if model_type == 'cg':
                f.write("group           beads subtract all nodes\n\n")

            f.write("neighbor        2.0 bin\n")
            f.write("neigh_modify    every 1 delay 0 check yes\n\n")
            
            f.write("# --- Switch to Soft Potential ---\n")
            if model_type == 'atomistic': f.write("kspace_style    none\n")
            f.write("pair_style      soft 1.0\n")
            f.write("pair_coeff      * * 0.0\n")
            f.write("variable        prefactor equal ramp(0,30)\n\n")
            
            if model_type == 'cg':
                f.write("# --- Stage A: Freeze All (Beads & Nodes) ---\n")
                f.write("fix             freeze_beads beads setforce 0 0 0\n")
                f.write("fix             freeze_nodes nodes setforce 0 0 0\n")
                f.write("fix             soft_push all adapt 1 pair soft a * * v_prefactor\n")
                f.write("min_style       cg\n")
                f.write("minimize        1.0e-4 1.0e-6 1000 10000\n")
                f.write("unfix           soft_push\n")
                f.write("unfix           freeze_beads\n")
                f.write("write_data      min_stage_A.data\n\n")

                f.write("# --- Stage B: Relax Beads (Nodes Fixed) ---\n")
                f.write("fix             soft_push all adapt 1 pair soft a * * v_prefactor\n")
                f.write("minimize        1.0e-4 1.0e-6 1000 10000\n")
                f.write("unfix           soft_push\n")
                f.write("unfix           freeze_nodes\n")
                f.write("write_data      min_stage_B.data\n\n")

                f.write("# --- Stage C: Relax All ---\n")
                f.write("fix             soft_push all adapt 1 pair soft a * * v_prefactor\n")
                f.write("minimize        1.0e-4 1.0e-6 10000 100000\n")
                f.write("unfix           soft_push\n")
                f.write("write_data      system_after_soft.data\n")
                f.write("write_restart   1.restart\n")
            else:
                f.write("# --- Freeze All Groups ---\n")
                f.write("if \"$(is_defined(group,si_atoms))\" then \"fix freeze_si si_atoms setforce 0 0 0\"\n")
                f.write("if \"$(is_defined(group,c_atoms))\"  then \"fix freeze_c  c_atoms  setforce 0 0 0\"\n")
                f.write("if \"$(is_defined(group,o_atoms))\"  then \"fix freeze_o  o_atoms  setforce 0 0 0\"\n")
                f.write("if \"$(is_defined(group,h_atoms))\"  then \"fix freeze_h  h_atoms  setforce 0 0 0\"\n")
                f.write("if \"$(is_defined(group,nodes))\"    then \"fix freeze_nodes nodes setforce 0 0 0\"\n\n")
                
                f.write("fix soft_push all adapt 1 pair soft a * * v_prefactor\n")
                f.write("min_style cg\nminimize 1.0e-4 1.0e-6 1000 10000\nunfix soft_push\n\n")
                
                f.write("if \"$(is_defined(fix,freeze_h))\" then \"unfix freeze_h\"\n")
                f.write("if \"$(is_defined(fix,freeze_c))\" then \"unfix freeze_c\"\n")
                f.write("fix soft_push all adapt 1 pair soft a * * v_prefactor\n")
                f.write("minimize 1.0e-4 1.0e-6 1000 10000\nunfix soft_push\n\n")
                
                f.write("if \"$(is_defined(fix,freeze_o))\" then \"unfix freeze_o\"\n")
                f.write("fix soft_push all adapt 1 pair soft a * * v_prefactor\n")
                f.write("minimize 1.0e-4 1.0e-6 1000 10000\nunfix soft_push\n\n")
                
                f.write("if \"$(is_defined(fix,freeze_si))\" then \"unfix freeze_si\"\n")
                f.write("fix soft_push all adapt 1 pair soft a * * v_prefactor\n")
                f.write("minimize 1.0e-4 1.0e-6 1000 10000\nunfix soft_push\n\n")
                
                f.write("if \"$(is_defined(fix,freeze_nodes))\" then \"unfix freeze_nodes\"\n")
                f.write("fix soft_push all adapt 1 pair soft a * * v_prefactor\n")
                f.write("minimize 1.0e-4 1.0e-6 1000 10000\nunfix soft_push\n\n")
                
                f.write("reset_timestep 0\ntimestep 1.0\n")
                f.write("fix soft_push all adapt 1 pair soft a * * v_prefactor\n")
                f.write("fix nve_limit all nve/limit 0.1\nrun 1000\nunfix nve_limit\nunfix soft_push\n")
                f.write("minimize 1.0e-4 1.0e-6 1000 10000\n")
                f.write("write_data system_after_soft.data\n")
                f.write("write_restart 1.restart\n")

        return script_path

    def write_parallel_production(self, settings_file="system.in.settings", model_type="atomistic"):
        """
        Parent function that generates the complete parallel minimization pipeline:
        1. minimize_2_parallel.in (Stage 2: Ramp)
        2. minimize_3_parallel.in (Stage 3: Tight Min + Equilibration)
        """
        if model_type == 'cg':
            self._write_cg_minimization_equil(settings_file)
            print(f"Generated parallel minimization scripts (CG Stages 2 & 3).")
            return

        self._write_stage2_ramp(settings_file, model_type)
        self._write_stage3_equilibration(settings_file, model_type)
        print(f"Generated parallel minimization scripts (Stages 2 & 3).")

    def _write_stage2_ramp(self, settings_file, model_type):
        """
        Stage 2: Parallel Ramp.
        METHODOLOGY: Set 1 (Slow 200k Step Ramp + Extended Cutoffs)
        Inputs: system_after_soft.data
        Outputs: system_ramped.data
        """
        script_path = os.path.join(self.sim_dir, "minimize_2_parallel.in")
        settings_path = os.path.relpath(os.path.join(self.chem_dir, settings_file), self.sim_dir).replace("\\", "/")
        
        with open(script_path, 'w') as f:
            f.write("# LAMMPS Stage 2: Parallel Ramp\n")
            f.write("# METHODOLOGY: Set 1 (Slow 200k Step Ramp + Extended Cutoffs)\n\n")
            
            f.write("units           real\n")
            f.write("atom_style      full\n")
            f.write("boundary        p p p\n")
            f.write("bond_style      harmonic\n")
            f.write("angle_style     harmonic\n")
            f.write("dihedral_style  harmonic\n")
            f.write("improper_style  cvff\n")
            f.write("pair_style      soft 1.0\n\n")
            
            f.write("# --- 1. Load Soft State ---\n")
            f.write("read_data       system_after_soft.data\n\n")
            
            f.write("# --- 2. CRITICAL SAFETY (From Set 1) ---\n")
            f.write("# Prevents \"Bond atoms missing\" and \"Neighbor list overflow\"\n")
            f.write("neigh_modify    one 10000\n")
            f.write("comm_modify     mode single cutoff 12.0\n\n")
            
            f.write("# --- 3. Soft Pre-Minimization ---\n")
            f.write("pair_style      soft 1.0\n")
            f.write("pair_coeff      * * 1.0\n")
            f.write("min_style       cg\n")
            f.write("minimize        1.0e-4 1.0e-6 1000 10000\n\n")
            
            f.write("# --- 4. Switch to Real Potential ---\n")
            f.write("pair_style      lj/cut/coul/long 10.0 10.0\n")
            f.write("kspace_style    pppm 1.0e-4\n")
            f.write(f"include         {settings_path}\n\n")
            
            f.write("# Enforce Set 1 Special Bonds\n")
            f.write("special_bonds   lj/coul 0.0 0.0 1.0\n\n")
            
            f.write("# --- 5. The Ramp (Set 1 Logic) ---\n")
            f.write("# Linearly scale epsilon/charges from 0.001 to 1.0\n")
            f.write("variable        scale equal \"ramp(0.001, 1.0)\"\n")
            f.write("timestep        1.0\n\n")
            
            f.write("fix             1 all adapt 1 pair lj/cut/coul/long epsilon * * v_scale\n")
            f.write("fix             fxnve all nve/limit 0.1\n")
            f.write("thermo          1000\n\n")
            
            # RUNTIME from experimental config
            run_steps = self._get_run_steps('atomistic')
            f.write(f"# RUNTIME: {run_steps} steps\n")
            f.write(f"run             {run_steps}\n\n")
            
            f.write("unfix           fxnve\n")
            f.write("unfix           1\n")
            f.write("kspace_modify   compute yes\n\n")
            
            f.write("write_data      system_ramped.data\n")

    def _write_stage3_equilibration(self, settings_file, model_type):
        """
        Stage 3: Parallel Equilibration.
        METHODOLOGY: Set 1 (Tight Min -> NVT -> NPT)
        Inputs: system_ramped.data
        Outputs: system_equilibrated.data
        """
        script_path = os.path.join(self.sim_dir, "minimize_3_parallel.in")
        settings_path = os.path.relpath(os.path.join(self.chem_dir, settings_file), self.sim_dir).replace("\\", "/")
        
        with open(script_path, 'w') as f:
            f.write("# LAMMPS Stage 3: Parallel Equilibration\n")
            f.write("# METHODOLOGY: Set 1 (Tight Min -> NVT -> NPT)\n\n")
            
            f.write("units           real\n")
            f.write("atom_style      full\n")
            f.write("boundary        p p p\n")
            f.write("bond_style      harmonic\n")
            f.write("angle_style     harmonic\n")
            f.write("dihedral_style  harmonic\n")
            f.write("improper_style  cvff\n")
            f.write("pair_style      lj/cut/coul/long 10.0 10.0\n\n")
            
            f.write("# --- 1. Load Ramped State ---\n")
            f.write("read_data       system_ramped.data\n\n")
            
            f.write("# --- 2. Safety Settings ---\n")
            f.write("# Keep these even in stage 3 to prevent random crashes\n")
            f.write("neigh_modify    one 10000\n\n")
            
            f.write("# --- 3. Define Potential ---\n")
            f.write("pair_style      lj/cut/coul/long 10.0 10.0\n")
            f.write("kspace_style    pppm 1.0e-4\n")
            f.write(f"include         {settings_path}\n")
            f.write("special_bonds   lj/coul 0.0 0.0 1.0\n\n")
            
            f.write("# --- 4. Tight Minimization (Set 1 Logic) ---\n")
            f.write("# High precision 1e-8/1e-10 tolerances\n")
            f.write("min_style       cg\n")
            f.write("minimize        1.0e-8 1.0e-10 10000000 100000000\n\n")
            
            f.write("write_data      system_minimized_final.data\n\n")
            
            f.write("# --- 5. Equilibration Loop (Set 1 Logic) ---\n")
            f.write("reset_timestep  0\n")
            f.write("variable        temp equal 300\n")
            f.write("velocity        all create ${temp} 12345\n\n")
            
            f.write("# NVT (1000 steps)\n")
            f.write("fix             1 all nvt temp ${temp} ${temp} 100.0\n")
            f.write("run             1000\n")
            f.write("unfix           1\n")
            f.write("write_data      after_nvt_real.data\n\n")
            
            f.write("# NPT (1000 steps)\n")
            f.write("fix             1 all npt temp ${temp} ${temp} 100.0 iso 1.0 1.0 1000.0\n")
            f.write("run             1000\n")
            f.write("unfix           1\n\n")
            
            f.write("write_data      system_equilibrated.data\n")
            f.write("print \"All Minimization Stages Complete.\"\n")

    def write_equilibration_sequence(self, settings_file="system.in.settings", model_type="atomistic"):
        """
        Generates a 10-step Equilibration Sequence (Atomistic).
        1-4: Annealing (1000K -> 300K)
        5-6: 300K Equilibration (1M steps)
        7-8: 373K Equilibration (1M steps)
        9-10: 800K Equilibration (1M steps)
        """
        if model_type == 'cg':
            self._write_cg_minimization_equil(settings_file)
            print(f"Generated CG minimization scripts (Stages 2 & 3).")
            return

    def _write_cg_minimization_equil(self, settings_file):
        """
        Generates Stage 2 & 3 scripts for CG model using Reference Logic (Harmonic Ramp).
        Stage 2: Harmonic Ramp (minimize_2_cg.in) - Switch to Harmonic for stability
        Stage 3: Equilibration (minimize_3_cg.in) - Switch back to FENE
        """
        # --- Stage 2: Harmonic Ramp Minimization ---
        script_path = os.path.join(self.sim_dir, "minimize_2_parallel.in")
        
        with open(script_path, 'w') as f:
            f.write("# LAMMPS Stage 2: CG Harmonic Ramp (Reference Logic)\n\n")
            
            # Read Restart from Stage 1 (Preserves state)
            f.write("read_restart    1.restart\n\n")
            
            f.write("neighbor        2.0 bin\n")
            f.write("neigh_modify    every 1 delay 0 check yes\n")
            f.write("comm_modify     mode single cutoff 5.0\n\n")

            f.write("# SWITCH TO HARMONIC for Robust Minimization\n")
            f.write("bond_style      harmonic\n")
            f.write("bond_coeff      1 466.1 0.97\n")
            if self.config.get('include_angles', True):
                f.write("angle_style     harmonic\n")
                f.write("angle_coeff     1 466.1 180.0\n\n") # Generic stiff angle

            f.write("# Soft to Real Potential Ramp\n")
            f.write("pair_style      soft 1.0\n")
            f.write("pair_coeff      * * 1.0\n")
            f.write("min_style       cg\n")
            f.write("minimize        1e-4 1e-6 1000 10000\n\n")

            # Switch to Real LJ
            pair_style = self.config.get('pair_style', 'attractive')
            pair_cutoff = 1.122462 if pair_style == 'repulsive' else 2.5
            f.write(f"pair_style      lj/cut {pair_cutoff}\n")
            f.write(f"pair_coeff      * * 1.0 1.0 {pair_cutoff}\n")
            
            # Ramp parameters from config
            cg_ramp = self.experimental.get('cg', {}).get('ramp', {})
            scale_min = cg_ramp.get('epsilon_scale_start', 0.001)
            scale_max = cg_ramp.get('epsilon_scale_end', 1.0)
            nve_limit = cg_ramp.get('nve_limit', 0.1)
            ramp_steps = cg_ramp.get('ramp_steps', 20000)
            
            f.write(f"variable        scale equal \"ramp({scale_min}, {scale_max})\"\n")
            f.write(f"fix             1 all adapt 1 pair lj/cut epsilon * * v_scale\n")
            f.write(f"fix             fxnve all nve/limit {nve_limit}\n")
            f.write("thermo          1000\n")
            f.write(f"run             {ramp_steps}\n")
            f.write("unfix           fxnve\n")
            f.write("unfix           1\n\n")
            
            f.write("write_restart   2.restart\n")
            f.write("write_data      system_ramped.data\n")
            
        # --- Stage 3: FENE Equilibration ---
        script_path = os.path.join(self.sim_dir, "minimize_3_parallel.in")
        
        with open(script_path, 'w') as f:
            f.write("# LAMMPS Stage 3: CG Equilibration (FENE Restore)\n\n")
            
            f.write("read_restart    2.restart\n\n")
            
            f.write("neighbor        2.0 bin\n")
            f.write("neigh_modify    every 1 delay 0 check yes\n\n")

            f.write("# --- Phase A: Harmonic Pre-Minimization ---\n")
            f.write("bond_style      harmonic\n")
            f.write("bond_coeff      1 466.1 0.97\n")
            if self.config.get('include_angles', True):
                f.write("angle_style     harmonic\n")
                f.write("angle_coeff     1 466.1 180.0\n")
            pair_style = self.config.get('pair_style', 'attractive')
            pair_cutoff = 1.122462 if pair_style == 'repulsive' else 2.5
            f.write(f"pair_style      lj/cut {pair_cutoff}\n")
            f.write(f"pair_coeff      * * 1.0 1.0 {pair_cutoff}\n")
            f.write("min_style       cg\n")
            f.write("minimize        1.0e-4 1.0e-6 1000 10000\n\n")

            # --- Phase B: Switch to FENE ---
            f.write("# --- Phase B: Switch to FENE ---\n")
            
            # Check config for angle removal (Default: True)
            # Only remove if they were included to begin with
            include_angles = self.config.get('include_angles', True)
            remove_angles = self.config.get('remove_cg_angles', True)
            
            if include_angles and remove_angles:
                f.write("delete_bonds    all angle 1-1 remove\n")
            
            f.write("bond_style      fene\n")
            f.write("special_bonds   fene\n")
            f.write("bond_coeff      1 30.0 1.5 1.0 1.0\n")
            f.write("minimize        1.0e-4 1.0e-6 1000 10000\n\n")
            
            f.write("# --- Phase C: Dynamics ---\n")
            # Get parameters from experimental config
            cg_dyn = self.experimental.get('cg', {}).get('dynamics', {})
            timestep = cg_dyn.get('timestep', 0.005)
            temp = cg_dyn.get('temperature', 1.0)
            tdamp = cg_dyn.get('tdamp', 1.0)
            thermo = cg_dyn.get('thermo_freq', 1000)
            run_steps = self._get_run_steps('cg')
            
            f.write("reset_timestep  0\n")
            f.write(f"timestep        {timestep}\n")
            f.write(f"velocity        all create {temp} 12345\n\n")
            f.write(f"fix             1 all nvt temp {temp} {temp} {tdamp}\n")
            f.write(f"thermo          {thermo}\n")
            f.write(f"run             {run_steps}\n")
            f.write("unfix           1\n\n")
            f.write("write_data      system_equilibrated.data\n")


        print(f"Generated parallel minimization scripts (CG Stages 2 & 3).")