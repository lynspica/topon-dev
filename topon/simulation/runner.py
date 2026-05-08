import os
import subprocess
import time
from pathlib import Path

class SimulationRunner:
    """
    Handles execution of LAMMPS simulation scripts.
    Supports local serial execution and MPI execution.
    """
    def __init__(self, sim_dir, executable="lmp", n_procs=1, use_mpi=False):
        self.sim_dir = Path(sim_dir)
        self.executable = executable
        self.n_procs = n_procs
        self.use_mpi = use_mpi
        
    def run_sequence(self, scripts, log_prefix="log"):
        """
        Runs a sequence of LAMMPS input scripts in order.
        Stops if any script fails.
        """
        if not self.sim_dir.exists():
            raise FileNotFoundError(f"Simulation directory not found: {self.sim_dir}")
            
        print(f"--- Starting Simulation Sequence in {self.sim_dir.name} ---")
        
        for i, script in enumerate(scripts):
            script_path = self.sim_dir / script
            if not script_path.exists():
                print(f"Error: Script not found: {script}")
                return False
                
            log_file = f"{log_prefix}.{script}.txt"
            print(f"Running: {script} -> {log_file} ... ", end="", flush=True)
            
            cmd = self._build_command(script, log_file)
            
            start_time = time.time()
            try:
                # Execute in the simulation directory
                result = subprocess.run(
                    cmd, 
                    cwd=str(self.sim_dir),
                    capture_output=True,
                    text=True
                )
                
                duration = time.time() - start_time
                
                if result.returncode == 0:
                    print(f"Success ({duration:.1f}s)")
                else:
                    print(f"FAILED (Exit Code: {result.returncode})")
                    print(f"Helper: Check {self.sim_dir / log_file} for details.")
                    print("STDERR Snippet:")
                    print('\n'.join(result.stderr.splitlines()[-5:]))
                    return False
                    
            except FileNotFoundError:
                print(f"FAILED (Executable '{self.executable}' not found)")
                return False
            except Exception as e:
                print(f"FAILED (Error: {e})")
                return False
                
        print("--- Simulation Sequence Complete ---")
        return True

    def _build_command(self, script, log_file):
        """Constructs the command line arguments."""
        cmd = []
        
        if self.use_mpi:
            cmd.extend(["mpirun", "-np", str(self.n_procs)])
            cmd.append(self.executable)
        else:
            cmd.append(self.executable)
            
        cmd.extend(["-in", script])
        cmd.extend(["-log", log_file])
        
        return cmd
