"""
C generator wrapper for Topon.

Handles invoking the C-based topology generator and SLURM script generation.
"""

import subprocess
from pathlib import Path
from typing import Optional, Union

from topon.config.schema import GeneratorConfig


def run_generator(
    config: GeneratorConfig,
    output_dir: Union[str, Path],
    exe_path: Optional[Union[str, Path]] = None,
) -> tuple[Path, Path]:
    """
    Run the C topology generator.
    
    Args:
        config: Generator configuration.
        output_dir: Directory for output files.
        exe_path: Path to generator executable (overrides config).
        
    Returns:
        Tuple of (nodes_file_path, edges_file_path).
        
    Raises:
        FileNotFoundError: If generator executable not found.
        RuntimeError: If generator fails.
    """
    exe = Path(exe_path or config.exe_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if not exe.exists():
        raise FileNotFoundError(f"Generator executable not found: {exe}")
    
    # Build command
    cmd = [
        str(exe),
        config.lattice_size,
        config.periodicity,
        str(config.max_functionality),
        str(config.max_trials),
        str(config.max_saves),
        str(config.degree_distribution),
        "0",  # extensive_logging
        config.lattice_type,
    ]
    
    print(f"Running generator: {' '.join(cmd)}")
    
    # Run generator
    result = subprocess.run(
        cmd,
        cwd=output_dir,
        capture_output=True,
        text=True,
    )
    
    if result.returncode != 0:
        print(f"Generator stdout: {result.stdout}")
        print(f"Generator stderr: {result.stderr}")
        raise RuntimeError(f"Generator failed with return code {result.returncode}")
    
    # Find output files
    nodes_files = list(output_dir.glob("output/*.nodes"))
    edges_files = list(output_dir.glob("output/*.edges"))
    
    if not nodes_files or not edges_files:
        raise RuntimeError("Generator did not produce output files")
    
    return nodes_files[0], edges_files[0]


def generate_slurm_script(
    config: GeneratorConfig,
    output_path: Union[str, Path],
    slurm_config: Optional[dict] = None,
) -> str:
    """
    Generate a SLURM batch script for running the generator on HPC.
    
    Args:
        config: Generator configuration.
        output_path: Path to write the SLURM script.
        slurm_config: SLURM-specific configuration (account, partition, etc.).
        
    Returns:
        Path to the generated script.
    """
    output_path = Path(output_path)
    
    # Default SLURM config
    slurm = slurm_config or {}
    account = slurm.get("account", "default_account")
    partition = slurm.get("partition", "short")
    nodes = slurm.get("nodes", 1)
    tasks = slurm.get("tasks_per_node", 1)
    time = slurm.get("time", "03:59:59")
    exe_path = slurm.get("generator_exe_path", "./generator.exe")
    module_loads = slurm.get("module_loads", [])
    
    # Build job name from lattice config
    job_name = f"gen-{config.lattice_size}-{config.lattice_type}"
    
    # Build command
    cmd = (
        f'"{exe_path}" {config.lattice_size} {config.periodicity} '
        f'{config.max_functionality} {config.max_trials} {config.max_saves} '
        f'"{config.degree_distribution}" 0 {config.lattice_type}'
    )
    
    # Build script
    script_lines = [
        "#!/bin/bash",
        f"#SBATCH -A {account}",
        f"#SBATCH -p {partition}",
        f"#SBATCH -N {nodes}",
        f"#SBATCH --ntasks-per-node={tasks}",
        f"#SBATCH -t {time}",
        f"#SBATCH --job-name={job_name}",
        "#SBATCH --export=ALL",
        "",
        "module purge all",
    ]
    
    for module in module_loads:
        script_lines.append(f"module load {module}")
    
    script_lines.extend([
        "",
        f'echo "Running generator with config: {config.lattice_size} {config.lattice_type}"',
        "",
        cmd,
        "",
        'echo "Generator complete"',
    ])
    
    script_content = "\n".join(script_lines)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(script_content)
    
    print(f"Generated SLURM script: {output_path}")
    
    return str(output_path)
