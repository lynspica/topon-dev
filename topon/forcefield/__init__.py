"""
Topon Force Field Module

Force field parameterization for LAMMPS.
"""

from pathlib import Path

# Path to bundled parameter files
FORCEFIELD_DIR = Path(__file__).parent
DREIDING_PARAM_FILE = FORCEFIELD_DIR / "DreidingX6parameters.txt"

def get_dreiding_params_path() -> Path:
    """Get path to bundled DREIDING parameters."""
    return DREIDING_PARAM_FILE

__all__ = [
    "DREIDING_PARAM_FILE",
    "get_dreiding_params_path",
]
