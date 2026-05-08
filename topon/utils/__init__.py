"""
Topon Utilities Module

Helper functions for polymer network generation.
"""

from pathlib import Path

# Module paths
UTILS_DIR = Path(__file__).parent

# Re-export commonly used items
from .network_config import PENDANT_GROUP_SMILES, DREIDING_PARAM_FILE
from .network_helpers import (
    calculate_entangled_kink,
    find_crossing_candidates,
    create_chain_mol,
    generate_chain_string,
    resolve_smiles,
    graft_side_chain,
    write_lammps_displacement_file,
    write_group_definitions_to_file,
    generate_approximate_side_chain_coords,
    generate_poss_coordinates,
)

__all__ = [
    "UTILS_DIR",
    "PENDANT_GROUP_SMILES",
    "DREIDING_PARAM_FILE",
    "calculate_entangled_kink",
    "find_crossing_candidates",
    "create_chain_mol",
    "generate_chain_string",
    "resolve_smiles",
    "graft_side_chain",
    "write_lammps_displacement_file",
    "write_group_definitions_to_file",
    "generate_approximate_side_chain_coords",
    "generate_poss_coordinates",
]
