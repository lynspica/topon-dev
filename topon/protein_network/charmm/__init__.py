"""CHARMM36m atomistic protein-network builder.

Forked verbatim from the legacy `topro` package (subprojects/protein_network/
topro) — the BFM topology stage was already mirrored into
`topon.protein_network.bfm`; this subpackage adds the *atomistic* chemistry
stage that wasn't migrated when the MARTINI port landed.

Public entry point: `topon.protein_network.charmm.build_systems` (CLI), or
`build_protein_system` + `add_water_and_ions` + the LAMMPS writers below
if you're driving it from Python.

Bundled CHARMM36m PRM/RTF/CMAP files live in `data/`. They are the
'updated_charges' variants that ship in the legacy `protein_generator`
fork; pass `--ff36` to `build_systems` to override with your own paths.
"""
from .charmm_ff import CHARMMForceField
from .builder import (
    Atom,
    build_protein_system,
    add_water_and_ions,
    compute_lattice_scale,
    find_cmap_crossterms,
)
from .lammps_writer import (
    find_angles,
    find_dihedrals,
    build_type_maps,
    write_lammps_data,
    write_lammps_settings,
    write_lammps_groups,
    write_lammps_input,
)
from .topology_io import (
    save_topology,
    load_topology,
    get_snapshot,
    list_snapshots,
)

__all__ = [
    "CHARMMForceField",
    "Atom",
    "build_protein_system",
    "add_water_and_ions",
    "compute_lattice_scale",
    "find_cmap_crossterms",
    "find_angles",
    "find_dihedrals",
    "build_type_maps",
    "write_lammps_data",
    "write_lammps_settings",
    "write_lammps_groups",
    "write_lammps_input",
    "save_topology",
    "load_topology",
    "get_snapshot",
    "list_snapshots",
]
