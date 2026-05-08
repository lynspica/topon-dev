"""Voxel-grid W (MARTINI water) bead packer.

Fills an orthogonal box with W beads on a regular grid, rejecting any voxel
whose centre lies within ``exclusion_radius_ang`` of a protein bead. Uses a
spatial cell list to keep the exclusion check linear in the number of beads
to place. Pack target density defaults to MARTINI 3 bulk water (~10 W/nm^3,
which corresponds to one W per ~0.464 nm grid step).

Each placed W is added to the system as its own molecule (distinct
``molecule_id``), matching how GROMACS .top files list water as ``W <count>``.

Dityrosine crosslinks are not implemented here -- they are emitted directly by
``builder.build_protein_system`` from the BFM snapshot's ``reactions`` list.
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

from .martini_ff import MartiniLibrary
from .system import Bead, System

DEFAULT_W_DENSITY_NM3: float = 10.0   # MARTINI 3 bulk water
DEFAULT_EXCLUSION_ANG: float = 4.0    # ~MARTINI bead diameter

# MARTINI 3 water variants. Mass = 18 g/mol per H2O * mapping ratio.
WATER_BEAD_TYPES: dict[str, int] = {
    "W":  4,   # regular: 1 bead = 4 H2O   (bulk water default)
    "SW": 3,   # small:   1 bead = 3 H2O   (confined water, e.g. pores)
    "TW": 2,   # tiny:    1 bead = 2 H2O   (very tight pockets)
}


def _grid_step_ang(density_w_per_nm3: float) -> float:
    """One W per voxel: voxel volume = 1 / density. Edge = volume^(1/3) (nm) * 10."""
    return (1.0 / density_w_per_nm3) ** (1.0 / 3.0) * 10.0


def _cell_list(positions: np.ndarray, box: np.ndarray, cell_size: float) -> dict[tuple[int, int, int], list[int]]:
    """Hash positions into spatial cells; cell coords wrap periodically."""
    cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    if positions.size == 0:
        return cells
    nx = max(1, int(math.ceil(box[0] / cell_size)))
    ny = max(1, int(math.ceil(box[1] / cell_size)))
    nz = max(1, int(math.ceil(box[2] / cell_size)))
    for i, p in enumerate(positions):
        cx = int(p[0] // cell_size) % nx
        cy = int(p[1] // cell_size) % ny
        cz = int(p[2] // cell_size) % nz
        cells[(cx, cy, cz)].append(i)
    return cells


def _too_close_to_any(
    p: np.ndarray,
    cells: dict[tuple[int, int, int], list[int]],
    positions: np.ndarray,
    box: np.ndarray,
    cell_size: float,
    excl_sq: float,
) -> bool:
    if positions.size == 0:
        return False
    nx = max(1, int(math.ceil(box[0] / cell_size)))
    ny = max(1, int(math.ceil(box[1] / cell_size)))
    nz = max(1, int(math.ceil(box[2] / cell_size)))
    cx = int(p[0] // cell_size) % nx
    cy = int(p[1] // cell_size) % ny
    cz = int(p[2] // cell_size) % nz
    # Inspect 3x3x3 neighbouring cells (with periodic wrap)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                key = ((cx + dx) % nx, (cy + dy) % ny, (cz + dz) % nz)
                idxs = cells.get(key)
                if not idxs:
                    continue
                others = positions[idxs]
                d = others - p
                d -= box * np.round(d / box)  # min-image
                if np.any(np.sum(d * d, axis=1) < excl_sq):
                    return True
    return False


def pack_water(
    sys_: System,
    library: MartiniLibrary,
    *,
    density_w_per_nm3: float = DEFAULT_W_DENSITY_NM3,
    exclusion_radius_ang: float = DEFAULT_EXCLUSION_ANG,
    seed: int = 42,
    max_beads: int | None = None,
    bead_type: str = "W",
) -> int:
    """Add water beads to `sys_` on a voxel grid, avoiding existing protein beads.

    `bead_type`: one of ``"W"`` (4 H2O/bead, default), ``"SW"`` (3 H2O/bead),
    ``"TW"`` (2 H2O/bead). Use SW/TW for confined-water studies where finer
    spatial resolution is needed than the regular 4-to-1 mapping. The number
    of waters PER BEAD is reported in ``WATER_BEAD_TYPES``.

    `density_w_per_nm3` is BEADS per nm^3, regardless of bead type. To preserve
    the actual H2O number-density when switching mappings, scale this value:
    ``density_SW = density_W * 4/3`` and ``density_TW = density_W * 4/2``.

    Returns the number of beads added. Mutates ``sys_.beads`` in place;
    each bead gets a fresh ``molecule_id`` (one molecule per water bead, matching
    GROMACS convention).
    """
    if bead_type not in WATER_BEAD_TYPES:
        raise ValueError(
            f"Unknown water bead type {bead_type!r}; supported: "
            f"{sorted(WATER_BEAD_TYPES)} (W=4 H2O/bead, SW=3, TW=2)."
        )
    box = np.array(sys_.box_dims_ang, dtype=float)
    if np.any(box <= 0):
        raise ValueError(f"system box not initialised: {sys_.box_dims_ang}")

    rng = np.random.default_rng(seed)
    grid = _grid_step_ang(density_w_per_nm3)
    n_grid = np.maximum(np.floor(box / grid).astype(int), 1)
    actual_grid = box / n_grid  # adjust to fit box exactly

    # Build cell list from existing protein beads
    if sys_.beads:
        protein_pos = np.array([b.position for b in sys_.beads], dtype=float)
    else:
        protein_pos = np.empty((0, 3), dtype=float)
    cell_size = max(exclusion_radius_ang, float(np.max(actual_grid)))
    cells = _cell_list(protein_pos, box, cell_size)
    excl_sq = exclusion_radius_ang ** 2

    next_atom_id = max((b.atom_id for b in sys_.beads), default=0) + 1
    next_mol_id = max((b.molecule_id for b in sys_.beads), default=0) + 1
    default_mass = 18.0 * WATER_BEAD_TYPES[bead_type]   # 18 g/mol per H2O * mapping
    w_mass = library.atomtypes[bead_type].mass if bead_type in library.atomtypes else default_mass

    placed = 0
    for ix in range(n_grid[0]):
        for iy in range(n_grid[1]):
            for iz in range(n_grid[2]):
                if max_beads is not None and placed >= max_beads:
                    return placed
                # voxel centre + small jitter (deterministic per voxel via rng)
                base = np.array([
                    (ix + 0.5) * actual_grid[0],
                    (iy + 0.5) * actual_grid[1],
                    (iz + 0.5) * actual_grid[2],
                ])
                jitter = (rng.random(3) - 0.5) * actual_grid * 0.1
                p = base + jitter
                p = p - box * np.floor(p / box)  # wrap into box
                if _too_close_to_any(p, cells, protein_pos, box, cell_size, excl_sq):
                    continue
                bead = Bead(
                    atom_id=next_atom_id,
                    bead_type=bead_type,
                    molecule_id=next_mol_id,
                    residue_idx=0,
                    residue_name=bead_type,
                    atom_name=bead_type,
                    charge=0.0,
                    mass=w_mass,
                    position=tuple(p.tolist()),
                )
                sys_.beads.append(bead)
                next_atom_id += 1
                next_mol_id += 1
                placed += 1
    return placed
