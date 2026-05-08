"""Pack NA / CL ions (MARTINI 3, both bead type TQ5) into a System.

The reference resilin systems in `Martini_Ahmet/top_files/*/system_wr.top` use
527 NA + 527 CL for charge neutrality at 0.15 M. This packer adds them at
random positions, avoiding existing protein and water beads, mirroring how
GROMACS handles ion insertion.
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

from .martini_ff import MartiniLibrary
from .system import Bead, System

ION_BEAD_TYPE: str = "TQ5"      # MARTINI 3: same bead type for NA and CL
ION_MASS: float = 36.0
NA_CHARGE: float = +1.0
CL_CHARGE: float = -1.0


def pack_ions(
    sys_: System,
    library: MartiniLibrary,
    *,
    n_na: int,
    n_cl: int,
    exclusion_radius_ang: float = 4.0,
    seed: int = 42,
) -> tuple[int, int]:
    """Add NA + CL ions (MARTINI 3 TQ5 bead) to `sys_`. Returns (n_na, n_cl) actually placed.

    Each ion is its own molecule (matches GROMACS [ molecules ] convention).
    Random rejection sampling: pick a random box position, accept if no
    existing bead is within `exclusion_radius_ang`.
    """
    box = np.array(sys_.box_dims_ang, dtype=float)
    if np.any(box <= 0):
        raise ValueError(f"system box not initialised: {sys_.box_dims_ang}")
    rng = np.random.default_rng(seed)

    # Existing positions (cell list for fast exclusion check).
    if sys_.beads:
        positions = np.array([b.position for b in sys_.beads], dtype=float)
    else:
        positions = np.empty((0, 3), dtype=float)

    cell_size = max(exclusion_radius_ang * 1.5, 5.0)
    nx = max(1, int(math.ceil(box[0] / cell_size)))
    ny = max(1, int(math.ceil(box[1] / cell_size)))
    nz = max(1, int(math.ceil(box[2] / cell_size)))
    cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for i, p in enumerate(positions):
        cx = int(p[0] // cell_size) % nx
        cy = int(p[1] // cell_size) % ny
        cz = int(p[2] // cell_size) % nz
        cells[(cx, cy, cz)].append(i)
    excl_sq = exclusion_radius_ang ** 2

    def _too_close(p: np.ndarray) -> bool:
        cx = int(p[0] // cell_size) % nx
        cy = int(p[1] // cell_size) % ny
        cz = int(p[2] // cell_size) % nz
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    key = ((cx + dx) % nx, (cy + dy) % ny, (cz + dz) % nz)
                    idxs = cells.get(key)
                    if not idxs: continue
                    others = positions[idxs]
                    d = others - p
                    d -= box * np.round(d / box)
                    if np.any(np.sum(d * d, axis=1) < excl_sq):
                        return True
        return False

    next_atom_id = max((b.atom_id for b in sys_.beads), default=0) + 1
    next_mol_id = max((b.molecule_id for b in sys_.beads), default=0) + 1
    mass = library.atomtypes[ION_BEAD_TYPE].mass if ION_BEAD_TYPE in library.atomtypes else ION_MASS

    placed_na = 0
    placed_cl = 0
    max_attempts_per_ion = 1000

    for kind, n_target, charge, atom_name in (
        ("NA", n_na, NA_CHARGE, "NA"),
        ("CL", n_cl, CL_CHARGE, "CL"),
    ):
        for _ in range(n_target):
            for _attempt in range(max_attempts_per_ion):
                p = rng.uniform(0.0, 1.0, size=3) * box
                if not _too_close(p):
                    break
            else:
                # could not place; skip
                continue
            sys_.beads.append(Bead(
                atom_id=next_atom_id,
                bead_type=ION_BEAD_TYPE,
                molecule_id=next_mol_id,
                residue_idx=0,
                residue_name="ION",
                atom_name=atom_name,
                charge=charge,
                mass=mass,
                position=tuple(p.tolist()),
            ))
            # update cell list with the new ion so subsequent ions avoid it
            cx = int(p[0] // cell_size) % nx
            cy = int(p[1] // cell_size) % ny
            cz = int(p[2] // cell_size) % nz
            cells[(cx, cy, cz)].append(len(positions))
            positions = np.vstack([positions, p[None, :]])
            next_atom_id += 1
            next_mol_id += 1
            if kind == "NA":
                placed_na += 1
            else:
                placed_cl += 1
    return placed_na, placed_cl
