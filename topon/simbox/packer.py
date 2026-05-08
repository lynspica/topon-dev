"""
Custom box packing algorithm for simbox.

Provides ``BoxPacker`` which places molecules at random positions and
orientations inside a periodic cubic box, using grid-based spatial
hashing for efficient overlap detection.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from topon.simbox.molecule import Molecule


# ---------------------------------------------------------------------------
# Rotation helpers
# ---------------------------------------------------------------------------
def _random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    """Uniform random rotation via Shoemake's quaternion method."""
    u1, u2, u3 = rng.random(3)
    q = np.array([
        np.sqrt(1 - u1) * np.sin(2 * np.pi * u2),
        np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
        np.sqrt(u1) * np.sin(2 * np.pi * u3),
        np.sqrt(u1) * np.cos(2 * np.pi * u3),
    ])
    # Quaternion → rotation matrix
    w, x, y, z = q
    return np.array([
        [1 - 2*(y*y + z*z),     2*(x*y - z*w),     2*(x*z + y*w)],
        [    2*(x*y + z*w), 1 - 2*(x*x + z*z),     2*(y*z - x*w)],
        [    2*(x*z - y*w),     2*(y*z + x*w), 1 - 2*(x*x + y*y)],
    ])


# ---------------------------------------------------------------------------
# Spatial hash grid for periodic overlap detection
# ---------------------------------------------------------------------------
class _SpatialGrid:
    """Cell-list for O(1) nearest-neighbour queries in a periodic box."""

    def __init__(self, box_lengths: np.ndarray, cell_size: float):
        self.box = box_lengths.copy()
        self.cell_size = cell_size
        self.n_cells = np.maximum(np.floor(self.box / cell_size).astype(int), 1)
        self.grid: dict[tuple, list[np.ndarray]] = defaultdict(list)

    def _cell_index(self, pos: np.ndarray) -> tuple[int, int, int]:
        wrapped = pos % self.box
        ci = np.minimum(
            (wrapped / self.cell_size).astype(int),
            self.n_cells - 1,
        )
        return tuple(ci)

    def _neighbour_offsets(self):
        """Yield 27 (dx, dy, dz) offsets covering the cell and its neighbours."""
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    yield np.array([dx, dy, dz])

    def insert(self, positions: np.ndarray) -> None:
        """Insert an array of positions (N, 3) into the grid."""
        for pos in positions:
            key = self._cell_index(pos)
            self.grid[key].append(pos.copy())

    def any_overlap(self, positions: np.ndarray, min_dist: float) -> bool:
        """Return True if *any* atom in *positions* is closer than *min_dist*
        to an existing atom (minimum-image convention)."""
        if min_dist <= 0.0:
            return False
        min_dist_sq = min_dist * min_dist
        half_box = 0.5 * self.box

        for pos in positions:
            base_cell = self._cell_index(pos)
            for offset in self._neighbour_offsets():
                nb_cell = tuple(
                    (np.array(base_cell) + offset) % self.n_cells
                )
                for other in self.grid.get(nb_cell, []):
                    delta = pos - other
                    # Minimum image convention
                    delta = delta - self.box * np.round(delta / self.box)
                    if np.dot(delta, delta) < min_dist_sq:
                        return True
        return False


# ---------------------------------------------------------------------------
# Placement record
# ---------------------------------------------------------------------------
@dataclass
class PlacedMolecule:
    """Bookkeeping for one placed molecule instance."""
    species_idx: int          # index into the species list
    instance_idx: int         # which copy of this species
    molecule: Molecule
    coordinates: np.ndarray   # (N_atoms, 3), already rotated & translated


# ---------------------------------------------------------------------------
# Packed system (output of packer)
# ---------------------------------------------------------------------------
@dataclass
class PackedBox:
    """Result of a packing run.

    Attributes
    ----------
    placements : list[PlacedMolecule]
        All placed molecule instances.
    box_lengths : ndarray
        (Lx, Ly, Lz) in Angstrom.
    density : float
        Target density used for box sizing (g/cm^3).
    """
    placements: list[PlacedMolecule] = field(default_factory=list)
    box_lengths: np.ndarray = field(default_factory=lambda: np.zeros(3))
    density: float = 0.0

    @property
    def total_atoms(self) -> int:
        return sum(p.molecule.num_atoms for p in self.placements)

    @property
    def total_molecules(self) -> int:
        return len(self.placements)


# ---------------------------------------------------------------------------
# Main packer
# ---------------------------------------------------------------------------
class BoxPacker:
    """Pack molecules into a periodic cubic (or orthorhombic) box.

    Parameters
    ----------
    density : float
        Target density in g/cm^3 (default 0.85).
    min_dist : float
        Minimum allowed interatomic distance in Angstrom (default 0.0).
        Lowering this allows minor overlaps which are resolved by soft-potential.
    seed : int | None
        Random seed for reproducibility.
    max_attempts : int
        Max random placement attempts per molecule before growing the box.
    growth_factor : float
        Factor by which to expand the box when placements fail.
    """

    AVOGADRO = 6.02214076e23

    def __init__(
        self,
        density: float = 0.85,
        min_dist: float = 0.0,
        seed: Optional[int] = None,
        max_attempts: int = 1000,
        growth_factor: float = 1.05,
    ):
        self.density = density
        self.min_dist = min_dist
        self.rng = np.random.default_rng(seed)
        self.max_attempts = max_attempts
        self.growth_factor = growth_factor

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    def pack(
        self,
        species: list[tuple[Molecule, int]],
        pre_placed: Optional[list[PlacedMolecule]] = None,
        box_lengths_override: Optional[np.ndarray] = None,
    ) -> PackedBox:
        """Pack molecules into a box.

        Parameters
        ----------
        species : list of (Molecule, count)
            Each tuple specifies a molecule type and how many copies to
            place.
        pre_placed : list of PlacedMolecule, optional
            Molecules that are already positioned (e.g. a polymer chain
            centered at the box midpoint).  Their atoms are registered in
            the spatial grid so new molecules avoid them, and they are
            prepended to the output placements list.
        box_lengths_override : ndarray, optional
            If given, use this box size instead of computing from density.

        Returns
        -------
        PackedBox
        """
        # --- 1. Compute box size from target density ---
        if box_lengths_override is not None:
            box_lengths = box_lengths_override.copy()
            L = box_lengths[0]
            # Compute mass for reporting
            total_mass_g = sum(mol.mw * count for mol, count in species)
            if pre_placed:
                total_mass_g += sum(pm.molecule.mw for pm in pre_placed)
        else:
            total_mass_g = sum(mol.mw * count for mol, count in species)
            if pre_placed:
                total_mass_g += sum(pm.molecule.mw for pm in pre_placed)
            total_mass_real_g = total_mass_g / self.AVOGADRO
            volume_cm3 = total_mass_real_g / self.density
            volume_A3 = volume_cm3 * 1e24
            L = volume_A3 ** (1.0 / 3.0)
            box_lengths = np.array([L, L, L])
            total_mass_g = total_mass_g  # keep g/mol for print

        print(f"[BoxPacker] Target density = {self.density:.3f} g/cm³")
        print(f"[BoxPacker] Initial box    = {L:.2f} Å  ({L:.2f} × {L:.2f} × {L:.2f})")

        # --- 2. Build insertion list (shuffled) ---
        insertion_list: list[tuple[int, int, Molecule]] = []
        sp_offset = len(pre_placed) if pre_placed else 0
        for sp_idx, (mol, count) in enumerate(species):
            for inst_idx in range(count):
                insertion_list.append((sp_idx + sp_offset, inst_idx, mol))
        self.rng.shuffle(insertion_list)

        # --- 3. Pack with overlap detection ---
        grid = _SpatialGrid(box_lengths, cell_size=self.min_dist)
        placements: list[PlacedMolecule] = []

        # Pre-register already-placed molecules (e.g. the chain)
        if pre_placed:
            for pm in pre_placed:
                grid.insert(pm.coordinates)
                placements.append(pm)

        n_total = len(insertion_list)

        for i, (sp_idx, inst_idx, mol) in enumerate(insertion_list):
            coords_template = mol.get_coordinates()
            # Centre the template at the origin
            coords_template -= coords_template.mean(axis=0)

            placed = False
            for attempt in range(self.max_attempts):
                # Random rotation + translation
                R = _random_rotation_matrix(self.rng)
                coords_rot = coords_template @ R.T
                shift = self.rng.random(3) * box_lengths
                coords_placed = coords_rot + shift
                # NOTE: do NOT wrap coords into the box.  Molecule atoms
                # must stay physically contiguous so that LAMMPS can later
                # assign correct image flags.  The _SpatialGrid overlap
                # check already uses minimum-image convention.

                if not grid.any_overlap(coords_placed, self.min_dist):
                    grid.insert(coords_placed)
                    placements.append(PlacedMolecule(
                        species_idx=sp_idx,
                        instance_idx=inst_idx,
                        molecule=mol,
                        coordinates=coords_placed,
                    ))
                    placed = True
                    break

            if not placed:
                # Grow box iteratively until molecule fits
                max_growth_rounds = 20
                for growth_round in range(max_growth_rounds):
                    box_lengths *= self.growth_factor
                    print(
                        f"[BoxPacker] Growing box to {box_lengths[0]:.2f} Å "
                        f"(molecule {i+1}/{n_total}, growth #{growth_round+1})"
                    )
                    grid = _SpatialGrid(box_lengths, cell_size=self.min_dist)
                    # Re-insert all previously placed molecules (no wrapping)
                    for pm in placements:
                        grid.insert(pm.coordinates)
                    # Retry this molecule
                    for attempt in range(self.max_attempts):
                        R = _random_rotation_matrix(self.rng)
                        coords_rot = coords_template @ R.T
                        shift = self.rng.random(3) * box_lengths
                        coords_placed = coords_rot + shift
                        if not grid.any_overlap(coords_placed, self.min_dist):
                            grid.insert(coords_placed)
                            placements.append(PlacedMolecule(
                                species_idx=sp_idx,
                                instance_idx=inst_idx,
                                molecule=mol,
                                coordinates=coords_placed,
                            ))
                            placed = True
                            break
                    if placed:
                        break

                if not placed:
                    raise RuntimeError(
                        f"Failed to place molecule {i+1}/{n_total} "
                        f"({mol.name} #{inst_idx}) after {max_growth_rounds} "
                        f"growth rounds. Try lowering the density or min_dist."
                    )

            if (i + 1) % max(1, n_total // 10) == 0 or i + 1 == n_total:
                print(f"[BoxPacker] Placed {i+1}/{n_total} molecules")

        print(
            f"[BoxPacker] Done - {len(placements)} molecules, "
            f"{sum(p.molecule.num_atoms for p in placements)} atoms, "
            f"box = {box_lengths[0]:.2f} × {box_lengths[1]:.2f} × {box_lengths[2]:.2f} Å"
        )

        return PackedBox(
            placements=placements,
            box_lengths=box_lengths,
            density=self.density,
        )
