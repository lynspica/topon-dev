"""
topon.simbox – Simulation Box Builder for Crosslink Reaction Simulations
========================================================================

Build simulation boxes from individual molecules (SMILES, PDB, or RDKit)
for use in LAMMPS crosslink reaction simulations.

Quick start::

    from topon.simbox import SimBox, MoleculeLibrary

    lib = MoleculeLibrary()
    epoxy = lib.epoxy_pdms(n_dms=2)      # ~500 g/mol
    amino = lib.amino_pdms(n_dms=8)      # ~850 g/mol
    poss  = lib.am0270_poss()             # ~1267 g/mol

    box = SimBox(density=0.85)
    box.add(epoxy, count=200)
    box.add(amino, count=100)
    box.add(poss, count=50)
    box.pack(seed=42)

    box.write("output/", forcefield="dreiding")
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from topon.simbox.molecule import Molecule
from topon.simbox.library import MoleculeLibrary

__all__ = ["SimBox", "Molecule", "MoleculeLibrary"]


class SimBox:
    """High-level facade for building a simulation box.

    Parameters
    ----------
    density : float
        Target packing density in g/cm^3 (default 0.85).
    min_dist : float
        Minimum inter-atomic distance (Angstrom) enforced during
        packing (default 2.0).
    temperature : float
        Target simulation temperature in K (default 300).
    pressure : float
        Target pressure in atm (default 1.0).
    """

    def __init__(
        self,
        density: float = 0.85,
        min_dist: float = 2.0,
        temperature: float = 300.0,
        pressure: float = 1.0,
    ):
        self.density = density
        self.min_dist = min_dist
        self.temperature = temperature
        self.pressure = pressure

        self._species: list[tuple[Molecule, int]] = []
        self._packed = None          # PackedBox after pack()
        self._system = None          # AssembledSystem after pack()

    # ------------------------------------------------------------------
    # Add molecules
    # ------------------------------------------------------------------
    def add(self, molecule: Molecule, count: int) -> None:
        """Register *count* copies of *molecule* for packing.

        Parameters
        ----------
        molecule : Molecule
            A molecule built via :class:`MoleculeLibrary`, ``Molecule.from_smiles``,
            ``Molecule.from_pdb``, etc.
        count : int
            Number of copies to place in the box.
        """
        if count < 1:
            raise ValueError("count must be >= 1")
        self._species.append((molecule, count))
        # Reset cached state
        self._packed = None
        self._system = None

    # ------------------------------------------------------------------
    # Pack
    # ------------------------------------------------------------------
    def pack(self, seed: Optional[int] = None, max_attempts: int = 1000) -> None:
        """Pack all registered molecules into a periodic box.

        Parameters
        ----------
        seed : int, optional
            Random seed for reproducibility.
        max_attempts : int
            Maximum random-placement attempts per molecule.
        """
        if not self._species:
            raise RuntimeError("No molecules added - call .add() first")

        from topon.simbox.packer import BoxPacker
        from topon.simbox.system import assemble

        packer = BoxPacker(
            density=self.density,
            min_dist=self.min_dist,
            seed=seed,
            max_attempts=max_attempts,
        )
        self._packed = packer.pack(self._species)
        self._system = assemble(self._packed)

    # ------------------------------------------------------------------
    # Write
    # ------------------------------------------------------------------
    def write(
        self,
        output_dir: str,
        forcefield: str = "dreiding",
        data_filename: str = "system.data",
    ) -> dict[str, str]:
        """Write LAMMPS data file, group definitions, and input scripts.

        Parameters
        ----------
        output_dir : str
            Target directory (created if it does not exist).
        forcefield : str
            Force field to use (currently only ``"dreiding"`` is supported).
        data_filename : str
            Name for the LAMMPS data file.

        Returns
        -------
        dict[str, str]
            Mapping of logical file names to their paths.
        """
        if self._system is None:
            raise RuntimeError("System not packed - call .pack() first")

        if forcefield.lower() != "dreiding":
            raise ValueError(
                f"Unsupported force field '{forcefield}'. "
                "Currently only 'dreiding' is supported."
            )

        from topon.simbox.writer import write_lammps
        from topon.simbox.inputs import write_inputs

        files: dict[str, str] = {}

        # Data file + groups
        written = write_lammps(
            self._system, output_dir, data_filename=data_filename,
        )
        files.update(written)

        # Input scripts
        input_files = write_inputs(
            self._system,
            output_dir,
            temperature=self.temperature,
            pressure=self.pressure,
            data_filename=data_filename,
        )
        files.update(input_files)

        return files

    # ------------------------------------------------------------------
    # Inspection helpers
    # ------------------------------------------------------------------
    @property
    def system(self):
        """The assembled system (available after ``.pack()``)."""
        return self._system

    @property
    def packed(self):
        """The raw packed box (available after ``.pack()``)."""
        return self._packed

    def summary(self) -> str:
        """Return a human-readable summary of the current state."""
        lines = ["SimBox Summary", "=" * 40]

        if not self._species:
            lines.append("  (no molecules added)")
        else:
            total_mols = sum(c for _, c in self._species)
            total_atoms = sum(m.num_atoms * c for m, c in self._species)
            lines.append(f"  Species:    {len(self._species)}")
            lines.append(f"  Molecules:  {total_mols}")
            lines.append(f"  Atoms:      {total_atoms}")
            lines.append(f"  Density:    {self.density} g/cm³")
            lines.append(f"  Min dist:   {self.min_dist} Å")
            lines.append("")
            for mol, count in self._species:
                lines.append(f"  {mol.name:30s}  ×{count:>5d}  (MW={mol.mw:.1f})")

        if self._system is not None:
            bx = self._system.box_lengths
            lines.append("")
            lines.append(f"  Box: {bx[0]:.2f} × {bx[1]:.2f} × {bx[2]:.2f} Å")
            n_sites = len(self._system.reactive_sites)
            lines.append(f"  Reactive sites: {n_sites}")

        return "\n".join(lines)

    def __repr__(self) -> str:
        n = sum(c for _, c in self._species) if self._species else 0
        packed = "packed" if self._system else "not packed"
        return f"SimBox(density={self.density}, molecules={n}, {packed})"
