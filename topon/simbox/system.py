"""
System assembly for simbox.

Merges individually packed molecules into a single unified RDKit Mol
with correct molecule-ID bookkeeping and a global reactive-site registry,
ready for the LAMMPS writer.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from topon.simbox.packer import PackedBox


# ---------------------------------------------------------------------------
# Reactive-site record
# ---------------------------------------------------------------------------
@dataclass
class ReactiveSiteEntry:
    """One reactive atom in the merged system."""
    global_atom_idx: int     # 0-based index in the merged Mol
    molecule_id: int         # 1-based LAMMPS molecule ID
    group_name: str          # e.g. "epoxide", "primary_amine"
    species_name: str        # e.g. "Epoxy-PDMS-n2"


# ---------------------------------------------------------------------------
# Assembled system
# ---------------------------------------------------------------------------
@dataclass
class AssembledSystem:
    """A fully assembled system ready for the LAMMPS writer.

    Attributes
    ----------
    mol : RDKit Mol
        Combined molecule (all atoms, explicit H, 3D coords).
    box_lengths : ndarray
        (Lx, Ly, Lz) in Angstrom.
    molecule_ids : list[int]
        Per-atom molecule ID (1-based, length = num_atoms).
    species_names : list[str]
        Per-molecule species name (length = num_molecules).
    reactive_sites : list[ReactiveSiteEntry]
        Global reactive-site registry.
    num_molecules : int
        Total number of molecule instances.
    """
    mol: object                                                # RDKit Mol
    box_lengths: np.ndarray
    molecule_ids: list[int] = field(default_factory=list)
    species_names: list[str] = field(default_factory=list)
    reactive_sites: list[ReactiveSiteEntry] = field(default_factory=list)
    num_molecules: int = 0


def assemble(packed: PackedBox) -> AssembledSystem:
    """Merge a *PackedBox* into a single :class:`AssembledSystem`.

    The function:
    1. Combines all individual RDKit mols into one via ``CombineMols``.
    2. Overwrites 3D coordinates with packed positions.
    3. Assigns per-atom molecule IDs (1-based).
    4. Builds a global reactive-site registry by offsetting each
       molecule's local reactive-site indices.

    Parameters
    ----------
    packed : PackedBox
        Output of :meth:`BoxPacker.pack`.

    Returns
    -------
    AssembledSystem
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    combined: Chem.RWMol | None = None
    atom_offset = 0
    molecule_ids: list[int] = []
    species_names: list[str] = []
    reactive_sites: list[ReactiveSiteEntry] = []

    for mol_id_0based, placement in enumerate(packed.placements):
        mol_id = mol_id_0based + 1  # 1-based for LAMMPS
        template_mol = placement.molecule.mol
        n_atoms = template_mol.GetNumAtoms()

        # Merge into combined Mol
        if combined is None:
            combined = Chem.RWMol(template_mol)
        else:
            combined = Chem.RWMol(Chem.CombineMols(combined.GetMol(), template_mol))

        # Write packed coordinates into the conformer
        conf = combined.GetConformer()
        for local_idx in range(n_atoms):
            global_idx = atom_offset + local_idx
            x, y, z = placement.coordinates[local_idx]
            conf.SetAtomPosition(global_idx, Point3D(float(x), float(y), float(z)))

        # Per-atom molecule IDs
        molecule_ids.extend([mol_id] * n_atoms)

        # Species tracking
        species_names.append(placement.molecule.name)

        # Reactive-site registry (offset local indices → global)
        for group_name, local_indices in placement.molecule.reactive_sites.items():
            for local_idx in local_indices:
                reactive_sites.append(ReactiveSiteEntry(
                    global_atom_idx=atom_offset + local_idx,
                    molecule_id=mol_id,
                    group_name=group_name,
                    species_name=placement.molecule.name,
                ))

        atom_offset += n_atoms

    if combined is None:
        raise ValueError("No molecules to assemble (empty PackedBox)")

    mol_final = combined.GetMol()

    print(
        f"[Assembler] Merged {len(packed.placements)} molecules -> "
        f"{mol_final.GetNumAtoms()} atoms, "
        f"{len(reactive_sites)} reactive sites"
    )

    return AssembledSystem(
        mol=mol_final,
        box_lengths=packed.box_lengths.copy(),
        molecule_ids=molecule_ids,
        species_names=species_names,
        reactive_sites=reactive_sites,
        num_molecules=len(packed.placements),
    )
