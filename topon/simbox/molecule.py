"""
Molecule definition and 3D conformer generation for simbox.

Provides the Molecule class with factory methods for creating molecules
from SMILES strings, PDB files, or existing RDKit mol objects.  Reactive
sites (epoxides, amines, etc.) are auto-detected via SMARTS patterns so
that downstream packing and writing stages know which atoms participate
in crosslinking reactions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

# ---------------------------------------------------------------------------
# SMARTS patterns for reactive-site auto-detection
# ---------------------------------------------------------------------------
REACTIVE_SMARTS: dict[str, str] = {
    "epoxide": "[C]1[O][C]1",                  # oxirane ring
    "primary_amine": "[NX3;H2;!$([NH2]C=O)]",  # -NH2 (not amide)
    "secondary_amine": "[NX3;H1]([#6])[#6]",   # >NH between two carbons
}


@dataclass
class Molecule:
    """
    A molecule with a 3D conformer and reactive-site annotations.

    Attributes
    ----------
    name : str
        Human-readable identifier (e.g. ``"Epoxy-PDMS"``).
    mol : rdkit.Chem.rdchem.Mol
        RDKit Mol with explicit H and an embedded 3D conformer.
    mw : float
        Exact molecular weight (g mol-1).
    smiles : str | None
        Canonical SMILES, if the molecule was built from one.
    reactive_sites : dict[str, list[int]]
        Mapping from reactive-group name to a list of heavy-atom indices
        that belong to that group (e.g. ``{"epoxide": [12, 14]}``).
    """

    name: str
    mol: object                        # RDKit Mol (lazy import avoids top-level dep)
    mw: float
    smiles: Optional[str] = None
    reactive_sites: dict[str, list[int]] = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------
    @classmethod
    def from_smiles(
        cls,
        name: str,
        smiles: str,
        reactive_smarts: Optional[dict[str, str]] = None,
    ) -> Molecule:
        """Build a *Molecule* from a SMILES string.

        The method adds explicit H, embeds a 3D conformer (ETKDGv3),
        and optimises geometry with MMFF (falling back to UFF).

        Parameters
        ----------
        name : str
            Identifier for this molecule.
        smiles : str
            Valid SMILES string.
        reactive_smarts : dict, optional
            Extra / override SMARTS patterns for reactive-site detection.
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        mol = Chem.AddHs(mol)

        # Embed 3D conformer
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            # Retry with random coordinates for difficult molecules
            params.useRandomCoords = True
            result = AllChem.EmbedMolecule(mol, params)
            if result == -1:
                raise RuntimeError(
                    f"Failed to generate 3D conformer for '{name}'"
                )

        # Geometry optimisation
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)

        mw = Descriptors.ExactMolWt(mol)
        sites = cls._detect_reactive_sites(mol, reactive_smarts)

        return cls(
            name=name, mol=mol, mw=mw,
            smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)),
            reactive_sites=sites,
        )

    @classmethod
    def from_pdb(
        cls,
        name: str,
        pdb_path: str,
        reactive_sites: Optional[dict[str, list[int]]] = None,
        reactive_smarts: Optional[dict[str, str]] = None,
    ) -> Molecule:
        """Build a *Molecule* from a PDB file.

        Parameters
        ----------
        name : str
            Identifier for this molecule.
        pdb_path : str
            Path to a .pdb file.
        reactive_sites : dict, optional
            Manually specified reactive sites (skip auto-detection).
        reactive_smarts : dict, optional
            Extra SMARTS patterns for auto-detection.
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors

        pdb_path = Path(pdb_path)
        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")

        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)
        if mol is None:
            raise ValueError(f"Failed to parse PDB file: {pdb_path}")

        # Add explicit H if none are present
        if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
            mol = Chem.AddHs(mol, addCoords=True)

        mw = Descriptors.ExactMolWt(mol)
        sites = reactive_sites or cls._detect_reactive_sites(mol, reactive_smarts)

        return cls(name=name, mol=mol, mw=mw, reactive_sites=sites)

    @classmethod
    def from_mol(
        cls,
        name: str,
        rdkit_mol,
        reactive_sites: Optional[dict[str, list[int]]] = None,
        reactive_smarts: Optional[dict[str, str]] = None,
    ) -> Molecule:
        """Wrap an existing RDKit Mol object.

        The molecule *must* already have explicit H and a 3D conformer.
        """
        from rdkit.Chem import Descriptors

        mw = Descriptors.ExactMolWt(rdkit_mol)
        sites = reactive_sites or cls._detect_reactive_sites(
            rdkit_mol, reactive_smarts
        )
        return cls(name=name, mol=rdkit_mol, mw=mw, reactive_sites=sites)

    # ------------------------------------------------------------------
    # Reactive-site detection
    # ------------------------------------------------------------------
    @staticmethod
    def _detect_reactive_sites(
        mol, custom_smarts: Optional[dict[str, str]] = None
    ) -> dict[str, list[int]]:
        """Return reactive sites found via SMARTS sub-structure matching."""
        from rdkit import Chem

        patterns: dict[str, str] = dict(REACTIVE_SMARTS)
        if custom_smarts:
            patterns.update(custom_smarts)

        sites: dict[str, list[int]] = {}
        for group_name, smarts in patterns.items():
            pat = Chem.MolFromSmarts(smarts)
            if pat is None:
                continue
            matches = mol.GetSubstructMatches(pat)
            if matches:
                atom_indices = sorted(set(idx for match in matches for idx in match))
                sites[group_name] = atom_indices
        return sites

    # ------------------------------------------------------------------
    # Coordinate helpers
    # ------------------------------------------------------------------
    def get_coordinates(self) -> np.ndarray:
        """Return atom positions as an (N, 3) array (Angstrom)."""
        conf = self.mol.GetConformer()
        return np.array(
            [list(conf.GetAtomPosition(i)) for i in range(self.mol.GetNumAtoms())]
        )

    def get_centroid(self) -> np.ndarray:
        """Centroid of the molecule (Angstrom)."""
        return self.get_coordinates().mean(axis=0)

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------
    @property
    def num_atoms(self) -> int:
        return self.mol.GetNumAtoms()

    @property
    def num_bonds(self) -> int:
        return self.mol.GetNumBonds()

    def __repr__(self) -> str:
        sites = list(self.reactive_sites.keys())
        return (
            f"Molecule(name='{self.name}', atoms={self.num_atoms}, "
            f"mw={self.mw:.1f}, reactive_sites={sites})"
        )
