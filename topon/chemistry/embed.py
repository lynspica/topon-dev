"""3D conformer embedding for polymer chain molecules.

The legacy :func:`topon.singlechain.workflow._assign_extended_linear_coords`
walks only the backbone heavy-atom path and leaves branch atoms (e.g. the
three methyls on each Si in PDMS) collapsed onto their parent atom.  At
DP=30 PDMS this produced an observed minimum pair-atom distance of
**0.436 Å** (Si-methyl C on top of a backbone O), which translates to
pair potentials of ~−10²⁷ kcal/mol and LAMMPS aborting at step 1 with
``Non-numeric box dimensions — simulation unstable``.

This module provides a replacement that uses RDKit's ETKDGv3 embedding
with ``useRandomCoords=True`` followed by MMFF94 (with UFF fallback).
It costs ~0.5 s per DP=30 chain and produces chemically sensible bond
lengths and angles from step 0.

The legacy placer is retained for backward compatibility; callers
should migrate to :func:`embed_with_etkdg` for any polymer with side
chains or stereochemistry.
"""

from __future__ import annotations

from typing import Optional


class EmbedFailedError(RuntimeError):
    """Raised when RDKit ETKDGv3 cannot place a 3D conformer."""


def embed_with_etkdg(mol, *, seed: int = 42, mmff_iters: int = 200,
                     uff_iters: int = 150) -> None:
    """Embed a 3D conformer with ETKDGv3 and relax with MMFF94 / UFF.

    Operates in place on the supplied RDKit ``Mol`` (which must already
    have explicit hydrogens via :func:`rdkit.Chem.AddHs`).

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule with explicit hydrogens.
    seed : int
        RNG seed for ETKDGv3.  Clamped to 31 bits internally because
        RDKit's ``randomSeed`` is a C ``long``.
    mmff_iters : int
        Max MMFF94 optimisation iterations (default 200).  Downstream
        LAMMPS NPT equilibration handles the residual strain, so spending
        many iterations here is wasted effort.  Set to 0 to skip.
    uff_iters : int
        Max UFF optimisation iterations (default 150) — used as a
        fallback when MMFF94 cannot parameterise the molecule (e.g.
        Si-containing chains).  Set to 0 to skip.

    Raises
    ------
    EmbedFailedError
        If ETKDGv3 cannot embed the molecule even with relaxed
        ``enforceChirality=False`` and ``ignoreSmoothingFailures=True``.
        Callers should treat this as a hard failure: the molecule's
        SMILES likely contains geometry RDKit cannot satisfy.
    """
    from rdkit.Chem import AllChem

    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed) & 0x7FFFFFFF
    params.useRandomCoords = True

    if AllChem.EmbedMolecule(mol, params) == -1:
        # Second chance: relax stereochemistry / distance geometry
        params.ignoreSmoothingFailures = True
        params.enforceChirality = False
        if AllChem.EmbedMolecule(mol, params) == -1:
            raise EmbedFailedError(
                "RDKit ETKDGv3 failed to embed a 3D conformer. "
                "The molecule's SMILES may contain geometry RDKit cannot "
                "satisfy (e.g. overconstrained distances or ring strain)."
            )

    mmff_ok = False
    if mmff_iters > 0:
        try:
            # Returns 0 on convergence, 1 if not converged, -1 if the FF
            # could not be parameterised (common for Si/transition metals).
            rc = AllChem.MMFFOptimizeMolecule(mol, maxIters=mmff_iters)
            mmff_ok = (rc != -1)
        except Exception:
            mmff_ok = False

    # UFF is a universal-fallback: only run it when MMFF failed to type the
    # molecule (e.g. Si containing).  For organics where MMFF succeeded,
    # UFF wastes several seconds per DP=30 chain (observed: 24 s extra
    # on a 554-atom Butyl oligomer — see solubility SPEC_QUESTIONS_V2.md).
    if not mmff_ok and uff_iters > 0:
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=uff_iters)
        except Exception:
            pass


def min_pairwise_distance(mol) -> float:
    """Return the minimum pairwise heavy+H distance (Å) over the conformer.

    Diagnostic helper — a post-embed value below ~0.9 Å signals collapsed
    branch atoms (the bug that motivated this module).
    """
    import numpy as np

    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    coords = np.empty((n, 3), dtype=float)
    for i in range(n):
        p = conf.GetAtomPosition(i)
        coords[i] = (p.x, p.y, p.z)
    diffs = coords[:, None, :] - coords[None, :, :]
    dists = (diffs * diffs).sum(-1) ** 0.5
    # Ignore the self-distance diagonal
    import numpy as _np
    _np.fill_diagonal(dists, _np.inf)
    return float(dists.min())
