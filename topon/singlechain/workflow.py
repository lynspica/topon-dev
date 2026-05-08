"""
Single-Chain in Solvent workflow for Topon.

Builds one atomistic polymer chain (DREIDING force field) with optional grafts
and copolymer sequence, then packs it into a periodic box with solvent molecules.

Supports **multi-solvent mixtures** via ``solvent_mixture`` parameter, and
always **centers the polymer chain** at the box center for clean analysis.

Usage::

    from topon.singlechain.workflow import run_workflow

    # Single solvent
    run_workflow(
        output_dir="chain_in_toluene",
        chain_smiles="[Si](C)(C)O",    # PDMS repeat unit
        dp=20,
        solvent_smiles="Cc1ccccc1",    # toluene
        n_solvent=200,
    )

    # Multi-solvent mixture (50/50 w/w iso-octane + toluene)
    run_workflow(
        output_dir="chain_in_fuelC",
        chain_smiles="[Si](C)(C)O",
        dp=20,
        solvent_mixture=[
            {"smiles": "CC(C)CC(C)(C)C", "weight_fraction": 0.5},
            {"smiles": "Cc1ccccc1",       "weight_fraction": 0.5},
        ],
    )

Or via CLI::

    topon chain --chain-smiles "[Si](C)(C)O" --dp 20 \\
                --solvent-smiles "Cc1ccccc1" --n-solvent 200

    topon chain --chain-smiles "[Si](C)(C)O" --dp 20 \\
                --solvent-mixture '[{"smiles":"CC(C)CC(C)(C)C","weight_fraction":0.5},{"smiles":"Cc1ccccc1","weight_fraction":0.5}]'
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Optional


def run_workflow(
    output_dir: str,
    chain_smiles: str,
    dp: int,
    solvent_smiles: Optional[str] = None,
    n_solvent: Optional[int] = None,
    solvent_mixture: Optional[list[dict]] = None,
    # Optional chain modifiers
    graft_density: float = 0.0,
    graft_smiles: Optional[str] = None,
    graft_dp: int = 5,
    copolymer_arrangement: Optional[str] = None,
    copolymer_composition: Optional[list[dict]] = None,
    # Box / packing
    density: float = 0.85,
    seed: int = 42,
    verbose: bool = True,
    # v2: embedder selector.  "linear" = legacy, buggy for branched chains;
    # "etkdg" = RDKit ETKDGv3 + MMFF94/UFF (recommended).
    embedder: str = "linear",
) -> dict:
    """
    Build a single polymer chain in solvent and write DREIDING LAMMPS files.

    Parameters
    ----------
    output_dir : str
        Directory for output files (created if absent).
    chain_smiles : str
        SMILES for the polymer repeat unit (e.g. ``"[Si](C)(C)O"`` for PDMS).
    dp : int
        Degree of polymerization (number of repeat units).
    solvent_smiles : str, optional
        SMILES for a single solvent molecule.  Ignored when ``solvent_mixture``
        is provided.  Defaults to toluene (``"Cc1ccccc1"``) if neither
        ``solvent_smiles`` nor ``solvent_mixture`` is given.
    n_solvent : int, optional
        Number of solvent molecules to pack (single-solvent mode only).  If
        ``None`` (default), the count is calculated automatically so the box
        edge is ≥ 1.5× the fully-extended backbone contour length.
    solvent_mixture : list[dict], optional
        Multi-component solvent specification.  Each entry is a dict with:

        - ``"smiles"`` — solvent SMILES string
        - ``"weight_fraction"`` — fraction of total solvent mass (values are
          normalised internally, so ``0.5 / 0.5`` is identical to ``1 / 1``).

        When provided, ``solvent_smiles`` and ``n_solvent`` are ignored; the
        total number of molecules per species is determined automatically from
        the target box size and weight fractions.
    graft_density : float
        Probability of a side-chain attachment at each backbone unit (0–1).
        Set 0 (default) to disable grafts.
    graft_smiles : str, optional
        SMILES for the graft repeat unit. Required when ``graft_density > 0``.
    graft_dp : int
        Number of repeat units per side chain (default 5).
    copolymer_arrangement : str, optional
        ``"block"``, ``"alternating"``, ``"random"``, or ``"gradient"``.
        ``None`` = homopolymer.
    copolymer_composition : list[dict], optional
        List of ``{"monomer": name_str, "fraction": float}`` entries.
        Monomer names are arbitrary labels that map to SMILES via the
        ``chain_smiles`` (primary, name ``"M0"``) and additional entries.
        Provide ``smiles`` key alongside ``monomer`` for copolymer components,
        e.g. ``[{"monomer": "M0", "smiles": "[Si](C)(C)O", "fraction": 0.5},
                 {"monomer": "M1", "smiles": "[Si](C)(CCC(F)(F)F)O", "fraction": 0.5}]``.
    density : float
        Target packing density in g/cm³ (default 0.85).
    seed : int
        Random seed.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        Paths to the written files plus metadata.
    """
    import networkx as nx
    import numpy as np

    random.seed(seed)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Resolve solvent specification: mixture > single > default
    # ------------------------------------------------------------------
    from topon.simbox.molecule import Molecule

    if solvent_mixture:
        solvent_species = _resolve_mixture(
            solvent_mixture, dp, chain_smiles, density, verbose,
        )
    else:
        if solvent_smiles is None:
            solvent_smiles = "Cc1ccccc1"  # default: toluene
        solvent_species = _resolve_single_solvent(
            solvent_smiles, n_solvent, dp, density, verbose,
        )

    # Total number of solvent molecules (for reporting)
    total_n_solvent = sum(count for _, count, _ in solvent_species)

    if verbose:
        print("=== Single-Chain in Solvent ===")
        print(f"  Chain SMILES : {chain_smiles}  dp={dp}")
        for mol_obj, count, label in solvent_species:
            print(f"  Solvent      : {label}  n={count}")
        print(f"  Output       : {output_dir}")

    # ------------------------------------------------------------------
    # 1. Build minimal 2-node graph representing a single linear chain
    # ------------------------------------------------------------------
    G = nx.MultiGraph()
    G.add_node(0, pos=(0.0, 0.0, 0.0), node_type="end")
    G.add_node(1, pos=(float(dp), 0.0, 0.0), node_type="end")

    edge_attrs: dict = {"dp": dp, "edge_type": "A"}

    # Copolymer sequence
    if copolymer_arrangement and copolymer_composition:
        from topon.chemistry.sequences import generate_monomer_sequence
        seq_cfg = {
            "arrangement": copolymer_arrangement,
            "composition": [
                {"monomer": e["monomer"], "fraction": e["fraction"]}
                for e in copolymer_composition
            ],
        }
        edge_attrs["monomer_sequence"] = generate_monomer_sequence(
            dp, seq_cfg, default_monomer="M0"
        )

    # Graft positions
    if graft_density > 0 and graft_smiles:
        positions = [i for i in range(dp) if random.random() < graft_density]
        edge_attrs["graft_positions"] = positions
        edge_attrs["graft_dp"] = graft_dp
        edge_attrs["graft_monomer"] = "G"
        if verbose:
            print(f"  Grafts       : {len(positions)} side chains  (density={graft_density})")

    G.add_edge(0, 1, key=0, **edge_attrs)

    # ------------------------------------------------------------------
    # 2. Build chemistry config dynamically
    # ------------------------------------------------------------------
    from topon.config.schema import (
        ChemistryConfig, NodeMoleculeConfig, EdgeChemistryConfig,
        MonomerConfig, ConnectionConfig,
    )

    # Primary chain monomer (name "M0")
    head = _guess_head(chain_smiles)
    tail = _guess_tail(chain_smiles)
    monomers: dict[str, MonomerConfig] = {
        "M0": MonomerConfig(smiles=chain_smiles, chain_head=head, chain_tail=tail),
    }

    # Additional copolymer monomers
    if copolymer_composition:
        for entry in copolymer_composition:
            mon_name = entry["monomer"]
            mon_smiles = entry.get("smiles", chain_smiles)
            if mon_name not in monomers:
                monomers[mon_name] = MonomerConfig(
                    smiles=mon_smiles,
                    chain_head=_guess_head(mon_smiles),
                    chain_tail=_guess_tail(mon_smiles),
                )

    # Graft monomer
    if graft_density > 0 and graft_smiles:
        monomers["G"] = MonomerConfig(
            smiles=graft_smiles,
            chain_head=_guess_head(graft_smiles),
            chain_tail=_guess_tail(graft_smiles),
        )

    chem_config = ChemistryConfig(
        model_type="atomistic",
        target_density=density,
        node_type_map={
            "end": NodeMoleculeConfig(molecule="[Si](C)(C)C", is_end_cap=True),
        },
        edge_type_map={"A": EdgeChemistryConfig(monomer="M0")},
        monomers=monomers,
        connection=ConnectionConfig(auto_bridge=True, default_bridge_atom="O"),
    )

    # ------------------------------------------------------------------
    # 3. Build the chain molecule with ChemistryBuilder + embed 3D
    # ------------------------------------------------------------------
    if verbose:
        print("Building chain molecule...")

    from topon.chemistry.builder import ChemistryBuilder
    from rdkit import Chem

    builder = ChemistryBuilder(G, dims=None, config=chem_config)
    chain_rwmol = builder.build()

    # Sanitize and add explicit H
    chain_mol_raw = chain_rwmol.GetMol() if hasattr(chain_rwmol, "GetMol") else chain_rwmol
    try:
        Chem.SanitizeMol(chain_mol_raw)
    except Exception:
        pass
    chain_mol = Chem.AddHs(chain_mol_raw)

    # Assign 3D coordinates.  Default is the legacy extended-linear
    # placement for backward compatibility with existing callers.  New
    # code should pass ``embedder="etkdg"`` for any branched polymer
    # (PDMS, PTFPMS, polyacrylates) — the legacy placer collapses branch
    # atoms onto their parent heavy atom.
    if embedder == "etkdg":
        from topon.chemistry.embed import embed_with_etkdg
        embed_with_etkdg(chain_mol, seed=seed)
    elif embedder == "linear":
        _assign_extended_linear_coords(chain_mol, bond_length=1.5)
    else:
        raise ValueError(
            f"Unknown embedder {embedder!r}; expected 'linear' or 'etkdg'."
        )

    if verbose:
        print(f"  Chain atoms: {chain_mol.GetNumAtoms()}")

    # ------------------------------------------------------------------
    # 4. Compute box, place chain at center, pack solvent around it
    # ------------------------------------------------------------------
    if verbose:
        solv_total = sum(c for _, c, _ in solvent_species)
        print(f"Packing {solv_total} solvent molecules ({len(solvent_species)} species)...")

    from topon.simbox.packer import BoxPacker, PackedBox, PlacedMolecule

    chain_mol_obj = Molecule.from_mol("chain", chain_mol)

    # Build the full species list (chain + solvents) for box sizing,
    # but we will NOT let the packer randomly place the chain.
    all_species: list[tuple[Molecule, int]] = [(chain_mol_obj, 1)]
    solvent_only_species: list[tuple[Molecule, int]] = []
    for mol_obj, count, _label in solvent_species:
        all_species.append((mol_obj, count))
        solvent_only_species.append((mol_obj, count))

    # --- Compute box size from total mass (same formula as BoxPacker) ---
    AVOGADRO = 6.02214076e23
    total_mass_g = sum(mol.mw * count for mol, count in all_species) / AVOGADRO
    volume_cm3 = total_mass_g / density
    volume_A3 = volume_cm3 * 1e24
    L = volume_A3 ** (1.0 / 3.0)
    box_lengths = np.array([L, L, L])

    # --- Place chain at box center (extended along x-axis) ---
    chain_coords = chain_mol_obj.get_coordinates()
    chain_centroid = chain_coords.mean(axis=0)
    box_center = box_lengths / 2.0
    shift = box_center - chain_centroid
    chain_coords_centered = chain_coords + shift

    chain_placement = PlacedMolecule(
        species_idx=0,
        instance_idx=0,
        molecule=chain_mol_obj,
        coordinates=chain_coords_centered,
    )

    if verbose:
        chain_span = np.ptp(chain_coords_centered[:, 0])  # extent along x
        print(f"  Chain extent along x: {chain_span:.1f} Å, box edge: {L:.1f} Å")
        print(f"  Chain centered at ({box_center[0]:.1f}, {box_center[1]:.1f}, {box_center[2]:.1f})")

    # --- Pack solvent around the pre-placed chain ---
    packer = BoxPacker(density=density, seed=seed)
    # Use the packer for solvent only, but pre-register chain atoms
    packed_box_solvent = packer.pack(solvent_only_species,
                                     pre_placed=[chain_placement],
                                     box_lengths_override=box_lengths)

    # The packed_box_solvent already includes the chain as the first placement
    packed_box = packed_box_solvent

    # Re-center the chain in case the BoxPacker scaled up the box lengths
    _center_chain_in_box(packed_box)

    if verbose:
        box = packed_box.box_lengths
        print(f"  Box length: {box[0]:.2f} Å  ({packed_box.total_atoms} atoms total)")
        print(f"  Chain centered at box midpoint.")


    # ------------------------------------------------------------------
    # 5. Assemble and write DREIDING LAMMPS files
    # ------------------------------------------------------------------
    if verbose:
        print("Writing LAMMPS files...")

    from topon.simbox.system import assemble
    from topon.simbox.writer import write_lammps
    from topon.simbox.inputs import write_inputs

    system = assemble(packed_box)

    written_data = write_lammps(system, str(output_path))
    written_inputs = write_inputs(system, str(output_path), include_crosslink=False)

    result = {
        "output_dir": str(output_path),
        "chain_atoms": chain_mol.GetNumAtoms(),
        "n_solvent": total_n_solvent,
        "solvent_species": [
            {"smiles": label, "count": count}
            for _, count, label in solvent_species
        ],
        "box_length_ang": float(packed_box.box_lengths[0]),
        **written_data,
        **written_inputs,
    }

    if verbose:
        print(f"=== Done. Files written to {output_dir} ===")

    return result


def _assign_extended_linear_coords(mol, bond_length: float = 1.5) -> None:
    """Assign 3D coordinates as an extended polymer chain.

    .. deprecated:: 0.2.0
       Leaves branch atoms collapsed onto their parent heavy atom — the
       pendant-placement step only displaces pendants *randomly* from
       the backbone, which commonly produces minimum pair-distances
       well below 1 Å (e.g. 0.44 Å Si-methyl-C on top of backbone O for
       PDMS DP=30).  Use :func:`topon.chemistry.embed.embed_with_etkdg`
       instead for any polymer with side chains or stereochemistry.

    Strategy (simple geometric placement):

    1. **Find the backbone** — longest path through the heavy-atom-only
       sub-graph (e.g. Si-O-Si-O for PDMS, C-C-C-C for PE).
    2. **Place backbone atoms along x-axis** — consecutive backbone atoms
       are separated by their standard bond distance along x.
    3. **Place pendant atoms** — every non-backbone atom is placed at its
       parent atom's position, then randomly displaced outward to the
       correct bond distance.

    This produces a fully extended chain where the backbone is linear along
    x and side groups / hydrogens radiate outward.  Bond angles are
    approximate; LAMMPS Stage 1 (soft push-off + minimization) corrects
    them to their harmonic equilibrium values.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule with explicit H (from ``Chem.AddHs``).
    bond_length : float
        Default bond length (Å) when no specific value is available.
    """
    import warnings
    warnings.warn(
        "_assign_extended_linear_coords places branch atoms on top of their "
        "parent (min pair distance can be < 0.5 Å for branched polymers). "
        "Use topon.chemistry.embed.embed_with_etkdg(mol, seed=...) instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    import numpy as np
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    n = mol.GetNumAtoms()
    rng = np.random.default_rng(42)

    # --- Step 1: Find backbone path ---
    backbone = _find_backbone_path(mol)
    backbone_set = set(backbone)

    # --- Step 2: Place backbone atoms along x-axis ---
    coords = np.zeros((n, 3))
    x_pos = 0.0

    for i, bb_idx in enumerate(backbone):
        if i == 0:
            coords[bb_idx] = [0.0, 0.0, 0.0]
        else:
            prev_idx = backbone[i - 1]
            d = _get_bond_length(mol, prev_idx, bb_idx, default=bond_length)
            x_pos += d
            coords[bb_idx] = [x_pos, 0.0, 0.0]

    # --- Step 3: Place pendant (non-backbone) atoms ---
    # BFS outward from backbone atoms through the molecular graph.
    # Each pendant atom is placed at its parent's position, then
    # displaced by a random unit vector × bond_distance.
    placed = set(backbone)
    queue: list[tuple[int, int]] = []  # (atom_idx, parent_idx)

    # Seed the queue with non-backbone neighbors of backbone atoms
    for bb_idx in backbone:
        atom = mol.GetAtomWithIdx(bb_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in placed:
                queue.append((nbr_idx, bb_idx))

    while queue:
        next_queue: list[tuple[int, int]] = []
        for atom_idx, parent_idx in queue:
            if atom_idx in placed:
                continue
            d = _get_bond_length(mol, parent_idx, atom_idx, default=bond_length)
            # Random direction (uniform on sphere)
            direction = rng.standard_normal(3)
            direction /= np.linalg.norm(direction)
            coords[atom_idx] = coords[parent_idx] + direction * d
            placed.add(atom_idx)

            # Enqueue this atom's unplaced neighbors
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in placed:
                    next_queue.append((nbr.GetIdx(), atom_idx))
        queue = next_queue

    # --- Center at origin ---
    centroid = coords.mean(axis=0)
    coords -= centroid

    # --- Write conformer ---
    conf = Chem.Conformer(n)
    for i in range(n):
        conf.SetAtomPosition(i, Point3D(
            float(coords[i, 0]), float(coords[i, 1]), float(coords[i, 2]),
        ))
    mol.AddConformer(conf, assignId=True)


# ---------------------------------------------------------------------------
# Standard bond lengths (Å) for common element pairs
# ---------------------------------------------------------------------------
_BOND_LENGTHS: dict[tuple[int, int, str], float] = {
    # (atomic_num_1, atomic_num_2, bond_type) → length in Å
    # C-C
    (6, 6, "SINGLE"):   1.54,
    (6, 6, "DOUBLE"):   1.34,
    (6, 6, "AROMATIC"): 1.40,
    (6, 6, "TRIPLE"):   1.20,
    # C-H
    (6, 1, "SINGLE"):   1.09,
    (1, 6, "SINGLE"):   1.09,
    # C-O
    (6, 8, "SINGLE"):   1.43,
    (6, 8, "DOUBLE"):   1.23,
    # C-N
    (6, 7, "SINGLE"):   1.47,
    (6, 7, "DOUBLE"):   1.29,
    (6, 7, "TRIPLE"):   1.16,
    # C-F
    (6, 9, "SINGLE"):   1.35,
    # C-Cl
    (6, 17, "SINGLE"):  1.77,
    # Si-O
    (14, 8, "SINGLE"):  1.65,
    (8, 14, "SINGLE"):  1.65,
    # Si-C
    (14, 6, "SINGLE"):  1.87,
    (6, 14, "SINGLE"):  1.87,
    # Si-H
    (14, 1, "SINGLE"):  1.48,
    (1, 14, "SINGLE"):  1.48,
    # O-H
    (8, 1, "SINGLE"):   0.96,
    (1, 8, "SINGLE"):   0.96,
    # N-H
    (7, 1, "SINGLE"):   1.01,
    (1, 7, "SINGLE"):   1.01,
}


def _get_bond_length(mol, idx1: int, idx2: int, default: float = 1.5) -> float:
    """Look up the standard bond length between two bonded atoms."""
    bond = mol.GetBondBetweenAtoms(idx1, idx2)
    if bond is None:
        return default
    z1 = mol.GetAtomWithIdx(idx1).GetAtomicNum()
    z2 = mol.GetAtomWithIdx(idx2).GetAtomicNum()
    btype = str(bond.GetBondType()).split(".")[-1]  # e.g. "SINGLE"
    key = (z1, z2, btype)
    if key in _BOND_LENGTHS:
        return _BOND_LENGTHS[key]
    key_rev = (z2, z1, btype)
    if key_rev in _BOND_LENGTHS:
        return _BOND_LENGTHS[key_rev]
    return default





def _find_backbone_path(mol) -> list[int]:
    """Find the backbone of a linear polymer as the longest path through
    the heavy-atom-only sub-graph.

    Returns a list of atom indices representing the backbone from one chain
    end to the other.
    """
    import networkx as nx

    # Build heavy-atom-only graph
    G = nx.Graph()
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if (mol.GetAtomWithIdx(a1).GetAtomicNum() > 1 and
                mol.GetAtomWithIdx(a2).GetAtomicNum() > 1):
            G.add_edge(a1, a2)

    if len(G) == 0:
        return list(range(mol.GetNumAtoms()))

    # Find leaf nodes (degree 1 in the heavy-atom graph = chain ends)
    leaves = [n for n in G.nodes() if G.degree(n) == 1]
    if len(leaves) < 2:
        leaves = list(G.nodes())

    # Find the pair of leaves with the longest shortest path (= graph diameter)
    best_path: list[int] = []
    for i, s in enumerate(leaves):
        for t in leaves[i + 1:]:
            try:
                path = nx.shortest_path(G, s, t)
                if len(path) > len(best_path):
                    best_path = path
            except nx.NetworkXNoPath:
                pass

    return best_path if best_path else list(G.nodes())


def _center_and_align(mol, backbone: list[int]) -> None:
    """Rotate so the backbone's principal axis aligns with x, then
    center the molecule at the origin.
    """
    import numpy as np
    from rdkit.Geometry import Point3D

    conf = mol.GetConformer()
    n = mol.GetNumAtoms()

    # Get all coordinates
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(n)])

    # Get backbone coordinates
    bb_coords = coords[backbone]

    # Principal axis of backbone (direction of greatest variance)
    bb_centered = bb_coords - bb_coords.mean(axis=0)
    if len(bb_centered) > 1:
        # SVD to find principal axis
        _, _, Vt = np.linalg.svd(bb_centered, full_matrices=False)
        principal_axis = Vt[0]  # first right singular vector
    else:
        principal_axis = np.array([1.0, 0.0, 0.0])

    # Build rotation matrix to align principal_axis with x-axis [1,0,0]
    R = _rotation_matrix_align(principal_axis, np.array([1.0, 0.0, 0.0]))

    # Apply rotation to all atoms
    centroid = coords.mean(axis=0)
    coords_centered = coords - centroid
    coords_rotated = (R @ coords_centered.T).T

    # Set new coordinates
    for i in range(n):
        conf.SetAtomPosition(i, Point3D(
            float(coords_rotated[i, 0]),
            float(coords_rotated[i, 1]),
            float(coords_rotated[i, 2]),
        ))


def _rotation_matrix_align(v_from: 'np.ndarray', v_to: 'np.ndarray') -> 'np.ndarray':
    """Compute the 3×3 rotation matrix that rotates vector *v_from* onto *v_to*.

    Uses Rodrigues' rotation formula.
    """
    import numpy as np

    a = v_from / np.linalg.norm(v_from)
    b = v_to / np.linalg.norm(v_to)

    cross = np.cross(a, b)
    dot = np.dot(a, b)
    sin_angle = np.linalg.norm(cross)

    if sin_angle < 1e-10:
        # Vectors are (anti-)parallel
        if dot > 0:
            return np.eye(3)
        else:
            # 180° rotation about any perpendicular axis
            perp = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
            perp = perp - np.dot(perp, a) * a
            perp /= np.linalg.norm(perp)
            return 2 * np.outer(perp, perp) - np.eye(3)

    # Skew-symmetric cross-product matrix
    K = np.array([
        [0, -cross[2], cross[1]],
        [cross[2], 0, -cross[0]],
        [-cross[1], cross[0], 0],
    ])
    R = np.eye(3) + K + K @ K * ((1 - dot) / (sin_angle ** 2))
    return R


def _assign_fallback_coords(mol, bond_length: float = 1.5) -> None:
    """Fallback coordinate assignment when ETKDGv3 fails.

    Places backbone atoms in an extended zigzag (tetrahedral angle) along the
    x-z plane, and pendant atoms at random angular offsets perpendicular to
    the local backbone direction.
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    n = mol.GetNumAtoms()
    rng = np.random.default_rng(42)
    conf = Chem.Conformer(n)

    # Simple zigzag: place atoms along x with tetrahedral zig-zag in z
    theta = np.radians(109.47 / 2.0)  # half tetrahedral angle
    for i in range(n):
        x = i * bond_length * np.cos(theta)
        z = bond_length * np.sin(theta) * (1 if i % 2 == 0 else -1)
        y = rng.uniform(-0.3, 0.3)  # small random y displacement
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

    mol.AddConformer(conf, assignId=True)


def _guess_head(smiles: str) -> str:
    """Guess the chain-head linking atom from a repeat-unit SMILES.

    Heuristic precedence:
    1. First bracketed heavy atom → e.g. ``[Si]`` → ``"Si"``
    2. First non-hydrogen heavy atom in the SMILES → ``"C"``, ``"N"``, etc.
    3. Fallback → ``"C"``

    This properly handles hydrocarbon monomers like ``CC``, ``CC(C)``,
    fluorinated monomers like ``C(F)(F)C``, and siloxanes like ``[Si](C)(C)O``.
    """
    import re
    # Try bracketed atom first (e.g. [Si], [N])
    m = re.search(r"\[([A-Z][a-z]?)", smiles)
    if m:
        return m.group(1)
    # Try first heavy element letter
    m = re.search(r"[A-Z]", smiles)
    if m:
        elem = m.group(0)
        # Check for two-letter elements at this position
        pos = m.start()
        if pos + 1 < len(smiles) and smiles[pos + 1].islower():
            elem = smiles[pos:pos + 2]
        return elem
    return "C"


def _guess_tail(smiles: str) -> str:
    """Guess the chain-tail linking atom from a repeat-unit SMILES.

    Heuristic:
    - If SMILES ends with ``O`` (e.g. siloxane ``[Si](C)(C)O``) → ``"O"``
    - If SMILES ends with a parenthesised group (e.g. ``CC(C)`` for propylene,
      ``CC(C#N)`` for acrylonitrile), the tail is the *backbone* atom before
      the parenthesis, typically ``"C"``.
    - Otherwise → same as head.
    """
    import re
    stripped = smiles.rstrip()
    # Siloxane-type ending: ...O
    if stripped.endswith("O") and not stripped.endswith(")"):
        return "O"
    # Ends with parenthesised side-group → backbone tail is the preceding atom
    if stripped.endswith(")"):
        # Walk backward to find the matching open paren
        depth = 0
        for i in range(len(stripped) - 1, -1, -1):
            if stripped[i] == ")":
                depth += 1
            elif stripped[i] == "(":
                depth -= 1
            if depth == 0:
                # The atom *before* this open paren is the backbone tail
                pre = stripped[:i]
                return _guess_head(pre[-2:]) if len(pre) >= 2 else _guess_head(pre)
    return _guess_head(smiles)


# ---------------------------------------------------------------------------
# Multi-solvent helpers
# ---------------------------------------------------------------------------

def _resolve_single_solvent(
    solvent_smiles: str,
    n_solvent: Optional[int],
    dp: int,
    density: float,
    verbose: bool,
) -> list[tuple]:
    """Resolve a single-solvent specification into a list of
    ``(Molecule, count, label)`` tuples.
    """
    from topon.simbox.molecule import Molecule

    mol_obj = Molecule.from_smiles("solvent", solvent_smiles)

    if n_solvent is None:
        # Estimate contour length including end-caps (approx 2 extra monomers)
        contour_ang = (dp + 2) * 2 * 1.65
        # Ensure box is at least 1.5x the chain length
        target_box = max(contour_ang * 1.5, 50.0)
        vol_A3 = target_box ** 3
        total_mass_g = density * vol_A3 * 1e-24
        n_solvent = max(1, int(total_mass_g * 6.022e23 / mol_obj.mw))
        if verbose:
            print(f"  Auto n_solvent: {n_solvent}  "
                  f"(target box {target_box:.0f} Å > contour {contour_ang:.0f} Å)")

    return [(mol_obj, n_solvent, solvent_smiles)]


def _resolve_mixture(
    solvent_mixture: list[dict],
    dp: int,
    chain_smiles: str,
    density: float,
    verbose: bool,
) -> list[tuple]:
    """Resolve a multi-solvent mixture into a list of
    ``(Molecule, count, label)`` tuples.

    Molecule counts are determined by:
    1. Computing the target box volume from the chain contour length.
    2. Allocating total solvent mass to each species by weight fraction.
    3. Converting mass → molecule count via molecular weight.
    """
    from topon.simbox.molecule import Molecule

    # Normalise weight fractions
    total_wf = sum(e["weight_fraction"] for e in solvent_mixture)
    if total_wf <= 0:
        raise ValueError("solvent_mixture weight fractions must sum to > 0")

    # Build Molecule objects and get MW
    mol_objects = []
    for entry in solvent_mixture:
        smi = entry["smiles"]
        mol_obj = Molecule.from_smiles(f"solv_{smi[:8]}", smi)
        mol_objects.append((mol_obj, entry["weight_fraction"] / total_wf, smi))

    # Target box size: including end-caps, and ensure >= 1.5x chain length
    contour_ang = (dp + 2) * 2 * 1.65
    target_box = max(contour_ang * 1.5, 50.0)
    vol_A3 = target_box ** 3
    total_mass_g = density * vol_A3 * 1e-24   # total mass of everything in box

    # Approximate: chain mass is small relative to solvent; allocate
    # ~95% of total mass to solvent (conservative estimate).
    solvent_mass_g = total_mass_g * 0.95

    result = []
    for mol_obj, wf_norm, smi in mol_objects:
        species_mass_g = solvent_mass_g * wf_norm
        n_mol = max(1, int(species_mass_g * 6.022e23 / mol_obj.mw))
        result.append((mol_obj, n_mol, smi))

    if verbose:
        print(f"  Auto mixture sizing (target box {target_box:.0f} Å):")
        for mol_obj, count, smi in result:
            print(f"    {smi}: {count} molecules (MW={mol_obj.mw:.1f})")

    return result


def _center_chain_in_box(packed_box) -> None:
    """Shift all coordinates so the chain's geometric center sits at the
    box midpoint.

    The chain is always the first placement (species_idx=0, instance_idx=0).
    All molecules (chain + solvent) are shifted by the same vector so their
    relative positions are preserved.
    """
    import numpy as np

    chain_placement = packed_box.placements[0]
    chain_centroid = chain_placement.coordinates.mean(axis=0)
    box_center = packed_box.box_lengths / 2.0
    shift = box_center - chain_centroid

    for pm in packed_box.placements:
        pm.coordinates = pm.coordinates + shift

