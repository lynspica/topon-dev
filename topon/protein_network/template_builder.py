"""Build a multi-chain System by replicating an ITP `ChainTemplate` per chain.

The builder offsets atom IDs by `chain_idx * template.n_atoms` so each chain's
bonded terms reference the right global atoms. No information loss vs polyply
-- every angle, dihedral, improper, exclusion is replicated exactly.
"""
from __future__ import annotations

import numpy as np

from . import bfm, sequence
from .itp_template import ChainTemplate
from .martini_ff import MartiniLibrary
from .system import (
    Angle, Bead, Bond, Constraint, Dihedral, Exclusion, System,
)
from .builder import (
    DEFAULT_BB_BB_TARGET_ANG, DEFAULT_SC_JITTER_ANG,
    DITYROSINE_BOND_LENGTH_NM, DITYROSINE_BOND_K_KJ,
    _auto_lattice_scale_ang, _interpolate_residue_positions,
)


def build_protein_system_from_template(
    snapshot: dict,
    template: ChainTemplate,
    library: MartiniLibrary,
    *,
    block_seq: str | None = None,
    crosslinker_letter: str = "Y",
    lattice_scale_ang: float | None = None,
    sc_jitter_ang: float = DEFAULT_SC_JITTER_ANG,
    seed: int = 42,
) -> System:
    """Build a System by replicating the ITP `template` per chain in the snapshot.

    Each chain gets its own copy of the template's atoms / bonds / angles /
    dihedrals / impropers / exclusions, with atom IDs offset by
    `chain_idx * template.n_atoms`. Coordinates come from min-image
    interpolation between BFM lattice nodes (same as `build_protein_system`).
    Dityrosine crosslinks from `snapshot["reactions"]` are added on top.
    Wrap-only writer convention: all crosslinks are kept; LAMMPS computes
    bond forces using min-image distance via the neighbor/ghost system.
    """
    rng = np.random.default_rng(seed)
    chains_flat: list[list[int]] = snapshot["chains"]
    Nx, Ny, Nz = snapshot["Nx"], snapshot["Ny"], snapshot["Nz"]
    n_chains = len(chains_flat)
    n_residues = template.n_residues
    n_nodes = len(chains_flat[0])
    n_segments_per_chain = n_nodes - 1

    if lattice_scale_ang is None:
        lattice_scale_ang = _auto_lattice_scale_ang(n_residues, n_segments_per_chain)
    box = np.array([Nx, Ny, Nz], dtype=float) * lattice_scale_ang

    if block_seq is None:
        block_seq, _ = template.block_sequence()
    node_to_res = sequence.get_node_residue_mapping(
        n_repeats=n_residues // len(block_seq),
        segs_per_block=n_segments_per_chain // (n_residues // len(block_seq)),
        block_seq=block_seq,
        crosslinker_letter=crosslinker_letter,
    )

    sys_ = System(box_dims_ang=tuple(box.tolist()))
    # global atom ID assigned sequentially across chains (1-based).
    next_global_id = 1
    sc4_atom_by_chain_residue: dict[tuple[int, int], int] = {}

    for ci, chain_flat in enumerate(chains_flat):
        molecule_id = ci + 1
        residue_positions = _interpolate_residue_positions(
            chain_flat, node_to_res, Nx, Ny, lattice_scale_ang, box
        )

        # Replicate atoms: for each ITP atom, place it at its residue's anchor
        # (BB at anchor, sidechains with deterministic jitter from BB).
        bb_pos_by_resnr: dict[int, np.ndarray] = {}
        local_to_global: dict[int, int] = {}
        for atom in template.atoms:
            if atom.atom_name == "BB":
                anchor = residue_positions.get(atom.resnr - 1, np.zeros(3))  # ITP resnr is 1-based
                pos = anchor.copy()
                bb_pos_by_resnr[atom.resnr] = pos
            else:
                bb = bb_pos_by_resnr.get(atom.resnr)
                if bb is None:
                    bb = residue_positions.get(atom.resnr - 1, np.zeros(3))
                jitter = rng.normal(0.0, sc_jitter_ang, size=3)
                pos = bb + jitter
            global_id = next_global_id
            local_to_global[atom.id] = global_id
            next_global_id += 1
            mass = library.atomtypes[atom.bead_type].mass if atom.bead_type in library.atomtypes else 72.0
            sys_.beads.append(Bead(
                atom_id=global_id, bead_type=atom.bead_type,
                molecule_id=molecule_id, residue_idx=atom.resnr,
                residue_name=atom.resname, atom_name=atom.atom_name,
                charge=atom.charge, mass=mass,
                position=tuple(pos.tolist()),
            ))
            if atom.atom_name == "SC4":
                sc4_atom_by_chain_residue[(ci, atom.resnr - 1)] = global_id

        # Replicate bonded terms (atom IDs offset via local_to_global).
        for b in template.bonds:
            sys_.bonds.append(Bond(
                a=local_to_global[b.i], b=local_to_global[b.j],
                funct=b.funct, length_nm=b.length_nm, k_kj=b.k_kj,
            ))
        for c in template.constraints:
            sys_.constraints.append(Constraint(
                a=local_to_global[c.i], b=local_to_global[c.j],
                length_nm=c.length_nm,
            ))
        for ang in template.angles:
            sys_.angles.append(Angle(
                a=local_to_global[ang.i], b=local_to_global[ang.j], c=local_to_global[ang.k],
                funct=ang.funct, angle_deg=ang.angle_deg, k_kj=ang.k_kj,
            ))
        for d in template.dihedrals_proper:
            sys_.dihedrals.append(Dihedral(
                a=local_to_global[d.i], b=local_to_global[d.j],
                c=local_to_global[d.k], d=local_to_global[d.l],
                funct=d.funct, angle_deg=d.angle_deg, k_kj=d.k_kj,
                mult=d.mult, is_improper=False,
            ))
        for d in template.dihedrals_improper:
            sys_.dihedrals.append(Dihedral(
                a=local_to_global[d.i], b=local_to_global[d.j],
                c=local_to_global[d.k], d=local_to_global[d.l],
                funct=d.funct, angle_deg=d.angle_deg, k_kj=d.k_kj,
                mult=1, is_improper=True,
            ))
        for excl in template.exclusions:
            ids = tuple(local_to_global[a] for a in excl.atoms if a in local_to_global)
            if len(ids) >= 2:
                sys_.exclusions.append(Exclusion(atoms=ids))

    # Collect dityrosine crosslinks (chain-pair, atom-pair) from BFM snapshot.
    crosslinks_pending: list[tuple[int, int, int, int]] = []  # (ci1, ci2, atom_a, atom_b)
    for rxn in snapshot.get("reactions", []):
        (ci1, ni1), (ci2, ni2) = rxn[0], rxn[1]
        r1 = node_to_res.get(ni1); r2 = node_to_res.get(ni2)
        if r1 is None or r2 is None: continue
        a = sc4_atom_by_chain_residue.get((ci1, r1))
        b = sc4_atom_by_chain_residue.get((ci2, r2))
        if a is None or b is None: continue
        crosslinks_pending.append((ci1, ci2, a, b))

    # Wrap-only convention (matches core topon): atoms are wrapped into
    # [0, box) at write time, no image flag column emitted. LAMMPS computes
    # bond forces via the neighbor/ghost system using min-image distance, so
    # a crosslink between two beads at the same lattice site is short
    # regardless of which chain wrap-counts placed them there. No per-chain
    # image shifts and no cycle-loser drops are needed.
    for (ci1, ci2, a, b) in crosslinks_pending:
        sys_.bonds.append(Bond(
            a=a, b=b, funct=1,
            length_nm=DITYROSINE_BOND_LENGTH_NM, k_kj=DITYROSINE_BOND_K_KJ,
            is_crosslink=True,
        ))

    return sys_
