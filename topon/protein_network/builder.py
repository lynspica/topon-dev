"""Build a MARTINI 3 protein system from a BFM snapshot + residue sequence.

Walks one BFM-lattice topology snapshot and a per-chain residue sequence and
emits a fully populated `System`: bead positions, intra- and inter-residue
bonded terms, constraints, exclusions. Sidechain beads are placed near their
backbone bead with deterministic jitter (seeded). Inter-chain dityrosine
crosslinks (SC4-SC4 bonds at TYR residues whose BFM lattice sites were merged)
are applied here, since the merger is recorded directly in the snapshot's
`reactions` list and the `chains` list (post-merger) shows the shared flat
index between two TYR positions.
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

from . import bfm, residues, sequence
from .martini_ff import MartiniLibrary
from .system import (
    Angle,
    Bead,
    Bond,
    Constraint,
    Dihedral,
    Exclusion,
    System,
)
import numpy as np  # noqa: E402  (re-import for type hint clarity)

# Default placement parameters.
# `lattice_scale_ang=None` triggers auto-scaling to give realistic BB-BB
# distances (MARTINI 3 default 0.36 nm = 3.6 A). Manual override is allowed.
DEFAULT_BB_BB_TARGET_ANG: float = 3.6     # MARTINI 3 backbone bond length, 0.36 nm
DEFAULT_SC_JITTER_ANG: float = 0.3        # sidechain bead offset (small, so initial constraint
                                           # bond lengths sit close to equilibrium ~ 2.7-3.0 A)
DITYROSINE_BOND_LENGTH_NM: float = 0.270  # SC4-SC4 dityrosine bond
DITYROSINE_BOND_K_KJ: float = 8000.0      # treat like a stiff backbone bond
CROSSLINK_SEP_ANG: float = 7.0            # separation of the two merged TYR residues at a
                                           # crosslink, so the rings form an adjacent dimer
                                           # instead of landing coincident (r~0). The SC4-SC4
                                           # bond relaxes them to equilibrium during minimisation.


def _auto_lattice_scale_ang(n_residues: int, n_segments_per_chain: int) -> float:
    """Auto-pick lattice_scale_ang so BB-BB equilibrium length is preserved.

    The BFM lattice has `n_segments_per_chain` segments per chain, but each
    segment may span many residues. Set scale so segment-length / residues-
    per-segment >= MARTINI BB-BB equilibrium (3.6 A).
    """
    if n_segments_per_chain <= 0:
        return DEFAULT_BB_BB_TARGET_ANG
    rps = max(1.0, (n_residues - 1) / n_segments_per_chain)
    return max(DEFAULT_BB_BB_TARGET_ANG, rps * DEFAULT_BB_BB_TARGET_ANG)


def _lattice_pos_ang(flat_idx: int, Nx: int, Ny: int, scale: float) -> np.ndarray:
    """Literal Cartesian position of a BFM lattice site (always inside [0, box))."""
    x, y, z = bfm.flat_to_xyz(flat_idx, Nx, Ny)
    return np.array([x * scale, y * scale, z * scale], dtype=float)


def _interpolate_residue_positions(
    chain_flat: list[int],
    node_to_res: dict[int, int],
    Nx: int,
    Ny: int,
    scale: float,
    box: np.ndarray,
) -> dict[int, np.ndarray]:
    """Walk a chain residue-by-residue; positions accumulate via min-image steps.

    Each chain starts at its first BFM node's literal `lattice_pos` and walks
    forward, using minimum-image displacement between consecutive anchors.
    Consecutive residues end up within ~scale/N_per_segment Angstroms of each
    other (=~ MARTINI bond length). Resulting unwrapped positions may fall
    outside [0, box); the writer assigns LAMMPS image flags via `floor(pos/box)`.

    This matches topro's `protein/builder.py:149-172` shape. For widely-spaced
    MARTINI lattices, two chains crosslinked at the same merged BFM node may
    end up far apart in their unwrapped frames (their walks reach the merge
    site from different sides of the box). The few resulting oversized
    crosslink bonds are dropped post-hoc by the caller; intra-chain bonds
    remain short via the min-image walk.
    """
    n_chain_nodes = len(chain_flat)
    node_positions = {ni: _lattice_pos_ang(fi, Nx, Ny, scale)
                       for ni, fi in enumerate(chain_flat)}
    residue_positions: dict[int, np.ndarray] = {}
    sorted_nodes = sorted(node_to_res.keys())
    for seg_i in range(len(sorted_nodes) - 1):
        ni_start = sorted_nodes[seg_i]
        ni_end = sorted_nodes[seg_i + 1]
        if ni_start >= n_chain_nodes or ni_end >= n_chain_nodes:
            continue
        r_start = node_to_res[ni_start]
        r_end = node_to_res[ni_end]
        # Walk from previous segment's endpoint (if any) so positions accumulate
        # consistently. For the very first segment, start at the literal lattice_pos.
        if r_start in residue_positions:
            p_start = residue_positions[r_start]
        else:
            p_start = node_positions[ni_start]
            residue_positions[r_start] = p_start
        p_end_lattice = node_positions[ni_end]
        # min-image step from p_start (which may be outside box) to p_end_lattice (in box):
        # bring p_end into the same image cell as p_start before computing diff.
        diff = p_end_lattice - p_start
        diff -= box * np.round(diff / box)
        n_seg_res = r_end - r_start
        for k in range(1, n_seg_res + 1):
            r = r_start + k
            if r not in residue_positions:
                frac = k / n_seg_res
                residue_positions[r] = p_start + frac * diff
    return residue_positions


def _build_inverse_mapping(node_to_res: dict[int, int]) -> list[tuple[int, int, int]]:
    """For every residue r, find (node_left, node_right, r_anchor_left, r_anchor_right)
    such that node_to_res[node_left] <= r <= node_to_res[node_right] and node_right
    is the smallest such anchor.

    Returns a list indexed by residue idx with (node_left, node_right, r_left, r_right).
    """
    sorted_nodes = sorted(node_to_res.keys())
    sorted_res = [node_to_res[n] for n in sorted_nodes]
    n_residues = sorted_res[-1] + 1
    out: list[tuple[int, int, int, int]] = []
    for r in range(n_residues):
        # binary search: find first anchor with res >= r
        right_i = next((i for i, rv in enumerate(sorted_res) if rv >= r), len(sorted_res) - 1)
        if right_i == 0:
            left_i = 0
            right_i = min(1, len(sorted_res) - 1)
        else:
            left_i = right_i - 1 if sorted_res[right_i] != r else right_i
            if left_i == right_i and right_i + 1 < len(sorted_nodes):
                right_i = right_i + 1
        out.append((sorted_nodes[left_i], sorted_nodes[right_i], sorted_res[left_i], sorted_res[right_i]))
    return out




def _patch_terminal_beads(
    raw_beads: list[tuple[str, str, float]],
    is_n_term: bool,
    is_c_term: bool,
) -> list[tuple[str, str, float]]:
    """Apply terminal patches by looking up the per-atom override in TERMINAL_PATCHES."""
    if not is_n_term and not is_c_term:
        return raw_beads
    overrides: dict[str, tuple[str, float]] = {}
    if is_n_term:
        for (name, bt, q) in residues.TERMINAL_PATCHES.get("N_term", []):
            overrides[name] = (bt, q)
    if is_c_term:
        for (name, bt, q) in residues.TERMINAL_PATCHES.get("C_term", []):
            overrides[name] = (bt, q)
    return [
        (name, overrides.get(name, (bt, q))[0], overrides.get(name, (bt, q))[1])
        for (name, bt, q) in raw_beads
    ]


def _select_bb_bond(prev_resname: str, this_resname: str) -> tuple[int, float, float]:
    """Look up the BB-BB bond record for a (prev, this) residue pair in the
    extracted backbone table. Falls back to the default (8000, 0.360) if absent."""
    for (a, b, funct, length, k) in residues.BACKBONE_BB_BONDS:
        if a == prev_resname and b == this_resname:
            return funct, length, k
    return 1, 0.360, 8000.0


def _select_bbbb_label(quad_resnames: tuple[str, str, str, str]) -> str:
    """Label for a BBBB dihedral based on which residues in the quadruplet are GLY.

    polyply emits 5 labels: BBBB (default), GGGX, GGXG, GXGG, XGGG depending on
    the position of GLY (G) versus a non-GLY (X) at each of the 4 backbone slots.
    """
    pattern = "".join("G" if r == "GLY" else "X" for r in quad_resnames)
    if pattern == "GGGG":
        # 4 GLYs in a row: prefer the GGGX bucket as a generic Gly-rich variant
        return "GGGX"
    if pattern in residues.BACKBONE_BBBB_DIHEDRALS:
        return pattern
    return "BBBB"


def _embed_residue_local(
    n_beads: int,
    edges: list[tuple[int, int, float]],
    seed: int,
    min_sep: float = 3.0,
    lr: float = 0.08,
    iters: int = 250,
) -> np.ndarray:
    """Embed one residue's beads in 3D from its intra bond/constraint edges.

    Bead 0 (BB) is pinned at the origin. ``edges`` are ``(i, j, target_len_ang)``
    local-index pairs (bonds AND ring constraints). Minimises sum of squared
    length errors plus a soft repulsion that keeps every NON-bonded pair at least
    ``min_sep`` Angstrom apart. This replaces the old ``bb + jitter`` placement,
    which piled all sidechain beads onto the backbone (down to ~0.02 nm apart) and
    was only tolerable in LAMMPS because ``special_bonds 0 0 0`` over-excluded the
    resulting intra-residue overlaps. Returns ``(n_beads, 3)`` offsets from BB.
    """
    if n_beads <= 1:
        return np.zeros((1, 3))
    rng = np.random.default_rng(seed)
    pos = np.zeros((n_beads, 3))
    for k in range(1, n_beads):
        pos[k] = np.array([0.0, 0.0, 3.0 * k]) + rng.normal(0.0, 0.3, 3)
    bonded = {(min(i, j), max(i, j)) for (i, j, _l) in edges}
    for _ in range(iters):
        grad = np.zeros((n_beads, 3))
        for (i, j, L) in edges:
            d = pos[i] - pos[j]
            r = float(np.linalg.norm(d)) + 1e-9
            g = (r - L) / r * d
            grad[i] += g
            grad[j] -= g
        for i in range(n_beads):
            for j in range(i + 1, n_beads):
                if (i, j) in bonded:
                    continue
                d = pos[i] - pos[j]
                r = float(np.linalg.norm(d)) + 1e-9
                if r < min_sep:
                    g = (min_sep - r) / r * d  # push i away from j
                    grad[i] -= g
                    grad[j] += g
        pos -= lr * grad
        pos[0] = 0.0
    return pos - pos[0]


def _orient_offsets(offsets: np.ndarray, tangent: np.ndarray, seed: int) -> np.ndarray:
    """Rotate residue-local offsets so the sidechain centroid points 'outward'
    (a per-residue direction perpendicular to the backbone tangent), reducing
    sidechain<->neighbour-backbone clashes. BB (offsets[0]==0) is unchanged."""
    if offsets.shape[0] <= 1:
        return offsets
    c = offsets[1:].mean(axis=0)
    nc = float(np.linalg.norm(c))
    if nc < 1e-6:
        return offsets
    axis_from = c / nc
    rng = np.random.default_rng(seed)
    t = tangent / (float(np.linalg.norm(tangent)) + 1e-9)
    v = rng.normal(0.0, 1.0, 3)
    outward = v - np.dot(v, t) * t
    no = float(np.linalg.norm(outward))
    if no < 1e-6:
        outward = np.array([1.0, 0.0, 0.0]); no = 1.0
    axis_to = outward / no
    vc = np.cross(axis_from, axis_to)
    s = float(np.linalg.norm(vc))
    cdot = float(np.dot(axis_from, axis_to))
    if s < 1e-9:
        R = np.eye(3) if cdot > 0 else -np.eye(3)
    else:
        k = vc / s
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        R = np.eye(3) + s * K + (1.0 - cdot) * (K @ K)
    return offsets @ R.T


def build_protein_system(
    snapshot: dict,
    sequence_3letter: list[str],
    library: MartiniLibrary,
    *,
    n_chains: int | None = None,
    block_seq: str | None = None,
    crosslinker_letter: str = "Y",
    lattice_scale_ang: float | None = None,
    sc_jitter_ang: float = DEFAULT_SC_JITTER_ANG,
    seed: int = 42,
) -> System:
    """Convert a BFM snapshot + residue sequence into a `System` of MARTINI beads.

    Parameters
    ----------
    snapshot
        One snapshot dict produced by `bfm.generate_topology` (or loaded from
        `topology_io`); has keys `chains`, `crosslinker_positions`, `Nx`,
        `Ny`, `Nz`, `reactions`, `conv`.
    sequence_3letter
        List of 3-letter residue codes per chain (e.g. ['GLY','GLY','ARG',...]).
        All chains share the same sequence in this first cut.
    library
        Loaded `MartiniLibrary` for bead mass/charge defaults.
    block_seq
        Optional one-letter repeat block; used to derive the crosslinker
        position within the block. Defaults to the resilin reference.

    Returns
    -------
    System with beads positioned in an orthogonal box of size
    ``(Nx, Ny, Nz) * lattice_scale_ang`` Angstroms.
    """
    rng = np.random.default_rng(seed)
    chains_flat: list[list[int]] = snapshot["chains"]
    Nx, Ny, Nz = snapshot["Nx"], snapshot["Ny"], snapshot["Nz"]
    if n_chains is None:
        n_chains = len(chains_flat)

    n_residues = len(sequence_3letter)
    n_nodes = len(chains_flat[0])
    n_segments_per_chain = n_nodes - 1
    if lattice_scale_ang is None:
        lattice_scale_ang = _auto_lattice_scale_ang(n_residues, n_segments_per_chain)

    box = np.array([Nx, Ny, Nz], dtype=float) * lattice_scale_ang
    # Build the BFM-node-to-residue mapping for this geometry; assumes every
    # chain has the same length (BFM enforces this).
    block = block_seq if block_seq is not None else residues.REFERENCE_BLOCK
    node_to_res = sequence.get_node_residue_mapping(
        n_repeats=n_residues // len(block),
        segs_per_block=(n_nodes - 1) // (n_residues // len(block)),
        block_seq=block,
        crosslinker_letter=crosslinker_letter,
    )
    inv = _build_inverse_mapping(node_to_res)

    # Crosslink dimer offsets: separate the two merged crosslinked TYR residues
    # so their rings don't land coincident (r~0). See template_builder for detail.
    crosslink_anchor_offset: dict[tuple[int, int], np.ndarray] = {}
    for rxn in snapshot.get("reactions", []):
        (ci1, ni1), (ci2, ni2) = rxn[0], rxn[1]
        r1 = node_to_res.get(ni1)
        r2 = node_to_res.get(ni2)
        if r1 is None or r2 is None:
            continue
        axrng = np.random.default_rng(seed * 991 + ni1 * 7 + ni2 * 13)
        axis = axrng.normal(0.0, 1.0, 3)
        axis = axis / (float(np.linalg.norm(axis)) + 1e-9)
        crosslink_anchor_offset[(ci1, r1)] = +axis * (CROSSLINK_SEP_ANG / 2.0)
        crosslink_anchor_offset[(ci2, r2)] = -axis * (CROSSLINK_SEP_ANG / 2.0)

    sys_ = System(box_dims_ang=tuple(box.tolist()))
    next_atom_id = 1
    # per-chain (chain_id, residue_idx, atom_name) -> atom_id
    name_to_atom: dict[tuple[int, int, str], int] = {}

    for ci, chain_flat in enumerate(chains_flat):
        molecule_id = ci + 1
        prev_bb_atom_id: int | None = None
        prev_prev_bb: int | None = None
        prev_prev_prev_bb: int | None = None
        prev_resname: str | None = None
        prev_prev_resname: str | None = None
        prev_prev_prev_resname: str | None = None

        # Topro-style placement: residues interpolated between adjacent BFM
        # nodes via min-image displacement. Positions may fall outside [0,box);
        # LAMMPS wraps them in on read_data and assigns image flags itself.
        residue_positions = _interpolate_residue_positions(
            chain_flat, node_to_res, Nx, Ny, lattice_scale_ang, box
        )

        for r_idx, resname in enumerate(sequence_3letter):
            res_def = residues.RESIDUES[resname]
            raw_beads = res_def["beads"]
            is_n_term = (r_idx == 0)
            is_c_term = (r_idx == n_residues - 1)
            patched = _patch_terminal_beads(raw_beads, is_n_term, is_c_term)

            anchor = (residue_positions.get(r_idx, np.zeros(3))
                      + crosslink_anchor_offset.get((ci, r_idx), np.zeros(3)))

            # proper sidechain geometry (embedded per residue, oriented outward),
            # not anchor+jitter (which piled SC beads on BB -> intra overlaps).
            _names = [bn for (bn, _bt, _q) in patched]
            _embed_order = (["BB"] if "BB" in _names else []) + [n for n in _names if n != "BB"]
            _ei = {n: i for i, n in enumerate(_embed_order)}
            _edges: list[tuple[int, int, float]] = []
            for _rec in res_def["intra_bonds"]:
                if _rec[0] in _ei and _rec[1] in _ei:
                    _edges.append((_ei[_rec[0]], _ei[_rec[1]], _rec[3] * 10.0))
            for _rec in res_def["intra_constraints"]:
                if _rec[0] in _ei and _rec[1] in _ei:
                    _edges.append((_ei[_rec[0]], _ei[_rec[1]], _rec[3] * 10.0))
            _offl = _embed_residue_local(len(_embed_order), _edges, seed=seed + r_idx)
            _tan = residue_positions.get(r_idx + 1, anchor) - residue_positions.get(r_idx - 1, anchor)
            if float(np.linalg.norm(_tan)) < 1e-6:
                _tan = np.array([0.0, 0.0, 1.0])
            _offl = _orient_offsets(_offl, _tan, seed=seed * 100003 + ci * 9176 + r_idx)
            _sc_off = {n: _offl[_ei[n]] for n in _embed_order}

            local_atom_ids: dict[str, int] = {}
            for atom_name, bead_type, charge in patched:
                if atom_name == "BB":
                    pos = anchor.copy()
                else:
                    pos = anchor + _sc_off.get(atom_name, np.zeros(3))
                bead = Bead(
                    atom_id=next_atom_id,
                    bead_type=bead_type,
                    molecule_id=molecule_id,
                    residue_idx=r_idx + 1,
                    residue_name=resname,
                    atom_name=atom_name,
                    charge=charge,
                    mass=library.get_mass(bead_type),
                    position=tuple(pos.tolist()),
                )
                sys_.beads.append(bead)
                local_atom_ids[atom_name] = next_atom_id
                name_to_atom[(ci, r_idx, atom_name)] = next_atom_id
                next_atom_id += 1

            # Intra-residue bonded terms
            for rec in res_def["intra_bonds"]:
                a_n, b_n, funct, length, k = rec
                if a_n in local_atom_ids and b_n in local_atom_ids:
                    sys_.bonds.append(Bond(local_atom_ids[a_n], local_atom_ids[b_n], funct, length, k))
            for rec in res_def["intra_constraints"]:
                a_n, b_n, _funct, length = rec
                if a_n in local_atom_ids and b_n in local_atom_ids:
                    sys_.constraints.append(Constraint(local_atom_ids[a_n], local_atom_ids[b_n], length))
            for rec in res_def["intra_angles"]:
                a_n, b_n, c_n, funct, angle, k = rec
                if all(n in local_atom_ids for n in (a_n, b_n, c_n)):
                    sys_.angles.append(Angle(local_atom_ids[a_n], local_atom_ids[b_n], local_atom_ids[c_n], funct, angle, k))
            for rec in res_def["intra_dihedrals_proper"]:
                a_n, b_n, c_n, d_n, funct, angle, k, mult = rec
                if all(n in local_atom_ids for n in (a_n, b_n, c_n, d_n)):
                    sys_.dihedrals.append(Dihedral(local_atom_ids[a_n], local_atom_ids[b_n], local_atom_ids[c_n], local_atom_ids[d_n], funct, angle, k, mult))
            for rec in res_def["intra_dihedrals_improper"]:
                a_n, b_n, c_n, d_n, funct, angle, k = rec
                if all(n in local_atom_ids for n in (a_n, b_n, c_n, d_n)):
                    sys_.dihedrals.append(Dihedral(local_atom_ids[a_n], local_atom_ids[b_n], local_atom_ids[c_n], local_atom_ids[d_n], funct, angle, k, mult=1, is_improper=True))
            for excl_atom_names in res_def["intra_exclusions"]:
                ids = tuple(local_atom_ids[n] for n in excl_atom_names if n in local_atom_ids)
                if len(ids) >= 2:
                    sys_.exclusions.append(Exclusion(atoms=ids))

            # Inter-residue: BB-BB bond + BBB angle + BBBB dihedral
            this_bb = local_atom_ids.get("BB")
            if this_bb is not None and prev_bb_atom_id is not None and prev_resname is not None:
                funct, length, k = _select_bb_bond(prev_resname, resname)
                sys_.bonds.append(Bond(prev_bb_atom_id, this_bb, funct, length, k))
            if (this_bb is not None and prev_bb_atom_id is not None
                    and prev_prev_bb is not None and residues.BACKBONE_BBB_ANGLE):
                funct, angle, k = residues.BACKBONE_BBB_ANGLE
                sys_.angles.append(Angle(prev_prev_bb, prev_bb_atom_id, this_bb, funct, angle, k))
            if (this_bb is not None and prev_bb_atom_id is not None
                    and prev_prev_bb is not None and prev_prev_prev_bb is not None
                    and prev_prev_resname is not None and prev_prev_prev_resname is not None
                    and prev_resname is not None):
                quad = (prev_prev_prev_resname, prev_prev_resname, prev_resname, resname)
                label = _select_bbbb_label(quad)
                for (funct, angle, k, mult) in residues.BACKBONE_BBBB_DIHEDRALS.get(label, []):
                    sys_.dihedrals.append(Dihedral(prev_prev_prev_bb, prev_prev_bb, prev_bb_atom_id, this_bb, funct, angle, k, mult))

            # advance backbone history
            prev_prev_prev_bb = prev_prev_bb
            prev_prev_bb = prev_bb_atom_id
            prev_bb_atom_id = this_bb
            prev_prev_prev_resname = prev_prev_resname
            prev_prev_resname = prev_resname
            prev_resname = resname

    # Apply dityrosine crosslinks recorded in the BFM snapshot.
    # Each reaction = [[ci1, ni1], [ci2, ni2]] where ni* are chain node indices.
    # All crosslinks are added here; the writer's priority-MST image-flag pass
    # (``lammps_writer._kruskal_image_flags_and_drop``) processes non-crosslink
    # bonds first and crosslinks last, so any winding-cycle back-edge in the
    # crosslink graph (≈ 0.1% of crosslinks in a typical BFM realisation) is
    # identified as a crosslink rather than a backbone bond and dropped from
    # the LAMMPS data file. Real backbone / sidechain bonds are tree edges by
    # construction and never drop.
    for rxn in snapshot.get("reactions", []):
        (ci1, ni1), (ci2, ni2) = rxn[0], rxn[1]
        r1 = node_to_res.get(ni1)
        r2 = node_to_res.get(ni2)
        if r1 is None or r2 is None:
            continue
        a = name_to_atom.get((ci1, r1, "SC4"))
        b = name_to_atom.get((ci2, r2, "SC4"))
        if a is None or b is None:
            continue
        sys_.bonds.append(
            Bond(a, b, funct=1, length_nm=DITYROSINE_BOND_LENGTH_NM,
                 k_kj=DITYROSINE_BOND_K_KJ, is_crosslink=True)
        )

    return sys_
