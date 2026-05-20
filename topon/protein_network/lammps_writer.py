"""LAMMPS data + settings + input-script writer for a MARTINI protein system.

Emits six files (the first three sibling to ``output_dir``, the three relaxation
stage scripts under ``output_dir/relaxation/``):

* ``<base>.data``                       LAMMPS data file (atom_style full):
                                         Atoms, Bonds, Angles, Dihedrals,
                                         Impropers + their type counts.
* ``<base>.in.settings``                pair_coeff / bond_coeff / angle_coeff /
                                         etc. lines covering every type id used
                                         in the data file.
* ``<base>.in.groups``                  group definitions for protein vs. water.
* ``relaxation/<base>_stage1.in``       soft-push overlap removal.
* ``relaxation/<base>_stage2.in``       LJ-epsilon ramp 0.001 -> 1.0 via
                                         nve/limit dynamics.
* ``relaxation/<base>_stage3.in``       tight CG min + brief NVT/NPT @ 310 K
                                         (uses ``units real``,
                                         ``lj/cut/coul/cut``, ``dielectric=15``,
                                         no PPPM).

LAMMPS units: ``real`` (Angstrom, kcal/mol). MARTINI bonded-term conversions
live in `martini_ff` and are exercised here.

GROMACS function-type mapping:

* bond funct 1 (harmonic)        ->  bond_style harmonic
* angle funct 2 (cosine sq.)     ->  angle_style cosine/squared
* angle funct 10 (restricted)    ->  angle_style cosine/squared (omits 1/sin^2 t,
                                                                 first-cut surrogate)
* dihedral funct 9 (proper, multi-term)  ->  dihedral_style charmm
                                              (one type per term; multiple
                                               types per quadruplet allowed)
* dihedral funct 2 (improper, harmonic)  ->  improper_style harmonic
                                              (placed in the Impropers section)
* `[ constraints ]`              ->  emitted as stiff harmonic bonds with
                                     K = CONSTRAINT_K_GMX = 10000 kJ/mol/nm^2
                                     (NOT 1e6 -- see the constant's rationale
                                      below; softened so mis-placed sidechains
                                      don't explode during NVT). This is the
                                      same K as the stiffest backbone bonds, so
                                      the constraint bonds are NOT the timestep
                                      limiter -- the dt~2 fs requirement comes
                                      from elsewhere (likely the r=0 dityrosine
                                      crosslinks before relaxation separates
                                      them; see internal/option_c_rigid/DESIGN.md).
                                     `fix shake` does NOT work on the TYR rings
                                     (over-constrained connected clusters);
                                     `fix rattle` or rigid bodies would be
                                     needed for true rigid constraints.

Reaction-field electrostatics are approximated with `dielectric 15.0` plus
`pair_style lj/cut/coul/cut`, mirroring the reference's `coulombtype = reaction-field
epsilon_r = 15`. The RF correction term is omitted (documented in the .in
header).

Image-flag assignment (added 2026-05; replaces the prior wrap-only convention):

  The Atoms section emits the full ``atom_style full`` row (10 columns:
  ``id mol type q x y z ix iy iz``). Image flags are computed by a
  priority-weighted MST (Kruskal) over the molecular bond graph using
  a 2-key sort:
    * priority 0 (non-crosslink bonds — backbone, sidechain, constraints):
      processed first, ordered by length within the priority.
    * priority 1 (crosslinks — dityrosine SC4-SC4): processed last,
      ordered by length within the priority.
  Image flags propagate from the spanning-tree root such that every
  tree-edge bond is minimum-image. Non-tree (cycle-closing) back-edges
  whose BFS-implied image-flag delta disagrees with the wrapped
  min-image delta are "winding-cycle bonds" — topologically forced
  long bonds around the periodic box, which no image-flag assignment
  can make MIC, and which break parallel-MPI ghost-shell construction
  in LAMMPS. They are dropped at write time.

  Why the priority key (2026-05-19 BB-BB-drop fix):
    The BFM topology generator merges two TYR residues onto the same
    lattice node at every dityrosine crosslink site. The two TYR/SC4
    beads belong to two different chains and are placed by two
    independent min-image chain walks — these walks may reach the
    shared lattice node from opposite sides of the periodic box, so
    the beads end up at the SAME wrapped position but with DIFFERENT
    natural (chain-walk) image flags. The crosslink bond, in wrapped
    MIC, is therefore length ≈ 0.05 A (essentially the per-axis
    coord_perturbation). With a pure length-sorted MST, every
    crosslink is added as a tree edge before any BB-BB bond. After
    the crosslinks merged chains into one giant component, BB-BB
    bonds (≈ 6.7 A at the BFM-projected scale) became the longest
    edges in chain-wraps-around-the-box cycles and got dropped.
    With the priority key, BB-BB bonds are now tree edges by
    construction and the crosslinks responsible for the cycle's
    winding are the ones dropped — matching the original design
    intent that crosslinks are the redundant, sacrificial elements.

  Why this is not the v33-v38 "phantom 240 A bonds" anti-pattern:
  v33-v38 assigned image flags per-chain by walking each chain
  independently and accumulating wrap counts. When two chains met at a
  dityrosine crosslink, their walk-accumulated images disagreed and
  the crosslink bond's image-flag-implied length was ~one box. This
  writer computes image flags over the GLOBAL bond graph using a
  spanning forest; tree edges are MIC by construction, regardless of
  which chain they belong to. The 240 A failure mode cannot recur.
"""
from __future__ import annotations

from collections import OrderedDict, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .martini_ff import (
    KJ_TO_KCAL,
    MartiniLibrary,
    gmx_angle_k_kjmol_to_lammps,
    gmx_bond_k_to_lammps,
    gmx_dihedral_k_to_lammps,
    gmx_lj_to_lammps,
)
from .system import System

CONSTRAINT_K_GMX: float = 10_000.0  # kJ/(mol*nm^2). Soft enough that mis-placed
# sidechain beads don't explode during NVT; users wanting true rigidity should
# add `fix shake` referencing the constraint bond_type ids in the input script.
HEADER_BANNER: str = (
    "# LAMMPS data file generated by topon.protein_network\n"
    "# MARTINI 3 protein network. Sequence- and BFM-driven.\n"
    "# Reaction-field electrostatics approximated as dielectric=15.0.\n"
)


def _bond_type_key(bond) -> tuple:
    if bond.is_crosslink:
        return ("crosslink", round(bond.length_nm, 4), round(bond.k_kj, 1))
    return (bond.funct, round(bond.length_nm, 4), round(bond.k_kj, 1))


def _constraint_type_key(c) -> tuple:
    return ("constraint", round(c.length_nm, 4))


def _angle_type_key(a) -> tuple:
    return (a.funct, round(a.angle_deg, 3), round(a.k_kj, 3))


def _dihedral_type_key(d) -> tuple:
    return (d.funct, round(d.angle_deg, 3), round(d.k_kj, 3), d.mult, d.is_improper)


def _assign_ids(records: list, key_fn) -> tuple[dict, list[int]]:
    """Return (key->type_id, per-record-type-id-list)."""
    type_ids: dict = OrderedDict()
    per_record: list[int] = []
    for r in records:
        k = key_fn(r)
        if k not in type_ids:
            type_ids[k] = len(type_ids) + 1
        per_record.append(type_ids[k])
    return type_ids, per_record


_BB_BEAD_TYPES: set[str] = {"Q5", "SP1", "P2", "SP2", "SP2a"}
_WATER_BEAD_TYPES: set[str] = {"W", "WF", "SW", "TW"}
_CROSSLINK_BEAD_TYPES: set[str] = {"TN6"}  # TYR SC4 = dityrosine site


# ----------------------------------------------------------------------
# Image-flag assignment via priority-weighted MST + winding-cycle drop
# ----------------------------------------------------------------------

def _kruskal_image_flags_and_drop(
    bead_positions: dict,
    all_bonds: list,
    Lx: float, Ly: float, Lz: float,
) -> tuple[dict, list[bool]]:
    """Assign per-atom image flags so every tree-edge of a priority-
    weighted MST over the bond graph is minimum-image; drop the back-
    edges whose unwrapped image-flag delta cannot be made MIC (cycles
    with non-zero winding around the periodic box).

    Args:
        bead_positions: mapping atom_id -> (x, y, z) WRAPPED into
            [0, Lx) x [0, Ly) x [0, Lz). Every atom that appears in
            ``all_bonds`` MUST have an entry here.
        all_bonds: list of bond records as 6-tuples
            ``(a_id, b_id, funct, length_nm, k_kj, is_crosslink)``.
            ``funct`` may be the int GROMACS function code or the
            string ``"constraint"`` for stiffened constraint bonds.
        Lx, Ly, Lz: box side lengths (A).

    Returns:
        ``(image_flags, keep_mask)`` where:
          * ``image_flags[atom_id] = (ix, iy, iz)`` for every atom.
            Atoms not reachable through bonds (e.g. water W beads with
            no bonds) get ``(0, 0, 0)``.
          * ``keep_mask`` is a list[bool] of length ``len(all_bonds)``;
            False = winding-cycle back-edge, drop.

    Sort priority (the structural fix for the "BB-BB drop" pathology):
        Bonds are sorted by a 2-key composite (priority, length):
          * priority 0 = non-crosslink bonds (backbone, sidechain,
            constraints) — chain-internal structural connectivity.
          * priority 1 = crosslinks (dityrosine SC4-SC4 between
            chains) — designed to be the redundant, sacrificial
            elements of the network.
        Within a priority, bonds are still ordered by length so that
        Kruskal's "longest back-edge gets dropped" property holds
        locally (e.g. inside a TYR ring triangle near a box boundary).

        Without the priority key, the BFM-merged dityrosine atoms
        place two SC4 beads at the same wrapped position (MIC distance
        ≈ 0.05 A), so crosslinks always sorted to the FRONT of the
        length list. After crosslinks merged chains into one giant
        component, BB-BB backbone bonds (~6.7 A at the BFM-projected
        scale, much longer than crosslinks) became the longest edges
        in every chain-wraps-around-the-box cycle and got dropped
        instead of the topologically responsible crosslink. The fix
        is to demote crosslinks to priority 1 so the longest edge in
        every cycle is now the crosslink itself.

    Why we still allow constraint drops:
        Even within a single TYR ring (SC1-SC2-SC3-SC4 triangles), if
        the 5 sidechain atoms straddle a box boundary in wrapped
        coords, NO image-flag assignment can satisfy all 5 constraint
        bonds simultaneously (one bond will always be ~half-box). That
        is a real topological winding inside the ring, separate from
        crosslink-induced winding, and dropping one constraint per
        offending ring is the principled response. Caller reports
        these counts.

    Hard invariant — real funct=1/9 bonds can NEVER drop:
        With the priority key, non-crosslink bonds form a forest
        (one tree per chain, plus the small TYR-ring cycles). The
        only non-crosslink edges that can land in a cycle are the
        intra-residue ring bonds (constraints). Any backbone /
        sidechain / inter-residue funct-int bond that drops is a
        topology corruption upstream; the caller raises AssertionError
        if this happens.

    Why an exact integer image-flag-delta criterion (not a magnitude
    threshold):
        10 A is project-specific magic that would be wrong for an
        atomistic CHARMM bond (~1.5 A) or a different lattice. The
        principled test is "after MST, does the bond's BFS-implied
        image-flag delta equal the wrapped min-image delta?" — pure
        integer arithmetic, no scale-dependent constants.
    """
    atom_ids = list(bead_positions.keys())
    id_to_idx = {aid: i for i, aid in enumerate(atom_ids)}
    box = np.array([Lx, Ly, Lz], dtype=float)

    def _mic_delta(a_id: int, b_id: int):
        """Return ``(d_mic, image_delta)`` for a bond a -> b.

        ``d_mic`` is the wrapped min-image displacement of B relative
        to A. ``image_delta`` is the integer (Δix, Δiy, Δiz) such that
        ``image_flags[B] = image_flags[A] + image_delta`` makes the
        bond minimum-image in unwrapped coords.

        Sign convention: raw difference ``d = B - A`` (wrapped); the
        number of boxes B needs shifting by is
        ``-round(d / L)``, hence ``image_delta = -round(d / L)``.
        """
        Ax, Ay, Az = bead_positions[a_id]
        Bx, By, Bz = bead_positions[b_id]
        d = np.array([Bx - Ax, By - Ay, Bz - Az])
        n_boxes = np.round(d / box).astype(int)
        d_mic = d - n_boxes * box
        return d_mic, -n_boxes

    # Pre-compute (length, image_delta) per bond. Sort key is
    # (priority, length): priority 0 for non-crosslinks (real bonds +
    # constraints) and priority 1 for crosslinks. This guarantees that
    # all non-crosslink bonds enter the MST as tree edges before any
    # crosslink is considered, so the only edges that can become
    # back-edges across chains are crosslinks themselves.
    bond_mic: list[tuple[np.ndarray, np.ndarray]] = []
    sortable: list[tuple[int, float, int]] = []
    for bidx, bond_tup in enumerate(all_bonds):
        a, b, _funct, _length_nm, _k_kj, is_xl = bond_tup
        d_mic, image_delta = _mic_delta(a, b)
        bond_mic.append((d_mic, image_delta))
        priority = 1 if is_xl else 0
        sortable.append((priority, float(np.linalg.norm(d_mic)), bidx))
    sortable.sort()

    # Kruskal MST with Union-Find
    n = len(atom_ids)
    parent = list(range(n))
    rank = [0] * n

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> bool:
        ra, rb = find(a), find(b)
        if ra == rb:
            return False
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        if rank[ra] == rank[rb]:
            rank[ra] += 1
        return True

    tree_edge_set: set[int] = set()
    backedge_bond_indices: list[int] = []
    for _priority, _length, bidx in sortable:
        a, b = all_bonds[bidx][0], all_bonds[bidx][1]
        ia, ib = id_to_idx[a], id_to_idx[b]
        if union(ia, ib):
            tree_edge_set.add(bidx)
        else:
            backedge_bond_indices.append(bidx)

    # BFS over each spanning-tree component to assign image flags
    tree_adj: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for bidx in tree_edge_set:
        a, b, *_ = all_bonds[bidx]
        tree_adj[a].append((b, bidx))
        tree_adj[b].append((a, bidx))

    image_flags: dict[int, tuple[int, int, int]] = {}
    visited_atoms: set[int] = set()
    for seed in atom_ids:
        if seed in visited_atoms:
            continue
        image_flags[seed] = (0, 0, 0)
        visited_atoms.add(seed)
        q = deque([seed])
        while q:
            u = q.popleft()
            for (v, bidx) in tree_adj[u]:
                if v in visited_atoms:
                    continue
                a_in_bond, b_in_bond, *_ = all_bonds[bidx]
                _, image_delta = bond_mic[bidx]
                if u == a_in_bond:
                    inc = image_delta
                else:
                    assert u == b_in_bond
                    inc = -image_delta
                u_img = np.array(image_flags[u], dtype=int)
                image_flags[v] = tuple(int(x) for x in (u_img + inc))
                visited_atoms.add(v)
                q.append(v)

    # Back-edge drop check: BFS-implied delta vs wrapped-MIC delta.
    # Drop bonds where the BFS-tree-implied image flags disagree with
    # the wrapped min-image delta (these are cycles with non-zero
    # winding around the box; topologically forced long bonds).
    keep: list[bool] = [True] * len(all_bonds)
    for bidx in backedge_bond_indices:
        ia = np.array(image_flags[all_bonds[bidx][0]], dtype=int)
        ib = np.array(image_flags[all_bonds[bidx][1]], dtype=int)
        bfs_delta = ib - ia
        _, mic_delta = bond_mic[bidx]
        mic_delta = np.asarray(mic_delta, dtype=int)
        if not np.array_equal(bfs_delta, mic_delta):
            keep[bidx] = False

    return image_flags, keep


def write_lammps(
    sys_: System,
    library: MartiniLibrary,
    output_dir: str | Path,
    *,
    base_name: str = "protein_network",
    include_input_script: bool = True,
    hierarchical_stage1: bool = False,
    coord_perturbation_ang: float = 0.05,
    coord_perturbation_seed: int = 7,
) -> dict[str, Path]:
    """Write the LAMMPS files for `sys_`. Returns paths keyed by artifact name
    (``data``, ``settings``, ``groups``, ``stage1``, ``stage2``, ``stage3``).

    If `hierarchical_stage1=True`, emits the core-topon-style progressive
    freeze/unfreeze stage 1 (mirrors `topon/writers/lammps_inputs.py:33` for
    CG networks), and adds per-class group definitions (bb/sc/crosslink/water)
    to the .in.groups file. Otherwise emits the flat topro-style stage 1.

    `coord_perturbation_ang`: per-axis Gaussian xyz jitter applied at write
    time (default 0.05 A). Mirrors core topon's pattern of breaking
    zero-distance degeneracies (e.g. dityrosine BB-BB pairs that BFM merges
    onto one lattice node) before LAMMPS sees them. Set to 0 to disable.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Type assignment ----------------------------------------------------------
    bead_types = sys_.bead_types_in_use()
    bead_type_id = {bt: i + 1 for i, bt in enumerate(bead_types)}

    # Treat constraints as additional bonds with ULTRA-stiff K
    constraint_bonds = [
        (c.a, c.b, "constraint", c.length_nm, CONSTRAINT_K_GMX, False)
        for c in sys_.constraints
    ]
    real_bonds = [
        (b.a, b.b, b.funct, b.length_nm, b.k_kj, b.is_crosslink)
        for b in sys_.bonds
    ]
    all_bonds = real_bonds + constraint_bonds

    bond_type_keys: OrderedDict = OrderedDict()
    bond_per_record: list[int] = []
    for (_a, _b, funct, length, k, is_crosslink) in all_bonds:
        if funct == "constraint":
            key = ("constraint", round(length, 4))
        elif is_crosslink:
            key = ("crosslink", round(length, 4), round(k, 1))
        else:
            key = (funct, round(length, 4), round(k, 1))
        if key not in bond_type_keys:
            bond_type_keys[key] = len(bond_type_keys) + 1
        bond_per_record.append(bond_type_keys[key])

    angle_type_keys, angle_per_record = _assign_ids(sys_.angles, _angle_type_key)

    propers = [d for d in sys_.dihedrals if not d.is_improper]
    impropers = [d for d in sys_.dihedrals if d.is_improper]
    dihedral_type_keys, dihedral_per_record = _assign_ids(propers, _dihedral_type_key)
    improper_type_keys, improper_per_record = _assign_ids(impropers, _dihedral_type_key)

    # Pre-compute wrapped bead positions + jitter (the same positions the
    # Atoms section will emit), then run priority-weighted MST over the
    # bond graph (non-crosslinks first, crosslinks last) to assign per-
    # atom image flags and identify winding-cycle back-edges that need
    # dropping. See _kruskal_image_flags_and_drop docstring for the
    # rationale; see module docstring for why this replaces the v33-v38
    # wrap-only convention.
    bx, by, bz = sys_.box_dims_ang
    _pert_rng = np.random.default_rng(coord_perturbation_seed)
    bead_positions: dict[int, tuple[float, float, float]] = {}
    for b in sys_.beads:
        if coord_perturbation_ang > 0.0:
            jx, jy, jz = _pert_rng.normal(0.0, coord_perturbation_ang, size=3)
        else:
            jx = jy = jz = 0.0
        x = (b.position[0] + jx) % bx
        y = (b.position[1] + jy) % by
        z = (b.position[2] + jz) % bz
        bead_positions[b.atom_id] = (x, y, z)

    image_flags, keep_mask = _kruskal_image_flags_and_drop(
        bead_positions, all_bonds, bx, by, bz,
    )
    kept_bonds = [bond for bond, k in zip(all_bonds, keep_mask) if k]
    kept_bond_per_record = [tid for tid, k in zip(bond_per_record, keep_mask) if k]
    n_dropped = sum(1 for k in keep_mask if not k)
    if n_dropped:
        n_xl_dropped = 0
        n_constraint_dropped = 0
        real_bond_drops: list[tuple[int, int, object]] = []
        for i, k in enumerate(keep_mask):
            if k:
                continue
            _a, _b, funct, _length, _kj, is_xl = all_bonds[i]
            if is_xl:
                n_xl_dropped += 1
            elif funct == "constraint":
                n_constraint_dropped += 1
            else:
                real_bond_drops.append((_a, _b, funct))

        # Hard invariant: with the priority-MST sort, real backbone /
        # sidechain bonds (funct=integer, is_crosslink=False) can NEVER
        # be back-edges in the spanning tree, so they must never drop.
        # If this assertion fires, something is corrupt upstream (e.g.
        # duplicate atom IDs in the System, or a chain whose bond graph
        # was made disconnected by a bug). The previous "warning only"
        # behaviour was masking the real BB-BB drop pathology that the
        # priority-MST fix addresses.
        assert not real_bond_drops, (
            f"BUG: {len(real_bond_drops)} real (non-crosslink, "
            f"non-constraint) bond(s) flagged as winding-cycle "
            f"back-edges by the priority-MST. First 5: "
            f"{real_bond_drops[:5]}. With the priority-MST sort key, "
            f"all non-crosslink bonds must be tree edges; this can "
            f"only happen if the System has structural corruption "
            f"(duplicate atom IDs, disconnected chains, etc.)."
        )

        print(
            f"[protein_network/lammps_writer] dropped {n_dropped} "
            f"winding-cycle bond(s) out of {len(all_bonds)} total "
            f"({100.0 * n_dropped / max(len(all_bonds), 1):.3f} %). "
            f"These are topologically forced long bonds in cycles with "
            f"non-zero winding around the periodic box and cannot be "
            f"made minimum-image by any image-flag assignment.\n"
            f"  Breakdown:  "
            f"crosslinks: {n_xl_dropped} (BFM crosslink graph has "
            f"winding cycles around the periodic box; the longest "
            f"crosslink in each cycle gets dropped), "
            f"constraints: {n_constraint_dropped} (TYR ring atoms "
            f"that straddle a box boundary in wrapped coords; "
            f"unavoidable when SC1-SC4 cluster crosses a face)"
        )

    # Data file ----------------------------------------------------------------
    data_path = out / f"{base_name}.data"
    with data_path.open("w", encoding="utf-8") as f:
        f.write(HEADER_BANNER)
        f.write("\n")
        f.write(f"{sys_.n_atoms()} atoms\n")
        f.write(f"{len(kept_bonds)} bonds\n")
        f.write(f"{len(sys_.angles)} angles\n")
        f.write(f"{len(propers)} dihedrals\n")
        f.write(f"{len(impropers)} impropers\n\n")
        f.write(f"{len(bead_types)} atom types\n")
        f.write(f"{len(bond_type_keys)} bond types\n")
        f.write(f"{max(1, len(angle_type_keys))} angle types\n")
        f.write(f"{max(1, len(dihedral_type_keys))} dihedral types\n")
        f.write(f"{max(1, len(improper_type_keys))} improper types\n\n")
        f.write(f"0.0 {bx:.6f} xlo xhi\n")
        f.write(f"0.0 {by:.6f} ylo yhi\n")
        f.write(f"0.0 {bz:.6f} zlo zhi\n\n")

        f.write("Masses\n\n")
        for bt, tid in bead_type_id.items():
            mass = library.atomtypes[bt].mass if bt in library.atomtypes else 72.0
            f.write(f"{tid} {mass:.4f}  # {bt}\n")
        f.write("\n")

        # Atoms section: 10 columns including image flags (id mol type q x y z
        # ix iy iz). Image flags come from the MST-assigned set computed
        # above. The xyz perturbation is the same one the previous wrap-only
        # writer used, applied during bead_positions build above.
        f.write("Atoms  # full\n\n")
        for b in sys_.beads:
            tid = bead_type_id[b.bead_type]
            x, y, z = bead_positions[b.atom_id]
            ix, iy, iz = image_flags.get(b.atom_id, (0, 0, 0))
            f.write(
                f"{b.atom_id} {b.molecule_id} {tid} {b.charge:.4f} "
                f"{x:.6f} {y:.6f} {z:.6f} {ix} {iy} {iz}  "
                f"# {b.residue_name}/{b.atom_name} {b.bead_type}\n"
            )
        f.write("\n")

        if kept_bonds:
            f.write("Bonds\n\n")
            for i, (a, b, funct, length, k, is_xl) in enumerate(kept_bonds, start=1):
                tid = kept_bond_per_record[i - 1]
                tag = "crosslink" if is_xl else ("constraint" if funct == "constraint" else "bond")
                f.write(f"{i} {tid} {a} {b}  # {tag}\n")
            f.write("\n")

        if sys_.angles:
            f.write("Angles\n\n")
            for i, ang in enumerate(sys_.angles, start=1):
                tid = angle_per_record[i - 1]
                f.write(f"{i} {tid} {ang.a} {ang.b} {ang.c}\n")
            f.write("\n")

        if propers:
            f.write("Dihedrals\n\n")
            for i, d in enumerate(propers, start=1):
                tid = dihedral_per_record[i - 1]
                f.write(f"{i} {tid} {d.a} {d.b} {d.c} {d.d}\n")
            f.write("\n")

        if impropers:
            f.write("Impropers\n\n")
            for i, d in enumerate(impropers, start=1):
                tid = improper_per_record[i - 1]
                f.write(f"{i} {tid} {d.a} {d.b} {d.c} {d.d}\n")
            f.write("\n")

    # Settings file ------------------------------------------------------------
    settings_path = out / f"{base_name}.in.settings"
    with settings_path.open("w", encoding="utf-8") as f:
        f.write("# pair_coeff: explicit MARTINI 3 nonbond_params (LAMMPS units real).\n")
        for i, bt_i in enumerate(bead_types, start=1):
            for j, bt_j in enumerate(bead_types[i - 1:], start=i):
                sigma_nm, eps_kj = library.get_lj_pair(bt_i, bt_j)
                sigma_a, eps_kcal = gmx_lj_to_lammps(sigma_nm, eps_kj)
                f.write(
                    f"pair_coeff {i} {j} {eps_kcal:.6f} {sigma_a:.6f}  "
                    f"# {bt_i}-{bt_j}\n"
                )
        f.write("\n# bond_coeff: harmonic for real bonds + constraint bonds (ULTRA stiff).\n")
        for key, tid in bond_type_keys.items():
            kind = key[0]
            if kind == "constraint":
                length_nm = key[1]
                k_l = gmx_bond_k_to_lammps(CONSTRAINT_K_GMX)
                f.write(
                    f"bond_coeff {tid} {k_l:.6f} {length_nm * 10:.6f}  "
                    f"# constraint (rigid-stiff)\n"
                )
            elif kind == "crosslink":
                length_nm, k_kj = key[1], key[2]
                k_l = gmx_bond_k_to_lammps(k_kj)
                f.write(
                    f"bond_coeff {tid} {k_l:.6f} {length_nm * 10:.6f}  "
                    f"# dityrosine crosslink\n"
                )
            else:
                funct, length_nm, k_kj = key
                k_l = gmx_bond_k_to_lammps(k_kj)
                f.write(
                    f"bond_coeff {tid} {k_l:.6f} {length_nm * 10:.6f}  "
                    f"# funct {funct} K={k_kj} kJ/mol/nm^2 r0={length_nm} nm\n"
                )

        if angle_type_keys:
            f.write("\n# angle_coeff: cosine/squared (covers funct 2 and approximates funct 10).\n")
            for key, tid in angle_type_keys.items():
                funct, angle_deg, k_kj = key
                k_l = gmx_angle_k_kjmol_to_lammps(k_kj)
                f.write(
                    f"angle_coeff {tid} {k_l:.6f} {angle_deg:.4f}  "
                    f"# funct {funct}\n"
                )

        if dihedral_type_keys:
            f.write("\n# dihedral_coeff: charmm style for funct 9 (proper, multi-term).\n")
            for key, tid in dihedral_type_keys.items():
                funct, angle_deg, k_kj, mult, _is_imp = key
                k_l = gmx_dihedral_k_to_lammps(k_kj)
                # CHARMM style: K (kcal/mol), n (mult), d (degrees), w=0 (no 1-4 weight override)
                f.write(
                    f"dihedral_coeff {tid} {k_l:.6f} {mult} {int(round(angle_deg))} 0.0  "
                    f"# funct {funct}\n"
                )

        if improper_type_keys:
            f.write("\n# improper_coeff: harmonic (funct 2 -- TYR ring planarity, etc.).\n")
            for key, tid in improper_type_keys.items():
                funct, angle_deg, k_kj, _mult, _is_imp = key
                k_l = gmx_angle_k_kjmol_to_lammps(k_kj)
                f.write(
                    f"improper_coeff {tid} {k_l:.6f} {angle_deg:.4f}  "
                    f"# funct {funct}\n"
                )

    # Groups file --------------------------------------------------------------
    groups_path = out / f"{base_name}.in.groups"
    with groups_path.open("w", encoding="utf-8") as f:
        protein_mols = sorted({b.molecule_id for b in sys_.beads if b.bead_type != "W"})
        water_mols = sorted({b.molecule_id for b in sys_.beads if b.bead_type == "W"})
        if protein_mols:
            lo, hi = protein_mols[0], protein_mols[-1]
            f.write(f"group protein molecule {lo}:{hi}\n")
        if water_mols:
            lo, hi = water_mols[0], water_mols[-1]
            f.write(f"group water molecule {lo}:{hi}\n")
        f.write(f"# {len(protein_mols)} protein chains, {len(water_mols)} W water beads\n")
        if hierarchical_stage1:
            # Per-class type groups for the hierarchical freeze/unfreeze protocol.
            # MARTINI 3 protein bead types fall into 4 classes:
            #   bb        = backbone (BB beads of all residues; types depend on residue)
            #   sc        = sidechains (excluding the TYR SC4 dityrosine site)
            #   crosslink = TN6 (TYR SC4) -- the merged BFM lattice site bead
            #   water     = W
            f.write("\n# Hierarchical-relaxation per-class groups (used by stage 1)\n")
            bb_ids = sorted(tid for bt, tid in bead_type_id.items() if bt in _BB_BEAD_TYPES)
            sc_ids = sorted(tid for bt, tid in bead_type_id.items()
                             if bt not in (_BB_BEAD_TYPES | _WATER_BEAD_TYPES | _CROSSLINK_BEAD_TYPES))
            xl_ids = sorted(tid for bt, tid in bead_type_id.items() if bt in _CROSSLINK_BEAD_TYPES)
            w_ids = sorted(tid for bt, tid in bead_type_id.items() if bt in _WATER_BEAD_TYPES)
            if bb_ids: f.write(f"group bb type {' '.join(map(str, bb_ids))}\n")
            if sc_ids: f.write(f"group sc type {' '.join(map(str, sc_ids))}\n")
            if xl_ids: f.write(f"group crosslink type {' '.join(map(str, xl_ids))}\n")
            if w_ids: f.write(f"group wbeads type {' '.join(map(str, w_ids))}\n")

    paths: dict[str, Path] = {
        "data": data_path,
        "settings": settings_path,
        "groups": groups_path,
    }

    if include_input_script:
        relax_dir = out / "relaxation"
        relax_dir.mkdir(parents=True, exist_ok=True)
        stage1_path = relax_dir / f"{base_name}_stage1.in"
        stage2_path = relax_dir / f"{base_name}_stage2.in"
        stage3_path = relax_dir / f"{base_name}_stage3.in"
        if hierarchical_stage1:
            stage1_path.write_text(_stage1_hierarchical(base_name, bool(w_ids if hierarchical_stage1 else False)), encoding="utf-8")
        else:
            stage1_path.write_text(_stage1_soft(base_name), encoding="utf-8")
        stage2_path.write_text(_stage2_ljramp(base_name), encoding="utf-8")
        stage3_path.write_text(_stage3_min_nvt_npt(base_name), encoding="utf-8")
        paths["stage1"] = stage1_path
        paths["stage2"] = stage2_path
        paths["stage3"] = stage3_path

    return paths


def _stage1_hierarchical(base_name: str, has_water: bool) -> str:
    """Hierarchical freeze/unfreeze soft-push, mirroring core topon
    `lammps_inputs.py:write_serial_soft_minimization` for CG networks.

    Adapted for MARTINI -- the key wrinkle is that BFM-derived dityrosine
    crosslinks place two SC4 atoms at exactly r=0 (same merged lattice node),
    so they MUST be free to separate during the first soft-push minimization
    or the LJ in stage 2 blows up. Therefore we only freeze the backbone in
    Stage A, then unfreeze for Stage B.

    Stages:
      A: freeze bb -> sc + crosslink + water relax (crosslinks separate to ~2.7A)
      B: unfreeze bb -> full soft minimisation
    """
    return f"""# Stage 1 -- Hierarchical soft-push (topon.protein_network, MARTINI 3).
# Mirrors core topon's progressive freeze/unfreeze pattern from
# topon/writers/lammps_inputs.py:write_serial_soft_minimization, adapted so
# the dityrosine crosslink atoms (TN6, initially at r=0 with their partner)
# can separate before LJ engagement in stage 2.
# Run from this relaxation/ directory.

{_COMMON_HEADER}
pair_style      lj/cut/coul/cut 12.0
pair_modify     shift yes
special_bonds   lj 0.0 1.0 1.0 coul 0.0 1.0 1.0
dielectric      15.0

read_data       ../{base_name}.data
include         ../{base_name}.in.settings
include         ../{base_name}.in.groups

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes one 10000
comm_modify     mode single cutoff 60.0

# Switch to soft potential
pair_style      soft 1.0
pair_coeff      * * 0.0
variable        prefactor equal ramp(0,60)

thermo          100
thermo_style    custom step pe ke etotal evdwl epair ebond eangle press vol

# === Stage A: freeze backbone; sc + crosslink (+ water if present) relax ===
# Backbone frozen so the BFM-derived chain skeleton is preserved while
# sidechains and dityrosine crosslink atoms (initially overlapping) find space.
fix             freeze_bb bb setforce 0 0 0
fix             soft_push all adapt 1 pair soft a * * v_prefactor
min_style       cg
minimize        1.0e-4 1.0e-6 5000 50000
unfix           soft_push
unfix           freeze_bb
write_data      min_stage_A.data

# === Stage B: full soft minimisation (everything free) ===
fix             soft_push all adapt 1 pair soft a * * v_prefactor
minimize        1.0e-4 1.0e-6 5000 50000
unfix           soft_push

# === Stage C: nve/limit dynamics + final min ===
# CG minimize gets stuck at the r=0 saddle of the soft potential (gradient=0
# at the cosine peak); brief NVE/limit dynamics nudges atoms off it. This is
# the same pattern core topon uses for atomistic relaxation
# (lammps_inputs.py:141-145).
reset_timestep  0
timestep        1.0
fix             soft_push all adapt 1 pair soft a * * v_prefactor
fix             nve_limit all nve/limit 0.1
run             1000
unfix           nve_limit
unfix           soft_push
fix             soft_push all adapt 1 pair soft a * * v_prefactor
minimize        1.0e-4 1.0e-6 5000 50000
unfix           soft_push

write_data      system_after_soft.data
write_restart   1.restart
print           "Stage 1 (hierarchical) done -> system_after_soft.data"
"""


_COMMON_HEADER = """units           real
boundary        p p p
atom_style      full
bond_style      harmonic
angle_style     cosine/squared
dihedral_style  charmm
improper_style  harmonic
"""


def _stage1_soft(base_name: str) -> str:
    return f"""# Stage 1 -- Soft-push overlap removal (topon.protein_network, MARTINI 3).
# Resolves overlaps from the BFM lattice -> Cartesian projection and the
# voxel-grid water packer. Adapted from the topro CHARMM stage 1 protocol.
# Run from this relaxation/ directory.

{_COMMON_HEADER}
pair_style      lj/cut/coul/cut 12.0
pair_modify     shift yes
special_bonds   lj 0.0 1.0 1.0 coul 0.0 1.0 1.0
dielectric      15.0

read_data       ../{base_name}.data
include         ../{base_name}.in.settings
include         ../{base_name}.in.groups

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes one 10000
comm_modify     mode single cutoff 60.0  # MARTINI bonds + sidechains can stretch under soft potential  # tolerate stretched IDP backbone bonds

# Switch to soft potential for overlap removal
pair_style      soft 1.0
pair_coeff      * * 0.0
variable        prefactor equal ramp(0,60)

thermo          100
thermo_style    custom step pe ke etotal evdwl ecoul epair ebond eangle edihed eimp press vol

# Stage A: ramped soft push under CG minimisation
fix             soft_push all adapt 1 pair soft a * * v_prefactor
min_style       cg
minimize        1.0e-4 1.0e-6 1000 10000
unfix           soft_push
write_data      min_stage_A.data

# Stage B: brief nve/limit dynamics under soft potential
reset_timestep  0
timestep        1.0
fix             soft_push all adapt 1 pair soft a * * v_prefactor
fix             nve_limit all nve/limit 0.1
run             1000
unfix           nve_limit
unfix           soft_push

# Stage C: final soft minimisation
fix             soft_push all adapt 1 pair soft a * * v_prefactor
minimize        1.0e-4 1.0e-6 1000 10000
unfix           soft_push

write_data      system_after_soft.data
write_restart   1.restart
print           "Stage 1 done -> system_after_soft.data"
"""


def _stage2_ljramp(base_name: str) -> str:
    return f"""# Stage 2 -- LJ epsilon ramp 0.001 -> 1.0 (topon.protein_network, MARTINI 3).
# Engages MARTINI 3 LJ interactions gradually under nve/limit dynamics.
# Adapted from the topro CHARMM stage 2 protocol; switched from
# lj/charmm/coul/long to MARTINI's lj/cut/coul/cut.
# Run from this relaxation/ directory.

{_COMMON_HEADER}
pair_style      soft 1.0

read_data       system_after_soft.data
include         ../{base_name}.in.groups

neigh_modify    one 10000
comm_modify     mode single cutoff 60.0  # MARTINI bonds + sidechains can stretch under soft potential

# Pre-minimise once more under soft potential
pair_coeff      * * 1.0
min_style       cg
minimize        1.0e-4 1.0e-6 1000 10000

# Switch to MARTINI lj/cut/coul/cut (supports `fix adapt` epsilon scaling)
pair_style      lj/cut/coul/cut 12.0
pair_modify     shift yes
special_bonds   lj 0.0 1.0 1.0 coul 0.0 1.0 1.0
dielectric      15.0
include         ../{base_name}.in.settings

variable        scale equal ramp(0.001,1.0)
timestep        1.0

thermo          1000
thermo_style    custom step pe ke etotal evdwl ecoul epair ebond eangle edihed eimp press vol temp

fix             1 all adapt 1 pair lj/cut/coul/cut epsilon * * v_scale scale yes
fix             fxnve all nve/limit 0.1
run             20000
unfix           fxnve
unfix           1

write_data      system_ramped.data
write_restart   2.restart
print           "Stage 2 done -> system_ramped.data"
"""


def _stage3_min_nvt_npt(base_name: str) -> str:
    return f"""# Stage 3 -- Tight min + short NVT + (optional) NPT.
# Faithful MARTINI 3 port of topro's CHARMM stage 3 (same structure, same
# 1 fs timestep, same 100.0 Tdamp). MARTINI-specific changes are noted inline:
#   * pair_style lj/cut/coul/cut    (vs lj/charmm/coul/long + kspace pppm)
#   * dielectric 15.0               (RF approximation; MARTINI default)
#   * special_bonds lj/coul 0 1 1   (nrexcl=1: exclude 1-2 only, matching the
#     reference high_pro.itp + GROMACS; was 0 0 0 which OVER-excluded 1-3/1-4 and
#     let the proper sidechain geometry collapse during dynamics)
#   * angle_style cosine/squared    (MARTINI 3 IDP backbone funct=10)
#   * no fix cmap                   (CMAP is CHARMM-only)
#   * comm_modify cutoff 60         (MARTINI bonds longer than atomistic)
# Run from this relaxation/ directory.

{_COMMON_HEADER}
pair_style      lj/cut/coul/cut 12.0    # MARTINI: was lj/charmm/coul/long + pppm
pair_modify     shift yes
special_bonds   lj 0.0 1.0 1.0 coul 0.0 1.0 1.0   # MARTINI: was special_bonds charmm
dielectric      15.0                              # MARTINI: RF approx (no equivalent in CHARMM)

read_data       system_ramped.data
include         ../{base_name}.in.groups
include         ../{base_name}.in.settings

neigh_modify    one 10000
comm_modify     mode single cutoff 60.0           # MARTINI: longer than CHARMM atomistic

thermo          100
thermo_style    custom step pe ke etotal evdwl ecoul epair ebond eangle edihed eimp press vol temp

# --- Tight minimisation (mirrors topro stage 3) ---
min_style       cg
minimize        1.0e-6 1.0e-8 100000 1000000
write_data      system_minimized_final.data

# --- Short NVT (5000 steps x 1 fs = 5 ps), T=310K is the MARTINI canonical ---
reset_timestep  0
variable        T equal 310.0                     # MARTINI: 310 K (topro CHARMM used 300 K)
velocity        all create ${{T}} 12345
timestep        1.0                               # same as topro CHARMM stage 3
fix             1 all nvt temp ${{T}} ${{T}} 100.0   # same Tdamp as topro
run             5000
unfix           1
write_data      after_nvt.data

# --- Short NPT at 1 atm (5000 steps x 1 fs = 5 ps) ---
# OFF by default: for dry / vacuum boxes (no W water packed), NPT collapses
# the cell because there is no solvent pressure to balance against.
# UNCOMMENT for hydrated systems (--water-density > 0) to let the box find
# its natural MARTINI density.
#
# fix             1 all npt temp ${{T}} ${{T}} 100.0 iso 1.0 1.0 1000.0
# run             5000
# unfix           1

write_data      ../system_equilibrated.data
print           "Stage 3 done -> ../system_equilibrated.data (NVT-only; uncomment NPT block for hydrated runs)."
"""
