"""6-neighbour cubic-lattice topology generator with stochastic crosslinking.

Forked from `legacy/subprojects/protein_network/topro/topro/bfm/{generator.py,
percolation.py}` (2026-01) without behavioural change. Single-file consolidation:
the previous `from .percolation import UnionFind` is replaced by an inlined
class. The matplotlib visualizer is intentionally not forked.

Algorithm
---------
1. Compute a cubic lattice large enough to hold n_chains at target_packing.
2. Place each chain as a self-avoiding random walk (SAW) with periodic
   boundary conditions and full inter-chain excluded volume.
3. Run optional Monte Carlo equilibration (end, kink/crankshaft, reptation moves).
4. Find all potential crosslink pairs: Y nodes on 6-adjacent lattice sites.
5. Shuffle candidates; apply crosslinks one by one, merging lattice sites.
6. After each crosslink, test for gel point (all chains in one connected component
   via Union-Find on the chain graph).
7. Save snapshots at the gel point and at n_extra_snapshots conversion levels beyond.

The "Y node" is a generic crosslinker position. For MARTINI protein networks
this corresponds to a TYR residue's lattice site (dityrosine bond formation);
nothing in this module is residue-aware.

Snapshot JSON format is preserved bit-for-bit so that JSONs produced by topro
are interchangeable with this module.
"""
from __future__ import annotations

import math
import warnings
from collections import defaultdict
from typing import Iterable

import numpy as np


# UnionFind (inlined from topro.bfm.percolation) ============================

class UnionFind:
    """Path-compressed, rank-based Union-Find."""

    def __init__(self, n: int) -> None:
        self.parent: list[int] = list(range(n))
        self.rank: list[int] = [0] * n

    def find(self, x: int) -> int:
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x: int, y: int) -> bool:
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return False
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1
        return True

    def n_components(self) -> int:
        return sum(1 for i in range(len(self.parent)) if self.find(i) == i)

    def components(self) -> dict[int, list[int]]:
        comps: dict[int, list[int]] = defaultdict(list)
        for i in range(len(self.parent)):
            comps[self.find(i)].append(i)
        return dict(comps)


def analyze_connectivity(snapshot: dict) -> dict:
    """Post-hoc cluster analysis of a snapshot dict."""
    reactions = snapshot["reactions"]
    n_chains = len(snapshot["chains"])
    uf = UnionFind(n_chains)
    for rxn in reactions:
        ci1, ci2 = rxn[0][0], rxn[1][0]
        if ci1 != ci2:
            uf.union(ci1, ci2)
    comps = uf.components()
    largest = max(len(v) for v in comps.values()) / n_chains if comps else 0.0
    return {
        "n_components": len(comps),
        "largest_cluster_frac": largest,
        "is_gel": len(comps) == 1,
        "conversion": snapshot["conv"],
    }


# Lattice helpers ===========================================================

def flat_to_xyz(flat_idx: int, Nx: int, Ny: int) -> tuple[int, int, int]:
    z = flat_idx // (Nx * Ny)
    rem = flat_idx % (Nx * Ny)
    y = rem // Nx
    x = rem % Nx
    return x, y, z


def xyz_to_flat(x: int, y: int, z: int, Nx: int, Ny: int) -> int:
    return z * Nx * Ny + y * Nx + x


def get_neighbors(x: int, y: int, z: int, Nx: int, Ny: int, Nz: int) -> list[int]:
    """Return the 6 lattice-adjacent flat indices (periodic)."""
    return [
        xyz_to_flat((x + 1) % Nx, y, z, Nx, Ny),
        xyz_to_flat((x - 1) % Nx, y, z, Nx, Ny),
        xyz_to_flat(x, (y + 1) % Ny, z, Nx, Ny),
        xyz_to_flat(x, (y - 1) % Ny, z, Nx, Ny),
        xyz_to_flat(x, y, (z + 1) % Nz, Nx, Ny),
        xyz_to_flat(x, y, (z - 1) % Nz, Nx, Ny),
    ]


def compute_lattice_size(n_chains: int, n_nodes_per_chain: int, target_packing: float = 0.45) -> int:
    """Return cubic lattice edge L with (n_chains*n_nodes)/L^3 ~= target_packing.

    Rounded up to the nearest odd integer (>= 9).
    """
    total = n_chains * n_nodes_per_chain
    L = math.ceil((total / target_packing) ** (1.0 / 3.0))
    if L % 2 == 0:
        L += 1
    return max(9, L)


# Chain placement ===========================================================

def _place_single_chain(occupied, Nx, Ny, Nz, n_nodes, rng, max_attempts=3000):
    N = Nx * Ny * Nz
    free = [i for i in range(N) if i not in occupied]
    if not free:
        return None
    for _ in range(max_attempts):
        start = free[int(rng.integers(len(free)))]
        chain = [start]
        chain_set = {start}
        failed = False
        for _ in range(n_nodes - 1):
            x, y, z = flat_to_xyz(chain[-1], Nx, Ny)
            nbs = get_neighbors(x, y, z, Nx, Ny, Nz)
            rng.shuffle(nbs)
            moved = False
            for nb in nbs:
                if nb not in chain_set and nb not in occupied:
                    chain.append(nb)
                    chain_set.add(nb)
                    moved = True
                    break
            if not moved:
                failed = True
                break
        if not failed:
            return chain
    return None


def place_chains(n_chains, n_nodes_per_chain, Nx, Ny, Nz, rng, max_restarts=300):
    """Place all chains as SAWs with full inter-chain excluded volume."""
    for restart in range(max_restarts):
        occupied: set[int] = set()
        chains: list[list[int]] = []
        success = True
        for _ in range(n_chains):
            chain = _place_single_chain(occupied, Nx, Ny, Nz, n_nodes_per_chain, rng)
            if chain is None:
                success = False
                break
            for fi in chain:
                occupied.add(fi)
            chains.append(chain)
        if success:
            return chains
    raise RuntimeError(
        f"Cannot place {n_chains} chains of length {n_nodes_per_chain} "
        f"on {Nx}x{Ny}x{Nz} lattice after {max_restarts} restarts. "
        "Try increasing target_packing, reducing n_chains, or using a larger lattice."
    )


# MC equilibration ==========================================================

def equilibrate(chains, Nx, Ny, Nz, n_steps, rng, report_interval=0):
    """Monte Carlo equilibration: end moves, kink/crankshaft moves, reptation."""
    n_chains = len(chains)
    occupied: set[int] = set()
    for chain in chains:
        for fi in chain:
            occupied.add(fi)
    accepted = 0
    for step in range(n_steps):
        ci = int(rng.integers(n_chains))
        chain = chains[ci]
        n = len(chain)
        move = int(rng.integers(3))

        if move == 0 and n >= 2:  # end move
            end_side = int(rng.integers(2))
            end_i, inner_i = (0, 1) if end_side == 0 else (n - 1, n - 2)
            old_fi = chain[end_i]
            ix, iy, iz = flat_to_xyz(chain[inner_i], Nx, Ny)
            candidates = [
                nb for nb in get_neighbors(ix, iy, iz, Nx, Ny, Nz)
                if nb not in occupied or nb == old_fi
            ]
            if candidates:
                new_fi = candidates[int(rng.integers(len(candidates)))]
                occupied.discard(old_fi)
                occupied.add(new_fi)
                chain[end_i] = new_fi
                accepted += 1

        elif move == 1 and n >= 3:  # kink / crankshaft
            i = int(rng.integers(1, n - 1))
            prev_fi = chain[i - 1]
            curr_fi = chain[i]
            next_fi = chain[i + 1]
            px, py, pz = flat_to_xyz(prev_fi, Nx, Ny)
            nx_, ny_, nz_ = flat_to_xyz(next_fi, Nx, Ny)
            prev_nbs = set(get_neighbors(px, py, pz, Nx, Ny, Nz))
            next_nbs = set(get_neighbors(nx_, ny_, nz_, Nx, Ny, Nz))
            valid = [
                v for v in (prev_nbs & next_nbs) - {prev_fi, next_fi}
                if v not in occupied or v == curr_fi
            ]
            if valid:
                new_fi = valid[int(rng.integers(len(valid)))]
                occupied.discard(curr_fi)
                occupied.add(new_fi)
                chain[i] = new_fi
                accepted += 1

        elif move == 2 and n >= 2:  # reptation
            direction = int(rng.integers(2))
            if direction == 0:
                remove_fi = chain[-1]
                grow_from = chain[0]
            else:
                remove_fi = chain[0]
                grow_from = chain[-1]
            gx, gy, gz = flat_to_xyz(grow_from, Nx, Ny)
            candidates = [
                nb for nb in get_neighbors(gx, gy, gz, Nx, Ny, Nz)
                if nb not in occupied or nb == remove_fi
            ]
            if candidates:
                new_fi = candidates[int(rng.integers(len(candidates)))]
                occupied.discard(remove_fi)
                occupied.add(new_fi)
                if direction == 0:
                    chain.pop()
                    chain.insert(0, new_fi)
                else:
                    chain.pop(0)
                    chain.append(new_fi)
                accepted += 1

        if report_interval and (step + 1) % report_interval == 0:
            rate = accepted / (step + 1)
            print(f"  MC step {step+1:>8,} / {n_steps:,}  acceptance={rate:.3f}")
    return accepted / n_steps if n_steps > 0 else 0.0


# Crosslinker (Y-node) positions ============================================

def get_y_positions(n_repeats: int, segs_per_block: int, y_offset_in_block: int = 0) -> list[int]:
    """Return the list of crosslinker (Y) chain-node indices.

    With segs_per_block=2, y_offset=0 -> [1, 3, 5, ..., 2*n_repeats-1].
    With segs_per_block=3, y_offset=1 -> [2, 5, 8, ..., 3*n_repeats-1].
    """
    return [k * segs_per_block + y_offset_in_block + 1 for k in range(n_repeats)]


def compute_chain_images(
    chain: list[int], Nx: int, Ny: int, Nz: int,
) -> list[tuple[int, int, int]]:
    """Track image-cell offset per node along a chain by detecting periodic-boundary crossings.

    Each step (i-1 -> i) is a 6-neighbour move on the periodic lattice. If the
    raw lattice-coordinate difference is +-1 the step stayed in the same image
    cell. If it's +-(N-1), the chain wrapped through a boundary, so the image
    flag changes by -+1 in that axis.

    Returns one (ix, iy, iz) tuple per chain node. Node 0 anchored at (0,0,0).

    Used by the projection layer (`template_builder` / `builder`) to write the
    correct LAMMPS image flag column for each bead so that intra-chain bonds
    are short under unwrapped distance. The BFM topology generator itself does
    not use this -- its crosslink discovery is purely lattice-adjacency, which
    is the physically meaningful proximity at the topology stage.
    """
    images: list[tuple[int, int, int]] = [(0, 0, 0)]
    for i in range(1, len(chain)):
        old_x, old_y, old_z = flat_to_xyz(chain[i - 1], Nx, Ny)
        new_x, new_y, new_z = flat_to_xyz(chain[i], Nx, Ny)
        ix, iy, iz = images[-1]
        dx = new_x - old_x
        if dx == -(Nx - 1): ix += 1     # wrapped Nx-1 -> 0 going +x
        elif dx == (Nx - 1): ix -= 1    # wrapped 0 -> Nx-1 going -x
        dy = new_y - old_y
        if dy == -(Ny - 1): iy += 1
        elif dy == (Ny - 1): iy -= 1
        dz = new_z - old_z
        if dz == -(Nz - 1): iz += 1
        elif dz == (Nz - 1): iz -= 1
        images.append((ix, iy, iz))
    return images


def find_crosslink_candidates_distance(
    chains: list[list[int]],
    y_positions: list[int],
    Nx: int, Ny: int, Nz: int,
    lattice_scale_ang: float,
    max_distance_ang: float,
    min_intrachain_sep: int = 2,
) -> list[tuple[tuple, tuple, float]]:
    """[Optional alt path] Find Y-Y pairs whose UNWRAPPED Cartesian distance is below `max_distance_ang`.

    NOTE: this is NOT the default discovery method anymore. The BFM topology
    layer is intentionally physics-unaware, so the default (`crosslink_method=
    "adjacent"`) uses lattice-adjacency only -- two Y at adjacent lattice
    sites are at min-image distance == 1 lattice unit, which is the
    physically meaningful proximity at the topology stage.

    The phantom-bond bug that motivated this function lives in the projection
    layer, not in topology. It is fixed in `template_builder` / `builder` via
    per-chain image-cell shifts (Option A) when writing LAMMPS coords.

    This function is kept available because it can be useful when comparing
    BFM-derived networks to a Cartesian-cutoff reference, but for normal use
    you do NOT want it: it conflates topology with projection and prematurely
    drops crosslinks that the projection layer can still resolve.

    Returns a list of (pair_a, pair_b, distance_ang) sorted by ascending distance.
    """
    box = np.array([Nx, Ny, Nz], dtype=float) * lattice_scale_ang
    y_pos_sorted = sorted(set(y_positions))
    y_rank = {v: i for i, v in enumerate(y_pos_sorted)}

    # Walk each chain, tracking image-cell offset. Y atoms get unwrapped
    # positions = (lattice_xyz + image_xyz) * scale -- exactly what LAMMPS
    # will see for bond force.
    y_records: list[tuple[int, int, np.ndarray]] = []
    for ci, chain in enumerate(chains):
        images = compute_chain_images(chain, Nx, Ny, Nz)
        for ni in y_positions:
            if ni < len(chain):
                x, y, z = flat_to_xyz(chain[ni], Nx, Ny)
                ix, iy, iz = images[ni]
                unwrapped = np.array([x + ix * Nx, y + iy * Ny, z + iz * Nz],
                                     dtype=float) * lattice_scale_ang
                y_records.append((ci, ni, unwrapped))

    candidates: list[tuple[tuple, tuple, float]] = []
    for i in range(len(y_records)):
        ci1, ni1, p1 = y_records[i]
        for j in range(i + 1, len(y_records)):
            ci2, ni2, p2 = y_records[j]
            if ci1 == ci2 and abs(y_rank[ni1] - y_rank[ni2]) < min_intrachain_sep:
                continue
            d = p1 - p2                                # raw unwrapped delta, NOT min-image
            dist = float(np.linalg.norm(d))
            if dist <= max_distance_ang:
                a = (ci1, ni1); b = (ci2, ni2)
                if a <= b: candidates.append((a, b, dist))
                else:      candidates.append((b, a, dist))
    candidates.sort(key=lambda x: x[2])
    return candidates


def find_crosslink_candidates(chains, y_positions, Nx, Ny, Nz, min_intrachain_sep=2):
    """Find Y-Y pairs at 6-adjacent lattice sites (default discovery method).

    Two Y atoms whose chain-walk landed them on lattice sites that are 6-adjacent
    on the periodic cubic lattice. This is the physically meaningful proximity
    test at the topology stage: BFM is physics-unaware, so distance-1 on the
    lattice IS the discovery rule.

    Image-flag mismatches (caused by chains reaching the same lattice site via
    different wrap counts) are a projection-stage concern handled in
    `template_builder` / `builder` via per-chain image shifts -- not here.

    `min_intrachain_sep` is the minimum Y-index gap allowed for intra-chain pairs
    (use a large value like 999 to disable intra-chain crosslinks).
    """
    y_pos_sorted = sorted(set(y_positions))
    y_rank = {v: i for i, v in enumerate(y_pos_sorted)}

    y_spatial: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for ci, chain in enumerate(chains):
        for ni in y_positions:
            if ni < len(chain):
                y_spatial[chain[ni]].append((ci, y_rank[ni]))

    candidates: list[tuple] = []
    seen: set[tuple] = set()
    for flat_idx, occupants in y_spatial.items():
        x, y, z = flat_to_xyz(flat_idx, Nx, Ny)
        for nb_flat in get_neighbors(x, y, z, Nx, Ny, Nz):
            if nb_flat not in y_spatial:
                continue
            for (ci1, yr1) in occupants:
                for (ci2, yr2) in y_spatial[nb_flat]:
                    if ci1 == ci2 and abs(yr1 - yr2) < min_intrachain_sep:
                        continue
                    a = (ci1, y_pos_sorted[yr1])
                    b = (ci2, y_pos_sorted[yr2])
                    pair = (min(a, b), max(a, b))
                    if pair not in seen:
                        seen.add(pair)
                        candidates.append(pair)
    return candidates


# Snapshot scheduling and crosslink loop ====================================

def _make_snapshot(label, conv, chains, y_positions, reactions, Nx, Ny, Nz,
                   reaction_distances_ang=None, gel_point_ang=None):
    out = {
        "label": label,
        "conv": round(conv, 6),
        "chains": [list(c) for c in chains],
        "crosslinker_positions": sorted(y_positions),
        "reactions": [
            [[r[0][0], r[0][1]], [r[1][0], r[1][1]]]
            for r in reactions
        ],
        "Nx": Nx,
        "Ny": Ny,
        "Nz": Nz,
    }
    if reaction_distances_ang is not None:
        out["reaction_distances_ang"] = [round(d, 4) for d in reaction_distances_ang]
    if gel_point_ang is not None:
        out["gel_point_distance_ang"] = round(gel_point_ang, 4)
    return out


def apply_crosslinks_distance_based(
    chains: list[list[int]],
    y_positions: list[int],
    Nx: int, Ny: int, Nz: int,
    *,
    lattice_scale_ang: float,
    max_distance_ang: float = 6.0,
    n_extra_snapshots: int = 4,
    snapshot_delta_conv: float = 0.05,
    min_intrachain_sep: int = 2,
    pre_gel_conversions: list[float] | None = None,
    rng=None,
) -> list[dict]:
    """Distance-based crosslink discovery + percolation tracking (v39 path).

    Same algorithmic structure as `apply_crosslinks_with_snapshots`, but
    candidate ORDERING is by ascending min-image Cartesian distance instead of
    random shuffle of lattice-adjacent pairs. Each TYR can crosslink at most
    once (one phenol ring -> one bond). Snapshots are saved at the gel point
    (first percolation) and at conversion increments thereafter; the
    `gel_point_distance_ang` is recorded so callers can know "the network gels
    once crosslinks shorter than X A are formed".

    For topro/CHARMM (small lattice_scale), this gives near-identical results
    to the legacy lattice-adjacent scheme. For MARTINI it gives the only
    physically meaningful answer.
    """
    if rng is None:
        rng = np.random.default_rng()

    n_chains = len(chains)
    total_y = n_chains * len(y_positions)

    candidates = find_crosslink_candidates_distance(
        chains, y_positions, Nx, Ny, Nz,
        lattice_scale_ang=lattice_scale_ang,
        max_distance_ang=max_distance_ang,
        min_intrachain_sep=min_intrachain_sep,
    )

    reacted_y: set[tuple[int, int]] = set()
    uf = UnionFind(n_chains)
    reactions: list[tuple] = []
    reaction_distances: list[float] = []
    snapshots: list[dict] = []

    gel_found = False
    gel_point_distance: float | None = None
    next_snap_conv = 0.0
    chains_work = [list(c) for c in chains]

    # Pre-gel snapshot targets (sorted, consumed in order as conv crosses each).
    pre_gel_targets = sorted(pre_gel_conversions) if pre_gel_conversions else []
    pre_gel_idx = 0

    for (a, b, dist) in candidates:
        ci1, ni1 = a; ci2, ni2 = b
        if (ci1, ni1) in reacted_y or (ci2, ni2) in reacted_y:
            continue
        reacted_y.add((ci1, ni1)); reacted_y.add((ci2, ni2))
        reactions.append((a, b)); reaction_distances.append(dist)
        merged_flat = chains_work[ci1][ni1]
        chains_work[ci2][ni2] = merged_flat
        if ci1 != ci2:
            uf.union(ci1, ci2)
        conv = len(reacted_y) / total_y

        # Pre-gel: emit a snapshot the first time conv crosses each target.
        while (pre_gel_idx < len(pre_gel_targets)
               and conv >= pre_gel_targets[pre_gel_idx]
               and not gel_found):
            target = pre_gel_targets[pre_gel_idx]
            label = f"pre_gel_conv{int(round(target*1000)):04d}"  # e.g. pre_gel_conv0250
            snapshots.append(_make_snapshot(
                label, conv, chains_work, y_positions, reactions,
                Nx, Ny, Nz, reaction_distances, None,
            ))
            pre_gel_idx += 1

        if not gel_found and uf.n_components() == 1:
            gel_found = True
            gel_point_distance = dist
            next_snap_conv = conv + snapshot_delta_conv
            snapshots.append(_make_snapshot(
                "gel_point", conv, chains_work, y_positions, reactions,
                Nx, Ny, Nz, reaction_distances, gel_point_distance,
            ))
        if gel_found and len(snapshots) <= n_extra_snapshots and conv >= next_snap_conv:
            label = f"post_gel_{len(snapshots)}"
            snapshots.append(_make_snapshot(
                label, conv, chains_work, y_positions, reactions,
                Nx, Ny, Nz, reaction_distances, gel_point_distance,
            ))
            next_snap_conv = conv + snapshot_delta_conv
        if len(snapshots) > n_extra_snapshots:
            break

    if not gel_found:
        conv = len(reacted_y) / total_y
        snapshots.append(_make_snapshot(
            "no_gel", conv, chains_work, y_positions, reactions,
            Nx, Ny, Nz, reaction_distances, None,
        ))
        warnings.warn(
            f"Gel point not reached at max_distance_ang={max_distance_ang} "
            f"(max conv = {conv:.3f}, n crosslinks = {len(reactions)}). "
            "Consider widening max_distance_ang, increasing n_chains, or "
            "lowering min_intrachain_sep."
        )
    return snapshots


def _lattice_image_delta(
    flat_a: int, flat_b: int, Nx: int, Ny: int, Nz: int,
) -> tuple[int, int, int]:
    """Integer image-flag delta along the BFM bond from flat_a to flat_b.

    Returns ``(dix, diy, diz)`` such that, in an unwrapped lattice
    frame, the position of ``flat_b`` equals the position of ``flat_a``
    plus the BFM bond vector. The wrapped delta is ``pos_b - pos_a``;
    the integer ``n_boxes = round(delta / N)`` says how many full
    periodic shifts the wrap step crossed; ``image_delta = -n_boxes``
    is the per-axis image-flag increment that makes the bond
    minimum-image when applied as ``image[b] = image[a] + image_delta``.

    For BFM lattice-adjacent pairs (`|diff| in {2, √5, √6, 3, √10}` and
    always small relative to N), this returns 0 along an axis unless
    the bond straddles that face (then it returns ±1).
    """
    xa, ya, za = flat_to_xyz(flat_a, Nx, Ny)
    xb, yb, zb = flat_to_xyz(flat_b, Nx, Ny)
    dx, dy, dz = xb - xa, yb - ya, zb - za
    nbx = int(round(dx / Nx))
    nby = int(round(dy / Ny))
    nbz = int(round(dz / Nz))
    return (-nbx, -nby, -nbz)


def _compute_chain_node_images(
    chain_flat: list[int], Nx: int, Ny: int, Nz: int,
) -> list[tuple[int, int, int]]:
    """Walk a chain forward, accumulating per-node image flags relative
    to the chain's own root (node 0 at image ``(0, 0, 0)``).

    Each consecutive (n_i, n_{i+1}) bond contributes its lattice
    image-delta to the accumulating image flag. The result is a list
    of length ``len(chain_flat)`` giving every chain node's image flag
    in the chain's own unwrapped lattice frame.
    """
    imgs: list[tuple[int, int, int]] = [(0, 0, 0)]
    for i in range(1, len(chain_flat)):
        d = _lattice_image_delta(chain_flat[i - 1], chain_flat[i], Nx, Ny, Nz)
        prev = imgs[-1]
        imgs.append((prev[0] + d[0], prev[1] + d[1], prev[2] + d[2]))
    return imgs


def apply_crosslinks_winding_safe(
    chains,
    y_positions,
    Nx,
    Ny,
    Nz,
    n_extra_snapshots=4,
    snapshot_delta_conv=0.05,
    min_intrachain_sep=2,
    rng=None,
    target_conversion: float | None = None,
):
    """Lattice-adjacent crosslinking with winding-cycle rejection.

    Same outer loop as `apply_crosslinks_with_snapshots` (random shuffle
    of lattice-adjacent TYR pairs, merge sites, snapshot at gel point +
    N post-gel snapshots), with one additional acceptance gate: a
    candidate crosslink is REJECTED if it would close a cycle in the
    chain-plus-crosslink bond graph whose accumulated lattice image-
    delta around the loop is non-zero (a winding cycle around the
    periodic box). Rejected candidates are skipped; the loop continues
    until ``target_conversion`` is reached or the candidate list is
    exhausted.

    The winding check is incremental: every (chain_idx, node_idx) atom
    carries an image flag in the unwrapped lattice frame of its
    union-find component. When merging two components via a crosslink,
    the entire component's image flags are rebased so that the new
    crosslink is minimum-image by construction. When closing a cycle
    within a single component, the image flags already give the
    BFS-implied delta along the existing spanning tree; that delta is
    compared to the crosslink's own lattice image-delta. Match =
    trivial cycle (accept); mismatch = winding cycle (reject).

    This is the BFM-level analogue of the writer's priority-MST
    drop logic, applied at reaction time so the dropped crosslinks
    are simply never formed, and replaced with non-winding candidates
    from the remaining pool until the target conversion is reached.

    Notes:
      * ``target_conversion`` defaults to ``None`` = "exhaust the
        candidate list, take whatever conversion the topology allows".
        Pass an explicit value (e.g. 0.25) to stop early once the
        target is met. The gel-point snapshot is still emitted on the
        first reaction that achieves single-component connectivity.
      * Statistics on rejections are written to the snapshot config
        so callers can audit how often the topology was winding-
        constrained.
    """
    if rng is None:
        rng = np.random.default_rng()

    n_chains = len(chains)
    total_y = n_chains * len(y_positions)

    candidates = find_crosslink_candidates(
        chains, y_positions, Nx, Ny, Nz, min_intrachain_sep
    )
    rng.shuffle(candidates)

    # Per-atom image flags in the unwrapped lattice frame. Start by
    # walking each chain independently (each chain is its own component).
    node_image: dict[tuple[int, int], tuple[int, int, int]] = {}
    for ci, chain_flat in enumerate(chains):
        imgs = _compute_chain_node_images(chain_flat, Nx, Ny, Nz)
        for ni, img in enumerate(imgs):
            node_image[(ci, ni)] = img

    # Union-find over CHAIN indices (each chain is internally connected
    # via backbone bonds; the only inter-chain bonds are crosslinks).
    uf = UnionFind(n_chains)
    # Reverse map from union-find root -> list of chain indices in that
    # component. Initialised lazily on first merge.
    component_members: dict[int, list[int]] = {ci: [ci] for ci in range(n_chains)}

    reacted_y: set[tuple[int, int]] = set()
    reactions: list[tuple] = []
    snapshots: list[dict] = []
    n_rejected_winding = 0

    gel_found = False
    next_snap_conv = 0.0
    chains_work = [list(c) for c in chains]

    for pair in candidates:
        (ci1, ni1), (ci2, ni2) = pair
        if (ci1, ni1) in reacted_y or (ci2, ni2) in reacted_y:
            continue

        # Crosslink's own lattice image-delta (a -> b convention).
        xl_delta = _lattice_image_delta(
            chains_work[ci1][ni1], chains_work[ci2][ni2], Nx, Ny, Nz,
        )

        ra, rb = uf.find(ci1), uf.find(ci2)
        if ra != rb:
            # Different components: merge. The crosslink edge becomes
            # a tree edge by construction. Shift component B's image
            # flags so that image[ni2] = image[ni1] + xl_delta.
            tgt = (
                node_image[(ci1, ni1)][0] + xl_delta[0],
                node_image[(ci1, ni1)][1] + xl_delta[1],
                node_image[(ci1, ni1)][2] + xl_delta[2],
            )
            cur = node_image[(ci2, ni2)]
            shift = (tgt[0] - cur[0], tgt[1] - cur[1], tgt[2] - cur[2])
            if shift != (0, 0, 0):
                for ci_shift in component_members[rb]:
                    for ni_shift in range(len(chains_work[ci_shift])):
                        old = node_image[(ci_shift, ni_shift)]
                        node_image[(ci_shift, ni_shift)] = (
                            old[0] + shift[0], old[1] + shift[1], old[2] + shift[2],
                        )
            # Merge component membership lists.
            merged_members = component_members[ra] + component_members[rb]
            del component_members[ra]
            del component_members[rb]
            uf.union(ci1, ci2)
            component_members[uf.find(ci1)] = merged_members
        else:
            # Same component: this crosslink closes a cycle. Compute
            # the BFS-implied image-delta along the existing tree and
            # compare to the crosslink's own lattice image-delta.
            ia = node_image[(ci1, ni1)]
            ib = node_image[(ci2, ni2)]
            tree_delta = (ib[0] - ia[0], ib[1] - ia[1], ib[2] - ia[2])
            if tree_delta != xl_delta:
                # Winding cycle — reject this crosslink. Try the next
                # candidate. The current TYRs remain unreacted and may
                # still pair with other partners later in the loop.
                n_rejected_winding += 1
                continue
            # Trivial cycle (same image-delta both ways): accept.
            # The crosslink is structurally redundant but doesn't wind.

        reacted_y.add((ci1, ni1))
        reacted_y.add((ci2, ni2))
        reactions.append(pair)
        merged_flat = chains_work[ci1][ni1]
        chains_work[ci2][ni2] = merged_flat
        conv = len(reacted_y) / total_y

        if not gel_found and uf.n_components() == 1:
            gel_found = True
            next_snap_conv = conv + snapshot_delta_conv
            snapshots.append(
                _make_snapshot("gel_point", conv, chains_work, y_positions,
                               reactions, Nx, Ny, Nz)
            )
        if gel_found and len(snapshots) <= n_extra_snapshots:
            if conv >= next_snap_conv:
                label = f"post_gel_{len(snapshots)}"
                snapshots.append(
                    _make_snapshot(label, conv, chains_work, y_positions,
                                   reactions, Nx, Ny, Nz)
                )
                next_snap_conv = conv + snapshot_delta_conv
        if target_conversion is not None and conv >= target_conversion:
            break
        if len(snapshots) > n_extra_snapshots:
            break

    if not gel_found:
        conv = len(reacted_y) / total_y
        snapshots.append(
            _make_snapshot("no_gel", conv, chains_work, y_positions,
                           reactions, Nx, Ny, Nz)
        )
        warnings.warn(
            f"Gel point not reached with winding-safe crosslinking "
            f"(max conv = {conv:.3f}, accepted = {len(reactions)}, "
            f"rejected_winding = {n_rejected_winding}). Consider more "
            f"equil_steps, larger n_chains, or smaller min_intrachain_sep."
        )

    # Annotate snapshots with rejection statistics so callers can audit.
    for snap in snapshots:
        snap["n_rejected_winding"] = n_rejected_winding

    return snapshots


def apply_crosslinks_with_snapshots(
    chains,
    y_positions,
    Nx,
    Ny,
    Nz,
    n_extra_snapshots=4,
    snapshot_delta_conv=0.05,
    min_intrachain_sep=2,
    rng=None,
):
    """Incrementally crosslink the network, saving snapshots at gel point and beyond."""
    if rng is None:
        rng = np.random.default_rng()

    n_chains = len(chains)
    total_y = n_chains * len(y_positions)

    candidates = find_crosslink_candidates(
        chains, y_positions, Nx, Ny, Nz, min_intrachain_sep
    )
    rng.shuffle(candidates)

    reacted_y: set[tuple[int, int]] = set()
    uf = UnionFind(n_chains)
    reactions: list[tuple] = []
    snapshots: list[dict] = []

    gel_found = False
    next_snap_conv = 0.0
    chains_work = [list(c) for c in chains]

    for pair in candidates:
        (ci1, ni1), (ci2, ni2) = pair
        if (ci1, ni1) in reacted_y or (ci2, ni2) in reacted_y:
            continue
        reacted_y.add((ci1, ni1))
        reacted_y.add((ci2, ni2))
        reactions.append(pair)
        merged_flat = chains_work[ci1][ni1]
        chains_work[ci2][ni2] = merged_flat
        if ci1 != ci2:
            uf.union(ci1, ci2)
        conv = len(reacted_y) / total_y
        if not gel_found and uf.n_components() == 1:
            gel_found = True
            next_snap_conv = conv + snapshot_delta_conv
            snapshots.append(
                _make_snapshot("gel_point", conv, chains_work, y_positions, reactions, Nx, Ny, Nz)
            )
        if gel_found and len(snapshots) <= n_extra_snapshots:
            if conv >= next_snap_conv:
                label = f"post_gel_{len(snapshots)}"
                snapshots.append(
                    _make_snapshot(label, conv, chains_work, y_positions, reactions, Nx, Ny, Nz)
                )
                next_snap_conv = conv + snapshot_delta_conv
        if len(snapshots) > n_extra_snapshots:
            break

    if not gel_found:
        conv = len(reacted_y) / total_y
        snapshots.append(
            _make_snapshot("no_gel", conv, chains_work, y_positions, reactions, Nx, Ny, Nz)
        )
        warnings.warn(
            f"Gel point not reached (max conv = {conv:.3f}). "
            "Consider more equil_steps, larger n_chains, or smaller min_intrachain_sep."
        )
    return snapshots


# Main entry point ==========================================================

def generate_topology(
    n_chains: int = 16,
    n_repeats: int = 12,
    segs_per_block: int = 2,
    y_offset_in_block: int | None = None,
    target_packing: float = 0.45,
    equil_steps: int = 100_000,
    n_extra_snapshots: int = 4,
    snapshot_delta_conv: float = 0.05,
    min_intrachain_sep: int = 2,
    seed: int | None = None,
    verbose: bool = True,
    *,
    lattice_scale_ang: float | None = None,
    max_crosslink_distance_ang: float | None = None,
    crosslink_method: str = "adjacent",
    pre_gel_conversions: list[float] | None = None,
) -> dict:
    """Generate a 6-neighbour cubic-lattice crosslinked network topology.

    Returns a dict ``{"config": {...}, "snapshots": [...]}`` whose JSON form is
    interchangeable with topro's `topo_*.json` topology files (preserves all
    snapshot keys and the snapshot label strings).

    ``pre_gel_conversions``: optional list of conversion fractions in [0,1)
    at which to also emit a snapshot *before* the gel point (only the
    distance-based crosslink method currently honours this). Useful when
    studying sub-percolated networks at controlled conversion. Snapshots
    are labelled ``pre_gel_conv{0250,0500,...}`` (1000x the value).
    """
    rng = np.random.default_rng(seed)

    if y_offset_in_block is None:
        y_offset_in_block = 0 if segs_per_block <= 2 else 1

    n_nodes = n_repeats * segs_per_block + 1
    y_positions = get_y_positions(n_repeats, segs_per_block, y_offset_in_block)

    end_nodes = {0, n_nodes - 1}
    bad = [y for y in y_positions if y in end_nodes]
    if bad:
        raise ValueError(
            f"Y positions {bad} collide with chain End nodes "
            f"(n_nodes={n_nodes}). Adjust y_offset_in_block."
        )

    L = compute_lattice_size(n_chains, n_nodes, target_packing)
    Nx = Ny = Nz = L
    actual_packing = n_chains * n_nodes / L ** 3

    if verbose:
        print(
            f"[BFM] {n_chains} chains x {n_nodes} nodes | "
            f"lattice {L}^3 | packing {actual_packing:.3f}"
        )
        print(f"[BFM] Y positions (per chain): {y_positions}")

    if verbose:
        print("[BFM] Placing chains ...")
    chains = place_chains(n_chains, n_nodes, Nx, Ny, Nz, rng)
    if verbose:
        print(f"[BFM]   {n_chains} chains placed.")

    if equil_steps > 0:
        if verbose:
            print(f"[BFM] Equilibrating ({equil_steps:,} MC steps) ...")
        acc = equilibrate(chains, Nx, Ny, Nz, equil_steps, rng,
                          report_interval=equil_steps // 5 if verbose else 0)
        if verbose:
            print(f"[BFM]   Acceptance rate: {acc:.3f}")

    if verbose:
        if crosslink_method == "none":
            method_label = "disabled (uncrosslinked melt)"
        elif crosslink_method == "distance":
            method_label = "distance-based (alt)"
        elif crosslink_method == "winding_safe":
            method_label = "lattice-adjacent + winding-cycle rejection"
        else:
            method_label = "lattice-adjacent"
        print(f"[BFM] Crosslinking ({method_label}) ...")

    if crosslink_method == "none":
        # Uncrosslinked melt: emit the equilibrated chains as a single conv=0
        # snapshot with no reactions. No Y node is ever merged onto another
        # chain's lattice site, so downstream builders see 4 (or n_chains)
        # independent molecules. This is the starting state for in-situ
        # crosslinking (e.g. LAMMPS `fix bond/react`), as opposed to the
        # a-priori crosslinked snapshots the other methods produce.
        snapshots = [_make_snapshot(
            "uncrosslinked", 0.0, chains, y_positions, [], Nx, Ny, Nz,
        )]
    elif crosslink_method == "distance":
        if lattice_scale_ang is None:
            lattice_scale_ang = 27.0  # MARTINI default; pass explicitly for atomistic CHARMM
            if verbose:
                print(f"[BFM]   lattice_scale_ang not given -> defaulting to {lattice_scale_ang} A "
                      f"(MARTINI scale; pass an explicit value for atomistic CHARMM)")
        if max_crosslink_distance_ang is None:
            max_crosslink_distance_ang = 1.05 * lattice_scale_ang
            if verbose:
                print(f"[BFM]   max_crosslink_distance_ang not given -> defaulting to "
                      f"{max_crosslink_distance_ang:.2f} A (= 1.05 x lattice_scale; captures the "
                      f"6 nearest-neighbour lattice sites, same as legacy lattice-adjacency)")
        snapshots = apply_crosslinks_distance_based(
            chains, y_positions, Nx, Ny, Nz,
            lattice_scale_ang=lattice_scale_ang,
            max_distance_ang=max_crosslink_distance_ang,
            n_extra_snapshots=n_extra_snapshots,
            snapshot_delta_conv=snapshot_delta_conv,
            min_intrachain_sep=min_intrachain_sep,
            pre_gel_conversions=pre_gel_conversions,
            rng=rng,
        )
    elif crosslink_method == "winding_safe":
        snapshots = apply_crosslinks_winding_safe(
            chains, y_positions, Nx, Ny, Nz,
            n_extra_snapshots=n_extra_snapshots,
            snapshot_delta_conv=snapshot_delta_conv,
            min_intrachain_sep=min_intrachain_sep,
            rng=rng,
        )
    else:
        snapshots = apply_crosslinks_with_snapshots(
            chains, y_positions, Nx, Ny, Nz,
            n_extra_snapshots=n_extra_snapshots,
            snapshot_delta_conv=snapshot_delta_conv,
            min_intrachain_sep=min_intrachain_sep,
            rng=rng,
        )

    if verbose:
        for snap in snapshots:
            print(f"[BFM]   {snap['label']:<22}  conv={snap['conv']:.4f}")

    config = {
        "n_chains": n_chains,
        "n_repeats": n_repeats,
        "segs_per_block": segs_per_block,
        "y_offset_in_block": y_offset_in_block,
        "n_nodes_per_chain": n_nodes,
        "y_positions": y_positions,
        "Nx": Nx, "Ny": Ny, "Nz": Nz,
        "target_packing": target_packing,
        "actual_packing": actual_packing,
        "equil_steps": equil_steps,
        "n_extra_snapshots": n_extra_snapshots,
        "snapshot_delta_conv": snapshot_delta_conv,
        "min_intrachain_sep": min_intrachain_sep,
        "seed": seed,
        "lattice_scale_ang": lattice_scale_ang,
        "max_crosslink_distance_ang": max_crosslink_distance_ang,
        "crosslink_method": crosslink_method,
    }
    return {"config": config, "snapshots": snapshots}
