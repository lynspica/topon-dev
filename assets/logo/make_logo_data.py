"""Build LAMMPS data for the topon logo: "topon" emerging as coloured beads in a
real topon-generated network slab.

Topology is genuinely topon's own generator (PythonTopologyGenerator, the
"strict sculpting" algorithm) on a rectangular SC slab -- not a hand-drawn
cartoon. Beads are laid along each network strand; a bead is ACCENT if it falls
inside a "topon" glyph outline, QUIET otherwise. The enclosed counters of the
two o's and the p are carved out as VACANCIES so the letters read as letters.

Variants (--variant):
  clean      carved counters, teal letters, full lattice background
  copolymer  + letters split A/B (teal / coral) -- signals the copolymer generator
  entangled  + REAL topon entanglements on the letter strands: topon's own
             find_crossing_candidates + select_entanglements pick the sites,
             realised with its KinkParams (two nearest parallel strands neck
             together, cross with `overshoot`, one over / one under by z_amp)
  organic    + background thinned by random junction vacancies (material, not graph paper)

Types:  1 quiet strand   2 accent-A strand   3 quiet junction   4 accent-A junction
        5 accent-B strand  6 accent-B junction  7 bead carrying a kink

Usage:  python make_logo_data.py --variant clean [-o out.data]
"""
from __future__ import annotations

import argparse
from collections import deque
from pathlib import Path
import sys

import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.path import Path as MplPath
from matplotlib.textpath import TextPath

sys.path.insert(0, r"C:/Users/ahmet/OneDrive - Northwestern University/DOW-Ahmet/topon")
from topon.topology.generator_python import PythonTopologyGenerator  # noqa: E402

# --- design knobs -------------------------------------------------------------
NX, NY, NZ = 84, 26, 2         # rectangular slab: wide banner, thin in z
BEADS_PER_STRAND = 4
TEXT = "topon"
FONT_FAMILY, FONT_WEIGHT = "DejaVu Sans", "bold"
TEXT_MARGIN_X = 0.07
COUNTER_RES = 8                # raster cells per lattice unit for counter finding
ORGANIC_VACANCY = 0.40         # fraction of background junctions vacated (organic)
N_ENTANGLEMENTS = 75           # kink sites requested from topon's entanglement stage
# Slight degree heterogeneity: the raw 84x26x2 SC lattice is a PERFECT 5-regular
# graph (NZ=2 collapses one z-bond), which renders as graph paper. topon's own
# strict-sculpting stage 4 thins edges to this target, giving a real degree
# spread (~5/4/3/2) instead of a uniform lattice. 10920 -> 9500 ~= 13% thinned.
EDGE_TARGET = 9500
SEED = 7

# NOTE: glyphs come from matplotlib's TextPath (exact vector outlines +
# contains_points), NOT PIL. PIL's FreeType is broken on this interpreter --
# ImageFont.truetype() loads any font and reports sane bboxes but rasterises
# .notdef ("tofu") boxes for every character, including its own bundled DejaVu.


class _Cfg:
    lattice_size = f"{NX}x{NY}x{NZ}"
    lattice_type = "SC"
    max_functionality = 6
    degree_distribution = f"e:{EDGE_TARGET}"
    periodicity = "111"


def build_text_contours():
    """TEXT scaled/centred onto the slab face, as a list of per-contour Paths.

    Returned as SEPARATE contours because matplotlib's compound-path
    `contains_points` fills glyph counters: an 'o' has two correctly-wound
    contours (outer CW, counter CCW) yet contains_point() is True at the hole's
    centre. Testing each contour and XOR-ing (even-odd rule) excludes the holes
    exactly, which is what makes the o/o/p counters carvable.
    """
    tp = TextPath((0, 0), TEXT, size=1.0,
                  prop=FontProperties(family=FONT_FAMILY, weight=FONT_WEIGHT))
    v = tp.vertices
    x0, x1 = v[:, 0].min(), v[:, 0].max()
    y0, y1 = v[:, 1].min(), v[:, 1].max()
    s = min(NX * (1 - 2 * TEXT_MARGIN_X) / (x1 - x0), NY * 0.72 / (y1 - y0))
    verts = (v - [x0, y0]) * s
    verts[:, 0] += (NX - (x1 - x0) * s) / 2
    verts[:, 1] += (NY - (y1 - y0) * s) / 2
    scaled = MplPath(verts, tp.codes)
    return [MplPath(p) for p in scaled.to_polygons(closed_only=False)]


def glyph_contains(contours, pts) -> np.ndarray:
    """Even-odd point-in-glyph over the contour list (holes excluded)."""
    pts = np.atleast_2d(pts)
    inside = np.zeros(len(pts), bool)
    for c in contours:
        inside ^= c.contains_points(pts)
    return inside


def enclosed_counters(contours):
    """Boolean grid of the glyphs' ENCLOSED counters -- the holes in o, o, p.

    Rasterise the glyph fill, flood-fill the exterior inward from the border
    through non-ink cells, and whatever non-ink is left unreached is enclosed by
    ink: exactly the counters. Open bays (under the n arch, around the t) reach
    the border and are correctly NOT counted.

    Returns (grid[row=y, col=x], res) with res cells per lattice unit.
    """
    res = COUNTER_RES
    W, H = NX * res, NY * res
    xs = (np.arange(W) + 0.5) / res
    ys = (np.arange(H) + 0.5) / res
    X, Y = np.meshgrid(xs, ys)
    ink = glyph_contains(contours, np.column_stack([X.ravel(), Y.ravel()])).reshape(H, W)

    outside = np.zeros_like(ink)
    dq = deque()
    for c in range(W):
        for r in (0, H - 1):
            if not ink[r, c] and not outside[r, c]:
                outside[r, c] = True; dq.append((r, c))
    for r in range(H):
        for c in (0, W - 1):
            if not ink[r, c] and not outside[r, c]:
                outside[r, c] = True; dq.append((r, c))
    while dq:
        r, c = dq.popleft()
        for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nr, nc = r + dr, c + dc
            if 0 <= nr < H and 0 <= nc < W and not ink[nr, nc] and not outside[nr, nc]:
                outside[nr, nc] = True; dq.append((nr, nc))
    return (~ink) & (~outside), res


def counter_blobs(counter: np.ndarray, res: int):
    """Connected components of the counter grid -> list of (cx, cy) centroids
    in slab coords, largest first."""
    seen = np.zeros_like(counter)
    blobs = []
    H, W = counter.shape
    for r0 in range(H):
        for c0 in range(W):
            if counter[r0, c0] and not seen[r0, c0]:
                dq = deque([(r0, c0)]); seen[r0, c0] = True; cells = []
                while dq:
                    r, c = dq.popleft(); cells.append((r, c))
                    for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        nr, nc = r + dr, c + dc
                        if 0 <= nr < H and 0 <= nc < W and counter[nr, nc] and not seen[nr, nc]:
                            seen[nr, nc] = True; dq.append((nr, nc))
                a = np.array(cells, dtype=float)
                blobs.append((len(cells), (a[:, 1].mean() + 0.5) / res, (a[:, 0].mean() + 0.5) / res))
    blobs.sort(reverse=True)
    return blobs


def in_grid(grid, res, xy):
    """Sample a raster grid at slab coords (vectorised)."""
    H, W = grid.shape
    c = np.clip((xy[:, 0] * res).astype(int), 0, W - 1)
    r = np.clip((xy[:, 1] * res).astype(int), 0, H - 1)
    return grid[r, c]


def apply_topon_kinks(G, pos, xyz, accent, edge_beads, n_target):
    """Entangle the LETTER strands using topon's own entanglement stage, then
    realise each entanglement as topon's kink geometry.

    This is the real thing, not a decorative link:
      * `find_crossing_candidates` pairs each strand with its NEAREST DISJOINT
        neighbour -- i.e. the closest roughly-parallel strand;
      * `select_entanglements` picks the kink sites (edge- and location-exclusive);
      * the pair is then realised with topon's own `KinkParams`: each strand is
        pulled toward its partner by (sep/2)(1+overshoot) so the two cross and
        OVERSHOOT past each other, one displaced +z_amp and the other -z_amp so
        they pass over/under, all under a Gaussian envelope of width `sigma`
        along the strand so the ends stay pinned at their junctions.

    Restricted to accent strands: candidate search is O(E^2) in pure Python, so
    the full 10.9k-edge lattice is out of reach -- and the letters are where the
    entanglements are wanted anyway.

    Returns (n_kinks, kink_bead_mask).
    """
    import networkx as nx
    from topon.assignment import entanglements as ent
    from topon.config.schema import EntanglementsConfig

    dims = np.array([NX, NY, NZ], float)
    acc_edges = [(u, v, chain) for (u, v, chain) in edge_beads
                 if accent[np.array(chain[1:-1]) - 1].mean() > 0.8]
    MG = nx.MultiGraph()
    for u, v, _ in acc_edges:
        MG.add_node(u, pos=tuple(pos[u])); MG.add_node(v, pos=tuple(pos[v]))
        MG.add_edge(u, v)
    print(f"[entangled] {len(acc_edges)} letter strands -> topon entanglement stage")

    cands = ent.find_crossing_candidates(MG, dims)
    kp = EntanglementsConfig(enabled=True).kink_params

    chain_of = {(min(u, v), max(u, v)): chain for u, v, chain in acc_edges}

    def midpoint(e):
        u, v = e[0], e[1]
        d = pos[v] - pos[u]
        d -= dims * np.round(d / dims)
        return pos[u] + 0.5 * d

    def sep_vec(c):
        d = midpoint(c[1]) - midpoint(c[0])
        return d - dims * np.round(d / dims)

    # (a) IN-PLANE pairs only. A pair stacked along z gets pulled along z --
    # which is the viewing axis -- so its kink collapses to a dot on screen
    # (that was the stray single beads). Keeping |dz|~0 means the pull is
    # in-plane and z_amp puts one strand over and one under: a visible X.
    inplane = [c for c in cands if abs(sep_vec(c)[2]) < 0.25]
    print(f"[entangled] {len(cands)} candidates -> {len(inplane)} in-plane "
          f"(rest are z-stacked and would be invisible edge-on)")
    if not inplane:
        return 0, np.zeros(len(xyz), bool)

    # (b) DELIBERATE placement: farthest-point sampling over candidate centres
    # so the kinks are spread evenly across the wordmark instead of clumping
    # where topon's random draw happens to land. (Candidate-finding and the
    # kink geometry below are still topon's; only the choice of WHICH sites is
    # logo-specific, because an even spread is a design requirement.)
    centres = np.array([0.5 * (midpoint(c[0]) + midpoint(c[1]))[:2] for c in inplane])
    order, d2 = [], np.full(len(inplane), np.inf)
    cur = int(np.argmin(((centres - centres.mean(0)) ** 2).sum(1)))
    for _ in range(min(len(inplane), 8 * n_target)):
        order.append(cur)
        d2 = np.minimum(d2, ((centres - centres[cur]) ** 2).sum(1))
        d2[cur] = -1.0
        cur = int(np.argmax(d2))
        if d2[cur] < 0:
            break

    kink_mask = np.zeros(len(xyz), bool)
    used, n_done = set(), 0
    for idx in order:
        if n_done >= n_target:
            break
        e1, e2 = inplane[idx]
        k1 = (min(e1[0], e1[1]), max(e1[0], e1[1]))
        k2 = (min(e2[0], e2[1]), max(e2[0], e2[1]))
        if k1 in used or k2 in used:          # edge exclusivity, as topon enforces
            continue
        c1, c2 = chain_of.get(k1), chain_of.get(k2)
        if c1 is None or c2 is None:
            continue
        d = sep_vec((e1, e2))
        sep = float(np.linalg.norm(d))
        if sep < 1e-6:
            continue
        dhat = d / sep
        pull = 0.5 * sep * (1.0 + kp.overshoot)     # cross AND overshoot past
        for chain, sign in ((c1, +1.0), (c2, -1.0)):
            for k, aid in enumerate(chain[1:-1]):    # interior beads; ends pinned
                t = (k + 1) / (len(chain) - 1)
                w = float(np.exp(-((t - 0.5) ** 2) / (2 * kp.sigma ** 2)))
                xyz[aid - 1] += sign * pull * w * dhat
                xyz[aid - 1][2] += sign * kp.z_amp * w
                if w > 0.3:
                    kink_mask[aid - 1] = True
        used.add(k1); used.add(k2)
        n_done += 1
    print(f"[entangled] {n_done} kinks realised, evenly spread "
          f"(overshoot={kp.overshoot}, z_amp={kp.z_amp}, sigma={kp.sigma})")
    return n_done, kink_mask


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="clean",
                    choices=["clean", "copolymer", "entangled", "organic", "organic_copoly"])
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()
    out = Path(args.out or f"topon_{args.variant}.data")
    rng = np.random.default_rng(SEED)

    # --- 1. topology from topon's own generator ---------------------------
    G = PythonTopologyGenerator(_Cfg()).generate(trials=1, max_saves=1)[0]
    pos = {n: np.array(G.nodes[n]["pos"], float) for n in G.nodes}
    print(f"[topology] {G.number_of_nodes()} junctions, {G.number_of_edges()} strands")

    global _CONTOURS
    _CONTOURS = build_text_contours()
    counter, res = enclosed_counters(_CONTOURS)
    blobs = counter_blobs(counter, res)
    print(f"[counters] {len(blobs)} enclosed voids (o, o, p): "
          + ", ".join(f"({x:.1f},{y:.1f})" for _, x, y in blobs[:4]))

    # --- 2. beads along the network ---------------------------------------
    xyz, is_junct, bonds, node_atom = [], [], [], {}
    for n in G.nodes:
        xyz.append(pos[n]); is_junct.append(True); node_atom[n] = len(xyz)

    edge_beads = []          # (u, v, [atom ids of interior beads])
    for u, v in G.edges:
        pu, pv = pos[u], pos[v]
        d = pv - pu
        if np.any(np.abs(d) > 1.5):
            continue
        chain = [node_atom[u]]
        for k in range(1, BEADS_PER_STRAND + 1):
            xyz.append(pu + (k / (BEADS_PER_STRAND + 1)) * d)
            is_junct.append(False); chain.append(len(xyz))
        chain.append(node_atom[v])
        bonds.extend(zip(chain[:-1], chain[1:]))
        edge_beads.append((u, v, chain))

    xyz = np.array(xyz); is_junct = np.array(is_junct)
    accent = glyph_contains(_CONTOURS, xyz[:, :2])
    print(f"[beads] {len(xyz)} beads, {len(bonds)} bonds, "
          f"{int(accent.sum())} accent ({100*accent.mean():.1f}%)")

    # --- 3. carve the counters as vacancies -------------------------------
    # Only QUIET beads sit in a counter (accent = in-glyph by definition), so
    # this removes exactly the background filling the o/o/p holes.
    keep = ~(in_grid(counter, res, xyz[:, :2]) & ~accent)
    print(f"[vacancy] carved {int((~keep).sum())} beads from the letter counters")

    # --- 4. variant: organic background ------------------------------------
    if args.variant in ("organic", "organic_copoly"):
        # Vacate background junctions, but only those whose strands carry no
        # accent bead -- the letters must stay intact.
        acc_nodes = set()
        for u, v, chain in edge_beads:
            if accent[np.array(chain) - 1].any():
                acc_nodes.add(u); acc_nodes.add(v)
        vac = {n for n in G.nodes
               if n not in acc_nodes and not accent[node_atom[n] - 1]
               and rng.random() < ORGANIC_VACANCY}
        drop = np.zeros(len(xyz), bool)
        for n in vac:
            drop[node_atom[n] - 1] = True
        for u, v, chain in edge_beads:
            if u in vac or v in vac:
                drop[np.array(chain[1:-1]) - 1] = True
        keep &= ~drop
        print(f"[organic] vacated {len(vac)} background junctions "
              f"({int(drop.sum())} beads)")

    # --- 5. types ----------------------------------------------------------
    types = np.where(is_junct, np.where(accent, 4, 3), np.where(accent, 2, 1))
    if args.variant in ("copolymer", "organic_copoly"):
        # A/B blocks along the chain contour: alternate by lattice cell in x,
        # so each letter shows both monomer species.
        b_block = (np.floor(xyz[:, 0] / 2.0).astype(int) % 2 == 1) & accent
        types = np.where(b_block, np.where(is_junct, 6, 5), types)
        print(f"[copolymer] {int(b_block.sum())} B-monomer accent beads")

    # --- 6. variant: real topon entanglements (kinks) -----------------------
    extra_xyz, extra_bonds = [], []
    if args.variant == "entangled":
        # No recolouring: a kink is a GEOMETRY feature and must read as one --
        # the strands neck together and cross over/under in teal like the rest
        # of the letter.
        apply_topon_kinks(G, pos, xyz, accent, edge_beads, N_ENTANGLEMENTS)

    # --- 7. reindex kept beads + emit -------------------------------------
    new_id = np.zeros(len(xyz), int)
    new_id[keep] = np.arange(1, int(keep.sum()) + 1)
    atoms = [(int(t), *p) for t, p in zip(types[keep], xyz[keep])]
    out_bonds = [(new_id[a - 1], new_id[b - 1]) for a, b in bonds
                 if keep[a - 1] and keep[b - 1]]
    off = len(atoms)
    for p in extra_xyz:
        atoms.append((7, *p))
    out_bonds.extend((a - len(xyz) + off, b - len(xyz) + off) for a, b in extra_bonds)

    n_types = 7
    pad = 1.0
    L = [f"topon logo [{args.variant}]: '{TEXT}' in a {NX}x{NY}x{NZ} SC network", "",
         f"{len(atoms)} atoms", f"{len(out_bonds)} bonds",
         f"{n_types} atom types", "1 bond types", "",
         f"{-pad:.3f} {NX+pad:.3f} xlo xhi",
         f"{-pad:.3f} {NY+pad:.3f} ylo yhi",
         f"{-pad-2:.3f} {NZ+pad+2:.3f} zlo zhi", "",
         "Masses", ""] + [f"{i} 1.0" for i in range(1, n_types + 1)] + \
        ["", "Atoms # full", ""]
    for i, (t, x, y, z) in enumerate(atoms, 1):
        L.append(f"{i} 1 {t} 0.0 {x:.4f} {y:.4f} {z:.4f}")
    L += ["", "Bonds", ""]
    for i, (a, b) in enumerate(out_bonds, 1):
        L.append(f"{i} 1 {a} {b}")
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"[OK] {len(atoms)} beads -> {out}")


if __name__ == "__main__":
    main()
