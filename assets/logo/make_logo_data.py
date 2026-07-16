"""Build a LAMMPS data file for the topon logo: the word "topon" emerging as
coloured beads inside a real topon-generated copolymer network slab.

Topology is genuinely topon's own generator (PythonTopologyGenerator, the
"strict sculpting" algorithm) on a rectangular SC slab -- not a hand-drawn
cartoon. Beads are then laid along each network strand, and a bead is painted
ACCENT if it falls inside a "topon" glyph rasterised over the slab face, QUIET
otherwise. Type ids carry that label into the render:

    type 1 = quiet strand bead      type 3 = quiet junction
    type 2 = accent strand bead     type 4 = accent junction

Usage:  python make_logo_data.py [out.data]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

sys.path.insert(0, r"C:/Users/ahmet/OneDrive - Northwestern University/DOW-Ahmet/topon")
from topon.topology.generator_python import PythonTopologyGenerator  # noqa: E402

# --- design knobs -------------------------------------------------------------
NX, NY, NZ = 84, 26, 2          # rectangular slab: wide banner, thin in z
BEADS_PER_STRAND = 4           # interior beads laid along each lattice edge
TEXT = "topon"
FONT_FAMILY, FONT_WEIGHT = "DejaVu Sans", "bold"
TEXT_MARGIN_X = 0.07           # fraction of slab width kept clear either side
SEED = 7

# NOTE: glyphs come from matplotlib's TextPath (exact vector outlines +
# contains_points), NOT PIL. PIL's FreeType is broken on this interpreter --
# ImageFont.truetype() loads any font and reports sane bboxes but rasterises
# .notdef ("tofu") boxes for every character, including its own bundled
# DejaVu. TextPath sidesteps the raster path entirely and gives an exact
# point-in-glyph test, which is what the bead labelling actually wants.


class _Cfg:
    lattice_size = f"{NX}x{NY}x{NZ}"
    lattice_type = "SC"
    max_functionality = 6
    degree_distribution = ""       # keep the full 6-regular lattice (dense, reads as a material)
    periodicity = "111"


def build_text_path():
    """TEXT as an exact vector outline, scaled/centred onto the slab face.
    Returns a matplotlib Path in slab coordinates (x in [0,NX], y in [0,NY])."""
    tp = TextPath((0, 0), TEXT, size=1.0,
                  prop=FontProperties(family=FONT_FAMILY, weight=FONT_WEIGHT))
    v = tp.vertices
    x0, x1 = v[:, 0].min(), v[:, 0].max()
    y0, y1 = v[:, 1].min(), v[:, 1].max()
    avail_w = NX * (1 - 2 * TEXT_MARGIN_X)
    avail_h = NY * 0.72
    s = min(avail_w / (x1 - x0), avail_h / (y1 - y0))       # uniform scale
    verts = (v - [x0, y0]) * s
    verts[:, 0] += (NX - (x1 - x0) * s) / 2                 # centre in x
    verts[:, 1] += (NY - (y1 - y0) * s) / 2                 # centre in y
    from matplotlib.path import Path as MplPath
    return MplPath(verts, tp.codes)


def label_beads(path, pts: np.ndarray) -> np.ndarray:
    """Boolean: is each (x, y) inside a glyph? Vectorised over all beads."""
    return path.contains_points(pts)


def main() -> None:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "topon_logo.data")

    # --- 1. topology from topon's own generator --------------------------
    gen = PythonTopologyGenerator(_Cfg())
    graphs = gen.generate(trials=1, max_saves=1)
    if not graphs:
        raise RuntimeError("generator returned no network")
    G = graphs[0]
    pos = {n: np.array(G.nodes[n]["pos"], dtype=float) for n in G.nodes}
    print(f"[topology] {G.number_of_nodes()} junctions, {G.number_of_edges()} strands "
          f"({NX}x{NY}x{NZ} SC slab)")

    path = build_text_path()
    print(f"[text] '{TEXT}' vector outline, {len(path.vertices)} vertices")

    box = np.array([NX, NY, NZ], dtype=float)
    xyz = []        # bead coords
    is_junct = []   # junction flag
    bonds = []      # (a, b) 1-based
    node_atom = {}

    # --- 2. junction beads ------------------------------------------------
    for n in G.nodes:
        p = pos[n]
        xyz.append(p)
        is_junct.append(True)
        node_atom[n] = len(xyz)

    # --- 3. strand beads along each edge (skip periodic wraps: a logo wants
    #        clean slab edges, not bonds shooting across the box) ----------
    n_skipped = 0
    for u, v in G.edges:
        pu, pv = pos[u], pos[v]
        d = pv - pu
        if np.any(np.abs(d) > 1.5):          # wrapping edge
            n_skipped += 1
            continue
        chain = [node_atom[u]]
        for k in range(1, BEADS_PER_STRAND + 1):
            xyz.append(pu + (k / (BEADS_PER_STRAND + 1)) * d)
            is_junct.append(False)
            chain.append(len(xyz))
        chain.append(node_atom[v])
        bonds.extend(zip(chain[:-1], chain[1:]))

    xyz = np.array(xyz)
    is_junct = np.array(is_junct)
    print(f"[beads] {len(xyz)} beads, {len(bonds)} bonds "
          f"({n_skipped} wrapping strands dropped)")

    # --- 4. label every bead against the glyph outline (one vector pass) ---
    accent = label_beads(path, xyz[:, :2])
    types = np.where(is_junct, np.where(accent, 4, 3), np.where(accent, 2, 1))
    atoms = [(int(t), p[0], p[1], p[2]) for t, p in zip(types, xyz)]
    print(f"[text] {int(accent.sum())} accent beads ({100*accent.mean():.1f}%)")

    # --- 4. LAMMPS data ---------------------------------------------------
    pad = 1.0
    L = [f"topon logo: '{TEXT}' in a {NX}x{NY}x{NZ} SC copolymer network", "",
         f"{len(atoms)} atoms", f"{len(bonds)} bonds", "4 atom types", "1 bond types", "",
         f"{-pad:.3f} {box[0]+pad:.3f} xlo xhi",
         f"{-pad:.3f} {box[1]+pad:.3f} ylo yhi",
         f"{-pad:.3f} {box[2]+pad:.3f} zlo zhi", "",
         "Masses", "", "1 1.0", "2 1.0", "3 1.0", "4 1.0", "", "Atoms # full", ""]
    for i, (t, x, y, z) in enumerate(atoms, 1):
        L.append(f"{i} 1 {t} 0.0 {x:.4f} {y:.4f} {z:.4f}")
    L += ["", "Bonds", ""]
    for i, (a, b) in enumerate(bonds, 1):
        L.append(f"{i} 1 {a} {b}")
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"[OK] -> {out}")


if __name__ == "__main__":
    main()
