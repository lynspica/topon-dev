"""High-quality arc movies: lattice -> non-ideal melt, densely sampled.

One builder for all three resolutions/knobs:
  cg      sculpt_250    homopolymer, single colour
  copoly  copoly_block  A/B copolymer sequence, two colours + dark junctions
  atom    atom_sculpt   DREIDING PDMS, element colours

The lattice->melt transition completes inside a single `minimize` under the
production decks, so a fixed-interval dump jump-cuts it. The movie decks
(`movie_cg.in`, `movie_atom.in`) inflate the lattice under a *ramped soft push*
with bounded dynamics (nve/limit) instead -- the transition is driven by excluded
volume, not by the bonds -- dumping every step. This script then selects ~35
frames evenly in a monotone PROGRESS metric so the pacing is even (no jump-cut),
renders them at high quality with a fixed camera, and assembles a seamless
boomerang GIF + an MP4.

Progress metric:
  disorder  bond-orientation disorder, 1 - <max axis component / |bond|>.
            0 on the axis-aligned lattice, rising through inflation AND coiling
            (bond LENGTH saturates while strands still coil). Used for CG/copoly.
  disp      mean minimum-image displacement from the lattice (frame 0). Robust
            for atomistic, where methyl H bonds give a nonzero disorder baseline.

Inputs (produced by the movie decks in <system>/04_Simulation/):
  infl.lammpstrj  soft-push inflation, dumped EVERY step
  therm.lammpstrj real-potential + thermal tail, dumped every ~15 steps
Coordinates only -- bond topology comes from the 03_Conformation .data file.

Run:  C:/v/ovito/Scripts/python.exe make_arc_movie.py [cg copoly atom]
"""
import os
import sys
from pathlib import Path

import numpy as np
import imageio.v2 as imageio
from PIL import Image
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import WrapPeriodicImagesModifier, LoadTrajectoryModifier

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from render_gallery import (node_ids, rebuild_bond_pbc, half_cell_offset, fresh,
                            element_palette, closest_entangled_pair,
                            A_RED, B_BLUE, JUNCT, ENT_GOLD, VIOLET, REST, TEAL)

OUT = HERE / "anim"
OUT.mkdir(parents=True, exist_ok=True)
FRAMESROOT = Path(os.environ.get("TOPON_GALLERY_FRAMES", HERE / "_frames"))

CONFIG = {
    "cg":     dict(system="sculpt_250",   paint="mono",     metric="disorder",
                   rad=0.11, bw=0.22),
    "copoly": dict(system="copoly_block", paint="ab",       metric="disorder",
                   rad=0.13, bw=0.26),
    # atomistic melt is high-entropy (many small atoms), so a GIF at the default
    # size/palette is heavy; trim it here (the MP4 stays full quality).
    "atom":   dict(system="atom_sculpt",  paint="elements", metric="disp",
                   rad=0.34, bw=0.0, gif_size=420, gif_colors=80),
    # Zoom into ONE entanglement as it relaxes: the two entangled chains in two
    # distinct colours, all grafts teal, the rest of the network grey. The crop
    # follows the pair's centroid every frame so it stays framed as it coils.
    "ent":    dict(system="entangled_grafted", paint="ent", metric="disorder",
                   rad=0.09, bw=0.09, zoom=True, keep=4.6),
}

RENDER = 660
GIF_SIZE = 460
GIF_COLORS = 96
GIF_MS = 91          # ~0.66x speed (was 60 ms/frame -> 91 ms is 1/0.66 slower)
MP4_SIZE = 660
FPS = 13             # ~0.66x speed (was 20)
N_TRANSITION = 30
N_LEADIN = 3
N_TAIL = 6


def progress(pipe, lt, topo, cell, kind):
    inv = np.linalg.inv(cell)
    p0 = np.array(pipe.compute(0).particles.positions)
    out = []
    for f in range(lt.source.num_frames):
        p = np.array(pipe.compute(f).particles.positions)
        if kind == "disorder":
            v = (p[topo[:, 0]] - p[topo[:, 1]]) @ inv.T
            v -= np.round(v)
            w = v @ cell.T
            n = np.linalg.norm(w, axis=1)
            g = n > 1e-9
            out.append(float(np.mean(1.0 - np.abs(w[g]).max(axis=1) / n[g])))
        else:  # disp: mean minimum-image displacement from the lattice
            d = (p - p0) @ inv.T
            d -= np.round(d)
            out.append(float(np.mean(np.linalg.norm(d @ cell.T, axis=1))))
    return np.array(out)


def select(prog):
    lo, hi = float(prog[0]), float(np.percentile(prog, 97))
    idx = [0]
    start = int(np.argmax(prog > lo + 0.02 * (hi - lo)))
    if start > 1:
        for t in np.linspace(1, start - 1, N_LEADIN, dtype=int):
            idx.append(int(t))
    for target in np.linspace(lo, hi, N_TRANSITION):
        idx.append(int(np.argmax(prog >= target)))
    tail0 = int(np.argmax(prog >= hi))
    for t in np.linspace(tail0, len(prog) - 1, N_TAIL, dtype=int):
        idx.append(int(t))
    seen, out = set(), []
    for i in idx:
        if i not in seen:
            seen.add(i)
            out.append(i)
    return out


def make_paint(cfg, ref, nodes, rad, pair=None):
    if cfg["paint"] == "elements":
        ec, er = element_palette(import_file(str(ref), atom_style="full").compute())

        def paint(frame, data):
            t = np.array(data.particles["Particle Type"])
            n = data.particles.count
            col = np.tile(np.array((0.36, 0.40, 0.45)), (n, 1))
            r = np.full(n, rad)
            for tid, c in ec.items():
                m = t == tid
                if m.any():
                    col[m], r[m] = c, er[tid] * rad / 0.42
            data.particles_.create_property("Color", data=col)
            data.particles_.create_property("Radius", data=r)
        return paint

    node_arr = np.fromiter(nodes, dtype=np.int64)

    if cfg["paint"] == "ent":
        # ONE entanglement: chain A gold, chain B violet, all grafts (type 3)
        # teal, everything else grey. Junctions larger but grey.
        a_arr = np.fromiter(pair[0], dtype=np.int64)
        b_arr = np.fromiter(pair[1], dtype=np.int64)

        def paint(frame, data):
            ident = np.array(data.particles["Particle Identifier"])
            t = np.array(data.particles["Particle Type"])
            n = data.particles.count
            col = np.tile(np.array(REST), (n, 1))
            r = np.full(n, rad)
            r[np.isin(ident, node_arr)] = rad * 2.0
            col[t == 3] = TEAL                     # side chains pop
            ma, mb = np.isin(ident, a_arr), np.isin(ident, b_arr)
            col[ma], col[mb] = ENT_GOLD, VIOLET
            r[ma | mb] = rad * 1.7                 # the two entangled strands
            data.particles_.create_property("Color", data=col)
            data.particles_.create_property("Radius", data=r)
        return paint

    def paint(frame, data):
        ident = np.array(data.particles["Particle Identifier"])
        t = np.array(data.particles["Particle Type"])
        n = data.particles.count
        isnode = np.isin(ident, node_arr)
        if cfg["paint"] == "ab":
            col = np.where((t == 1)[:, None], np.array(A_RED), np.array(B_BLUE))
            col = col.astype(float)
            col[isnode] = JUNCT                    # match the copolymer stills
        else:                                      # mono
            col = np.tile(np.array(B_BLUE), (n, 1))
        r = np.full(n, rad)
        r[isnode] = rad * 2.4
        data.particles_.create_property("Color", data=col)
        data.particles_.create_property("Radius", data=r)
    return paint


def build(tag):
    cfg = CONFIG[tag]
    ref = fresh(cfg["system"], "03_Conformation/system_relaxed.data")
    sim = ref.parent.parent / "04_Simulation"
    nodes = node_ids(ref.parent.parent)

    # arc = inflation, then the thermal tail if the deck produced one. For the
    # atomistic deck the ramped soft push alone already carries the arc (bonds
    # 0.44 -> 1.12 A AND continued coiling, mean displacement rising to 3.5 A),
    # so it has no separate thermal file.
    arc = sim / "movie_arc.lammpstrj"
    data = (sim / "infl.lammpstrj").read_bytes()
    tail = sim / "therm.lammpstrj"
    if tail.exists():
        data += tail.read_bytes()
    arc.write_bytes(data)

    pipe = import_file(str(ref), atom_style="full")
    lt = LoadTrajectoryModifier()
    lt.source.load(str(arc))
    pipe.modifiers.append(lt)

    probe = pipe.compute(0)
    cell = np.array(probe.cell[:, :3])
    topo = np.array(probe.particles.bonds.topology)
    prog = progress(pipe, lt, topo, cell, cfg["metric"])
    sel = select(prog)
    gaps = np.abs(np.diff(prog[sel]))
    print(f"[{tag}] selected {len(sel)} of {lt.source.num_frames}; "
          f"{cfg['metric']} {prog[sel[0]]:.3f} -> {prog[sel[-1]]:.3f}; "
          f"largest gap {gaps.max():.3f}", flush=True)

    # Entanglement zoom: pick the pair on the lattice, then FOLLOW it -- recentre
    # on the two chains' centroid every frame so it stays framed as it coils.
    pair = None
    if cfg.get("zoom"):
        dmin, _, pr = closest_entangled_pair(probe, nodes)
        ident0 = np.array(probe.particles["Particle Identifier"])
        pa = {int(ident0[i]) for i in pr[0][1]}
        pb = {int(ident0[i]) for i in pr[1][1]}
        pair = (pa, pb)
        focus = np.fromiter(pa | pb, dtype=np.int64)
        inv = np.linalg.inv(cell)
        boxc = 0.5 * (cell[0] + cell[1] + cell[2])
        print(f"[{tag}] focus pair {dmin:.2f} sigma apart, {len(focus)} beads",
              flush=True)

        def follow(frame, data, foc=focus, inv=inv, cell=cell, c=boxc):
            ident = np.array(data.particles["Particle Identifier"])
            p = np.array(data.particles.positions)
            fp = p[np.isin(ident, foc)]
            d = (fp - fp[0]) @ inv.T           # unwrap the focus set (may straddle)
            d -= np.round(d)
            centroid = fp[0] + (d @ cell.T).mean(0)
            data.particles_.positions_[...] = p + (c - centroid)
        pipe.modifiers.append(follow)
    else:
        off = half_cell_offset(probe, nodes)
        if np.any(off):
            pipe.modifiers.append(
                lambda f, d, s=off: d.particles_.positions_.__setitem__(
                    Ellipsis, d.particles.positions + s))

    pipe.modifiers.append(WrapPeriodicImagesModifier())
    pipe.modifiers.append(rebuild_bond_pbc)
    pipe.modifiers.append(make_paint(cfg, ref, nodes, cfg["rad"], pair))

    if cfg.get("zoom"):
        keep = cfg.get("keep", 4.6)
        boxc = 0.5 * (cell[0] + cell[1] + cell[2])

        def crop(frame, data, c=boxc, r=keep):
            d = np.linalg.norm(np.array(data.particles.positions) - c, axis=1)
            data.particles_.delete_elements((d > r).astype(np.int8))
        pipe.modifiers.append(crop)

    pipe.add_to_scene()

    d0 = pipe.compute(0)
    if d0.particles.bonds is not None and cfg["bw"] > 0:
        d0.particles.bonds.vis.width = cfg["bw"]
        d0.particles.bonds.vis.use_particle_colors = True
    if cfg.get("zoom"):
        d0.cell.vis.enabled = False            # the crop overruns the cell

    vp = Viewport(type=Viewport.Type.Perspective)
    vp.camera_dir = (-0.55, -0.75, -0.45)
    if cfg.get("zoom"):
        # Frame the FIXED crop sphere (radius keep at box centre), not the
        # per-frame atoms -- zoom_all would frame frame 0's sparse near-lattice
        # crop tightly and let the fuller melt crops drift in scale. Manual camera
        # keeps the pair the same size in every frame.
        keep = cfg.get("keep", 4.6)
        c = 0.5 * (cell[0] + cell[1] + cell[2])
        u = np.array((-0.55, -0.75, -0.45))
        u = u / np.linalg.norm(u)
        D = keep * 4.0
        vp.camera_pos = tuple(c - u * D)
        vp.camera_dir = tuple(u)
        vp.fov = float(2.0 * np.arctan(keep * 1.12 / D))
    else:
        vp.zoom_all(size=(RENDER, RENDER))
    rr = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_samples=16,
                         shadows=True)
    fdir = FRAMESROOT / f"{tag}_movie"
    fdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for k, f in enumerate(sel):
        out = fdir / f"{k:03d}.png"
        vp.render_image(filename=str(out), size=(RENDER, RENDER), frame=f,
                        background=(1, 1, 1), renderer=rr)
        paths.append(out)
    pipe.remove_from_scene()

    imgs = [Image.open(p).convert("RGB") for p in paths]

    def boomerang(seq):
        return seq + seq[-2:0:-1] if len(seq) > 2 else seq

    gsize = cfg.get("gif_size", GIF_SIZE)
    gcolors = cfg.get("gif_colors", GIF_COLORS)
    g = boomerang([im.resize((gsize, gsize), Image.LANCZOS) for im in imgs])
    palsrc = g[len(g) // 2].quantize(colors=gcolors, method=Image.FASTOCTREE,
                                     dither=Image.NONE)
    pal = [im.quantize(palette=palsrc, dither=Image.NONE) for im in g]
    gif = OUT / f"{tag}_arc.gif"
    pal[0].save(gif, save_all=True, append_images=pal[1:], loop=0,
                duration=GIF_MS, optimize=True, disposal=2)

    m = boomerang([im.resize((MP4_SIZE, MP4_SIZE), Image.LANCZOS) for im in imgs])
    mp4 = OUT / f"{tag}_arc.mp4"
    w = imageio.get_writer(mp4, fps=FPS, codec="libx264", macro_block_size=8,
                           ffmpeg_params=["-crf", "28", "-preset", "veryslow",
                                          "-pix_fmt", "yuv420p"])
    for im in m:
        w.append_data(np.asarray(im))
    w.close()
    print(f"[{tag}] GIF {len(g)}f {gif.stat().st_size/1024:.0f} KB | "
          f"MP4 {len(m)}f {mp4.stat().st_size/1024:.0f} KB", flush=True)


if __name__ == "__main__":
    for tag in (sys.argv[1:] or ["cg", "copoly", "atom"]):
        build(tag)
    print("arc movies done")
