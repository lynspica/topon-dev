"""Entanglement arc, TWO synchronized panels in one animation:

  LEFT  the full network, with ONE entanglement highlighted -- its two chains in
        two distinct colours (gold, violet), all grafts teal, the rest grey. A
        box is drawn around the entangled region.
  RIGHT that same region zoomed in (in a box), the crop following the pair.

Both panels play the same MD arc -- lattice -> minimised -> equilibrated -- so you
see where the entanglement sits in the whole network and, beside it, how it and
its side chains relax. Seamless boomerang loop, ~0.66x speed.

Trajectory: entangled_grafted, produced by movie_ent.in (reads the pristine
03_Conformation lattice so frame 0 is the tight 0.39-sigma crossing the stills
show, then inflates under a ramped soft push).

Run:  C:/v/ovito/Scripts/python.exe make_ent_movie.py
"""
import os
import sys
from pathlib import Path

import numpy as np
import imageio.v2 as imageio
from PIL import Image, ImageDraw
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import WrapPeriodicImagesModifier, LoadTrajectoryModifier

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from render_gallery import (node_ids, rebuild_bond_pbc, half_cell_offset, fresh,
                            closest_entangled_pair,
                            ENT_GOLD, VIOLET, REST, TEAL)

OUT = HERE / "anim"
OUT.mkdir(parents=True, exist_ok=True)
FRAMES = Path(os.environ.get("TOPON_GALLERY_FRAMES", HERE / "_frames")) / "ent_movie"
FRAMES.mkdir(parents=True, exist_ok=True)

NAME = "entangled_grafted"
PANEL = 560            # render size of each panel
KEEP = 4.6             # zoom crop radius (sigma)
RAD_FULL, BW_FULL = 0.11, 0.13
RAD_ZOOM, BW_ZOOM = 0.11, 0.13
GIF_MS = 91            # ~0.66x speed
FPS = 13
N_TRANSITION, N_LEADIN, N_TAIL = 26, 3, 8
CAM = np.array((-0.55, -0.75, -0.45))


def disorder(pipe, lt, topo, cell):
    inv = np.linalg.inv(cell)
    out = []
    for f in range(lt.source.num_frames):
        p = np.array(pipe.compute(f).particles.positions)
        v = (p[topo[:, 0]] - p[topo[:, 1]]) @ inv.T
        v -= np.round(v)
        w = v @ cell.T
        n = np.linalg.norm(w, axis=1)
        g = n > 1e-9
        out.append(float(np.mean(1.0 - np.abs(w[g]).max(axis=1) / n[g])))
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


# A paler grey for the whole network, so the ONE highlighted entanglement reads.
FAINT = (0.74, 0.76, 0.79)


def focus_grafts(probe, nodes, focus_ids):
    """Graft (type-3) bead IDs that hang off the two focus chains.

    ONLY these are coloured -- colouring every graft in the network made the one
    entanglement look like many highlighted spots. Walk each type-3 branch out
    from any focus backbone bead."""
    ident = np.array(probe.particles["Particle Identifier"])
    t = np.array(probe.particles["Particle Type"])
    topo = np.array(probe.particles.bonds.topology)
    tb = {int(ident[i]): int(t[i]) for i in range(len(ident))}
    adj = {}
    for a, b in topo:
        adj.setdefault(int(ident[a]), []).append(int(ident[b]))
        adj.setdefault(int(ident[b]), []).append(int(ident[a]))
    out = set()
    for bid in focus_ids:
        for nb in adj.get(bid, []):
            if tb.get(nb) == 3 and nb not in out:
                out.add(nb)
                stack = [nb]
                while stack:
                    u = stack.pop()
                    for v in adj.get(u, []):
                        if tb.get(v) == 3 and v not in out:
                            out.add(v)
                            stack.append(v)
    return out


def paint_fn(nodes, pa, pb, fgrafts, rad, pair_scale, dim_rest=False):
    """Colour ONLY the one entanglement: chain A gold, chain B violet, their own
    side chains teal. Everything else -- all other strands, all other grafts,
    junctions -- one faint grey. That is what keeps it a *single* entanglement."""
    a_arr = np.fromiter(pa, dtype=np.int64)
    b_arr = np.fromiter(pb, dtype=np.int64)
    g_arr = np.fromiter(fgrafts, dtype=np.int64)

    def paint(frame, data):
        ident = np.array(data.particles["Particle Identifier"])
        n = data.particles.count
        col = np.tile(np.array(FAINT), (n, 1))
        r = np.full(n, rad * (0.7 if dim_rest else 1.0))
        mg = np.isin(ident, g_arr)
        col[mg], r[mg] = TEAL, rad * 1.1          # the entanglement's side chains
        ma, mb = np.isin(ident, a_arr), np.isin(ident, b_arr)
        col[ma], col[mb] = ENT_GOLD, VIOLET
        r[ma | mb] = rad * pair_scale
        data.particles_.create_property("Color", data=col)
        data.particles_.create_property("Radius", data=r)
    return paint


def project(P, vp, size):
    """Perspective-project a world point to pixel (x, y). Camera params are read
    from the viewport after zoom_all; aspect is 1 (square panels)."""
    eye = np.array(vp.camera_pos)
    fwd = np.array(vp.camera_dir, float); fwd /= np.linalg.norm(fwd)
    up = np.array(vp.camera_up, float)
    right = np.cross(fwd, up); right /= np.linalg.norm(right)
    tup = np.cross(right, fwd)
    v = np.array(P) - eye
    depth = v @ fwd
    half = np.tan(vp.fov / 2)
    xc = (v @ right) / depth / half
    yc = (v @ tup) / depth / half
    return (xc * 0.5 + 0.5) * size, (0.5 - yc * 0.5) * size


def build_pipe(ref, arc):
    pipe = import_file(str(ref), atom_style="full")
    lt = LoadTrajectoryModifier()
    lt.source.load(str(arc))
    pipe.modifiers.append(lt)
    return pipe, lt


def main():
    ref = fresh(NAME, "03_Conformation/system_relaxed.data")
    sim = ref.parent.parent / "04_Simulation"
    nodes = node_ids(ref.parent.parent)
    arc = sim / "movie_arc.lammpstrj"
    data = (sim / "infl.lammpstrj").read_bytes()
    if (sim / "therm.lammpstrj").exists():
        data += (sim / "therm.lammpstrj").read_bytes()
    arc.write_bytes(data)

    # progress + frame selection, and the focus pair (from frame 0)
    pipe0, lt0 = build_pipe(ref, arc)
    probe = pipe0.compute(0)
    cell = np.array(probe.cell[:, :3])
    topo = np.array(probe.particles.bonds.topology)
    prog = disorder(pipe0, lt0, topo, cell)
    sel = select(prog)
    dmin, _, pr = closest_entangled_pair(probe, nodes)
    ident0 = np.array(probe.particles["Particle Identifier"])
    pa = {int(ident0[i]) for i in pr[0][1]}
    pb = {int(ident0[i]) for i in pr[1][1]}
    focus = np.fromiter(pa | pb, dtype=np.int64)
    fgrafts = focus_grafts(probe, nodes, pa | pb)
    print(f"pair {dmin:.2f} sigma, {len(focus)} beads, {len(fgrafts)} focus grafts; "
          f"selected {len(sel)} of {lt0.source.num_frames}", flush=True)
    inv = np.linalg.inv(cell)
    boxc = 0.5 * (cell[0] + cell[1] + cell[2])

    # Both panels FOLLOW the pair (recentre its centroid to the box centre each
    # frame), so it is always whole and central -- the still-lattice pair sits on
    # a periodic face, which otherwise splits it across two box faces and makes the
    # locator jump. Full = whole network around the centred pair; zoom = a crop.
    def follow(frame, data, foc=focus, inv=inv, cell=cell, c=boxc):
        ident = np.array(data.particles["Particle Identifier"])
        p = np.array(data.particles.positions)
        fp = p[np.isin(ident, foc)]
        d = (fp - fp[0]) @ inv.T
        d -= np.round(d)
        centroid = fp[0] + (d @ cell.T).mean(0)
        data.particles_.positions_[...] = p + (c - centroid)

    rr = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_samples=14,
                         shadows=True)

    # ---- FULL pipeline: whole network, the ONE entanglement highlighted ---
    full, ltf = build_pipe(ref, arc)
    full.modifiers.append(follow)
    full.modifiers.append(WrapPeriodicImagesModifier())
    full.modifiers.append(rebuild_bond_pbc)
    full.modifiers.append(paint_fn(nodes, pa, pb, fgrafts, RAD_FULL,
                                   pair_scale=2.6))
    full.add_to_scene()
    df0 = full.compute(0)
    df0.particles.bonds.vis.width = BW_FULL
    df0.particles.bonds.vis.use_particle_colors = True
    vpf = Viewport(type=Viewport.Type.Perspective)
    vpf.camera_dir = tuple(CAM)
    vpf.zoom_all(size=(PANEL, PANEL))
    # The pair is at the box centre in every frame, so the locator is fixed.
    loc_px, loc_py = project(boxc, vpf, PANEL)

    # ---- ZOOM pipeline: crop that follows the pair ------------------------
    zoom, ltz = build_pipe(ref, arc)
    zoom.modifiers.append(follow)
    zoom.modifiers.append(WrapPeriodicImagesModifier())
    zoom.modifiers.append(rebuild_bond_pbc)
    zoom.modifiers.append(paint_fn(nodes, pa, pb, fgrafts, RAD_ZOOM,
                                   pair_scale=2.0, dim_rest=True))

    def crop(frame, data, c=boxc, r=KEEP):
        d = np.linalg.norm(np.array(data.particles.positions) - c, axis=1)
        data.particles_.delete_elements((d > r).astype(np.int8))
    zoom.modifiers.append(crop)
    zoom.add_to_scene()
    dz0 = zoom.compute(0)
    dz0.particles.bonds.vis.width = BW_ZOOM
    dz0.particles.bonds.vis.use_particle_colors = True
    dz0.cell.vis.enabled = False
    vpz = Viewport(type=Viewport.Type.Perspective)
    u = CAM / np.linalg.norm(CAM)
    D = KEEP * 4.0
    vpz.camera_pos = tuple(boxc - u * D)
    vpz.camera_dir = tuple(u)
    vpz.fov = float(2.0 * np.arctan(KEEP * 1.12 / D))

    # ---- render both panels per selected frame, composite -----------------
    LOC = (200, 130, 20)          # locator box colour (matches the gold pair)
    gap = 18
    paths = []
    for k, f in enumerate(sel):
        pf = FRAMES / f"full_{k:03d}.png"
        pz = FRAMES / f"zoom_{k:03d}.png"
        vpf.render_image(filename=str(pf), size=(PANEL, PANEL), frame=f,
                         background=(1, 1, 1), renderer=rr)
        vpz.render_image(filename=str(pz), size=(PANEL, PANEL), frame=f,
                         background=(1, 1, 1), renderer=rr)
        full_im = Image.open(pf).convert("RGB")
        zoom_im = Image.open(pz).convert("RGB")

        # locator box on the full panel, at the pair's projected position
        df = full.compute(f)
        ident = np.array(df.particles["Particle Identifier"])
        p = np.array(df.particles.positions)
        fp = p[np.isin(ident, focus)]
        d = (fp - fp[0]) @ inv.T
        d -= np.round(d)
        C = fp[0] + (d @ cell.T).mean(0)
        px, py = project(C, vpf, PANEL)
        r = 78
        ImageDraw.Draw(full_im).rectangle([px - r, py - r, px + r, py + r],
                                          outline=LOC, width=4)
        # matching box on the zoom panel
        ImageDraw.Draw(zoom_im).rectangle([2, 2, PANEL - 3, PANEL - 3],
                                          outline=LOC, width=5)

        canvas = Image.new("RGB", (PANEL * 2 + gap, PANEL), "white")
        canvas.paste(full_im, (0, 0))
        canvas.paste(zoom_im, (PANEL + gap, 0))
        out = FRAMES / f"c_{k:03d}.png"
        canvas.save(out)
        paths.append(out)
    full.remove_from_scene()
    zoom.remove_from_scene()

    imgs = [Image.open(p).convert("RGB") for p in paths]

    def boomerang(seq):
        return seq + seq[-2:0:-1] if len(seq) > 2 else seq

    # GIF (subsample width for size) + MP4 (full)
    W = 900
    def rs(im, w):
        return im.resize((w, int(im.height * w / im.width)), Image.LANCZOS)
    g = boomerang([rs(im, 760) for im in imgs])
    palsrc = g[len(g) // 2].quantize(colors=96, method=Image.FASTOCTREE,
                                     dither=Image.NONE)
    pal = [im.quantize(palette=palsrc, dither=Image.NONE) for im in g]
    gif = OUT / "ent_arc.gif"
    pal[0].save(gif, save_all=True, append_images=pal[1:], loop=0,
                duration=GIF_MS, optimize=True, disposal=2)

    m = boomerang([rs(im, 1120) for im in imgs])
    mp4 = OUT / "ent_arc.mp4"
    w = imageio.get_writer(mp4, fps=FPS, codec="libx264", macro_block_size=8,
                           ffmpeg_params=["-crf", "28", "-preset", "veryslow",
                                          "-pix_fmt", "yuv420p"])
    for im in m:
        w.append_data(np.asarray(im))
    w.close()
    print(f"[ent] GIF {len(g)}f {gif.stat().st_size/1024:.0f} KB | "
          f"MP4 {len(m)}f {mp4.stat().st_size/1024:.0f} KB", flush=True)


if __name__ == "__main__":
    main()
