"""'topon sculpts any lattice' — SC / BCC / FCC pruned side by side.

Three panels, each a real strict-sculpting run replayed from its own
`move_history`: SC (6-coordinate), BCC (8), FCC (12), all thinned to the SAME
mean degree 4.0. Shows that the target degree distribution is reached regardless
of the starting lattice.

Node radius is scaled per lattice by its minimum-image nearest-neighbour distance,
so the three read at a comparable visual density (FCC's neighbours are 0.707 apart
vs SC's 1.0 -- using one radius makes FCC a solid blob).

Colours are per-lattice (not per-degree): BCC/FCC reach degree 8/12, beyond the
degree palette, and `write_data` clamps the atom type at 6.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import imageio.v2 as imageio
from PIL import Image, ImageDraw, ImageFont
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import WrapPeriodicImagesModifier

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from render_gallery import rebuild_bond_pbc, hq_renderer, render_hq

ROOT = Path(os.environ.get("TOPON_LAT3_FRAMES", "C:/tmp/lat3"))
OUT = HERE / "anim"; OUT.mkdir(parents=True, exist_ok=True)
TMP = ROOT / "render"; TMP.mkdir(parents=True, exist_ok=True)

LATS = [("SC", "simple cubic", 6, (0.20, 0.46, 0.70)),
        ("BCC", "body-centred", 8, (0.10, 0.62, 0.55)),
        ("FCC", "face-centred", 12, (0.70, 0.45, 0.20))]
PANEL = 520
HOLD = 5


def nn_distance(probe):
    p = np.array(probe.particles.positions); cell = np.array(probe.cell[:, :3])
    topo = np.array(probe.particles.bonds.topology)
    f = (p[topo[:, 0]] - p[topo[:, 1]]) @ np.linalg.inv(cell).T
    f -= np.round(f)
    return float(np.median(np.linalg.norm(f @ cell.T, axis=1)))


def main():
    metas = {lat: json.loads((ROOT / lat / "meta.json").read_text())
             for lat, _, _, _ in LATS}
    files = {lat: sorted((ROOT / lat).glob("f[0-9]*.data")) for lat, _, _, _ in LATS}
    nframes = min(len(v) for v in files.values())
    print(f"[lat3] {nframes} frames per lattice", flush=True)

    rr = hq_renderer()
    vps, radii = {}, {}
    composites = []
    try:
        fbold = ImageFont.truetype("arialbd.ttf", 24)
        freg = ImageFont.truetype("arial.ttf", 17)
    except Exception:
        fbold = freg = ImageFont.load_default()

    for k in range(nframes):
        panels = []
        for lat, name, coord, col in LATS:
            dpath = files[lat][k]
            pipe = import_file(str(dpath), atom_style="full")
            probe = pipe.compute()
            if lat not in radii:
                nn = nn_distance(probe)
                radii[lat] = (0.15 * nn, 0.055 * nn)
            R, BW = radii[lat]
            nedges = probe.particles.bonds.count if probe.particles.bonds else 0
            N = probe.particles.count
            mean_deg = 2 * nedges / N

            pipe.modifiers.append(WrapPeriodicImagesModifier())
            pipe.modifiers.append(rebuild_bond_pbc)

            def paint(frame, data, c=col, r=R):
                n = data.particles.count
                data.particles_.create_property(
                    "Color", data=np.tile(np.array(c), (n, 1)))
                data.particles_.create_property("Radius", data=np.full(n, r))
            pipe.modifiers.append(paint)
            pipe.add_to_scene()
            d = pipe.compute()
            if d.particles.bonds is not None:
                d.particles.bonds.vis.width = BW
                d.particles.bonds.vis.use_particle_colors = True
            if lat not in vps:
                vp = Viewport(type=Viewport.Type.Perspective)
                vp.camera_dir = (-0.55, -0.75, -0.45)
                vp.zoom_all(size=(PANEL, PANEL))
                vps[lat] = vp
            ppath = TMP / f"{lat}_{k:03d}.png"
            render_hq(vps[lat], ppath, (PANEL, PANEL), renderer=rr)
            pipe.remove_from_scene()
            panels.append((lat, name, coord, ppath, nedges, mean_deg,
                           metas[lat]["base_edges"]))

        pad = 52
        canvas = Image.new("RGB", (PANEL * 3, PANEL + pad), "white")
        dr = ImageDraw.Draw(canvas)
        for i, (lat, name, coord, ppath, ne, md, be) in enumerate(panels):
            canvas.paste(Image.open(ppath).convert("RGB"), (i * PANEL, pad))
            dr.text((i * PANEL + 14, 6), lat, font=fbold, fill=(30, 30, 30))
            dr.text((i * PANEL + 66, 11), f"{name} · {coord}-coordinate",
                    font=freg, fill=(110, 110, 110))
            dr.text((i * PANEL + 14, 30),
                    f"{ne} / {be} edges     mean degree {md:.2f}",
                    font=freg, fill=(60, 60, 60))
        composites.append(canvas)
        print(f"  frame {k}: " +
              "  ".join(f"{p[0]} {p[5]:.2f}" for p in panels), flush=True)

    fwd = [composites[0]] * (HOLD - 1) + composites + [composites[-1]] * (HOLD - 1)

    def rs(im, w):
        return im.resize((w, int(im.height * w / im.width)), Image.LANCZOS)

    gseq = [rs(im, 900) for im in fwd]
    durs = [91] * len(gseq)
    rew = composites[-2:0:-1][::2]
    gseq += [rs(im, 900) for im in rew]; durs += [45] * len(rew)
    palsrc = gseq[len(gseq) // 3].quantize(colors=96, method=Image.FASTOCTREE,
                                           dither=Image.NONE)
    pal = [im.quantize(palette=palsrc, dither=Image.NONE) for im in gseq]
    gif = OUT / "lattices_arc.gif"
    pal[0].save(gif, save_all=True, append_images=pal[1:], loop=0,
                duration=durs, optimize=True, disposal=2)
    mseq = [rs(im, 1320) for im in fwd] + [rs(im, 1320) for im in rew]
    mp4 = OUT / "lattices_arc.mp4"
    w = imageio.get_writer(mp4, fps=14, codec="libx264", macro_block_size=8,
                           ffmpeg_params=["-crf", "27", "-preset", "veryslow",
                                          "-pix_fmt", "yuv420p"])
    for im in mseq:
        w.append_data(np.asarray(im))
    w.close()
    print(f"[lat3] GIF {len(gseq)}f {gif.stat().st_size/1024:.0f} KB | "
          f"MP4 {mp4.stat().st_size/1024:.0f} KB", flush=True)


if __name__ == "__main__":
    main()
