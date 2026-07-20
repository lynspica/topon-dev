"""'How topon prunes a lattice to a target degree distribution' — animation.

LEFT: the network, nodes coloured by current degree (deep blue = 6, warming as
they lose edges), rendered in the gallery OVITO style. RIGHT: a live degree
histogram, bars coloured to match, with the TARGET distribution as an outline.

Driven by the REAL generator: it replays the exact `move_history` the
strict-sculpting algorithm recorded, so the edges vanish in the true order the
algorithm removed them. Forward play (pruning) then a fast rewind, looped.

Run:  C:/v/ovito/Scripts/python.exe make_sculpt_movie.py
"""
"""Render stage (OVITO python) -- reads the .data frames + meta.json that
sculpt_frames.build_all() wrote in the topon python. No networkx here: each node's
degree is its atom type - 1."""
import collections
import json
import os
import sys
from pathlib import Path

import numpy as np
import imageio.v2 as imageio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import WrapPeriodicImagesModifier

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from render_gallery import rebuild_bond_pbc, hq_renderer, render_hq

# degree -> colour (must match sculpt_frames.DEG_COLOUR)
DEG_COLOUR = {
    6: (0.16, 0.28, 0.52), 5: (0.20, 0.46, 0.70), 4: (0.10, 0.62, 0.55),
    3: (0.65, 0.68, 0.20), 2: (0.86, 0.52, 0.18), 1: (0.78, 0.25, 0.20),
    0: (0.80, 0.80, 0.83),
}

OUT = HERE / "anim"
OUT.mkdir(parents=True, exist_ok=True)
# where sculpt_frames.build_all() wrote the .data frames + meta.json
FRAMEDIR = Path(os.environ.get("TOPON_SCULPT_FRAMES", "C:/tmp/sculpt_frames"))
TMP = FRAMEDIR / "render"
TMP.mkdir(parents=True, exist_ok=True)

PANEL = 560
NODE_R, BOND_W = 0.17, 0.07
HOLD = 5


def hist_panel(deg_counts, base_counts, target_counts, removed, total, mean_deg,
               out):
    """Degree histogram, bars coloured by degree, target as an outline."""
    fig, ax = plt.subplots(figsize=(PANEL / 100, PANEL / 100), dpi=100)
    degs = list(range(0, 7))
    cur = [deg_counts.get(d, 0) for d in degs]
    tgt = [target_counts.get(d, 0) for d in degs]
    cols = [DEG_COLOUR[d] for d in degs]
    ax.bar(degs, cur, width=0.8, color=cols, edgecolor="white", linewidth=0.6,
           zorder=2)
    # target outline
    ax.bar(degs, tgt, width=0.8, facecolor="none", edgecolor=(0.2, 0.2, 0.2),
           linewidth=1.4, linestyle=(0, (4, 2)), zorder=3)
    ax.set_xlabel("node degree (functionality)", fontsize=11)
    ax.set_ylabel("number of nodes", fontsize=11)
    ax.set_xticks(degs)
    ax.set_ylim(0, max(base_counts.values()) * 1.08)
    ax.set_xlim(-0.6, 6.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=10)
    ax.text(0.02, 0.97, f"edges removed  {removed} / {total}",
            transform=ax.transAxes, va="top", fontsize=11.5, fontweight="bold")
    ax.text(0.02, 0.90, f"mean degree  {mean_deg:.2f}",
            transform=ax.transAxes, va="top", fontsize=11)
    ax.text(0.98, 0.90, "– – target", transform=ax.transAxes, va="top",
            ha="right", fontsize=10, color=(0.2, 0.2, 0.2))
    fig.tight_layout(pad=0.6)
    fig.savefig(out, facecolor="white")
    plt.close(fig)


def main():
    meta = json.loads((FRAMEDIR / "meta.json").read_text())
    N = meta["n_nodes"]
    total = meta["total"]
    at = meta["at"]
    base_counts = {int(k): v for k, v in meta["base_counts"].items()}
    target_counts = {int(k): v for k, v in meta["target_counts"].items()}
    dfiles = sorted(FRAMEDIR.glob("f[0-9]*.data"))
    print(f"[sculpt] {N} nodes, {meta['base_edges']} -> {meta['final_edges']} "
          f"edges over {total} removals; {len(dfiles)} frames", flush=True)

    rr = hq_renderer()
    vp = None
    composites = []
    for k, dpath in enumerate(dfiles):
        pipe = import_file(str(dpath), atom_style="full")
        probe = pipe.compute()
        types = np.array(probe.particles["Particle Type"])
        dcount = collections.Counter(int(t) - 1 for t in types)   # degree = type-1
        nedges = probe.particles.bonds.count if probe.particles.bonds else 0
        mean_deg = 2 * nedges / N
        rmv = at[k] if k < len(at) else total

        pipe.modifiers.append(WrapPeriodicImagesModifier())
        pipe.modifiers.append(rebuild_bond_pbc)

        def paint(frame, data):
            t = np.array(data.particles["Particle Type"])
            n = data.particles.count
            col = np.zeros((n, 3))
            for tid in range(1, 8):
                m = t == tid
                if m.any():
                    col[m] = DEG_COLOUR[tid - 1]
            data.particles_.create_property("Color", data=col)
            data.particles_.create_property("Radius", data=np.full(n, NODE_R))
        pipe.modifiers.append(paint)
        pipe.add_to_scene()
        d = pipe.compute()
        if d.particles.bonds is not None:
            d.particles.bonds.vis.width = BOND_W
            d.particles.bonds.vis.use_particle_colors = True
        if vp is None:
            vp = Viewport(type=Viewport.Type.Perspective)
            vp.camera_dir = (-0.55, -0.75, -0.45)
            vp.zoom_all(size=(PANEL, PANEL))       # fix camera on the full lattice
        npath = TMP / f"net_{k:03d}.png"
        render_hq(vp, npath, (PANEL, PANEL), renderer=rr)
        pipe.remove_from_scene()

        hpath = TMP / f"hist_{k:03d}.png"
        hist_panel(dcount, base_counts, target_counts, rmv, total, mean_deg, hpath)

        net = Image.open(npath).convert("RGB")
        his = Image.open(hpath).convert("RGB").resize((PANEL, PANEL), Image.LANCZOS)
        canvas = Image.new("RGB", (PANEL * 2 + 16, PANEL), "white")
        canvas.paste(net, (0, 0)); canvas.paste(his, (PANEL + 16, 0))
        composites.append(canvas)
        print(f"  frame {k}: {nedges} edges, mean deg {mean_deg:.2f}", flush=True)

    # forward (prune) + fast rewind, with holds
    fwd = [composites[0]] * (HOLD - 1) + composites + [composites[-1]] * (HOLD - 1)

    def rs(im, w):
        return im.resize((w, int(im.height * w / im.width)), Image.LANCZOS)

    gseq = [rs(im, 760) for im in fwd]
    durs = [91] * len(gseq)
    # fast rewind (skip endpoints, 2x speed)
    rew = composites[-2:0:-1][::2]
    gseq += [rs(im, 760) for im in rew]
    durs += [45] * len(rew)

    palsrc = gseq[len(gseq) // 3].quantize(colors=128, method=Image.FASTOCTREE,
                                           dither=Image.NONE)
    pal = [im.quantize(palette=palsrc, dither=Image.NONE) for im in gseq]
    gif = OUT / "sculpt_arc.gif"
    pal[0].save(gif, save_all=True, append_images=pal[1:], loop=0,
                duration=durs, optimize=True, disposal=2)

    mseq = [rs(im, 1120) for im in fwd] + [rs(im, 1120) for im in rew]
    mp4 = OUT / "sculpt_arc.mp4"
    w = imageio.get_writer(mp4, fps=14, codec="libx264", macro_block_size=8,
                           ffmpeg_params=["-crf", "27", "-preset", "veryslow",
                                          "-pix_fmt", "yuv420p"])
    for im in mseq:
        w.append_data(np.asarray(im))
    w.close()
    print(f"[sculpt] GIF {len(gseq)}f {gif.stat().st_size/1024:.0f} KB | "
          f"MP4 {mp4.stat().st_size/1024:.0f} KB", flush=True)


if __name__ == "__main__":
    main()
