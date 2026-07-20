"""'How BFM reaches its gel point' — animation (render stage, OVITO python).

LEFT: SAW chains on the lattice, nodes coloured by connected component (the giant
cluster bold teal, other clusters muted), crosslinks in red accumulating. RIGHT: a
gelation curve — largest-component fraction vs conversion — with a moving marker
and the gel point (first percolation to one component) flagged.

Reads f{k}.data + meta.json from bfm_frames.build_all() (topon python).
"""
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

GIANT = (0.06, 0.52, 0.46)          # the percolating cluster
XLINK = (0.78, 0.20, 0.17)          # crosslink bonds (red)
BACKBONE = (0.80, 0.82, 0.85)       # backbone bonds (faint)
# muted palette for the non-giant clusters (cycled by rank)
CLUSTERS = [(0.45, 0.52, 0.62), (0.60, 0.55, 0.42), (0.52, 0.45, 0.58),
            (0.42, 0.58, 0.60), (0.62, 0.48, 0.45), (0.50, 0.58, 0.48)]

FRAMEDIR = Path(os.environ.get("TOPON_BFM_FRAMES", "C:/tmp/bfm_frames"))
OUT = HERE / "anim"; OUT.mkdir(parents=True, exist_ok=True)
TMP = FRAMEDIR / "render"; TMP.mkdir(parents=True, exist_ok=True)
PANEL = 560
NODE_R, BW = 0.20, 0.085
HOLD = 5


def curve_panel(conv_s, frac_s, cur_conv, cur_frac, gel_conv, out,
                equil=None):
    """Right-hand panel. During phase A (`equil` = (step, total)) the chains are
    still moving and no crosslinks exist, so the gelation curve is drawn faint as
    'what comes next'; in phase B it is live with a marker."""
    fig, ax = plt.subplots(figsize=(PANEL / 100, PANEL / 100), dpi=100)
    faint = equil is not None
    ax.plot(conv_s, frac_s, color=(0.5, 0.5, 0.55),
            lw=1.4, alpha=0.25 if faint else 1.0, zorder=1)
    ax.fill_between(conv_s, frac_s, color=GIANT,
                    alpha=0.03 if faint else 0.10, zorder=0)
    if gel_conv is not None:
        ax.axvline(gel_conv, color=XLINK, lw=1.4, ls=(0, (4, 2)),
                   alpha=0.25 if faint else 1.0, zorder=2)
        ax.text(gel_conv, 1.02, "gel point", color=XLINK, ha="center",
                fontsize=10.5, fontweight="bold", alpha=0.25 if faint else 1.0)
    if not faint:
        ax.plot([cur_conv], [cur_frac], "o", ms=11, color=GIANT,
                markeredgecolor="white", markeredgewidth=1.4, zorder=4)
    ax.set_xlabel("conversion  (fraction of crosslinkers reacted)", fontsize=11)
    ax.set_ylabel("largest connected cluster  (fraction of chains)", fontsize=11)
    ax.set_xlim(0, max(conv_s) * 1.02); ax.set_ylim(0, 1.08)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=10)
    if faint:
        step, total = equil
        ax.text(0.02, 0.97, "1 · chains equilibrate", transform=ax.transAxes,
                va="top", fontsize=11.5, fontweight="bold")
        ax.text(0.02, 0.90, f"Monte-Carlo moves   {step:,} / {total:,}",
                transform=ax.transAxes, va="top", fontsize=11,
                color=(0.35, 0.35, 0.40))
        ax.text(0.02, 0.83, "end · kink/crankshaft · reptation",
                transform=ax.transAxes, va="top", fontsize=10,
                color=(0.55, 0.55, 0.60))
        ax.text(0.5, 0.45, "no crosslinks yet", transform=ax.transAxes,
                ha="center", fontsize=12, color=(0.65, 0.65, 0.70))
    else:
        state = "GEL — one spanning network" if cur_frac > 0.999 else \
                ("percolating" if cur_frac > 0.5 else "sol — disconnected clusters")
        ax.text(0.02, 0.97, f"2 · crosslink   conversion {cur_conv:.2f}",
                transform=ax.transAxes, va="top", fontsize=11.5,
                fontweight="bold")
        ax.text(0.02, 0.90, state, transform=ax.transAxes, va="top", fontsize=11,
                color=(GIANT if cur_frac > 0.5 else (0.4, 0.4, 0.45)))
    fig.tight_layout(pad=0.6)
    fig.savefig(out, facecolor="white"); plt.close(fig)


def main():
    meta = json.loads((FRAMEDIR / "meta.json").read_text())
    gel_conv = meta["gel_conv"]
    conv_s, frac_s = meta["conv_series"], meta["frac_series"]
    fconv, ffrac = meta["frame_conv"], meta["frame_frac"]
    n_equil = meta.get("n_equil_frames", 0)
    equil_total = 120000
    dfiles = sorted(FRAMEDIR.glob("f[0-9]*.data"))
    print(f"[bfm] {meta['n_chains']} chains, {meta['Nx']}^3, "
          f"{n_equil} equilibration + {len(dfiles)-n_equil} crosslink frames, "
          f"{meta['n_crosslinks']} crosslinks, gel@{gel_conv}", flush=True)

    rr = hq_renderer()
    vp = None
    composites = []
    for k, dpath in enumerate(dfiles):
        pipe = import_file(str(dpath), atom_style="full")
        pipe.modifiers.append(WrapPeriodicImagesModifier())
        pipe.modifiers.append(rebuild_bond_pbc)

        def paint(frame, data):
            mol = np.array(data.particles["Molecule Identifier"])
            n = data.particles.count
            col = np.zeros((n, 3))
            col[mol == 1] = GIANT                       # rank 0 -> mol 1 = giant
            for r in range(2, mol.max() + 1):
                m = mol == r
                if m.any():
                    col[m] = CLUSTERS[(r - 2) % len(CLUSTERS)]
            data.particles_.create_property("Color", data=col)
            data.particles_.create_property("Radius", data=np.full(n, NODE_R))
            # bond colours: backbone faint, crosslinks red
            if data.particles.bonds is not None:
                bt = np.array(data.particles.bonds["Bond Type"])
                bcol = np.tile(np.array(BACKBONE), (len(bt), 1))
                bcol[bt == 2] = XLINK
                data.particles_.bonds_.create_property("Color", data=bcol)
        pipe.modifiers.append(paint)
        pipe.add_to_scene()
        d = pipe.compute()
        if d.particles.bonds is not None:
            d.particles.bonds.vis.width = BW
        if vp is None:
            vp = Viewport(type=Viewport.Type.Perspective)
            vp.camera_dir = (-0.55, -0.75, -0.45)
            vp.zoom_all(size=(PANEL, PANEL))
        npath = TMP / f"net_{k:03d}.png"
        render_hq(vp, npath, (PANEL, PANEL), renderer=rr)
        pipe.remove_from_scene()

        hpath = TMP / f"curve_{k:03d}.png"
        if k < n_equil:                      # phase A: chains still moving
            step = int(round(equil_total * k / max(1, n_equil - 1)))
            curve_panel(conv_s, frac_s, 0.0, 0.0, gel_conv, hpath,
                        equil=(step, equil_total))
            label = f"equilibrating {step:,} MC steps"
        else:                                # phase B: crosslinking
            j = k - n_equil
            curve_panel(conv_s, frac_s, fconv[j], ffrac[j], gel_conv, hpath)
            label = f"conv {fconv[j]:.2f}, largest {ffrac[j]:.2f}"

        net = Image.open(npath).convert("RGB")
        cur = Image.open(hpath).convert("RGB").resize((PANEL, PANEL), Image.LANCZOS)
        canvas = Image.new("RGB", (PANEL * 2 + 16, PANEL), "white")
        canvas.paste(net, (0, 0)); canvas.paste(cur, (PANEL + 16, 0))
        composites.append(canvas)
        print(f"  frame {k}: {label}", flush=True)

    fwd = [composites[0]] * (HOLD - 1) + composites + [composites[-1]] * (HOLD - 1)

    def rs(im, w):
        return im.resize((w, int(im.height * w / im.width)), Image.LANCZOS)

    gseq = [rs(im, 760) for im in fwd]
    durs = [91] * len(gseq)
    rew = composites[-2:0:-1][::2]
    gseq += [rs(im, 760) for im in rew]; durs += [45] * len(rew)
    palsrc = gseq[2 * len(gseq) // 3].quantize(colors=128, method=Image.FASTOCTREE,
                                               dither=Image.NONE)
    pal = [im.quantize(palette=palsrc, dither=Image.NONE) for im in gseq]
    gif = OUT / "bfm_arc.gif"
    pal[0].save(gif, save_all=True, append_images=pal[1:], loop=0,
                duration=durs, optimize=True, disposal=2)
    mseq = [rs(im, 1120) for im in fwd] + [rs(im, 1120) for im in rew]
    mp4 = OUT / "bfm_arc.mp4"
    w = imageio.get_writer(mp4, fps=14, codec="libx264", macro_block_size=8,
                           ffmpeg_params=["-crf", "27", "-preset", "veryslow",
                                          "-pix_fmt", "yuv420p"])
    for im in mseq:
        w.append_data(np.asarray(im))
    w.close()
    print(f"[bfm] GIF {len(gseq)}f {gif.stat().st_size/1024:.0f} KB | "
          f"MP4 {mp4.stat().st_size/1024:.0f} KB", flush=True)


if __name__ == "__main__":
    main()
