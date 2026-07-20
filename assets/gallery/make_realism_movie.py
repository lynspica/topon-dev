"""'Making the junctions realistic' — render stage (OVITO python).

LEFT: the crosslinker point set + its Gaussian-matched strands. RIGHT: the live
junction g(r), with the crystalline-vs-liquid threshold marked and a running
peak-height readout.

Higher-quality render than the earlier demos: 2x supersampling (the single
biggest win -- true anti-aliasing on the sphere/cylinder edges), more AO samples,
and tuned direct light so the spheres read as solid rather than flat.
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
from PIL import Image, ImageDraw, ImageFont
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import WrapPeriodicImagesModifier

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from render_gallery import rebuild_bond_pbc

FRAMEDIR = Path(os.environ.get("TOPON_REALISM_FRAMES", "C:/tmp/realism_frames"))
OUT = HERE / "anim"; OUT.mkdir(parents=True, exist_ok=True)
TMP = FRAMEDIR / "render"; TMP.mkdir(parents=True, exist_ok=True)

PANEL = 560
SS = 2                       # supersampling factor
NODE = (0.16, 0.30, 0.52)    # junction
EDGE = (0.72, 0.76, 0.82)    # strand (pale: must not bury the ordering)
HOT = (0.80, 0.24, 0.18)     # crystalline
COOL = (0.06, 0.52, 0.46)    # liquid-like
GIF_MS = 150          # slow: each stage needs reading time

STAGE_TEXT = {
    "A": ("1 · one lattice", "SC only — every junction identical"),
    "B": ("2 · mix SC + BCC + FCC", "more neighbour shells… still crystalline"),
    "C": ("3 · add site jitter", "junctions become liquid-like"),
}


def gr_panel(r, g, peak, sq0, stage, jitter, n, out):
    fig, ax = plt.subplots(figsize=(PANEL / 100, PANEL / 100), dpi=100)
    liquid = peak < 1.5
    col = COOL if liquid else HOT
    ax.axhline(1.0, color=(0.6, 0.6, 0.65), lw=1.0, ls=(0, (3, 3)), zorder=1)
    ax.fill_between(r, g, 1.0, where=(np.array(g) > 1), color=col, alpha=0.13,
                    zorder=1)
    ax.plot(r, g, color=col, lw=2.0, zorder=3)
    ax.set_xlabel("junction–junction distance  r  (spacing units)", fontsize=11)
    ax.set_ylabel("g(r)", fontsize=11)
    ax.set_xlim(min(r), max(r)); ax.set_ylim(0, 7.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=10)
    title, sub = STAGE_TEXT[stage]
    ax.text(0.02, 0.97, title, transform=ax.transAxes, va="top", fontsize=12.5,
            fontweight="bold")
    ax.text(0.02, 0.905, sub, transform=ax.transAxes, va="top", fontsize=10.5,
            color=(0.35, 0.35, 0.4))
    ax.text(0.98, 0.97, f"g(r) peak   {peak:.2f}", transform=ax.transAxes,
            va="top", ha="right", fontsize=12, fontweight="bold", color=col)
    ax.text(0.98, 0.905,
            "liquid-like" if liquid else "crystalline",
            transform=ax.transAxes, va="top", ha="right", fontsize=10.5, color=col)
    ax.text(0.98, 0.845, f"S(q→0) {sq0:.3f}", transform=ax.transAxes, va="top",
            ha="right", fontsize=9.5, color=(0.5, 0.5, 0.55))
    if stage == "C" and jitter > 0:
        ax.text(0.02, 0.83, f"jitter  {jitter:.2f} × spacing",
                transform=ax.transAxes, va="top", fontsize=10.5,
                color=(0.35, 0.35, 0.4))
    ax.text(0.02, 0.06, f"{n} junctions", transform=ax.transAxes, fontsize=9.5,
            color=(0.55, 0.55, 0.6))
    fig.tight_layout(pad=0.6)
    fig.savefig(out, facecolor="white"); plt.close(fig)


def main():
    meta = json.loads((FRAMEDIR / "meta.json").read_text())
    frames = meta["frames"]
    dfiles = sorted(FRAMEDIR.glob("f[0-9]*.data"))
    print(f"[realism] {len(dfiles)} frames, L={meta['L']:.2f}", flush=True)

    rr = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_samples=26,
                         ambient_occlusion_brightness=0.85,
                         antialiasing=True, antialiasing_samples=10,
                         direct_light=True, direct_light_intensity=0.95,
                         shadows=True)
    vp = None
    composites = []
    try:
        fbold = ImageFont.truetype("arialbd.ttf", 22)
    except Exception:
        fbold = ImageFont.load_default()

    for k, (dpath, fm) in enumerate(zip(dfiles, frames)):
        pipe = import_file(str(dpath), atom_style="full")
        pipe.modifiers.append(WrapPeriodicImagesModifier())
        pipe.modifiers.append(rebuild_bond_pbc)

        def paint(frame, data):
            n = data.particles.count
            data.particles_.create_property(
                "Color", data=np.tile(np.array(NODE), (n, 1)))
            data.particles_.create_property("Radius", data=np.full(n, 0.15))
            # Strands drawn THIN and PALE: heavy strands turn this into a
            # hairball and bury the junction ordering, which is the subject --
            # but pale ones show the network without competing with it.
            if data.particles.bonds is not None:
                nb = data.particles.bonds.count
                data.particles_.bonds_.create_property(
                    "Color", data=np.tile(np.array(EDGE), (nb, 1)))
        pipe.modifiers.append(paint)
        pipe.add_to_scene()
        d = pipe.compute()
        if d.particles.bonds is not None:
            d.particles.bonds.vis.width = 0.020
        if vp is None:
            vp = Viewport(type=Viewport.Type.Perspective)
            # Near-axis: on a perfect lattice the sites line up into clean rows,
            # so jitter destroying that order is unmistakable. A 3/4 view muddles
            # rows together and the transition reads far weaker.
            vp.camera_dir = (-0.06, -1.0, -0.05)
            vp.zoom_all(size=(PANEL, PANEL))
        big = TMP / f"net_{k:03d}.png"
        vp.render_image(filename=str(big), size=(PANEL * SS, PANEL * SS),
                        background=(1, 1, 1), renderer=rr)     # supersample
        pipe.remove_from_scene()

        gpath = TMP / f"gr_{k:03d}.png"
        gr_panel(fm["r"], fm["g"], fm["peak"], fm["sq0"], fm["stage"],
                 fm["jitter"], fm["n"], gpath)

        net = Image.open(big).convert("RGB").resize((PANEL, PANEL), Image.LANCZOS)
        gr = Image.open(gpath).convert("RGB").resize((PANEL, PANEL), Image.LANCZOS)
        canvas = Image.new("RGB", (PANEL * 2 + 16, PANEL), "white")
        canvas.paste(net, (0, 0)); canvas.paste(gr, (PANEL + 16, 0))
        composites.append(canvas)
        print(f"  {k:03d} {fm['stage']} peak {fm['peak']:.2f}", flush=True)

    def rs(im, w):
        return im.resize((w, int(im.height * w / im.width)), Image.LANCZOS)

    seq = composites + composites[-2:0:-1][::2]        # forward + fast rewind
    durs = [GIF_MS] * len(composites) + [70] * len(composites[-2:0:-1][::2])
    g = [rs(im, 760) for im in seq]
    palsrc = g[len(g) // 2].quantize(colors=112, method=Image.FASTOCTREE,
                                     dither=Image.NONE)
    pal = [im.quantize(palette=palsrc, dither=Image.NONE) for im in g]
    gif = OUT / "realism_arc.gif"
    pal[0].save(gif, save_all=True, append_images=pal[1:], loop=0,
                duration=durs, optimize=True, disposal=2)
    m = [rs(im, 1200) for im in seq]
    mp4 = OUT / "realism_arc.mp4"
    w = imageio.get_writer(mp4, fps=7, codec="libx264", macro_block_size=8,
                           ffmpeg_params=["-crf", "26", "-preset", "veryslow",
                                          "-pix_fmt", "yuv420p"])
    for im in m:
        w.append_data(np.asarray(im))
    w.close()
    print(f"[realism] GIF {len(g)}f {gif.stat().st_size/1024:.0f} KB | "
          f"MP4 {mp4.stat().st_size/1024:.0f} KB", flush=True)


if __name__ == "__main__":
    main()
