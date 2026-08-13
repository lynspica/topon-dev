"""Render the shell gallery with OVITO, in the house style.

    "C:/v/ovito/Scripts/python.exe" tests/workflows/render_shells.py
    "C:/v/ovito/Scripts/python.exe" tests/workflows/render_shells.py --lattice FCC

One panel per shell, combined into a strip, all at one scale. Each panel
carries a routed chain wound around a partner from that shell.

What the strip does *not* show is a reach that grows with the shell number.
Measured on SC: 7.8, 19.1, 18.9, 16.6 sigma for shells 1 to 4. Shell 1 is
genuinely compact, but beyond it the extent is dominated by the length of the
chains themselves, about 13 sigma, rather than by the gap between them, which
is 0, 13, 18 and 22 sigma. The shells differ in *which* partner is reached,
and only the first is visibly closer.

The background network is drawn thin and quiet: it is 8217 of the 8438 beads,
and at full size it buries the two chains the figure is about.

Colours follow the gallery skill -- quiet grey for the background, teal for
the routed chain -- with a warm red for its partner, since the message here is
a pair and needs two foreground colours rather than one.

Note the type convention differs from the copolymer renders: junctions are
type 4 here, not type 3.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier
from ovito.vis import TachyonRenderer, Viewport

GREY = np.array([0.788, 0.808, 0.831])     # #c9ced4, the background network
TEAL = np.array([0.055, 0.486, 0.420])     # #0e7c6b, the routed chain
RED = np.array([0.690, 0.227, 0.180])      # #B03A2E, its partner
DARK = np.array([0.173, 0.243, 0.314])     # #2C3E50, crosslink junctions

# Type ids as written by entangle_shellshow.py.
T_REST, T_ROUTED, T_PARTNER, T_JUNCTION = 1, 2, 3, 4


def subject_of(data_file):
    """Centre and half-extent of the designed pair."""
    pipe = import_file(str(data_file), atom_style="full")
    data = pipe.compute()
    t = data.particles["Particle Type"].array
    pos = data.particles.positions.array
    sub = pos[(t == T_ROUTED) | (t == T_PARTNER)]
    c = sub.mean(axis=0)
    return c, float(np.abs(sub - c).max())


def render_one(data_file, out_png, size=900, span=None):
    pipe = import_file(str(data_file), atom_style="full")
    pipe.modifiers.append(CreateBondsModifier(cutoff=1.6))

    def paint(frame, data):
        t = data.particles["Particle Type"].array
        colors = np.tile(GREY, (data.particles.count, 1))
        colors[t == T_ROUTED] = TEAL
        colors[t == T_PARTNER] = RED
        colors[t == T_JUNCTION] = DARK
        data.particles_.create_property("Color", data=colors)

    def set_radii(frame, data):
        t = data.particles["Particle Type"].array
        # The background is 97 percent of the beads. At the usual 0.28 it
        # buries the two chains the figure exists to show, so it is drawn
        # thin and the designed pair keeps the house size.
        radii = np.full(data.particles.count, 0.10)
        radii[t == T_ROUTED] = 0.28
        radii[t == T_PARTNER] = 0.28
        radii[t == T_JUNCTION] = 0.30
        data.particles_.create_property("Radius", data=radii)

    pipe.modifiers.append(paint)
    pipe.modifiers.append(set_radii)
    pipe.add_to_scene()

    data = pipe.compute()
    data.cell.vis.enabled = False
    if data.particles.bonds is not None:
        data.particles.bonds.vis.width = 0.16

    # House camera direction, but framed on the designed pair.
    #
    # zoom_all() fits the whole cell, and at that scale the two chains the
    # figure is about occupy a tenth of the frame and read as a smudge. The
    # 3/4 view is kept; only the framing changes.
    centre, reach = subject_of(data_file)
    d = np.array([-0.55, -0.75, -0.45], float)
    d = d / np.linalg.norm(d)
    vp = Viewport(type=Viewport.Type.Perspective)
    vp.camera_dir = tuple(d)
    vp.fov = np.deg2rad(35.0)
    # One scale for every panel in a strip.
    #
    # Framing each panel on its own pair makes each legible and the strip
    # meaningless: the message is that an outer shell reaches further, and
    # that is exactly what a per-panel zoom cancels out. The caller passes the
    # largest reach in the set so every panel is at the same scale.
    use = span if span is not None else max(reach * 2.6, 6.0)
    vp.camera_pos = tuple(centre - d * use / np.tan(vp.fov / 2.0))
    vp.render_image(
        filename=str(out_png), size=(size, size),
        background=(1.0, 1.0, 1.0),
        renderer=TachyonRenderer(ambient_occlusion=True,
                                 ambient_occlusion_samples=12,
                                 shadows=True))
    pipe.remove_from_scene()
    return out_png


def strip(panels, labels, out_png, title):
    """Combine panels left to right with a title and per-panel captions."""
    from PIL import Image, ImageDraw, ImageFont

    ims = [Image.open(p) for p in panels]
    w, h = ims[0].size
    pad, head, foot = 12, 64, 52
    W = len(ims) * w + (len(ims) + 1) * pad
    H = head + h + foot + pad
    sheet = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(sheet)

    def font(sz):
        for name in ("segoeui.ttf", "arial.ttf", "DejaVuSans.ttf"):
            try:
                return ImageFont.truetype(name, sz)
            except OSError:
                continue
        return ImageFont.load_default()

    d.text((pad, 18), title, fill="black", font=font(30))
    for i, (im, lab) in enumerate(zip(ims, labels)):
        x = pad + i * (w + pad)
        sheet.paste(im, (x, head))
        d.text((x + 6, head + h + 8), lab, fill="black", font=font(22))
    sheet.save(out_png)
    return out_png


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lattice", default="SC")
    ap.add_argument("--orient", default="any")
    ap.add_argument("--shells", type=int, nargs="*", default=[1, 2, 3, 4])
    ap.add_argument("--gallery", default=None)
    ap.add_argument("--size", type=int, default=900)
    args = ap.parse_args()

    root = Path(args.gallery) if args.gallery else (
        Path(__file__).resolve().parents[2]
        / "tests/output/entangle_steps/shell_gallery")
    out = root / "renders"
    out.mkdir(parents=True, exist_ok=True)

    gaps = {"SC": {1: "0.00", 2: "1.00", 3: "1.41", 4: "1.73"},
            "BCC": {1: "0.00", 2: "0.71", 3: "0.82", 4: "0.87"},
            "FCC": {1: "0.00", 2: "0.53", 3: "0.60", 4: "0.64"},
            "MIX": {1: "0.00", 2: "0.45", 3: "0.57", 4: "0.67"}}

    # Measure every pair first, so the strip can share one scale.
    found = []
    for sh in args.shells:
        src = root / f"{args.lattice}_{args.orient}_shell{sh}.data"
        if not src.exists():
            print(f"  shell {sh}: no file, skipping", flush=True)
            continue
        _c, reach = subject_of(src)
        found.append((sh, src, reach))
    if not found:
        print("  nothing to render")
        return 1
    span = max(max(r for _s, _f, r in found) * 2.6, 6.0)
    print(f"  common scale: reach "
          + ", ".join(f"shell {s} {r:.1f}" for s, _f, r in found)
          + f" sigma -> span {span:.1f}", flush=True)

    panels, labels = [], []
    for sh, src, _reach in found:
        png = out / f"{args.lattice}_{args.orient}_shell{sh}.png"
        print(f"  rendering shell {sh} ...", flush=True)
        render_one(src, png, args.size, span)
        panels.append(png)
        g = gaps.get(args.lattice, {}).get(sh, "?")
        labels.append(f"shell {sh}   gap {g} lattice units"
                      + ("   (share a crosslink)" if sh == 1 else ""))

    if not panels:
        print("  nothing rendered")
        return 1
    final = out / f"{args.lattice}_{args.orient}_shells.png"
    strip(panels, labels,
          final,
          f"{args.lattice}, {args.orient} partners: the same chain (teal) "
          f"wound twice around a partner (red) drawn from each shell")
    print(f"  {final}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
