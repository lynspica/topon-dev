"""Render entangle_density runs: chains cold-to-hot by delivered count.

    "C:/v/ovito/Scripts/python.exe" tests/workflows/render_entangled.py
    "C:/v/ovito/Scripts/python.exe" tests/workflows/render_entangled.py \
        relaxed_ctrl3c relaxed_grad5

Needs OVITO's python, no topon imports: everything comes from the
``designed_pairs.json`` each run writes next to its ``designed.data`` --
which chains form each designed pair, which shell, and how many
entanglements Z1+ measured on it. Chains are drawn solid, blue for e=0
through orange and red to dark red for e>=3, entangled chains slightly
fatter, nothing else styled. Arguments are run directory names under
``tests/output/entangle_steps``; with none, the five verification systems.
"""
import json
import sys
from pathlib import Path

import numpy as np

from ovito.io import import_file
from ovito.modifiers import DeleteSelectedModifier
from ovito.vis import TachyonRenderer, Viewport

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "tests/output/entangle_steps"

DEFAULT = [
    ("relaxed_grad5",   "gradient z^3, shells 1:0.5/2:0.5"),
    ("relaxed_hot5",    "hotspot x8 (centre sphere), shells 1:0.5/2:0.5"),
    ("relaxed_ctrl3c",  "four-shell 20/50/25/5, one-run controller"),
    ("relaxed_ctrl4c",  "uniform over four shells, one-run controller"),
    ("relaxed_ctrl5c",  "weighted outward 5/15/30/50, one-run controller"),
]

HEAT = {0: (0.25, 0.42, 0.78),
        1: (0.89, 0.50, 0.30),
        2: (0.80, 0.14, 0.11),
        3: (0.45, 0.03, 0.05)}


def render_one(root, out_png, size=900):
    dump = json.loads((root / "04_Simulation" /
                       "designed_pairs.json").read_text())
    pipe = import_file(str(root / "04_Simulation" / "designed.data"),
                       atom_style="full")
    d0 = pipe.compute()
    ids = d0.particles["Particle Identifier"].array
    row = {int(a): i for i, a in enumerate(ids)}
    n = d0.particles.count

    dens = np.zeros(n, int)
    for p in dump["pairs"]:
        if p["delivered"]:
            atoms = [row[a] for a in p["atoms_a"] + p["atoms_b"] if a in row]
            dens[atoms] = np.maximum(dens[atoms],
                                     min(int(p["delivered"]), 3))

    def style(frame, dd):
        c = np.empty((n, 3))
        r = np.full(n, 0.16)
        for lvl, col in HEAT.items():
            m = dens == lvl
            c[m] = col
            if lvl > 0:
                r[m] = 0.26
        dd.particles_.create_property("Color", data=c)
        dd.particles_.create_property("Radius", data=r)

    def drop_bonds(frame, dd):
        if dd.particles.bonds is not None and dd.particles.bonds.count:
            dd.particles_.bonds_.create_property(
                "Selection",
                data=np.ones(dd.particles.bonds.count, dtype=np.int8))

    pipe.modifiers.append(style)
    pipe.modifiers.append(drop_bonds)
    pipe.modifiers.append(DeleteSelectedModifier(operate_on={"bonds"}))
    pipe.add_to_scene()

    d2 = pipe.compute()
    d2.cell.vis.enabled = True
    cell = np.asarray(d2.cell[...])
    corners = np.array([cell[:, 3] + i * cell[:, 0] + j * cell[:, 1]
                        + k * cell[:, 2]
                        for i in (0, 1) for j in (0, 1) for k in (0, 1)])
    lo, hi = corners.min(0), corners.max(0)
    centre = 0.5 * (lo + hi)
    radius = 0.5 * float(np.linalg.norm(hi - lo)) + 1

    dv = np.array([-0.55, -0.75, -0.45])
    dv /= np.linalg.norm(dv)
    vp = Viewport(type=Viewport.Type.Perspective, camera_dir=tuple(dv))
    vp.fov = np.deg2rad(35.0)
    vp.camera_pos = tuple(centre - dv * radius / np.sin(vp.fov / 2.0))
    vp.render_image(filename=str(out_png), size=(size, size),
                    background=(1.0, 1.0, 1.0),
                    renderer=TachyonRenderer(ambient_occlusion=False,
                                             shadows=False))
    pipe.remove_from_scene()

    got = dump["delivered"]
    order = sorted(int(k) for k in got)
    return ("delivered " + "/".join(f"{got[str(s)]:.2f}".lstrip("0")
                                    for s in order)
            + "  (shells " + "/".join(str(s) for s in order) + ")")


def main():
    names = sys.argv[1:]
    panels = ([(t, lab) for t, lab in DEFAULT] if not names
              else [(t, t) for t in names])
    items = []
    for tag, label in panels:
        run = OUT / tag
        if not (run / "04_Simulation" / "designed_pairs.json").exists():
            print(f"  {tag}: no designed_pairs.json, skipped")
            continue
        png = OUT / f"heat_{tag}.png"
        line = render_one(run, png)
        print(f"  {tag}: {line}", flush=True)
        items.append((png, label, line))
    if not items:
        return 1

    from PIL import Image, ImageDraw, ImageFont
    ims = [(Image.open(p), lab, ln) for p, lab, ln in items]
    w, h = ims[0][0].size
    cols = min(3, len(ims))
    rows = (len(ims) + 1 + cols - 1) // cols       # +1 for the legend cell
    pad, head, cap = 14, 64, 66
    sheet = Image.new("RGB", (cols * w + (cols + 1) * pad,
                              head + rows * (h + cap) + pad), "white")
    dr = ImageDraw.Draw(sheet)

    def font(sz):
        for nm in ("segoeui.ttf", "arial.ttf", "DejaVuSans.ttf"):
            try:
                return ImageFont.truetype(nm, sz)
            except OSError:
                continue
        return ImageFont.load_default()

    dr.text((pad, 16), "Chains coloured cold to hot by their delivered "
            "entanglement count", fill="black", font=font(30))
    for i, (im, lab, ln) in enumerate(ims):
        r, c = divmod(i, cols)
        x, y = pad + c * (w + pad), head + r * (h + cap)
        sheet.paste(im, (x, y))
        dr.text((x + 6, y + h + 4), lab, fill="black", font=font(24))
        dr.text((x + 6, y + h + 34), ln, fill=(90, 90, 90), font=font(21))

    r, c = divmod(len(ims), cols)
    x, y = pad + c * (w + pad), head + r * (h + cap)
    dr.text((x + 40, y + 80), "chain colour", fill="black", font=font(30))
    for k, (col, txt) in enumerate([(HEAT[0], "e = 0"), (HEAT[1], "e = 1"),
                                    (HEAT[2], "e = 2"),
                                    (HEAT[3], "e = 3 or more")]):
        yy = y + 160 + 60 * k
        dr.ellipse([x + 40, yy, x + 78, yy + 38],
                   fill=tuple(int(255 * v) for v in col))
        dr.text((x + 98, yy + 4), txt, fill="black", font=font(24))
    out = OUT / "chains_by_e.png"
    sheet.save(out)
    print(f"\n  {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
