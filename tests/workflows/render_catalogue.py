"""A catalogue of render styles for the same system, to choose from.

    "C:/v/ovito/Scripts/python.exe" tests/workflows/render_catalogue.py

Six treatments of one shell-2 configuration. They differ in what they hide,
which is the only real decision in a figure like this: the designed pair is
221 beads out of 8438, so anything that leaves the background at full weight
buries the subject, and anything that removes it entirely loses the context
that makes the shell meaningful.

    A  network in context   background thin and pale, the house 3/4 view
    B  pair alone           background hidden, nothing but the two chains
    C  local neighbourhood  only beads within a radius of the pair
    D  fat pair             background thin, designed chains oversized
    E  ghost cage           background as thin sticks, no beads at all
    F  close crop           tight on the winding itself

Every panel uses the same camera direction so they can be compared, and each
is labelled with what it drops.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from ovito.io import import_file
from ovito.modifiers import DeleteSelectedModifier
from ovito.vis import TachyonRenderer, Viewport

GREY = np.array([0.788, 0.808, 0.831])
PALE = np.array([0.882, 0.898, 0.914])
TEAL = np.array([0.055, 0.486, 0.420])
RED = np.array([0.690, 0.227, 0.180])
DARK = np.array([0.173, 0.243, 0.314])

T_REST, T_ROUTED, T_PARTNER, T_JUNCTION = 1, 2, 3, 4

# Transparency, and which junctions anchor the designed chains.
#
# Two things the first catalogue got wrong. A pale grey background still reads
# as solid, so it competes with the subject instead of receding; real
# transparency lets the design sit in front of the network rather than among
# it. And the local views hid the crosslinks the routed chain actually runs
# between, which is the thing that makes a winding legible as a design rather
# than a squiggle -- those junctions are now found and drawn large in the
# chain's own colour.

STYLES = {
    "D1_plain": dict(
        label="D1  fat pair, opaque network",
        note="the original D, for reference",
        rest_r=0.08, pair_r=0.44, junc_r=0.26, bond=0.20,
        alpha=0.0, anchors=False, keep="all", zoom=2.4),
    "D2_transparent": dict(
        label="D2  fat pair, network at 70 percent",
        note="background recedes, design sits in front",
        rest_r=0.08, pair_r=0.44, junc_r=0.26, bond=0.20,
        alpha=0.7, anchors=False, keep="all", zoom=2.4),
    "D3_anchors": dict(
        label="D3  transparent, anchors marked",
        note="the four crosslinks the two chains run between",
        rest_r=0.08, pair_r=0.44, junc_r=0.22, bond=0.20,
        alpha=0.75, anchors=True, keep="all", zoom=2.4),
    "D4_anchors_near": dict(
        label="D4  anchors, local view",
        note="within 14 sigma, anchors kept whatever the distance",
        rest_r=0.12, pair_r=0.44, junc_r=0.24, bond=0.20,
        alpha=0.55, anchors=True, keep="near", zoom=1.6),
    "D5_ghost": dict(
        label="D5  near-invisible network",
        note="background at 90 percent, anchors marked",
        rest_r=0.06, pair_r=0.46, junc_r=0.20, bond=0.18,
        alpha=0.9, anchors=True, keep="all", zoom=2.4),
    "D6_crop": dict(
        label="D6  close crop with anchors",
        note="tight, but the anchoring crosslinks stay in frame",
        rest_r=0.09, pair_r=0.46, junc_r=0.22, bond=0.09,
        alpha=0.55, anchors=True, keep="near", near=14.0, zoom=1.15),
    "D7_dense": dict(
        label="D7  close crop, dense lattice",
        note="tighter cutoff and thinner strands for a packed network",
        rest_r=0.05, pair_r=0.46, junc_r=0.20, bond=0.07,
        alpha=0.8, anchors=True, keep="near", near=7.0, zoom=1.15),
}


def reach_of(data_file):
    """Half-extent of the designed pair, for a shared scale."""
    d = import_file(str(data_file), atom_style="full").compute()
    t = d.particles["Particle Type"].array
    p = d.particles.positions.array
    sub = p[(t == T_ROUTED) | (t == T_PARTNER)]
    return float(np.abs(sub - sub.mean(axis=0)).max())


def render(data_file, out_png, st, size=900, span=None,
           frame_box=False):
    # The data file's own bonds, not distance-generated ones.
    #
    # CreateBondsModifier turned 71680 real bonds into 1492823 by joining any
    # two beads within the cutoff, bonded or not. At melt-like spacing that is
    # most of the neighbours, which is why the network rendered as solid white
    # pipework and buried every design inside it.
    pipe = import_file(str(data_file), atom_style="full")

    probe = pipe.compute()
    t0 = probe.particles["Particle Type"].array
    p0 = probe.particles.positions.array
    sub = p0[(t0 == T_ROUTED) | (t0 == T_PARTNER)]
    centre = sub.mean(axis=0)
    reach = float(np.abs(sub - centre).max())

    # How much neighbourhood to keep. A dense lattice needs less: on MIX
    # a 14 sigma crop still holds enough strands to bury the design
    # entirely, where on SC it is the right amount of context.
    near = np.linalg.norm(p0 - centre, axis=1) < st.get("near", 14.0)

    # The crosslinks the two designed chains are anchored to.
    #
    # A junction is type 4 whichever chain it belongs to, so the ones that
    # matter here are found by proximity: a type-4 bead within a bond length
    # of a designed bead is an endpoint of that chain. Without these the local
    # views show a winding with no visible attachment, which is the one thing
    # that makes it a design rather than a squiggle.
    from scipy.spatial import cKDTree
    junc = np.flatnonzero(t0 == T_JUNCTION)
    anchor_r, anchor_p = np.zeros(len(p0), bool), np.zeros(len(p0), bool)
    if len(junc):
        tree = cKDTree(p0[junc])
        for mask, out in ((t0 == T_ROUTED, anchor_r),
                          (t0 == T_PARTNER, anchor_p)):
            # One bond length, not two: at 1.6 a neighbouring chain's
            # junction is also caught and the routed chain appeared to
            # have three endpoints.
            hit = tree.query_ball_point(p0[mask], r=1.15)
            for group in hit:
                for j in group:
                    out[junc[j]] = True

    def paint(frame, data):
        t = data.particles["Particle Type"].array
        c = np.tile(PALE, (data.particles.count, 1))
        c[t == T_ROUTED] = TEAL
        c[t == T_PARTNER] = RED
        c[t == T_JUNCTION] = DARK
        if st["anchors"]:
            # An anchor takes its chain's colour, darkened, so it reads as
            # that chain's endpoint rather than as network furniture.
            c[anchor_r] = np.array([0.02, 0.24, 0.20])   # deep teal
            c[anchor_p] = np.array([0.36, 0.08, 0.06])   # deep red
        data.particles_.create_property("Color", data=c)

    def fade(frame, data):
        if not st["alpha"]:
            return
        t = data.particles["Particle Type"].array
        a = np.full(data.particles.count, float(st["alpha"]))
        a[(t == T_ROUTED) | (t == T_PARTNER)] = 0.0
        if st["anchors"]:
            a[anchor_r] = 0.0
            a[anchor_p] = 0.0
        data.particles_.create_property("Transparency", data=a)

    def radii(frame, data):
        t = data.particles["Particle Type"].array
        r = np.full(data.particles.count, st["rest_r"])
        r[t == T_ROUTED] = st["pair_r"]
        r[t == T_PARTNER] = st["pair_r"]
        r[t == T_JUNCTION] = st["junc_r"]
        if st["anchors"]:
            # Larger than the chain beads by a clear margin. At 0.62 against
            # a 0.44 chain they disappeared into the winding.
            r[anchor_r] = 1.15
            r[anchor_p] = 1.15
        data.particles_.create_property("Radius", data=r)

    def cull(frame, data):
        # Delete, do not shrink.
        #
        # A radius of zero hides the bead and leaves its bonds, which then
        # render as full-width tubes: four panels of the first catalogue were
        # dominated by white pipework that was supposed to be hidden. Removing
        # the particles removes their bonds with them.
        t = data.particles["Particle Type"].array
        keep = (t == T_ROUTED) | (t == T_PARTNER)
        if st["anchors"]:
            keep = keep | anchor_r | anchor_p
        if st["keep"] == "pair":
            sel = ~keep
        elif st["keep"] == "near":
            sel = ~keep & ~near
        else:
            return
        data.particles_.create_property("Selection",
                                        data=sel.astype(np.int8))

    pipe.modifiers.append(paint)
    pipe.modifiers.append(fade)
    pipe.modifiers.append(radii)
    if st["keep"] != "all":
        pipe.modifiers.append(cull)
        pipe.modifiers.append(DeleteSelectedModifier())
    pipe.add_to_scene()

    def cut_wrapped(frame, data):
        # Drop bonds that span the box.
        #
        # Chains are written unwrapped, so a bond from a chain end to a
        # junction that sits on the far side of the cell is drawn as a line
        # straight across the picture. Those are the long dark streaks; they
        # are an artefact of the coordinates, not a feature of the network.
        b = data.particles.bonds
        if b is None or b.count == 0:
            return
        topo = b.topology.array
        pos = data.particles.positions.array
        v = pos[topo[:, 1]] - pos[topo[:, 0]]
        long = np.linalg.norm(v, axis=1) > 2.5
        if long.any():
            data.particles_.bonds_.create_property(
                "Selection", data=long.astype(np.int8))

    pipe.modifiers.append(cut_wrapped)
    pipe.modifiers.append(DeleteSelectedModifier(operate_on={"bonds"}))

    data = pipe.compute()
    data.cell.vis.enabled = False
    if data.particles.bonds is not None:
        data.particles.bonds.vis.width = st["bond"]

    d = np.array([-0.55, -0.75, -0.45], float)
    d = d / np.linalg.norm(d)
    vp = Viewport(type=Viewport.Type.Perspective)
    vp.camera_dir = tuple(d)
    vp.fov = np.deg2rad(35.0)
    # A caller comparing panels passes one span for all of them: per-panel
    # framing makes shell 1 look huge and the rest small, which is a
    # difference in zoom being read as a difference in the design.
    if frame_box:
        # Frame the cell, not the pair.
        #
        # A network overview is about the lattice, and framing it on two
        # chains gave four panels at four different scales: SC and MIX came
        # out as specks while FCC overflowed its frame. The box is the subject
        # there, so its own extent sets the distance.
        # From the particles, not the cell. Chains are written unwrapped and
        # run well outside the periodic box, so a cell-framed camera left
        # them overflowing: FCC filled its panel while SC and MIX were
        # specks. A robust extent is the middle of the distribution, which
        # ignores the few chains that wander furthest.
        lo = np.percentile(p0, 2, axis=0)
        hi = np.percentile(p0, 98, axis=0)
        centre = 0.5 * (lo + hi)
        span = float(np.max(hi - lo)) * 0.62
    else:
        span = span if span is not None else max(reach * st["zoom"], 5.0)
    vp.camera_pos = tuple(centre - d * span / np.tan(vp.fov / 2.0))
    vp.render_image(filename=str(out_png), size=(size, size),
                    background=(1.0, 1.0, 1.0),
                    renderer=TachyonRenderer(ambient_occlusion=True,
                                             ambient_occlusion_samples=12,
                                             shadows=True))
    pipe.remove_from_scene()
    return out_png


def sheet(items, out_png, title, cols=3):
    from PIL import Image, ImageDraw, ImageFont

    ims = [(Image.open(p), lab, note) for p, lab, note in items]
    w, h = ims[0][0].size
    pad, head, cap = 14, 70, 62
    rows = (len(ims) + cols - 1) // cols
    W = cols * w + (cols + 1) * pad
    H = head + rows * (h + cap) + pad
    out = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(out)

    def font(sz):
        for n in ("segoeui.ttf", "arial.ttf", "DejaVuSans.ttf"):
            try:
                return ImageFont.truetype(n, sz)
            except OSError:
                continue
        return ImageFont.load_default()

    d.text((pad, 20), title, fill="black", font=font(32))
    for i, (im, lab, note) in enumerate(ims):
        r, c = divmod(i, cols)
        x = pad + c * (w + pad)
        y = head + r * (h + cap)
        out.paste(im, (x, y))
        d.text((x + 6, y + h + 6), lab, fill="black", font=font(26))
        d.text((x + 6, y + h + 34), note, fill=(90, 90, 90), font=font(20))
    out.save(out_png)
    return out_png


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lattice", default="SC")
    ap.add_argument("--orient", default="any")
    ap.add_argument("--shell", type=int, default=2)
    ap.add_argument("--size", type=int, default=900)
    args = ap.parse_args()

    root = (Path(__file__).resolve().parents[2]
            / "tests/output/entangle_steps/shell_gallery")
    src = root / f"{args.lattice}_{args.orient}_shell{args.shell}.data"
    if not src.exists():
        print(f"  no such file: {src}")
        return 1
    out = root / "catalogue"
    out.mkdir(parents=True, exist_ok=True)

    items = []
    for key, st in STYLES.items():
        png = out / f"{key}.png"
        print(f"  {st['label']} ...", flush=True)
        render(src, png, st, args.size)
        items.append((png, st["label"], st["note"]))

    final = out / f"catalogue_{args.lattice}_shell{args.shell}.png"
    sheet(items, final,
          f"Render styles, {args.lattice} shell {args.shell}: "
          f"teal chain wound twice around the red partner")
    print(f"\n  {final}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
