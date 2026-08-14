"""Plain lattice catalogue: the networks themselves, nothing entangled.

    "C:/v/ovito/Scripts/python.exe" tests/workflows/render_lattices.py

Seven panels on one sheet: the four regular lattices (SC, Diamond, BCC, FCC)
and three mixtures at different SC:BCC:FCC ratios. Chains are drawn straight
so what is on show is the connectivity, and every panel is framed on its own
box.

Default OVITO look on purpose: no colour overrides, so the panels match what
the user sees opening the same file in the OVITO GUI. The one styling change
is the display radius (0.22 instead of OVITO's 0.5, which fuses the beads
into a solid block). Bonds are the file's own; bonds that span the periodic
cell are dropped rather than drawn across the picture.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

# Imports live inside build() and render(): building needs topon (and rdkit),
# rendering needs ovito, and no one interpreter here has both. The script runs
# twice -- --build-only under the project python, --render-only under OVITO's.
OUT = ROOT / "tests/output/entangle_steps"

# "max functionality", not "functionality": at dims 3 the sculpted SC and
# FCC graphs realise mean degrees 3.6 and 7.0, so an unqualified
# "functionality 4/8" caption contradicts the built graph. BCC and Diamond
# do come out uniform.
PANELS = [
    ("SC",      None,             "regular, max functionality 4"),
    ("DIAMOND", None,             "regular, max functionality 4"),
    ("BCC",     None,             "regular, max functionality 8"),
    ("FCC",     None,             "regular, max functionality 8"),
    ("MIX", (0.8, 0.1, 0.1),    "mix 80/10/10 SC-rich"),
    ("MIX", (0.34, 0.33, 0.33), "mix 34/33/33 even"),
    ("MIX", (0.1, 0.1, 0.8),    "mix 10/10/80 FCC-rich"),
]


def build(lat, mix, dims, dp, coil, root):
    from tests.workflows.entangle_all import CASES
    from tests.workflows.entangle_steps import (BOND, LATTICE, build_network,
                                                geometry, write_system)
    from topon.conformation.paths import straight

    spec = dict(LATTICE)
    if lat == "DIAMOND":
        # Not in the shared CASES table: that table is what entangle_all
        # sweeps, and adding Diamond there would silently widen every sweep.
        spec.update(dict(lattice="DIAMOND", max_func=4,
                         degree_dist="0:0,1:0", mix={}))
    else:
        spec.update(CASES[lat])
    spec["dims"] = (dims,) * 3
    if mix is not None:
        spec["mix"] = {"SC": mix[0], "BCC": mix[1], "FCC": mix[2]}
        spec["lattice"] = "MIX"
        # Forbid only the unusable degrees (0 and 1); everything else is
        # left to the sites. The MIX case's degree list demands degree 7 and 8
        # nodes, and a mix with few BCC/FCC sites has almost
        # none that can carry them: the rejection sampler then
        # never converges, and both this build and the lattice
        # catalogue hung on it for half an hour or more.
        spec["degree_dist"] = "0:0,1:0"
    graph = build_network(spec)
    geo = geometry(graph, dp=dp, bond=BOND, coil=coil)
    paths = {k: straight(c0, c1, dp + 1)
             for k, (c0, c1) in geo["chords"].items()}
    shutil.rmtree(root, ignore_errors=True)
    _n, node_atom, chain_atoms = write_system(graph, geo, paths, root)
    f = root / "01_Topology" / "system.data"
    if not f.exists():
        f = next(root.rglob("system.data"))

    # write_system writes zero coordinates; the conform stage this build
    # skips is what applies them. Left as zeros, every atom sits at the
    # origin: Tachyon grinds on thousands of coincident spheres and the
    # percentile camera degenerates to a point, which is why the catalogue
    # render ran for forty minutes without producing a panel.
    from tests.workflows.entangle_relaxed import rewrite_coords
    from tests.workflows.entangle_steps import chain_ids
    xyz = {}
    for k in sorted(geo["chords"]):
        for aid, q in zip(chain_ids(k, node_atom, chain_atoms, geo["ends"]),
                          paths[k]):
            xyz[aid] = q
    # Wrap every bead into the cell. The straight paths run in each chain's
    # own minimum image, so a junction shared by a wrapped and an unwrapped
    # chain was written at whichever image came last and the other chains'
    # bonds to it were torn -- a quarter to a half of the chains per panel
    # rendered as stubs. Wrapping collapses all images of a junction to the
    # one in-cell position; boundary-crossing chains then show as two stubs
    # at opposite faces, which is what cut_wrapped's bond cut is for.
    L = np.asarray(geo["L"], float)
    xyz = {a: q - L * np.floor(q / L) for a, q in xyz.items()}
    rewrite_coords(f, f, xyz)

    # The raw file carries write_system's placeholder box (+-50); the render
    # draws the cell, so the header has to be the real periodic box or the
    # contour dwarfs the network.
    txt = f.read_text().splitlines()
    for j, ln in enumerate(txt):
        s = ln.strip()
        if s.endswith("xlo xhi"):
            txt[j] = f"0.0 {geo['L'][0]:.6f} xlo xhi"
        elif s.endswith("ylo yhi"):
            txt[j] = f"0.0 {geo['L'][1]:.6f} ylo yhi"
        elif s.endswith("zlo zhi"):
            txt[j] = f"0.0 {geo['L'][2]:.6f} zlo zhi"
    f.write_text("\n".join(txt) + "\n")
    return f, graph.number_of_edges()


def render(data_file, out_png, size):
    from ovito.io import import_file
    from ovito.modifiers import DeleteSelectedModifier
    from ovito.vis import TachyonRenderer, Viewport

    pipe = import_file(str(data_file), atom_style="full")

    def cut_wrapped(frame, data):
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
    pipe.add_to_scene()

    data = pipe.compute()
    # Colours and the cell wireframe stay OVITO's defaults. Only the display
    # radius moves: at the default 0.5 the beads fuse into a solid block and
    # the topology is the one thing the picture does not show.
    data.particles.vis.radius = 0.22

    # Frame the whole network plus the simulation box. Fitting a sphere
    # around the full extent keeps every junction and the cell contour in
    # view; percentile framing cropped the panels mid-network, which hid
    # exactly the thing on show.
    cell = np.asarray(data.cell[...])
    corners = np.array([cell[:, 3] + i * cell[:, 0] + j * cell[:, 1]
                        + k * cell[:, 2]
                        for i in (0, 1) for j in (0, 1) for k in (0, 1)])
    pos = np.vstack([data.particles.positions.array, corners])
    lo, hi = pos.min(axis=0), pos.max(axis=0)
    centre = 0.5 * (lo + hi)
    radius = 0.5 * float(np.linalg.norm(hi - lo)) + 0.6

    d = np.array([-0.55, -0.75, -0.45], float)
    d = d / np.linalg.norm(d)
    vp = Viewport(type=Viewport.Type.Perspective)
    vp.camera_dir = tuple(d)
    vp.fov = np.deg2rad(35.0)
    vp.camera_pos = tuple(centre - d * radius / np.sin(vp.fov / 2.0))
    vp.render_image(filename=str(out_png), size=(size, size),
                    background=(1.0, 1.0, 1.0),
                    renderer=TachyonRenderer(ambient_occlusion=False,
                                             shadows=False))
    pipe.remove_from_scene()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dims", type=int, default=3,
                    help="cells per side. 3 keeps the render legible and "
                         "fast; the lattices repeat, so more cells add beads "
                         "without adding information.")
    ap.add_argument("--dp", type=int, default=20)
    ap.add_argument("--coil", type=float, default=8.0)
    ap.add_argument("--size", type=int, default=850)
    ap.add_argument("--build-only", action="store_true")
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()

    out = OUT / "lattice_catalogue"
    out.mkdir(parents=True, exist_ok=True)

    items = []
    for i, (lat, mix, note) in enumerate(PANELS):
        tag = lat if mix is None else f"MIX_{int(mix[0]*100)}"
        root = out / f"build_{i}_{tag}"
        marker = root / "chains.txt"
        if not args.render_only and not marker.exists():
            print(f"  {tag}: building ...", flush=True)
            f, n = build(lat, mix, args.dims, args.dp, args.coil, root)
            marker.write_text(f"{f}\n{n}\n")
        if args.build_only:
            continue
        f, n = marker.read_text().splitlines()
        png = out / f"{i}_{tag}.png"
        print(f"  {tag}: rendering {n} chains ...", flush=True)
        render(f, png, args.size)
        items.append((png, tag if mix is None else
                      f"MIX {int(mix[0]*100)}/{int(mix[1]*100)}/"
                      f"{int(mix[2]*100)}", f"{note}, {n} chains"))
    if args.build_only:
        print("  built; now run --render-only under OVITO python")
        return 0

    from PIL import Image, ImageDraw, ImageFont
    ims = [(Image.open(p), lab, nt) for p, lab, nt in items]
    w, h = ims[0][0].size
    cols = 4
    rows = (len(ims) + cols - 1) // cols
    pad, head, cap = 14, 70, 58
    W = cols * w + (cols + 1) * pad
    H = head + rows * (h + cap) + pad
    sheet = Image.new("RGB", (W, H), "white")
    dr = ImageDraw.Draw(sheet)

    def font(sz):
        for nm in ("segoeui.ttf", "arial.ttf", "DejaVuSans.ttf"):
            try:
                return ImageFont.truetype(nm, sz)
            except OSError:
                continue
        return ImageFont.load_default()

    dr.text((pad, 20), "Lattice catalogue: regular networks and three "
            "mixtures (chains drawn straight)",
            fill="black", font=font(30))
    for i, (im, lab, nt) in enumerate(ims):
        r, c = divmod(i, cols)
        x = pad + c * (w + pad)
        y = head + r * (h + cap)
        sheet.paste(im, (x, y))
        dr.text((x + 6, y + h + 6), lab, fill="black", font=font(26))
        dr.text((x + 6, y + h + 34), nt, fill=(90, 90, 90), font=font(19))
    final = out / "lattice_catalogue.png"
    sheet.save(final)
    print(f"\n  {final}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
