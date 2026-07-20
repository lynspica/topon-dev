"""Render the topon gallery.

THE STORY IS THE LATTICE. topon places every bead/atom on a lattice: chains are
drawn taut between junctions, each edge of the graph strung with `dp` beads. That
state is not physical: CG bonds sit at 0.20 sigma against a 0.97 equilibrium
(4.9x short) and atomistic at a 0.45 A median against a 1.09 A C-H (2.8x short).
MD is what makes it real -- the chains coil out and the lattice becomes a melt.
The gallery shows that arc for both resolutions, then the design knobs on the
lattice state.

Sections
  1. lattice -> pushed-off/ramped -> equilibrated, for CG and atomistic
  2. copolymer sequences (lattice, CG beads)
  3. entanglements + grafts (full network + zoom)

PROVENANCE (deliberately mixed -- see assets/gallery/README.md):
  * The CG arc is generated FRESH by gen_systems.py and run through LAMMPS.
    The copolymer and entanglement panels are also fresh but show only the
    03_Conformation lattice state -- no MD is run for them.
  * The v21_cg_combined golden is NOT used: 190 of its bonds sit at 7.45 sigma
    (= sqrt(3) x its 4.2996 spacing) and 10 at 10.53, with the bead-ends pinned
    to lattice nodes. Fresh output shows nothing like it. Note the golden LOADS a
    210-edge sculpted graph rather than generating a 375-edge lattice, and its
    config no longer exists, so it cannot be re-run for a controlled comparison.
  * The atomistic arc uses the v21_atomistic_combined golden, which is clean
    (0 stretched bonds after minimisation) and already carries the full MD arc;
    regenerating it means a ~10M-iteration PPPM minimisation.

Run with:  C:/v/ovito/Scripts/python.exe render_gallery.py [panel ...]
"""
import os
import re
import sys
from pathlib import Path

import numpy as np
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import WrapPeriodicImagesModifier

HERE = Path(__file__).resolve().parent      # assets/gallery
REPO = HERE.parents[1]
OUT = HERE                                  # renders sit beside this script
GOLD = REPO / "tests" / "output"
# Generated systems are big and transient, so they live outside the repo.
# gen_systems.py writes them here; override with TOPON_GALLERY_SYSTEMS.
FRESH = Path(os.environ.get("TOPON_GALLERY_SYSTEMS", HERE / "systems"))

# House semantic colours
A_RED = (0.690, 0.227, 0.180)      # #B03A2E  monomer A / backbone
B_BLUE = (0.227, 0.463, 0.690)     # #3A76B0  monomer B
JUNCT = (0.145, 0.196, 0.251)      # #252F40  lattice junction
TEAL = (0.055, 0.486, 0.420)       # #0e7c6b  graft / side chain
ENT_GOLD = (0.902, 0.639, 0.153)   # #E6A327  entangled chain 1
VIOLET = (0.478, 0.271, 0.522)     # #7A4585  entangled chain 2
REST = (0.55, 0.58, 0.63)          # #8C949F  the whole rest of the network
QUIET = (0.788, 0.808, 0.831)      # #c9ced4

SI = (0.94, 0.78, 0.31)
O_ = (0.83, 0.28, 0.22)
C_ = (0.36, 0.40, 0.45)
N_ = (0.23, 0.46, 0.69)
H_ = (0.88, 0.89, 0.91)
F_ = (0.42, 0.72, 0.38)

# Type IDs do NOT mean the same element across systems, so element colours are
# keyed off each file's own mass table.
BY_MASS = [(28.0855, SI, 0.62), (18.9984, F_, 0.40), (15.9994, O_, 0.46),
           (14.0067, N_, 0.44), (12.0110, C_, 0.42), (1.0079, H_, 0.24)]


def element_palette(probe):
    cmap, rmap = {}, {}
    for t in probe.particles.particle_types.types:
        best = min(BY_MASS, key=lambda e: abs(e[0] - t.mass))
        cmap[t.id], rmap[t.id] = ((best[1], best[2]) if abs(best[0] - t.mass) < 0.5
                                  else (C_, 0.42))
    return cmap, rmap


def node_ids(sysdir):
    """Junction atom IDs, straight from topon's own LAMMPS group file.

    Necessary, not cosmetic: in the copolymer systems the junctions and the A
    monomers are BOTH particle type 1, so types alone cannot separate them.
    """
    g = (sysdir / "02_Chemistry" / "system.groups").read_text()
    m = re.search(r"group nodes id ([\d\s]+)", g)
    return set(int(x) for x in m.group(1).split())


def fresh(name, stage):
    return FRESH / name / name / stage


# rad = bead radius, FIXED ACROSS EACH ARC ROW so the viewer sees the real
# change (taut lattice -> coiled melt) rather than a radius trick. Bond width is
# deliberately >= the 0.198 bead spacing so a strand renders as a solid rod, not
# a dotted line.
CG_R, CG_BW = 0.11, 0.22
AT_R = 0.32

# The whole gallery is built on a HETEROGENEOUS network, not an ideal lattice:
# topon's strict-sculpting stage thins the 375-edge 5x5x5 SC lattice to 250 edges
# (degree_distribution = "e:250"), giving junction functionality 2-6 with a mean
# of 4.0 -- a tetrafunctional crosslinked network, the standard elastomer
# connectivity -- instead of every node at 6. That is the point of topon, and a
# perfect grid hides it.
PANELS = {
    # ---- 1. the arc (sculpt_250: mean functionality 4.0) ----------------
    # mono=True: the whole homopolymer network is ONE colour (junctions by size,
    # not colour). Topological colouring is reserved for the entanglement panel;
    # chemistry colouring (copolymer A/B, atomistic elements) is kept elsewhere.
    "cg_lattice": dict(src=fresh("sculpt_250", "03_Conformation/system_relaxed.data"),
                       mode="cg", sysdir=fresh("sculpt_250", ""), rad=CG_R, bw=CG_BW,
                       asbuilt=True, mono=True),
    # End of stages 1-2, NOT system_after_soft. Stage 1 is a complete NO-OP:
    # every `minimize` in the CG deck stops after 1 iteration because
    # lammps_inputs.py hardcodes etol=1e-4 while LAMMPS tests |dE|/(|E1|+|E2|),
    # which at E~3.9e6 already scores ~4e-7. Max displacement from as-built:
    # 2e-4 sigma. The bonds reach 0.96 via the `run 20000` NVE ramp in stage 2,
    # so this is "pushed off and ramped", not "minimised" (it ends at T=18.4).
    "cg_minimised": dict(src=fresh("sculpt_250", "04_Simulation/system_ramped.data"),
                         mode="cg", sysdir=fresh("sculpt_250", ""), rad=CG_R, bw=CG_BW,
                         mono=True),
    "cg_equilibrated": dict(src=fresh("sculpt_250", "04_Simulation/system_equilibrated.data"),
                            mode="cg", sysdir=fresh("sculpt_250", ""), rad=CG_R, bw=CG_BW,
                            mono=True),
    # Atomistic row is ALSO sculpted now (3x3x3 PDMS, 81 edges -> 54, mean f 4.0),
    # generated fresh and run through the atomistic MD arc locally, so the whole
    # gallery sits on heterogeneous networks. element_palette keys off the mass
    # table, which here is 1=Si_ (dummy, m=1) / 2=Si / 3=C / 4=O / 5=H -- a
    # different order from the golden, handled automatically.
    "atom_lattice": dict(
        src=fresh("atom_sculpt", "03_Conformation/system_relaxed.data"),
        mode="elements", rad=AT_R, bw=0.22),
    "atom_minimised": dict(
        src=fresh("atom_sculpt", "04_Simulation/system_minimized_final.data"),
        mode="elements", rad=AT_R, bw=0.22),
    "atom_equilibrated": dict(
        src=fresh("atom_sculpt", "04_Simulation/system_equilibrated.data"),
        mode="elements", rad=AT_R, bw=0.22),
    # ---- 2. copolymer sequences (sculpted 4x4x4, lattice state) ----------
    # Also sculpted (192 edges -> 128, mean functionality 4.0), so these are
    # heterogeneous networks too -- just coarse enough that each strand's A/B
    # pattern still reads.
    #
    # NO gradient PANEL. `arrangement: "gradient"` IGNORES the requested
    # composition: sequences.py scales the weight window by n, so for n=2 the
    # ranges [0,0.5] and [0.5,1] never overlap and nothing blends -- every
    # composition comes out 50:50 (ask A=0.1, get A=0.50). At equal fractions and
    # even dp that is byte-identical to `block`, i.e. the same picture twice.
    **{f"copoly_{k}": dict(src=fresh(f"copoly_{k}", "03_Conformation/system_relaxed.data"),
                           mode="cg", sysdir=fresh(f"copoly_{k}", ""),
                           rad=CG_R, bw=CG_BW, asbuilt=True)
       for k in ("block", "random", "alternating")},
    # ---- 3. entanglements + grafts (sculpted 5x5x5, e:250) --------------
    # ONE entanglement is highlighted: its two chains in two distinct colours,
    # the entire rest of the network (backbone, other strands, grafts, junctions)
    # in a single third colour. Same scheme in the full-network view and the zoom.
    "ent_full": dict(src=fresh("entangled_grafted", "03_Conformation/system_relaxed.data"),
                     mode="cg", sysdir=fresh("entangled_grafted", ""),
                     rad=CG_R, bw=CG_BW, pair_scheme=True, asbuilt=True),
    "ent_zoom": dict(src=fresh("entangled_grafted", "03_Conformation/system_relaxed.data"),
                     mode="cg", sysdir=fresh("entangled_grafted", ""),
                     rad=0.075, bw=0.11, zoom=True, keep=4.0, pair_scheme=True, asbuilt=True),
}


def hq_renderer():
    """The gallery's high-quality Tachyon settings.

    Pair with SUPERSAMPLE: render at SUPERSAMPLE x the final size and downscale
    with LANCZOS. That is the single biggest quality win -- true anti-aliasing on
    every sphere and cylinder silhouette, which Tachyon's own AA only partly
    gives. Costs ~5x the render time per frame (0.4s -> 2.3s at 560px), which is
    irrelevant for a 30-frame demo.
    """
    return TachyonRenderer(
        ambient_occlusion=True, ambient_occlusion_samples=26,
        ambient_occlusion_brightness=0.85,
        antialiasing=True, antialiasing_samples=10,
        direct_light=True, direct_light_intensity=0.95,
        shadows=True)


SUPERSAMPLE = 2


def render_hq(vp, path, size, renderer=None, **kw):
    """Render supersampled, then downscale to `size` (a (w, h) tuple)."""
    from PIL import Image
    w, h = size
    vp.render_image(filename=str(path), size=(w * SUPERSAMPLE, h * SUPERSAMPLE),
                    background=(1, 1, 1), renderer=renderer or hq_renderer(),
                    **kw)
    Image.open(path).convert("RGB").resize((w, h), Image.LANCZOS).save(path)


def median_bond(data):
    topo = np.array(data.particles.bonds.topology)
    p = np.array(data.particles.positions)
    cell = np.array(data.cell[:, :3])
    f = (p[topo[:, 0]] - p[topo[:, 1]]) @ np.linalg.inv(cell).T
    f -= np.round(f)
    return float(np.median(np.linalg.norm(f @ cell.T, axis=1)))


def rebuild_bond_pbc(frame, data):
    """Recompute each bond's PBC shift from the minimum image.

    OVITO fixes these when the file is READ; any later position change (wrapping)
    leaves them describing the old geometry and every boundary-crossing bond gets
    drawn as a ray straight across the box.
    """
    if data.particles.bonds is None or "Periodic Image" not in data.particles.bonds:
        return
    p = np.array(data.particles.positions)
    cell = np.array(data.cell[:, :3])
    topo = np.array(data.particles.bonds.topology)
    f = (p[topo[:, 1]] - p[topo[:, 0]]) @ np.linalg.inv(cell).T
    shift = -np.round(f)
    shift[:, ~np.array(data.cell.pbc, dtype=bool)] = 0
    data.particles_.bonds_["Periodic Image_"][...] = shift.astype(np.int32)


def half_cell_offset(probe, nodes):
    """Shift that moves the lattice planes off the box faces.

    THIS IS WHY THE LATTICE PANELS LOOKED "DOTTED". The as-built lattice puts
    nodes at x=0 -- exactly ON the periodic face. The conformation stage then
    jitters every bead by +-1e-4, so wrapping throws the beads that land at
    -1e-4 across to x~L while their neighbours stay at x~0. Measured on
    cg_basic: 2 761 of 7 625 beads sit within 0.01 of a face and 1 596 of 7 875
    bonds (20%) end up spanning the box, each drawn as two stubs leaving
    opposite faces. Whole strands render as dashes -- and no bead radius or bond
    width fixes it, because the bonds genuinely aren't there.

    Shifting by half a lattice spacing puts every plane in the interior. The
    picture is identical, just translated.
    """
    p = np.array(probe.particles.positions)
    ident = np.array(probe.particles["Particle Identifier"])
    q = p[np.isin(ident, list(nodes))]
    if len(q) < 2:
        return np.zeros(3)
    cell = np.array(probe.cell[:, :3])
    inv = np.linalg.inv(cell)
    # nearest-neighbour distance among junctions = the lattice spacing
    d = q[:, None, :] - q[None, :, :]
    f = d @ inv.T
    f -= np.round(f)
    dist = np.linalg.norm(f @ cell.T, axis=2)
    np.fill_diagonal(dist, np.inf)
    a = float(np.median(dist.min(axis=1)))
    return np.full(3, 0.5 * a)


def walk_chains(data, nodes):
    """Every chain between two junctions, with how far it bows off its own axis.

    topon exports no entanglement metadata, so a kink is found the way the eye
    finds it. An un-entangled chain is drawn dead straight along its lattice edge
    (bow = 0.00); an entangled one is pushed toward its partner and bows ~3 sigma.
    The split is unambiguous -- in entangled_grafted, 351 of 375 chains bow 0.00
    and exactly 24 bow > 1.0, which is 2 x the 12 entanglements the config asked
    for: BOTH partners of every pair get kinked.

    Returns [(bow, [atom indices], unwrapped points)].
    """
    p = np.array(data.particles.positions)
    cell = np.array(data.cell[:, :3])
    inv = np.linalg.inv(cell)
    ident = np.array(data.particles["Particle Identifier"])
    topo = np.array(data.particles.bonds.topology)
    isnode = np.isin(ident, list(nodes))
    adj = [[] for _ in range(data.particles.count)]
    for a, b in topo:
        adj[a].append(b)
        adj[b].append(a)

    def mic(v):
        f = v @ inv.T
        return (f - np.round(f)) @ cell.T

    # Strip the grafts before walking. A backbone bead carrying a side chain has
    # degree 3, and a walk that refuses degree-3 beads stops dead there --
    # splitting one chain into two fragments that then look like two chains
    # "0.40 sigma apart" when they are literally the same beads. Side chains are
    # dead-end branches, so repeatedly peeling degree-1 non-junction beads
    # removes them entirely and leaves a clean degree-2 backbone.
    deg = np.array([len(a) for a in adj])
    pruned = np.zeros(data.particles.count, dtype=bool)
    while True:
        leaves = [i for i in range(data.particles.count)
                  if not pruned[i] and not isnode[i] and deg[i] == 1]
        if not leaves:
            break
        for i in leaves:
            pruned[i] = True
            for n in adj[i]:
                if not pruned[n]:
                    deg[n] -= 1
        for i in leaves:
            deg[i] = 0

    out, seen = [], set()
    for start in np.where(isnode)[0]:
        for first in adj[start]:
            if isnode[first] or first in seen or pruned[first]:
                continue
            chain, prev, cur = [first], start, first
            while len(chain) < 400:
                nxt = [n for n in adj[cur]
                       if n != prev and not isnode[n] and not pruned[n]]
                if not nxt:
                    break
                prev, cur = cur, nxt[0]
                chain.append(cur)
            if len(chain) < 5:
                continue
            seen.update(chain)
            pts = [p[chain[0]]]
            for k in range(1, len(chain)):
                pts.append(pts[-1] + mic(p[chain[k]] - p[chain[k - 1]]))
            pts = np.array(pts)
            axis = pts[-1] - pts[0]
            n = np.linalg.norm(axis)
            if n < 1e-6:
                continue
            u = axis / n
            rel = pts - pts[0]
            bow = float(np.linalg.norm(rel - np.outer(rel @ u, u), axis=1).max())
            out.append((bow, chain, pts))
    return out


def closest_entangled_pair(data, nodes, thresh=1.0):
    """The two kinked chains that neck closest together, and their contact point.

    That contact is what an entanglement IS here, so it is what the zoom frames.
    """
    cell = np.array(data.cell[:, :3])
    inv = np.linalg.inv(cell)

    def mic(v):
        f = v @ inv.T
        return (f - np.round(f)) @ cell.T

    kinked = sorted((c for c in walk_chains(data, nodes) if c[0] > thresh),
                    key=lambda c: -c[0])[:12]
    best = (1e9, None, None, None)
    for i in range(len(kinked)):
        for j in range(i + 1, len(kinked)):
            a, b = kinked[i][2], kinked[j][2]
            for x in a[::2]:
                d = np.linalg.norm(mic(b[::2] - x), axis=1)
                k = int(d.argmin())
                if d[k] < best[0]:
                    best = (float(d[k]), x, b[::2][k], (kinked[i], kinked[j]))
    dmin, xa, xb, pair = best
    contact = xa + 0.5 * mic(xb - xa)      # midpoint of the closest approach
    return dmin, contact, pair


def render(name, size=(900, 900)):
    spec = PANELS[name]
    src = Path(spec["src"])
    if not src.exists():
        print(f"[SKIP] {name}: missing {src}")
        return
    mode = spec["mode"]
    pipe = import_file(str(src), atom_style="full")
    probe = pipe.compute()
    mb = median_bond(probe)
    elem_c, elem_r = element_palette(probe) if mode == "elements" else ({}, {})
    nodes = node_ids(Path(spec["sysdir"])) if mode == "cg" else set()

    # ONE entanglement highlighted: its two chains distinct, the rest uniform.
    # Both the full-network and zoom views use the same focus pair.
    pair_a = pair_b = contact = None
    if spec.get("pair_scheme"):
        dmin, contact, pair = closest_entangled_pair(probe, nodes)
        pair_a, pair_b = pair[0][1], pair[1][1]
        print(f"    focus entanglement: 2 chains {dmin:.2f} sigma apart", flush=True)

    # Move the lattice planes off the box faces before any wrapping, or 20% of
    # the bonds render as stubs and the strands come out dashed. ONLY for the
    # as-built lattice states -- after MD the junctions have left their sites, so
    # there are no planes on the faces and the shift would be meaningless.
    if mode == "cg" and spec.get("asbuilt"):
        off = half_cell_offset(probe, nodes)
        if np.any(off):
            def deboundary(frame, data, s=off):
                data.particles_.positions_[...] = data.particles.positions + s
            pipe.modifiers.append(deboundary)

    cell = np.array(probe.cell[:, :3])
    if spec.get("zoom"):
        # Slide the contact to the box centre BEFORE wrapping, so the whole
        # neighbourhood ends up contiguous. Cropping on minimum-image distance
        # instead drags in beads from the far side of the box, which render as a
        # detached island floating next to the kink.
        shift = 0.5 * (cell[0] + cell[1] + cell[2]) - contact

        def recentre(frame, data, s=shift):
            data.particles_.positions_[...] = data.particles.positions + s
        pipe.modifiers.append(recentre)

    pipe.modifiers.append(WrapPeriodicImagesModifier())
    pipe.modifiers.append(rebuild_bond_pbc)

    def paint(frame, data):
        t = np.array(data.particles["Particle Type"])
        ident = np.array(data.particles["Particle Identifier"])
        n = data.particles.count
        col = np.tile(np.array(C_), (n, 1))
        rad = np.full(n, spec["rad"])
        if mode == "elements":
            for tid, c in elem_c.items():
                m = t == tid
                if m.any():
                    col[m], rad[m] = c, elem_r[tid] * spec["rad"] / 0.42
        elif pair_a is not None:
            # ONE entanglement: two chains in two distinct colours, the ENTIRE
            # rest of the network (backbone, other strands, grafts, junctions) in
            # a single third colour. Junctions stay larger so they read as
            # landmarks, but they carry no separate colour.
            isnode = np.isin(ident, list(nodes))
            col[:] = REST
            rad[isnode] = spec["rad"] * 2.4
            for idx, c in ((pair_a, ENT_GOLD), (pair_b, VIOLET)):
                m = np.zeros(n, bool)
                m[np.array(idx)[np.array(idx) < n]] = True
                col[m] = c
                rad[m] = spec["rad"] * 1.9
        elif spec.get("mono"):
            # Homopolymer arc: one colour for the whole network. Junctions are
            # larger, not a different colour.
            isnode = np.isin(ident, list(nodes))
            col[:] = B_BLUE
            rad[isnode] = spec["rad"] * 2.4
        else:
            isnode = np.isin(ident, list(nodes))
            # Among non-junction beads the type IS the design knob: for a
            # copolymer 1/2 = A/B; with grafts 2/3 = backbone/side chain.
            col[~isnode & (t == 1)] = A_RED
            col[~isnode & (t == 2)] = A_RED if 3 in set(t) else B_BLUE
            col[~isnode & (t == 3)] = TEAL
            col[isnode] = JUNCT
            rad[isnode] = spec["rad"] * 2.4
        data.particles_.create_property("Color", data=col)
        data.particles_.create_property("Radius", data=rad)

    pipe.modifiers.append(paint)

    if spec.get("zoom"):
        keep = spec.get("keep", 4.0)
        cen = 0.5 * (cell[0] + cell[1] + cell[2])

        def crop(frame, data, c=cen, r=keep):
            d = np.linalg.norm(np.array(data.particles.positions) - c, axis=1)
            data.particles_.delete_elements((d > r).astype(np.int8))
        pipe.modifiers.append(crop)

    pipe.add_to_scene()
    d = pipe.compute()
    if d.particles.bonds is not None:
        d.particles.bonds.vis.width = spec["bw"]
        d.particles.bonds.vis.use_particle_colors = True
    if spec.get("zoom"):
        d.cell.vis.enabled = False

    vp = Viewport(type=Viewport.Type.Perspective)
    vp.camera_dir = (-0.55, -0.75, -0.45)
    vp.zoom_all(size=size)
    r = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_samples=12,
                        shadows=True)
    out = OUT / f"{name}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    vp.render_image(filename=str(out), size=size, background=(1, 1, 1), renderer=r)
    pipe.remove_from_scene()
    box = np.diag(np.array(probe.cell[:, :3]))[0]
    print(f"[OK] {name}: {d.particles.count} shown, box {box:.1f}, "
          f"median bond {mb:.2f} -> {out.name}", flush=True)


if __name__ == "__main__":
    for nm in (sys.argv[1:] or list(PANELS)):
        render(nm)
    print("gallery done")
