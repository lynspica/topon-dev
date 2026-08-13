"""Showcase sheets: what the design stage can be asked for.

    "C:/v/ovito/Scripts/python.exe" tests/workflows/render_showcase.py --style D6

Four sheets, each answering one question a reader would ask:

    shells     which neighbour is entangled, and what "neighbour" means
    lattices   whether it works away from the simple cubic case
    counts     whether the number of windings is a knob
    orient     whether parallel and perpendicular partners are distinguishable

Every panel in a sheet shares one camera and one scale, so what differs
between them is the design rather than the framing.

Colours: teal is the routed chain, red its partner, and the large dark spheres
are the crosslinks each chain is anchored to -- without those a winding reads
as a squiggle rather than a design between two fixed points.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from render_catalogue import (STYLES, reach_of, render,  # noqa: E402
                              sheet)

ROOT = (Path(__file__).resolve().parents[2]
        / "tests/output/entangle_steps/shell_gallery")

GAPS = {"SC": {1: "0.00", 2: "1.00", 3: "1.41", 4: "1.73"},
        "BCC": {1: "0.00", 2: "0.71", 3: "0.82", 4: "0.87"},
        "FCC": {1: "0.00", 2: "0.53", 3: "0.60", 4: "0.64"},
        "MIX": {1: "0.00", 2: "0.45", 3: "0.57", 4: "0.67"}}


def common_span(paths, style):
    """One scale for a sheet, from its largest pair."""
    r = [reach_of(p) for p in paths if p.exists()]
    return max(max(r) * STYLES[style]["zoom"], 5.0) if r else None


def sheet_shells(style, out, size):
    span = common_span([ROOT / f'SC_any_shell{s}.data'
                        for s in (1, 2, 3, 4)], style)
    items = []
    for sh in (1, 2, 3, 4):
        src = ROOT / f"SC_any_shell{sh}.data"
        if not src.exists():
            continue
        png = out / f"showcase_shell{sh}.png"
        render(src, png, STYLES[style], size, span)
        note = ("they share a crosslink" if sh == 1
                else f"{GAPS['SC'][sh]} lattice units apart")
        items.append((png, f"shell {sh}", note))
    return items, ("Which neighbour is entangled. A shell is every chain at "
                   "the same closest approach; on SC those are 0, 1, "
                   "sqrt(2), sqrt(3).")


def sheet_lattices(style, out, size):
    span = common_span([ROOT / f'{l}_any_shell2.data'
                        for l in ('SC', 'BCC', 'FCC', 'MIX')], style)
    items = []
    for lat in ("SC", "BCC", "FCC", "MIX"):
        src = ROOT / f"{lat}_any_shell2.data"
        if not src.exists():
            continue
        png = out / f"showcase_{lat}.png"
        render(src, png, STYLES[style], size, span)
        items.append((png, lat, f"shell 2 at {GAPS[lat][2]} lattice units"))
    return items, ("The same request on four lattices. SC shells are widely "
                   "spaced; BCC, FCC and MIX are bunched within a few "
                   "hundredths.")


def sheet_counts(style, out, size):
    span = common_span([ROOT / f'count{c}_SC_shell2.data'
                        for c in (1, 2)], style)
    items = []
    for c in (1, 2):
        src = ROOT / f"count{c}_SC_shell2.data"
        if not src.exists():
            continue
        png = out / f"showcase_count{c}.png"
        render(src, png, STYLES[style], size, span)
        items.append((png, f"{c} entanglement" + ("s" if c > 1 else ""),
                      "same pair, same shell, different winding"))
    return items, ("How many. The count is set by how far round the partner "
                   "the chain is taken, and is verified after minimisation "
                   "rather than assumed.")


def sheet_orient(style, out, size):
    span = common_span([ROOT / f'SC_{o}_shell2.data'
                        for o in ('parallel', 'perpendicular')], style)
    items = []
    for o in ("parallel", "perpendicular"):
        src = ROOT / f"SC_{o}_shell2.data"
        if not src.exists():
            continue
        png = out / f"showcase_{o}.png"
        render(src, png, STYLES[style], size, span)
        items.append((png, o, "shell 2, 1.00 lattice units either way"))
    return items, ("A shell is a distance, not an orientation. On SC shell 2 "
                   "holds 740 parallel pairs and 1898 perpendicular, and "
                   "either can be asked for.")


def sheet_lattice_shells(style, out, size, lat):
    """One lattice, all four shells: the shell control is not SC-only."""
    span = common_span([ROOT / f"{lat}_any_shell{s}.data"
                        for s in (1, 2, 3, 4)], style)
    items = []
    for sh in (1, 2, 3, 4):
        src = ROOT / f"{lat}_any_shell{sh}.data"
        if not src.exists():
            continue
        png = out / f"showcase_{lat}_shell{sh}.png"
        render(src, png, STYLES[style], size, span)
        note = ("they share a crosslink" if sh == 1
                else f"{GAPS[lat][sh]} lattice units apart")
        items.append((png, f"shell {sh}", note))
    n = {"SC": 106, "BCC": 512, "FCC": 896, "MIX": 354}[lat]
    return items, (f"{lat}, {n} chains: the same four shells, asked for and "
                   f"built the same way as on the simple cubic case.")


def sheet_networks(style, out, size):
    """The lattices themselves, so the shell spacings have a picture."""
    from render_catalogue import STYLES as _S
    wide = dict(_S["D2_transparent"])
    wide.update(rest_r=0.16, pair_r=0.16, junc_r=0.42, bond=0.14,
                alpha=0.25, anchors=False, keep="all", zoom=7.0)
    items = []
    for lat in ("SC", "BCC", "FCC", "MIX"):
        src = ROOT / f"{lat}_any_shell2.data"
        if not src.exists():
            continue
        png = out / f"showcase_net_{lat}.png"
        render(src, png, wide, size, frame_box=True)
        n = {"SC": 106, "BCC": 512, "FCC": 896, "MIX": 354}[lat]
        f = {"SC": "functionality 4",
             "BCC": "functionality 8",
             "FCC": "functionality 8",
             "MIX": "mixed SC/BCC/FCC, functionality 2 to 8"}[lat]
        items.append((png, lat, f"{n} chains, {f}"))
    return items, ("The four networks. Shell spacing follows from how "
                   "densely the lattice packs its strands: SC 1.00, BCC "
                   "0.71, FCC 0.53, MIX 0.45 lattice units to the second "
                   "shell.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--style", default="D6", choices=list(STYLES))
    ap.add_argument("--size", type=int, default=1000)
    ap.add_argument("--only", default=None)
    args = ap.parse_args()

    out = ROOT / "showcase"
    out.mkdir(parents=True, exist_ok=True)

    builders = {"shells": sheet_shells, "lattices": sheet_lattices,
                "counts": sheet_counts, "orient": sheet_orient,
                "networks": sheet_networks}
    for lat in ("SC", "BCC", "FCC", "MIX"):
        builders[f"shells_{lat}"] = (
            lambda st, o, sz, _l=lat: sheet_lattice_shells(st, o, sz, _l))
    if args.only:
        builders = {args.only: builders[args.only]}

    for name, fn in builders.items():
        print(f"  {name} ...", flush=True)
        items, caption = fn(args.style, out, args.size)
        if not items:
            print(f"    nothing to render for {name}")
            continue
        final = out / f"showcase_{name}_{args.style}.png"
        sheet(items, final, caption, cols=len(items))
        print(f"    {final}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
