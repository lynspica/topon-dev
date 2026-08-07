"""Entanglement paths on the real mixed lattice, across pairs and regimes.

    python tests/workflows/waypoint_gallery.py
    python tests/workflows/waypoint_gallery.py --bands 1 2 3 --pairs 2

Every panel is a real chain pair from the MIX 4x4x4 network, drawn through
waypoints placed by hand along the chain. Rows are neighbour regimes, how
far apart the two strands are; columns are what was asked for at each site.

Nothing here refuses a request. The paths are built from the points given,
measured, and reported: linking number against the same pair drawn straight,
and closest approach.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.entanglement.braid import far_closed_linking  # noqa: E402
from topon.conformation.entanglement.waypoints import (  # noqa: E402
    Site,
    entangled_pair,
)
from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    DP,
    build_network,
    geometry,
    separation_bands,
)

OUT = ROOT / "tests/output/entangle_steps"
A_COLOUR, B_COLOUR = "#e6194b", "#4363d8"

# What to ask for in each column. The point of the sweep is that these are
# choices, not the output of a search.
RECIPES = [
    ("one at the middle", [Site(0.5, 1)]),
    ("two turns, one site", [Site(0.5, 2)]),
    ("three sites", [Site(0.25, 1), Site(0.5, 1), Site(0.75, 1)]),
    ("four sites", [Site(0.2, 1), Site(0.4, 1), Site(0.6, 1), Site(0.8, 1)]),
]


def straight(a0, a1, n):
    t = np.linspace(0.0, 1.0, n)[:, None]
    return a0 + t * (a1 - a0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bands", type=int, nargs="*", default=[1, 3, 5])
    ap.add_argument("--pairs", type=int, default=1,
                    help="pairs to draw from each band")
    ap.add_argument("--reach", type=float, default=0.45)
    ap.add_argument("--beads", type=int, default=DP + 2)
    args = ap.parse_args()

    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    graph = build_network()
    geo = geometry(graph, bond=BOND)
    bands = separation_bands(geo)
    chords = geo["chords"]

    picks = []
    for b in args.bands:
        if b > len(bands):
            continue
        for gap_u, ka, kb in bands[b - 1][:args.pairs]:
            picks.append((b, gap_u, ka, kb))
    if not picks:
        raise SystemExit("no pairs in the requested bands")

    rows, cols = len(picks), len(RECIPES)
    titles = []
    results = []
    for b, gap_u, ka, kb in picks:
        a0, a1 = chords[ka]
        b0, b1 = chords[kb]
        base = far_closed_linking(straight(a0, a1, 60), straight(b0, b1, 60))
        for label, spec in RECIPES:
            pa, pb, info = entangled_pair(a0, a1, b0, b1, spec,
                                          n_beads=args.beads,
                                          reach=args.reach)
            lk = far_closed_linking(pa, pb) - base
            sep = float(np.linalg.norm(pa[:, None, :] - pb[None, :, :],
                                       axis=-1).min())
            bond = float(np.linalg.norm(np.diff(pa, axis=0), axis=1).max())
            asked = sum(s.turns for s in spec)
            results.append((b, ka, kb, label, spec, pa, pb, info,
                            asked, lk, sep, bond))
            titles.append(f"band {b}, {ka}-{kb} &middot; {label}<br>"
                          f"<sub>asked {asked}, got {lk:+.1f}, "
                          f"clearance {sep:.1f} sigma</sub>")

    fig = make_subplots(rows=rows, cols=cols,
                        specs=[[{"type": "scene"}] * cols] * rows,
                        subplot_titles=titles,
                        horizontal_spacing=0.01, vertical_spacing=0.06)

    print(f"  {'band':>5} {'pair':>10} {'recipe':>20} {'asked':>6} "
          f"{'got':>6} {'clearance':>10} {'max bond':>9}")
    for i, (b, ka, kb, label, spec, pa, pb, info,
            asked, lk, sep, bond) in enumerate(results):
        r, c = i // cols + 1, i % cols + 1
        first = (i == 0)
        for p, colour, name in ((pa, A_COLOUR, "chain A"),
                                (pb, B_COLOUR, "chain B")):
            fig.add_trace(go.Scatter3d(
                x=p[:, 0], y=p[:, 1], z=p[:, 2], mode="lines",
                line=dict(color=colour, width=5), name=name,
                legendgroup=name, showlegend=first, hoverinfo="skip"),
                row=r, col=c)
        # Where each site was asked for.
        m = np.array([d["mid"] for d in info])
        fig.add_trace(go.Scatter3d(
            x=m[:, 0], y=m[:, 1], z=m[:, 2], mode="markers",
            marker=dict(size=5, color="#f58231", symbol="diamond"),
            name="site", legendgroup="site", showlegend=first,
            hovertext=[f"at {d['at']:.2f}, {d['turns']} turn(s)"
                       for d in info], hoverinfo="text"), row=r, col=c)
        # Junctions.
        for p in (pa, pb):
            fig.add_trace(go.Scatter3d(
                x=p[[0, -1], 0], y=p[[0, -1], 1], z=p[[0, -1], 2],
                mode="markers",
                marker=dict(size=4, color="#222", symbol="square"),
                name="crosslink", legendgroup="crosslink",
                showlegend=(first and p is pa), hoverinfo="skip"),
                row=r, col=c)

        both = np.vstack([pa, pb])
        cen = both.mean(axis=0)
        rad = 0.62 * float(np.abs(both - cen).max())
        key = "scene" if i == 0 else f"scene{i + 1}"
        fig.update_layout(**{key: dict(
            xaxis=dict(range=[cen[0] - rad, cen[0] + rad], title="",
                       showticklabels=False),
            yaxis=dict(range=[cen[1] - rad, cen[1] + rad], title="",
                       showticklabels=False),
            zaxis=dict(range=[cen[2] - rad, cen[2] + rad], title="",
                       showticklabels=False),
            aspectmode="cube")})
        print(f"  {b:5d} {f'{ka}-{kb}':>10} {label:>20} {asked:6d} "
              f"{lk:6.1f} {sep:10.2f} {bond:9.2f}")

    fig.update_layout(
        height=420 * rows, width=380 * cols,
        title=("Entanglement paths through chosen waypoints, MIX 4x4x4. "
               "Rows are neighbour regimes, columns are what was asked for. "
               "Drag any panel to rotate."),
        template="plotly_white", margin=dict(l=0, r=0, t=110, b=0),
        legend=dict(orientation="h", y=-0.03))

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "waypoint_gallery.html"
    fig.write_html(str(path), include_plotlyjs="inline")
    print(f"\n  wrote {path.name}  ({path.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
