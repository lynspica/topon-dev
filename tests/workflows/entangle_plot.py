"""Interactive view of one entanglement, at the scale you can actually read.

    python tests/workflows/entangle_plot.py                  # 1 to 4 windings
    python tests/workflows/entangle_plot.py --windings 3     # just one
    python tests/workflows/entangle_plot.py --context        # add neighbours

Two chains, drawn as their real bead paths, filling the frame. Earlier views
put 120 braids in one box and nothing could be followed; this draws the pair
and leaves the rest out unless asked for.

Each requested winding count gets its own panel, so what changes with the
count is visible side by side: the braid gets longer, and the clearance does
not, which is the property the allocator now guarantees.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    ContactRequest,
    build_network,
    build_with_braids,
    geometry,
    min_separation,
    separation_bands,
)

OUT = ROOT / "tests/output/entangle_steps"

A_COLOUR = "#e6194b"
B_COLOUR = "#4363d8"


def braid_mask(path, contact, half_span):
    """Which beads lie inside the braid, in the contact's own frame."""
    u = (path - contact.origin) @ contact.axis
    return np.abs(u) <= half_span


def add_pair(fig, row, col, paths, ka, kb, alloc, context=None):
    import plotly.graph_objects as go

    a = alloc.accepted[0]
    for chain, colour, name in ((ka, A_COLOUR, "chain A"),
                                (kb, B_COLOUR, "chain B")):
        p = paths[chain]
        fig.add_trace(go.Scatter3d(
            x=p[:, 0], y=p[:, 1], z=p[:, 2], mode="lines+markers",
            line=dict(color=colour, width=6),
            marker=dict(size=2.5, color=colour),
            name=f"{name} ({chain})", legendgroup=name,
            showlegend=(row == 1 and col == 1),
            hovertemplate=f"{name} bead %{{pointNumber}}<extra></extra>"),
            row=row, col=col)

        # The braid itself, thicker, so the entangled stretch is obvious
        # against the straight run that leads into it.
        m = braid_mask(p, a.contact, a.half_span)
        if m.any():
            fig.add_trace(go.Scatter3d(
                x=p[m, 0], y=p[m, 1], z=p[m, 2], mode="lines",
                line=dict(color=colour, width=14), opacity=0.55,
                showlegend=False, hoverinfo="skip"), row=row, col=col)

    for chain in (ka, kb):
        p = paths[chain]
        fig.add_trace(go.Scatter3d(
            x=p[[0, -1], 0], y=p[[0, -1], 1], z=p[[0, -1], 2],
            mode="markers", marker=dict(size=6, color="#222", symbol="square"),
            name="crosslinks", legendgroup="crosslinks",
            showlegend=(row == 1 and col == 1 and chain == ka),
            hoverinfo="skip"), row=row, col=col)

    if context:
        xs, ys, zs = [], [], []
        for p in context:
            xs += list(p[:, 0]) + [None]
            ys += list(p[:, 1]) + [None]
            zs += list(p[:, 2]) + [None]
        fig.add_trace(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="lines",
            line=dict(color="rgba(150,150,150,0.35)", width=2),
            name="other chains", legendgroup="other",
            showlegend=(row == 1 and col == 1), hoverinfo="skip"),
            row=row, col=col)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--windings", type=int, nargs="*", default=[1, 2, 3, 4])
    ap.add_argument("--shell", type=int, default=1)
    ap.add_argument("--bond", type=float, default=BOND)
    ap.add_argument("--context", action="store_true",
                    help="draw the chains that pass nearby, in grey")
    ap.add_argument("--clearance", type=float, default=1.0)
    args = ap.parse_args()

    from plotly.subplots import make_subplots

    graph = build_network()
    probe = geometry(graph, bond=args.bond)
    band = separation_bands(probe)[min(args.shell, 6) - 1]

    panels, pair = [], None
    for want in args.windings:
        for gap_u, ka, kb in ([pair] if pair else band):
            geo, alloc, paths = build_with_braids(
                graph, [ContactRequest(ka, kb, windings=want)], args.bond)
            if alloc.accepted:
                pair = (gap_u, ka, kb)
                break
        else:
            print(f"  e={want}: no pair could be built")
            continue
        a = alloc.accepted[0]
        sep = min_separation(paths[ka], paths[kb])
        ctx = None
        if args.context:
            mid = 0.5 * (paths[ka][0] + paths[ka][-1])
            near = sorted(
                (k for k in paths if k not in (ka, kb)),
                key=lambda k: np.linalg.norm(paths[k].mean(axis=0) - mid))[:14]
            ctx = [paths[k] for k in near]
        panels.append((want, a, paths, ka, kb, sep, ctx))
        print(f"  e={want}: granted {a.windings}, clearance {sep:.2f} sigma, "
              f"pitch {a.shape.pitch:.2f}, span {2 * a.half_span:.1f}")

    if not panels:
        raise SystemExit("nothing to draw")

    n = len(panels)
    fig = make_subplots(
        rows=1, cols=n, specs=[[{"type": "scene"}] * n],
        subplot_titles=[f"asked {w}, granted {a.windings}"
                        f"<br><sub>clearance {s:.2f} sigma</sub>"
                        for w, a, _, _, _, s, _ in panels])

    for i, (want, a, paths, ka, kb, sep, ctx) in enumerate(panels, start=1):
        add_pair(fig, 1, i, paths, ka, kb,
                 type("A", (), {"accepted": [a]})(), ctx)
        # Frame each panel on the braid, not the whole box, or the feature
        # we are looking at is a few pixels in the middle of a large cube.
        both = np.vstack([paths[ka], paths[kb]])
        c = both.mean(axis=0)
        r = 0.6 * float(np.abs(both - c).max())
        scene = dict(
            xaxis=dict(range=[c[0] - r, c[0] + r], title=""),
            yaxis=dict(range=[c[1] - r, c[1] + r], title=""),
            zaxis=dict(range=[c[2] - r, c[2] + r], title=""),
            aspectmode="cube")
        fig.update_layout(**{f"scene{'' if i == 1 else i}": scene})

    fig.update_layout(
        title=("One entanglement between two chains, at 1 to 4 windings. "
               "Drag any panel to rotate."),
        template="plotly_white", margin=dict(l=0, r=0, t=90, b=0),
        legend=dict(orientation="h", y=-0.02))

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "entanglement_curves.html"
    fig.write_html(str(path), include_plotlyjs="inline")
    print(f"  wrote {path.name}  ({path.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
