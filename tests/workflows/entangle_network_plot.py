"""The entangled pair shown inside the real mixed-lattice network.

    python tests/workflows/entangle_network_plot.py
    python tests/workflows/entangle_network_plot.py --sites 3

Left panel: the whole network, every chain, with the entangled pair picked
out in colour. Right panel: the same pair on its own, framed on the braid.

The pair is chosen for the length of its *close run*, not for its minimum
gap. Several entanglements between one pair need the two chains to stay
alongside each other over a long stretch; a pair that merely touches at one
point has room for exactly one site, however close that point is.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.entanglement.braid import (  # noqa: E402
    closest_approach,
    feasible_window,
    gap_at,
)
from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    ContactRequest,
    build_network,
    build_with_braids,
    geometry,
    min_separation,
)

OUT = ROOT / "tests/output/entangle_steps"
A_COLOUR, B_COLOUR = "#e6194b", "#4363d8"


def rank_pairs(geo, want_sites, span_guess=7.8):
    """Pairs ranked by how many sites their close run could hold.

    Two quantities matter and they pull against each other. The run is how
    much chord the pair spends near each other, and sets how many sites fit.
    The gap is how far apart they are there, and sets whether a braid can
    reach across at all. Ranking on the run alone finds chains that are
    parallel and far apart; ranking on the gap alone finds chains that touch
    at a point and diverge. This wants both.
    """
    ch, s = geo["chords"], geo["scale"]
    ids = sorted(ch)
    out = []
    for i, ka in enumerate(ids):
        a0, a1 = ch[ka]
        mid_a, L = 0.5 * (a0 + a1), np.linalg.norm(a1 - a0)
        for kb in ids[i + 1:]:
            b0, b1 = ch[kb]
            if np.linalg.norm(mid_a - 0.5 * (b0 + b1)) > 1.6 * s:
                continue
            t, _ = closest_approach(a0, a1, b0, b1)
            gap, _ = gap_at(a0, a1, b0, b1, t)
            if not (1e-6 < gap < 0.35 * s):
                continue
            lo, hi = feasible_window(a0, a1, b0, b1)
            run = (hi - lo) * L
            # Sites the run could hold, then closeness as the tie-break.
            out.append((min(want_sites, int(run // span_guess)), -gap,
                        run, gap, ka, kb))
    out.sort(reverse=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sites", type=int, default=3)
    ap.add_argument("--bond", type=float, default=BOND)
    ap.add_argument("--tries", type=int, default=25)
    ap.add_argument("--tolerance", type=float, default=1.0,
                    help="how far the pair may drift apart and still "
                         "host a site, as a fraction of their closest gap")
    ap.add_argument("--no-verify", action="store_true",
                    help="skip the linking-number check on each "
                         "placement; use when Z1+ is the arbiter")
    args = ap.parse_args()

    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    graph = build_network()
    probe = geometry(graph, bond=args.bond)
    ranked = rank_pairs(probe, args.sites)
    if not ranked:
        raise SystemExit("no candidate pairs")
    print(f"  {len(ranked)} candidate pairs; best close runs:")
    for fit, _, run, gap, ka, kb in ranked[:5]:
        print(f"    {ka}-{kb}: run {run:.1f} sigma, gap {gap:.1f}, "
              f"room for {fit} site(s)")

    best = None
    for fit, _, run, gap, ka, kb in ranked[:args.tries]:
        reqs = [ContactRequest(ka, kb, windings=1, priority=-i)
                for i in range(args.sites)]
        geo, alloc, paths = build_with_braids(graph, reqs, args.bond,
                                              verify=not args.no_verify,
                                              tolerance=args.tolerance)
        n = len(alloc.accepted)
        if best is None or n > best[0]:
            best = (n, geo, alloc, paths, ka, kb, run, gap)
        if n >= args.sites:
            break

    n, geo, alloc, paths, ka, kb, run, gap = best
    sep = min_separation(paths[ka], paths[kb])
    c0, c1 = geo["chords"][ka]
    chord = c1 - c0
    at = sorted(float((x.contact.origin - c0) @ chord / (chord @ chord))
                for x in alloc.accepted)
    print()
    print(f"  chosen {ka}-{kb}: {n} of {args.sites} sites granted, "
          f"at {', '.join(f'{v:.2f}' for v in at)} along the chain")
    print(f"  close run {run:.1f} sigma, gap {gap:.1f}, clearance {sep:.2f}")
    reasons = {}
    for r in alloc.rejected:
        reasons[r.reason] = reasons.get(r.reason, 0) + 1
    for why, k in reasons.items():
        print(f"    {k} refused: {why}")

    fig = make_subplots(
        rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
        subplot_titles=(f"the whole network ({len(paths)} chains)",
                        f"the pair, {n} site(s), clearance {sep:.2f} sigma"))

    # Every other chain, as one trace: hundreds of traces makes the legend
    # useless and the file large.
    xs, ys, zs = [], [], []
    for k, p in paths.items():
        if k in (ka, kb):
            continue
        xs += list(p[:, 0]) + [None]
        ys += list(p[:, 1]) + [None]
        zs += list(p[:, 2]) + [None]
    fig.add_trace(go.Scatter3d(
        x=xs, y=ys, z=zs, mode="lines", name="network",
        line=dict(color="rgba(140,140,140,0.30)", width=1.5),
        hoverinfo="skip"), row=1, col=1)

    jx = np.array(list(geo["pos"].values()))
    fig.add_trace(go.Scatter3d(
        x=jx[:, 0], y=jx[:, 1], z=jx[:, 2], mode="markers", name="crosslinks",
        marker=dict(size=2, color="#222", opacity=0.5), hoverinfo="skip"),
        row=1, col=1)

    for col in (1, 2):
        for chain, colour, label in ((ka, A_COLOUR, "chain A"),
                                     (kb, B_COLOUR, "chain B")):
            p = paths[chain]
            fig.add_trace(go.Scatter3d(
                x=p[:, 0], y=p[:, 1], z=p[:, 2], mode="lines",
                line=dict(color=colour, width=7 if col == 2 else 5),
                name=f"{label} ({chain})", legendgroup=label,
                showlegend=(col == 1),
                hovertemplate=f"{label} bead %{{pointNumber}}<extra></extra>"),
                row=1, col=col)
            if col == 2:
                fig.add_trace(go.Scatter3d(
                    x=p[:, 0], y=p[:, 1], z=p[:, 2], mode="markers",
                    marker=dict(size=2.5, color=colour), showlegend=False,
                    hoverinfo="skip"), row=1, col=2)

        for i, x in enumerate(alloc.accepted):
            o = x.contact.origin
            fig.add_trace(go.Scatter3d(
                x=[o[0]], y=[o[1]], z=[o[2]], mode="markers",
                marker=dict(size=7, color="#f58231", symbol="diamond"),
                name="entanglement site", legendgroup="site",
                showlegend=(col == 1 and i == 0),
                hovertemplate=(f"site {i + 1} of {n}<br>"
                               f"at {at[i]:.2f} along the chain<br>"
                               f"gap {x.contact.gap:.2f}<extra></extra>")),
                row=1, col=col)

    L = geo["L"]
    both = np.vstack([paths[ka], paths[kb]])
    c = both.mean(axis=0)
    r = 0.62 * float(np.abs(both - c).max())
    fig.update_layout(
        scene=dict(xaxis=dict(range=[0, L[0]], title="x / sigma"),
                   yaxis=dict(range=[0, L[1]], title="y / sigma"),
                   zaxis=dict(range=[0, L[2]], title="z / sigma"),
                   aspectmode="cube"),
        scene2=dict(xaxis=dict(range=[c[0] - r, c[0] + r], title=""),
                    yaxis=dict(range=[c[1] - r, c[1] + r], title=""),
                    zaxis=dict(range=[c[2] - r, c[2] + r], title=""),
                    aspectmode="cube"),
        title=(f"MIX 4x4x4, {graph.number_of_nodes()} crosslinks, "
               f"{len(paths)} chains. Drag either panel to rotate."),
        template="plotly_white", margin=dict(l=0, r=0, t=80, b=0),
        legend=dict(orientation="h", y=-0.02))

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "entanglement_in_network.html"
    fig.write_html(str(path), include_plotlyjs="inline")
    print(f"  wrote {path.name}  ({path.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
