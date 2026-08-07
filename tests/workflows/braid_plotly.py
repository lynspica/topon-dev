"""Interactive 3-D view of a CG network's chain paths, entangled and linear.

Run directly::

    python tests/workflows/braid_plotly.py                       # 4x4x4 MIX
    python tests/workflows/braid_plotly.py --dims 5x5x5
    python tests/workflows/braid_plotly.py --only-entangled      # declutter

Writes a self-contained HTML file. Open it in a browser and drag to rotate:
whether two strands hook or merely pass near each other is a question about
depth, and no fixed projection answers it. Rotating does.

What is drawn:

* every chain is its real bead path between two crosslinks, built by the same
  code the writer uses -- not a schematic;
* entangled pairs are drawn in colour, one hue per braid, so the two partners
  of a contact share a colour and can be followed around each other;
* unentangled chains are drawn thin and grey, as the linear reference;
* crosslinks are marked, since a braid near a junction means something
  different from one at mid-chain.

Traces are grouped in the legend, so clicking hides all the grey chains at
once and leaves the entanglements alone in the box.
"""
from __future__ import annotations

import argparse
import random
import sys
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.entanglement import (  # noqa: E402
    BraidShape,
    ContactRequest,
    allocate_contacts,
    closest_approach,
    compose_chain_path,
    far_closed_linking,
    gap_at,
    min_separation,
)
from topon.topology.generator_python import PythonTopologyGenerator  # noqa: E402

OUT = ROOT / "tests/output/braid_plotly"

# Distinct hues for braids. Deliberately not a continuous colormap: adjacent
# braids must be told apart, which ranks separation above ordering.
PALETTE = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4", "#42d4f4",
    "#f032e6", "#bfef45", "#fabed4", "#469990", "#dcbeff", "#9a6324",
    "#800000", "#aaffc3", "#808000", "#000075", "#a9a9a9", "#ffe119",
]


class _Cfg:
    def __init__(self, lattice, dims, max_func, degree_dist, mix, cutoff):
        self.lattice_type = lattice
        self.lattice_size = dims
        self.max_functionality = max_func
        self.degree_distribution = degree_dist
        self.periodicity = "111"
        self.mix_fractions = mix
        self.mix_cutoff = cutoff


def build(lattice, dims, max_func, degree_dist, mix, cutoff, seed):
    random.seed(seed)
    gen = PythonTopologyGenerator(
        _Cfg(lattice, dims, max_func, degree_dist, mix, cutoff))
    with redirect_stdout(StringIO()):
        graphs = gen.generate(trials=6000, max_saves=1, time_limit=180)
    if not graphs:
        raise SystemExit(f"no {lattice} graph for {dims}")
    return graphs[0]


def density_scale(graph, dp, density):
    """Lattice units -> sigma, from bead count and target density.

    Taken from physics rather than chosen, because it sets the junction
    spacing and therefore whether any pair of strands is close enough to
    braid at all.
    """
    n_beads = graph.number_of_edges() * dp + graph.number_of_nodes()
    cells = float(np.prod(np.asarray(graph.graph["box"], float)))
    return ((n_beads / density) / cells) ** (1.0 / 3.0)


def chain_chords(graph, scale):
    """One chord per chain, unwrapped across the periodic boundary.

    A chain whose junctions sit on opposite faces is physically short: it
    crosses the boundary. Built from the raw positions it would be drawn as
    a line straight across the system, which is both wrong and unreadable.
    """
    pos = {n: np.asarray(d["pos"], float) * scale
           for n, d in graph.nodes(data=True)}
    L = np.asarray(graph.graph["box"], float) * scale

    chords, ends, wraps = {}, {}, set()
    for k, (u, v) in enumerate(sorted(graph.edges())):
        a = pos[u]
        raw = pos[v] - a
        mic = raw - L * np.round(raw / L)
        if not np.allclose(mic, raw):
            wraps.add(k)
        chords[k] = (a, a + mic)
        ends[k] = (u, v)
    return chords, pos, ends, wraps


def candidates(chords, max_gap, limit):
    ids = sorted(chords)
    found = []
    for i, ka in enumerate(ids):
        a0, a1 = chords[ka]
        mid_a = 0.5 * (a0 + a1)
        for kb in ids[i + 1:]:
            b0, b1 = chords[kb]
            if np.linalg.norm(mid_a - 0.5 * (b0 + b1)) > 3.0 * max_gap:
                continue
            s, _ = closest_approach(a0, a1, b0, b1)
            gap, _ = gap_at(a0, a1, b0, b1, s)
            if 1e-6 < gap <= max_gap:
                found.append((gap, ka, kb))
    found.sort()
    return [ContactRequest(ka, kb, windings=1, priority=-g)
            for g, ka, kb in found[:limit]]


def figure(paths, alloc, junctions, box, wraps, only_entangled, title):
    import plotly.graph_objects as go

    braid_of = {}
    for i, a in enumerate(alloc.accepted):
        for chain in (a.request.chain_a, a.request.chain_b):
            braid_of.setdefault(chain, []).append(i)

    fig = go.Figure()

    if not only_entangled:
        # One trace for all plain chains: hundreds of separate traces make
        # the legend useless and the file large.
        xs, ys, zs = [], [], []
        for k, p in paths.items():
            if k in braid_of:
                continue
            xs += list(p[:, 0]) + [None]
            ys += list(p[:, 1]) + [None]
            zs += list(p[:, 2]) + [None]
        if xs:
            fig.add_trace(go.Scatter3d(
                x=xs, y=ys, z=zs, mode="lines", name="linear chains",
                line=dict(color="rgba(150,150,150,0.45)", width=1.5),
                hoverinfo="skip"))

    for i, a in enumerate(alloc.accepted):
        col = PALETTE[i % len(PALETTE)]
        r = a.request
        for chain in (r.chain_a, r.chain_b):
            p = paths[chain]
            fig.add_trace(go.Scatter3d(
                x=p[:, 0], y=p[:, 1], z=p[:, 2], mode="lines",
                name=f"braid {r.chain_a}-{r.chain_b}",
                legendgroup=f"b{i}",
                showlegend=(chain == r.chain_a),
                line=dict(color=col, width=5),
                hovertemplate=(f"chain {chain}<br>braid {r.chain_a}-{r.chain_b}"
                               f"<br>e={a.windings}<extra></extra>")))
        # Mark the contact so a braid can be found without hunting.
        o = a.contact.origin
        fig.add_trace(go.Scatter3d(
            x=[o[0]], y=[o[1]], z=[o[2]], mode="markers",
            marker=dict(size=4, color=col, symbol="diamond"),
            legendgroup=f"b{i}", showlegend=False,
            hovertemplate=(f"contact {r.chain_a}-{r.chain_b}"
                           f"<br>gap {a.contact.gap:.2f}"
                           f"<br>e={a.windings}<extra></extra>")))

    jx = np.array([p for p in junctions.values()])
    fig.add_trace(go.Scatter3d(
        x=jx[:, 0], y=jx[:, 1], z=jx[:, 2], mode="markers",
        name="crosslinks", marker=dict(size=2.5, color="#222", opacity=0.55),
        hoverinfo="skip"))

    L = np.asarray(box, float)
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis=dict(title="x / sigma", range=[0, L[0]]),
            yaxis=dict(title="y / sigma", range=[0, L[1]]),
            zaxis=dict(title="z / sigma", range=[0, L[2]]),
            aspectmode="cube"),
        legend=dict(itemsizing="constant", font=dict(size=10)),
        margin=dict(l=0, r=0, t=44, b=0),
        template="plotly_white")
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lattice", default="MIX")
    ap.add_argument("--dims", default="4x4x4")
    ap.add_argument("--mix", default="0.2,0.4,0.4")
    ap.add_argument("--mix-cutoff", type=float, default=1.0)
    ap.add_argument("--max-func", type=int, default=12)
    ap.add_argument("--degree-dist", default="0:0,1:0")
    ap.add_argument("--dp", type=int, default=200)
    ap.add_argument("--density", type=float, default=0.30)
    ap.add_argument("--gap-factor", type=float, default=0.6)
    ap.add_argument("--max-requests", type=int, default=400)
    ap.add_argument("--beads", type=int, default=0,
                    help="beads drawn per chain; 0 uses dp")
    ap.add_argument("--only-entangled", action="store_true")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    dims = tuple(int(v) for v in args.dims.lower().split("x"))
    mix = dict(zip(("SC", "BCC", "FCC"),
                   (float(v) for v in args.mix.split(","))))

    graph = build(args.lattice, dims, args.max_func, args.degree_dist,
                  mix, args.mix_cutoff, args.seed)
    scale = density_scale(graph, args.dp, args.density)
    chords, junctions, ends, wraps = chain_chords(graph, scale)

    reqs = candidates(chords, args.gap_factor * scale, args.max_requests)
    alloc = allocate_contacts(reqs, chords, BraidShape())

    n_draw = args.beads or (args.dp + 2)
    paths = {k: compose_chain_path(k, alloc, chords, n_draw) for k in chords}

    exact = wrong = 0
    clear = []
    for a in alloc.accepted:
        r = a.request
        lk = far_closed_linking(paths[r.chain_a], paths[r.chain_b], a.contact)
        clear.append(min_separation(paths[r.chain_a], paths[r.chain_b]))
        # Magnitude only: handedness follows the chord orientations, so both
        # signs occur on a real lattice and only the count is prescribed.
        if round(abs(lk)) == a.windings:
            exact += 1
        else:
            wrong += 1

    multi = sum(1 for v in alloc.partners.values() if len(v) > 1)
    tag = f"{args.lattice}_{dims[0]}x{dims[1]}x{dims[2]}"
    title = (f"{tag} CG network - {graph.number_of_nodes()} crosslinks, "
             f"{len(chords)} chains, {len(alloc.accepted)} braids "
             f"({exact} exact), {multi} chains with 2+ partners")

    print(f"{tag}: {graph.number_of_nodes()} junctions, {len(chords)} chains "
          f"({len(wraps)} wrap), spacing {scale:.1f} sigma")
    print(f"  {len(reqs)} candidates -> {len(alloc.accepted)} braids, "
          f"{exact} exact, {wrong} wrong, {multi} multi-partner chains")
    if clear:
        print(f"  clearance {min(clear):.2f} to {max(clear):.2f} sigma")
    reasons = {}
    for rj in alloc.rejected:
        reasons[rj.reason] = reasons.get(rj.reason, 0) + 1
    for why, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"    {n:4d} refused: {why}")

    fig = figure(paths, alloc, junctions,
                 np.asarray(graph.graph["box"], float) * scale,
                 wraps, args.only_entangled, title)

    OUT.mkdir(parents=True, exist_ok=True)
    suffix = "_entangled_only" if args.only_entangled else ""
    path = OUT / f"{tag}{suffix}.html"
    fig.write_html(str(path), include_plotlyjs="inline")
    print(f"  wrote {path.name}  ({path.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
