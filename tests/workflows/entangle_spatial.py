"""Spatially heterogeneous entanglement, built with the kink, run through MD.

    python tests/workflows/entangle_spatial.py --bias clusters
    python tests/workflows/entangle_spatial.py --bias region --stages 2
    python tests/workflows/entangle_spatial.py --bias uniform      # control

Uses the pipeline's own entanglement path: `select_entanglements` draws the
pairs with a spatial bias, and `calculate_entangled_kink` builds the kinked
backbone. Chains that were not selected are drawn straight, so the picture
shows the network and the entanglements in it rather than the entanglements
alone.

Spatial control needs headroom. At 0.1 entanglements per chain a region bias
gives about eightfold enrichment in the target volume; by 1.2 per chain the
draw has taken nearly every candidate and there is nothing left to bias.
Keep the density well below the size of the candidate pool.

Writes an interactive Plotly view and runs the LAMMPS stages.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.assignment.entanglements import (  # noqa: E402
    find_crossing_candidates,
    select_entanglements,
)
from topon.config.schema import EntanglementsConfig  # noqa: E402
from topon.utils.network_helpers import calculate_entangled_kink  # noqa: E402
from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    DP,
    LATTICE,
    OUT,
    build_network,
    conform_and_script,
    geometry,
    report_bonds,
    resample,
    run_md,
    write_system,
)

BIAS = {
    "uniform": {},
    "region": {"center": [0.5, 0.5, 0.5], "radius": 0.28, "strength": 50.0},
    "anti_region": {"center": [0.5, 0.5, 0.5], "radius": 0.28,
                    "strength": 50.0},
    "gradient": {"axis": "x", "strength": 4.0},
    "clusters": {"centers": [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]],
                 "sigma": 0.13, "strength": 60.0},
}


def select(graph, kind, per_chain, seed):
    """Draw entangled pairs with a spatial bias, on a MultiGraph copy."""
    import random

    G = nx.MultiGraph()
    G.add_nodes_from(graph.nodes(data=True))
    for u, v in graph.edges():
        G.add_edge(u, v)
    dims = np.asarray(graph.graph["box"], float)

    cands = find_crossing_candidates(G, dims)
    random.seed(seed)
    cfg = EntanglementsConfig(enabled=True, avg_crosslinks_per_chain=per_chain,
                              placement_bias_kind=kind,
                              placement_bias_params=BIAS[kind])
    sel = select_entanglements(G, cfg, dims, candidates=list(cands),
                               num_chains=G.number_of_edges())
    return G, dims, sel, len(cands)


def kinked_paths(graph, geo, sel, dims, dp, params):
    """Kinked backbones for the entangled chains, straight ones for the rest.

    The kink is aimed at the partner's midpoint under the minimum image, the
    same way `generate_cg_combined.py` does it, so a pair whose crosslinks
    sit on opposite faces still bulges toward each other rather than across
    the box.
    """
    ch, ends = geo["chords"], geo["ends"]
    scale = geo["scale"]
    key_of = {frozenset(v): k for k, v in ends.items()}

    partner = {}
    for e1, e2, count in sel:
        k1 = key_of.get(frozenset((e1[0], e1[1])))
        k2 = key_of.get(frozenset((e2[0], e2[1])))
        if k1 is None or k2 is None:
            continue
        partner[k1] = (k2, count)
        partner[k2] = (k1, count)

    pos = {n: np.asarray(d["pos"], float) for n, d in graph.nodes(data=True)}
    paths, sites = {}, []
    for k, (c0, c1) in ch.items():
        if k not in partner:
            paths[k] = resample(np.stack([c0, c1]), dp + 2)
            continue
        other, count = partner[k]
        u, v = ends[k]
        pu, pv = pos[u], pos[v]
        mic = (pv - pu) - dims * np.round((pv - pu) / dims)
        mid = pu + 0.5 * mic

        ou, ov = ends[other]
        qu, qv = pos[ou], pos[ov]
        pmic = (qv - qu) - dims * np.round((qv - qu) / dims)
        pmid = qu + 0.5 * pmic
        d = pmid - mid
        d -= dims * np.round(d / dims)
        orient = d if np.linalg.norm(d) > 1e-6 else None

        raw = calculate_entangled_kink(pu, pu + mic, dp + 2, params,
                                       orientation_vec=orient,
                                       num_entanglements=float(count))
        p = np.array([raw[i] for i in sorted(raw)], float) * scale
        paths[k] = p
        sites.append((k, other, (mid + 0.5 * d) * scale, count))
    return paths, partner, sites


def plot(name, geo, paths, partner, sites, path, dims):
    import plotly.graph_objects as go

    fig = go.Figure()
    xs, ys, zs = [], [], []
    for k, p in paths.items():
        if k in partner:
            continue
        xs += list(p[:, 0]) + [None]
        ys += list(p[:, 1]) + [None]
        zs += list(p[:, 2]) + [None]
    if xs:
        fig.add_trace(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="lines", name="linear chains",
            line=dict(color="rgba(150,150,150,0.28)", width=1),
            hoverinfo="skip"))

    ex, ey, ez = [], [], []
    for k in partner:
        p = paths[k]
        ex += list(p[:, 0]) + [None]
        ey += list(p[:, 1]) + [None]
        ez += list(p[:, 2]) + [None]
    if ex:
        fig.add_trace(go.Scatter3d(
            x=ex, y=ey, z=ez, mode="lines", name="entangled chains",
            line=dict(color="#e6194b", width=3), hoverinfo="skip"))

    if sites:
        S = np.array([s[2] for s in sites])
        cnt = [s[3] for s in sites]
        fig.add_trace(go.Scatter3d(
            x=S[:, 0], y=S[:, 1], z=S[:, 2], mode="markers",
            name="entanglement sites",
            marker=dict(size=5, color=cnt, colorscale="YlOrRd",
                        cmin=1, showscale=True,
                        colorbar=dict(title="e", thickness=12, len=0.5)),
            hovertemplate="e=%{marker.color}<extra></extra>"))

    jx = np.array(list(geo["pos"].values()))
    fig.add_trace(go.Scatter3d(
        x=jx[:, 0], y=jx[:, 1], z=jx[:, 2], mode="markers", name="crosslinks",
        marker=dict(size=2, color="#222", opacity=0.45), hoverinfo="skip"))

    L = geo["L"]
    fig.update_layout(
        title=(f"{name} bias — {len(sites)} entanglement sites on "
               f"{len(paths)} chains. Grey is the linear network, red the "
               f"kinked chains. Drag to rotate."),
        scene=dict(xaxis=dict(range=[0, L[0]], title="x / σ"),
                   yaxis=dict(range=[0, L[1]], title="y / σ"),
                   zaxis=dict(range=[0, L[2]], title="z / σ"),
                   aspectmode="cube"),
        template="plotly_white", margin=dict(l=0, r=0, t=60, b=0))
    fig.write_html(str(path), include_plotlyjs="inline")
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bias", default="clusters", choices=sorted(BIAS))
    ap.add_argument("--per-chain", type=float, default=0.20,
                    help="entanglements per chain; spatial bias needs this "
                         "well below saturation of the candidate pool")
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--density", type=float, default=0.30)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--no-md", action="store_true")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    graph = build_network()
    geo = geometry(graph, dp=args.dp, density=args.density)
    G, dims, sel, n_cand = select(graph, args.bias, args.per_chain, args.seed)
    paths, partner, sites = kinked_paths(
        graph, geo, sel, dims, args.dp,
        {"overshoot": 0.2, "z_amp": 0.5, "sigma": 0.15})

    L = geo["L"]
    frac = np.array([s[2] for s in sites]) / L if sites else np.zeros((0, 3))
    r = np.linalg.norm(frac - 0.5, axis=1) if len(frac) else np.zeros(0)
    bonds = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                            for p in paths.values()])
    print(f"  {graph.number_of_nodes()} crosslinks, "
          f"{graph.number_of_edges()} chains, DP {args.dp}")
    print(f"  box {L[0]:.1f} sigma, density {geo['density']:.3f}")
    print(f"  {n_cand} candidates, {len(sel)} kinks, "
          f"{len(partner)} chains entangled")
    if len(r):
        print(f"  sites: {int((r < 0.28).sum())} in the central quarter, "
              f"mean |r-centre| {r.mean():.3f}")
    print(f"  bonds {bonds.min():.3f} to {bonds.max():.3f}")

    root = OUT / f"spatial_{args.bias}"
    n_atoms, node_atom, chain_atoms = write_system(graph, geo, paths, root)
    sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                             protocol="hardcore")
    if not args.no_md:
        print(f"\n  --- LAMMPS, stages 1 to {args.stages} ---")
        run_md(sim, args.stages)
        print()
        report_bonds(root)

    out = OUT / "spatial"
    out.mkdir(parents=True, exist_ok=True)
    p = plot(args.bias, geo, paths, partner, sites,
             out / f"{args.bias}.html", dims)
    print(f"\n  wrote {p.name} ({p.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
