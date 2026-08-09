"""Designed entanglements on every lattice, end to end, verified.

    python tests/workflows/entangle_all.py                    # build + MD + Z1+
    python tests/workflows/entangle_all.py --no-md            # geometry only
    python tests/workflows/entangle_all.py --lattices SC MIX

Builds each lattice at a coil where chains actually interpenetrate, draws
them as random walks, winds pairs together where their paths already meet,
runs the three-stage protocol, and measures every prescribed pair on its own
with Z1+.

This supersedes the chord-based construction in the step scripts. That one
sites an entanglement between two chords and sends both chains to the
midpoint, which costs contour in proportion to chord separation and so
refuses on uniform lattices. This one rotates the two chains' separation
where they are already adjacent, which costs about 0.4 sigma per winding.

Writes an interactive Plotly view and a markdown summary per run.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.entanglement.braid import far_closed_linking  # noqa: E402
from topon.conformation.entanglement.contact import (  # noqa: E402
    find_contacts,
    wind_all,
)
from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    DP,
    LATTICE,
    OUT,
    build_network,
    chain_ids,
    conform_and_script,
    geometry,
    report_bonds,
    run_md,
    run_z1,
    write_system,
    z1_export,
)

# Contour over chord. Below about 4 the chains are too extended to meet at
# all; from 6 they interpenetrate and contacts are everywhere. See the
# 2026-08-08 journal entry -- the low value used by the step scripts is what
# made every earlier limit look like a property of the lattice.
COIL = 6.0

CASES = {
    "SC":  dict(lattice="SC", max_func=4, degree_dist="0:0,1:0", mix={}),
    "BCC": dict(lattice="BCC", max_func=8, degree_dist="0:0,1:0", mix={}),
    "FCC": dict(lattice="FCC", max_func=8, degree_dist="0:0,1:0", mix={}),
    "MIX": dict(lattice="MIX", max_func=8, mix={"SC": 0.2, "BCC": 0.4, "FCC": 0.4},
                degree_dist="0:0,1:0,2:15,3:40,4:55,5:35,6:15,7:7,8:3"),
}


def random_walk(a0, a1, n_bonds, bond, rng):
    """Unit-bond walk from a0 to a1: a chain shape, at exactly the contour
    its beads need.

    Each step is drawn from the cone of directions that still leaves the far
    end reachable, so the walk closes on its junction without any bond being
    rescaled afterwards.
    """
    pts = [np.asarray(a0, float)]
    end = np.asarray(a1, float)
    for k in range(n_bonds - 1):
        d = end - pts[-1]
        r = float(np.linalg.norm(d))
        rem = n_bonds - k - 1
        cmin = ((r * r + bond * bond - (rem * bond) ** 2) / (2 * r * bond)
                if r > 1e-12 else -1.0)
        cmin = min(1.0, max(-1.0, cmin))
        c = rng.uniform(cmin, 1.0)
        s = np.sqrt(max(0.0, 1.0 - c * c))
        dh = d / r if r > 1e-12 else np.array([0.0, 0.0, 1.0])
        t = np.cross(dh, [0, 0, 1.0] if abs(dh[2]) < 0.9 else [1.0, 0, 0])
        t /= np.linalg.norm(t)
        u = np.cross(dh, t)
        phi = rng.uniform(0.0, 2.0 * np.pi)
        pts.append(pts[-1] + bond * (c * dh + s * (np.cos(phi) * t
                                                   + np.sin(phi) * u)))
    pts.append(end)
    return np.array(pts)


def build(name, coil, dp, bond, want, turns, half, seed):
    spec = dict(LATTICE)
    spec.update(CASES[name])
    graph = build_network(spec)
    geo = geometry(graph, dp=dp, bond=bond, coil=coil)
    ch, ends = geo["chords"], geo["ends"]
    rng = np.random.default_rng(seed)

    paths = {k: random_walk(c0, c1, dp + 1, bond, rng)
             for k, (c0, c1) in ch.items()}

    # Chains sharing a crosslink touch there by construction; a winding on a
    # junction is not an entanglement between the two chains.
    shared = {frozenset((a, b)) for a in ch for b in ch
              if a < b and set(ends[a]) & set(ends[b])}
    found = find_contacts(paths, box=geo["L"], max_sep=2.5, margin=6,
                          exclude=shared)
    wound, applied = wind_all(paths, found[:want], turns=turns, half=half,
                              box=geo["L"])
    return graph, geo, paths, wound, applied, len(found)


def measure(paths, wound, applied, box):
    """Linking before and after, on the built paths."""
    rows = []
    for c in applied:
        pa0, pb0 = paths[c.chain_a], paths[c.chain_b]
        pa1, pb1 = wound[c.chain_a], wound[c.chain_b]
        d = pa0.mean(0) - pb0.mean(0)
        shift = box * np.round(d / box)
        before = far_closed_linking(pa0, pb0 + shift)
        after = far_closed_linking(pa1, pb1 + shift)
        sep = float(np.linalg.norm(pa1[:, None, :] - (pb1 + shift)[None, :, :],
                                   axis=-1).min())
        rows.append(dict(pair=(c.chain_a, c.chain_b), sep_found=c.sep,
                         before=before, after=after, delta=after - before,
                         sep_built=sep))
    return rows


def plot(name, geo, wound, applied, path):
    import plotly.graph_objects as go

    fig = go.Figure()
    xs, ys, zs = [], [], []
    used = {k for c in applied for k in c.pair}
    for k, p in wound.items():
        if k in used:
            continue
        xs += list(p[:, 0]) + [None]
        ys += list(p[:, 1]) + [None]
        zs += list(p[:, 2]) + [None]
    if xs:
        fig.add_trace(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="lines", name="other chains",
            line=dict(color="rgba(150,150,150,0.25)", width=1),
            hoverinfo="skip"))

    palette = ["#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
               "#42d4f4", "#f032e6", "#bfef45", "#469990", "#9a6324"]
    for i, c in enumerate(applied):
        col = palette[i % len(palette)]
        for k in c.pair:
            p = wound[k]
            fig.add_trace(go.Scatter3d(
                x=p[:, 0], y=p[:, 1], z=p[:, 2], mode="lines",
                line=dict(color=col, width=4),
                name=f"{c.chain_a}-{c.chain_b}",
                legendgroup=f"e{i}", showlegend=(k == c.chain_a),
                hovertemplate=f"chain {k}<extra></extra>"))
        m = 0.5 * (wound[c.chain_a][c.i_a] + wound[c.chain_b][c.i_b])
        fig.add_trace(go.Scatter3d(
            x=[m[0]], y=[m[1]], z=[m[2]], mode="markers",
            marker=dict(size=5, color=col, symbol="diamond"),
            legendgroup=f"e{i}", showlegend=False,
            hovertemplate=(f"{c.chain_a}-{c.chain_b}<br>"
                           f"paths were {c.sep:.2f} sigma apart"
                           f"<extra></extra>")))

    jx = np.array(list(geo["pos"].values()))
    fig.add_trace(go.Scatter3d(
        x=jx[:, 0], y=jx[:, 1], z=jx[:, 2], mode="markers", name="crosslinks",
        marker=dict(size=2, color="#222", opacity=0.5), hoverinfo="skip"))

    L = geo["L"]
    fig.update_layout(
        title=(f"{name}: {len(applied)} designed entanglements, "
               f"coil {COIL:g}, {len(wound)} chains. Drag to rotate."),
        scene=dict(xaxis=dict(range=[0, L[0]], title="x / sigma"),
                   yaxis=dict(range=[0, L[1]], title="y / sigma"),
                   zaxis=dict(range=[0, L[2]], title="z / sigma"),
                   aspectmode="cube"),
        template="plotly_white", margin=dict(l=0, r=0, t=60, b=0),
        legend=dict(font=dict(size=9)))
    fig.write_html(str(path), include_plotlyjs="inline")
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lattices", nargs="*", default=list(CASES))
    ap.add_argument("--coil", type=float, default=COIL)
    ap.add_argument("--want", type=int, default=10)
    ap.add_argument("--turns", type=int, default=1)
    ap.add_argument("--half", type=int, default=10)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--bond", type=float, default=BOND)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--no-md", action="store_true")
    ap.add_argument("--stages", type=int, default=3)
    args = ap.parse_args()

    out = OUT / "all"
    out.mkdir(parents=True, exist_ok=True)
    summary = []

    for name in args.lattices:
        print(f"\n{'=' * 66}\n{name}\n{'=' * 66}")
        graph, geo, paths, wound, applied, n_found = build(
            name, args.coil, args.dp, args.bond, args.want, args.turns,
            args.half, args.seed)
        L = geo["L"]
        rows = measure(paths, wound, applied, L)
        exact = sum(1 for r in rows if abs(abs(r["delta"]) - args.turns) < 0.5)
        bonds = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                                for p in wound.values()])
        print(f"  {graph.number_of_nodes()} crosslinks, "
              f"{graph.number_of_edges()} chains, DP {args.dp}")
        print(f"  box {L[0]:.1f} sigma, density {geo['density']:.4f}, "
              f"coil {args.coil:g}")
        print(f"  {n_found} contacts available, {len(applied)} wound")
        print(f"  bonds {bonds.min():.3f} to {bonds.max():.3f}")
        print(f"  linking as built: {exact} of {len(rows)} exact")

        root = OUT / f"all_{name}"
        n_atoms, node_atom, chain_atoms = write_system(graph, geo, wound, root)
        seq = {k: chain_ids(k, node_atom, chain_atoms, geo["ends"])
               for k in {x for c in applied for x in c.pair}}
        sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                                 protocol="hardcore")

        z1_dir = OUT / f"all_{name}_z1"
        if z1_dir.exists():
            for f in z1_dir.glob("*.Z1"):
                f.unlink()
        z1_dir.mkdir(parents=True, exist_ok=True)
        z1_export(root / "03_Conformation/system_conformed.data",
                  [seq[c.chain_a] for c in applied[:1]]
                  + [seq[c.chain_b] for c in applied[:1]],
                  z1_dir / "_probe.Z1")
        for i, c in enumerate(applied):
            z1_export(root / "03_Conformation/system_conformed.data",
                      [seq[c.chain_a], seq[c.chain_b]],
                      z1_dir / f"built_{i:02d}.Z1")

        if not args.no_md:
            print(f"\n  --- LAMMPS, {args.stages} stage(s) ---")
            run_md(sim, args.stages)
            final = {1: "system_after_soft.data", 2: "system_ramped.data",
                     3: "system_equilibrated.data"}[args.stages]
            f = root / "04_Simulation" / final
            if f.exists():
                for i, c in enumerate(applied):
                    z1_export(f, [seq[c.chain_a], seq[c.chain_b]],
                              z1_dir / f"final_{i:02d}.Z1")
            print()
            report_bonds(root)

        z = run_z1(z1_dir) or {}
        hit_b = hit_f = 0
        for i in range(len(applied)):
            vb = z.get(f"built_{i:02d}")
            vf = z.get(f"final_{i:02d}")
            if vb and all(v == args.turns for v in vb):
                hit_b += 1
            if vf and all(v == args.turns for v in vf):
                hit_f += 1
        n = len(applied)
        print(f"\n  Z1+ as built:  {hit_b} of {n} exact")
        if not args.no_md:
            print(f"  Z1+ after MD:  {hit_f} of {n} exact")

        html = plot(name, geo, wound, applied, out / f"{name}.html")
        print(f"  wrote {html.name}")
        summary.append(dict(name=name, nodes=graph.number_of_nodes(),
                            chains=graph.number_of_edges(), box=float(L[0]),
                            density=geo["density"], found=n_found, wound=n,
                            exact_lk=exact, z1_built=hit_b, z1_final=hit_f,
                            bmin=float(bonds.min()), bmax=float(bonds.max())))

    md = out / "RESULTS.md"
    with open(md, "w", encoding="utf-8") as f:
        f.write("# Designed entanglements on every lattice\n\n")
        f.write(f"Coil {args.coil:g}, DP {args.dp}, bond {args.bond}, "
                f"{args.turns} winding per pair, seed {args.seed}. "
                f"Hard-core three-stage protocol.\n\n")
        f.write("Each pair measured on its own with Z1+, the rest of the "
                "network removed, so every count is between those two chains "
                "and nothing else.\n\n")
        f.write("| lattice | crosslinks | chains | box | density | contacts "
                "found | wound | Z1+ as built | Z1+ after MD | bonds |\n")
        f.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for s in summary:
            f.write(f"| {s['name']} | {s['nodes']} | {s['chains']} | "
                    f"{s['box']:.0f} σ | {s['density']:.3f} | {s['found']} | "
                    f"{s['wound']} | **{s['z1_built']} / {s['wound']}** | "
                    f"**{s['z1_final']} / {s['wound']}** | "
                    f"{s['bmin']:.2f}–{s['bmax']:.2f} |\n")
        f.write("\n## How the winding is made\n\n"
                "The two chains are not moved. Where their paths already come "
                "close, their separation is rotated about the local axis "
                "through `2*pi*turns` across a short window. The midpoint of "
                "the pair is unchanged everywhere, both ends of the window "
                "are whole turns so the path outside is untouched, and the "
                "cost is the circumference of a circle of radius `sep/2` — "
                "about 0.4 σ for one winding, against the tens of σ the "
                "chord-based construction needed.\n\n"
                "## Why the coil matters\n\n"
                "A chain carrying only 1.8× its chord in contour is nearly "
                "extended and meets nothing, so there is nowhere to put a "
                "winding. From about 6×, chains interpenetrate and contacts "
                "are everywhere. Every limit reported before this — SC "
                "failing, outer shells failing — came from running at 1.8.\n")
    print(f"\nwrote {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
