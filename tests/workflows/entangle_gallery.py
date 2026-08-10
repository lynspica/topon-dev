"""A gallery of designed entanglement cases, each verified and drawn.

    python tests/workflows/entangle_gallery.py
    python tests/workflows/entangle_gallery.py --cases pair2 ring
    python tests/workflows/entangle_gallery.py --rounds 6

Each case names a small topology, searches for chain paths that deliver it,
measures the result with Z1+, and writes an interactive view. The markers in
the view are not where the sites were aimed -- they are the entanglement
points Z1+ found, read back out of its shortest-path output, so what is drawn
is what was measured.

The cases vary one thing at a time: how many entanglements on a pair, how far
apart the partners are, how many partners one chain carries, whether the
requested topology is a chain or a closed ring, and which lattice underneath.
"""
from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.paths import bridging_walk  # noqa: E402
from tests.workflows.entangle_all import CASES as LATTICES  # noqa: E402
from tests.workflows.entangle_design import measure_one, route_one  # noqa: E402
from tests.workflows.entangle_search import Wish, _both  # noqa: E402
from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    DP,
    LATTICE,
    OUT,
    build_network,
    geometry,
)

PALETTE = ["#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
           "#42d4f4", "#f032e6", "#9a6324"]


def kink_points(sp_file, chains_of_interest):
    """Entanglement coordinates from Z1+'s shortest-path output.

    Returns ``[(x, y, z, chain, partner), ...]`` for the chains asked about.
    These are measured positions, not the waypoints the search proposed --
    which is the point of drawing them.
    """
    rows = [ln.split() for ln in Path(sp_file).read_text().splitlines()
            if ln.strip()]
    i = 0
    n = int(float(rows[i][0]))
    i += 2
    out = []
    for c in range(1, n + 1):
        m = int(float(rows[i][0]))
        i += 1
        for _ in range(m):
            r = rows[i]
            i += 1
            if len(r) >= 7 and int(float(r[4])) == 1 and c in chains_of_interest:
                out.append((float(r[0]), float(r[1]), float(r[2]),
                            c, int(float(r[5]))))
    return out


def plan_for(name, keys, ch, ends, L):
    """The named topology, as [(routed chain, {partner: count}), ...]."""

    def order(a):
        m = 0.5 * (ch[a][0] + ch[a][1])
        return [b for _, b in sorted(
            (float(np.linalg.norm(
                (lambda v: v - L * np.round(v / L))(
                    0.5 * (ch[b][0] + ch[b][1]) - m))), b)
            for b in keys if b != a and not set(ends[a]) & set(ends[b]))]

    A = keys[0]
    near, mid, far = order(A)[2], order(A)[len(keys) // 3], order(A)[len(keys) // 2]

    if name == "pair1":
        return [(A, {mid: 1})], "one entanglement, mid-range partner"
    if name == "pair2":
        return [(A, {mid: 2})], "two entanglements on the same pair"
    if name == "pair3":
        return [(A, {mid: 3})], "three entanglements on the same pair"
    if name == "near":
        return [(A, {near: 1})], "one entanglement, nearest available partner"
    if name == "far":
        return [(A, {far: 1})], "one entanglement, a distant partner"
    if name == "hub":
        o = order(A)
        return ([(A, {o[len(keys) // 3]: 1, o[len(keys) // 5]: 1,
                      o[len(keys) // 8]: 1})],
                "one chain carrying three different partners")
    if name == "ring":
        # Walk a cycle through near neighbours rather than picking four
        # chains by index. Chains spread across the box cannot reach one
        # another: taking every 26th chain gave a ring where no link formed
        # at all, because each pair was most of a box apart.
        a = keys[0]
        ring = [a]
        for _ in range(3):
            nxt = next(x for x in order(ring[-1])
                       if x not in ring and all(
                           not set(ends[x]) & set(ends[y]) for y in ring))
            ring.append(nxt)
        w, x, y, z = ring
        return ([(w, {x: 1}), (x, {y: 1}), (y, {z: 1}), (z, {w: 1})],
                "a closed ring: A-B, B-C, C-D, D-A")
    if name == "mixed":
        o = order(A)
        return ([(A, {o[len(keys) // 3]: 2}), (keys[40], {o[len(keys) // 4]: 1})],
                "two chains, different counts, sharing a neighbourhood")
    raise SystemExit(f"unknown case {name}")


def draw(name, note, geo, paths, keys, plan, idx, kinks, path_out):
    import plotly.graph_objects as go

    L = geo["L"]
    routed = {r for r, _ in plan}
    targets = {t for _, w in plan for t in w}
    involved = routed | targets

    fig = go.Figure()
    xs, ys, zs = [], [], []
    for k, p in paths.items():
        if k in involved:
            continue
        xs += list(p[:, 0]) + [None]
        ys += list(p[:, 1]) + [None]
        zs += list(p[:, 2]) + [None]
    fig.add_trace(go.Scatter3d(
        x=xs, y=ys, z=zs, mode="lines", name=f"{len(paths) - len(involved)} "
        f"other chains", line=dict(color="rgba(150,150,150,0.16)", width=1),
        hoverinfo="skip"))

    for i, k in enumerate(sorted(involved)):
        p = paths[k]
        role = "routed" if k in routed else "target"
        fig.add_trace(go.Scatter3d(
            x=p[:, 0], y=p[:, 1], z=p[:, 2], mode="lines",
            line=dict(color=PALETTE[i % len(PALETTE)],
                      width=6 if k in routed else 4,
                      dash=None if k in routed else "dot"),
            name=f"chain {k} ({role})",
            hovertemplate=f"chain {k}, {role}<extra></extra>"))

    if kinks:
        K = np.array([[x, y, z] for x, y, z, _, _ in kinks])
        lab = [f"chain {c} threaded by {p}" for _, _, _, c, p in kinks]
        fig.add_trace(go.Scatter3d(
            x=K[:, 0], y=K[:, 1], z=K[:, 2], mode="markers",
            marker=dict(size=6, color="#000", symbol="x"),
            name="entanglement points (measured)",
            hovertext=lab, hoverinfo="text"))

    fig.update_layout(
        title=f"{name} — {note}",
        scene=dict(xaxis=dict(title="x / σ"), yaxis=dict(title="y / σ"),
                   zaxis=dict(title="z / σ"), aspectmode="data"),
        template="plotly_white", margin=dict(l=0, r=0, t=54, b=0),
        legend=dict(font=dict(size=10)))
    fig.write_html(str(path_out), include_plotlyjs="inline")


def run_case(name, lattice, rounds, per_round, dp, seed, out_dir):
    spec = dict(LATTICE)
    spec.update(LATTICES[lattice])
    spec["dims"] = (4, 4, 4)
    graph = build_network(spec)
    geo = geometry(graph, dp=dp, density=0.85)
    ch, ends, L = geo["chords"], geo["ends"], geo["L"]
    keys = sorted(ch)
    idx = {k: i + 1 for i, k in enumerate(keys)}

    rng = np.random.default_rng(seed)
    paths = {k: bridging_walk(c0, c1, dp + 1, BOND, rng)
             for k, (c0, c1) in ch.items()}
    plan, note = plan_for(name, keys, ch, ends, L)

    work = OUT / "gallery_work"
    base = measure_one(paths, keys, geo, work, "base")
    for _ in range(2):
        todo = [(r, w) for r, w in plan
                if any(_both(base, idx[r], idx[t]) != n for t, n in w.items())]
        if not todo:
            break
        for routed, want in todo:
            w = {idx[t]: n for t, n in want.items()}
            got0 = base.get(idx[routed], collections.Counter())
            wish = Wish(chain=idx[routed], want=w, penalty=1.0,
                        baseline=sum(v for q, v in got0.items() if q not in w))
            best = route_one(paths, keys, geo, routed, wish, rounds,
                             per_round, rng, work)
            if best is not None:
                paths[routed] = best[4]
        base = measure_one(paths, keys, geo, work, "step")

    rows = []
    for routed, want in plan:
        for t, n in want.items():
            rows.append((routed, t, n, _both(base, idx[routed], idx[t])))
    involved = {idx[r] for r, _ in plan} | {idx[t] for _, w in plan for t in w}
    sp = work / "SP_step.dat"
    kinks = kink_points(sp, involved) if sp.exists() else []
    bonds = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                            for p in paths.values()])

    out_dir.mkdir(parents=True, exist_ok=True)
    draw(f"{name} on {lattice}", note, geo, paths, keys, plan, idx, kinks,
         out_dir / f"{name}_{lattice}.html")
    return rows, note, bonds, len(keys)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases", nargs="*",
                    default=["pair1", "pair2", "pair3", "near", "far",
                             "hub", "ring", "mixed"])
    ap.add_argument("--lattice", default="SC", choices=sorted(LATTICES))
    ap.add_argument("--also", nargs="*", default=["MIX"],
                    help="repeat the first case on these lattices too")
    ap.add_argument("--rounds", type=int, default=5)
    ap.add_argument("--per-round", type=int, default=16)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    out_dir = OUT / "gallery"
    summary = []
    jobs = [(c, args.lattice) for c in args.cases]
    jobs += [(args.cases[0], lat) for lat in args.also]

    for name, lat in jobs:
        print(f"\n=== {name} on {lat} ===")
        rows, note, bonds, n_chains = run_case(
            name, lat, args.rounds, args.per_round, args.dp, args.seed,
            out_dir)
        ok = sum(1 for _, _, w, g in rows if w == g)
        for a, b, w, g in rows:
            print(f"  {a}-{b}: asked {w}, got {g}" + ("   ok" if w == g else ""))
        print(f"  {ok} of {len(rows)} exact, bonds "
              f"{bonds.min():.2f} to {bonds.max():.2f}")
        summary.append((name, lat, note, rows, ok, bonds, n_chains))

    md = out_dir / "GALLERY.md"
    with open(md, "w", encoding="utf-8") as f:
        f.write("# Designed entanglement: worked examples\n\n")
        f.write("Each case names a topology, searches for chain paths that "
                "deliver it, and verifies the result with Z1+ on the built "
                "configuration. The markers in each view are the entanglement "
                "points Z1+ found, not the waypoints the search proposed.\n\n")
        f.write("| case | lattice | what it asks for | delivered | bonds |\n")
        f.write("|---|---|---|---|---|\n")
        for name, lat, note, rows, ok, bonds, n in summary:
            detail = ", ".join(f"{a}-{b} x{w}" for a, b, w, _ in rows)
            f.write(f"| [{name}]({name}_{lat}.html) | {lat} | {note}"
                    f" ({detail}) | **{ok} of {len(rows)}** | "
                    f"{bonds.min():.2f}–{bonds.max():.2f} |\n")
        f.write("\n## Reading a view\n\n"
                "Solid coloured chains are the ones that were routed; dotted "
                "ones are their requested partners. Everything else is grey "
                "and faint. Black crosses are measured entanglement points — "
                "hover one to see which chain is threaded by which.\n\n"
                "All cases are DP 80 at melt density, so every chain in the "
                "box carries several entanglements of its own. Only the "
                "coloured ones were asked for.\n")
    print(f"\nwrote {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
