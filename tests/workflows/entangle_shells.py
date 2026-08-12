"""Can a designed entanglement be placed in an outer neighbour shell?

    python tests/workflows/entangle_shells.py
    python tests/workflows/entangle_shells.py --shells 1 2 3 4 --per-shell 4

The old kink bulges a chain sideways toward a neighbour, so it can only reach
whatever is already adjacent. `compute_shell_weights` records what that costs:
band 1 delivered 5 of 7, band 2 delivered 2 of 16, band 3 delivered 0 of 16.
Weighting the outer shells up did not help, because the reach was the limit,
not the weighting.

That measurement is about the old construction and does not have to hold for
the new one. `construct` routes a chain out to its partner and back rather than
bulging it, and a chain carries far more contour than its chord needs, so an
outer shell should be reachable. This measures whether it is.

One pair per shell at a time, each routed on its own into the same relaxed
melt, so the shells are compared under identical conditions and nothing any
other design did can be mistaken for a shell effect.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import networkx as nx

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.assignment.entanglements import (  # noqa: E402
    chain_distances,
    neighbour_shells,
)
from topon.conformation.paths import Clearance, bridging_walk  # noqa: E402
from tests.workflows.entangle_all import CASES  # noqa: E402
from tests.workflows.entangle_relaxed import (  # noqa: E402
    construct_exact,
    measure_pairs,
    paths_from,
    rewrite_coords,
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
    run_md,
    scale_for_design,
    write_system,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shells", type=int, nargs="*", default=[1, 2, 3, 4])
    ap.add_argument("--per-shell", type=int, default=3)
    ap.add_argument("--want", type=int, default=2)
    ap.add_argument("--dims", type=int, default=4)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--ring", type=float, default=2.0)
    ap.add_argument("--margin", type=float, default=2.0)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--ranks", type=int, default=4)
    ap.add_argument("--phases", type=int, default=4)
    ap.add_argument("--clearance", type=float, default=0.9)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--survive", action="store_true",
                    help="relax each built pair on its own and measure it "
                         "again, which is the like-for-like comparison "
                         "against the old kink's shell figures")
    ap.add_argument("--tag", default="_shells")
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (args.dims,) * 3
    graph = build_network(spec)
    dims = np.asarray(graph.graph["box"], float)

    G = nx.MultiGraph()
    G.add_nodes_from(graph.nodes(data=True))
    for u, v in graph.edges():
        G.add_edge(u, v)

    dist = chain_distances(G, dims)
    shells = neighbour_shells(G, dims, max_shell=max(args.shells),
                              distances=dist)

    # Chain keys here are (u, v, key); geometry indexes chains by their sorted
    # edge position. Map one to the other so a shell lookup can name a pair the
    # conformation code understands.
    edges = sorted(graph.edges())
    idx_of = {frozenset(e): i for i, e in enumerate(edges)}
    # `dist` is keyed by the MultiGraph's (u, v, key) triples; record the gap
    # by chain index so it can be reported next to the result.
    gap_of = {}
    for (ca, cb), r in dist.items():
        ia = idx_of.get(frozenset((ca[0], ca[1])))
        ib = idx_of.get(frozenset((cb[0], cb[1])))
        if ia is not None and ib is not None:
            gap_of[(min(ia, ib), max(ia, ib))] = r

    # Pick pairs shell by shell, each from a different routed chain so no chain
    # is asked to serve two shells at once.
    rng = np.random.default_rng(args.seed)
    picks, used = {}, set()
    for s in args.shells:
        got = []
        for chain, by in sorted(shells.items()):
            if s not in by:
                continue
            a = idx_of.get(frozenset((chain[0], chain[1])))
            if a is None or a in used:
                continue
            for other in by[s]:
                b = idx_of.get(frozenset((other[0], other[1])))
                if b is None or b in used or b == a:
                    continue
                got.append((a, b))
                used.update((a, b))
                break
            if len(got) >= args.per_shell:
                break
        picks[s] = got

    flat = [(s, a, b) for s, ps in picks.items() for a, b in ps]
    if not flat:
        print("  no pairs found in the requested shells")
        return 1
    print("  pairs picked: "
          + "; ".join(f"shell {s}: {len(ps)}" for s, ps in picks.items()))

    # Size the box for the longest route in the set, so an outer-shell pair is
    # not judged on a box that could never hold its detour.
    pairs_g = [(a, b) for _s, a, b in flat]
    sc = scale_for_design(graph, pairs_g, dp=args.dp, bond=BOND,
                          radius=args.ring, margin=args.margin,
                          site_span=(0.3, 0.7))
    geo = geometry(graph, dp=args.dp, scale=sc)
    keys = sorted(geo["chords"])
    n_beads = graph.number_of_edges() * args.dp + graph.number_of_nodes()
    print(f"  box {geo['L'][0]:.1f} sigma, density "
          f"{n_beads / float(np.prod(geo['L'])):.3f}")

    rng0 = np.random.default_rng(args.seed)
    paths0 = {k: bridging_walk(c0, c1, args.dp + 1, BOND, rng0)
              for k, (c0, c1) in geo["chords"].items()}
    root = OUT / f"relaxed{args.tag}"
    shutil.rmtree(root, ignore_errors=True)
    _n, node_atom, chain_atoms = write_system(graph, geo, paths0, root)
    seq = {k: chain_ids(k, node_atom, chain_atoms, geo["ends"]) for k in keys}
    sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                             protocol="hardcore")
    print("  --- relaxing a plain melt ---")
    run_md(sim, args.stages)
    after = {1: "system_after_soft.data", 2: "system_ramped.data",
             3: "system_equilibrated.data"}[args.stages]
    relaxed = root / "04_Simulation" / after
    if not relaxed.exists():
        print("  relaxation produced no output")
        return 1
    box, paths, xyz0 = paths_from(relaxed, keys, seq)

    work = OUT / f"shells_work_{os.getpid()}"
    print(f"\n  {'shell':>6} {'pair':>10} {'gap':>7} {'built':>7} "
          f"{'wanted':>7}")
    rows = []
    for s, a, b in flat:
        mine = set(seq[a])
        avoid = Clearance(np.array([xyz0[i] for i in sorted(xyz0)
                                    if i not in mine]), box, args.clearance)
        try:
            p, got = construct_exact(
                paths, a, b, args.want, box, avoid, seq, relaxed, work,
                (0.3, 0.7), args.ring, args.dp,
                ranks=args.ranks, phases=args.phases)
        except ValueError as e:
            print(f"  {s:>6} {f'{a}-{b}':>10} {'':>7} {str(e)[:34]}")
            rows.append((s, None, a, b, None))
            continue
        if p is None:
            print(f"  {s:>6} {f'{a}-{b}':>10} {'':>7} {'nothing built':>7}")
            rows.append((s, None, a, b, None))
            continue
        gap = gap_of.get((min(a, b), max(a, b)), float("nan"))
        print(f"  {s:>6} {f'{a}-{b}':>10} {gap:7.2f} "
              f"{got:>7} {args.want:>7}"
              + ("   ok" if got == args.want else ""))
        rows.append((s, got, a, b, p))

    # Survival, one pair at a time.
    #
    # The figures this is being compared against were measured after the full
    # protocol, so "built as designed" is not the same claim and must not be
    # presented as though it were. Each pair is relaxed in its own copy of the
    # melt, alone: designed windings interfere with each other, and that is a
    # separate effect that would otherwise be charged to the shell.
    survived = {}
    if args.survive:
        cand = OUT / f"shells_cand_{os.getpid()}"
        (cand / "03_Conformation").mkdir(parents=True, exist_ok=True)
        for sub in ("02_Chemistry", "04_Simulation"):
            (cand / sub).mkdir(parents=True, exist_ok=True)
            for f in (root / sub).glob("*"):
                if f.is_file() and (sub == "02_Chemistry"
                                    or f.suffix == ".in"):
                    (cand / sub / f.name).write_text(f.read_text())

        print(f"\n  --- relaxing each pair on its own ---")
        for s, got, a, b, p in rows:
            if p is None:
                continue
            rewrite_coords(relaxed, cand / "03_Conformation" /
                           "system_relaxed.data", dict(zip(seq[a], p)))
            for stale in ("system_after_soft.data", "system_ramped.data",
                          "system_equilibrated.data"):
                q = cand / "04_Simulation" / stale
                if q.exists():
                    q.unlink()
            try:
                run_md(cand / "04_Simulation", args.stages)
            except Exception:
                continue
            out = cand / "04_Simulation" / after
            if not out.exists():
                continue
            pair = (min(a, b), max(a, b))
            survived[pair] = measure_pairs(out, seq, [pair], work)[pair]

    print(f"\n  {'shell':>6} {'asked':>7} {'built exactly':>14} "
          f"{'built nothing':>14}"
          + (f" {'survived exactly':>17}" if args.survive else ""))
    for s in args.shells:
        mine = [(g, a, b) for t, g, a, b, _p in rows if t == s]
        if not mine:
            continue
        exact = sum(1 for g, _a, _b in mine if g == args.want)
        none = sum(1 for g, _a, _b in mine if not g)
        line = f"  {s:>6} {len(mine):>7} {exact:>14} {none:>14}"
        if args.survive:
            kept = sum(1 for _g, a, b in mine
                       if survived.get((min(a, b), max(a, b))) == args.want)
            line += f" {kept:>17}"
        print(line)
    if not args.survive:
        print("\n  built as designed only. The old-kink figures this is "
              "compared against were measured after the full protocol, so "
              "pass --survive for a like-for-like number.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
