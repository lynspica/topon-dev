"""What do you set to get N entanglements between two named chains?

    python tests/workflows/entangle_calibrate.py
    python tests/workflows/entangle_calibrate.py --spans 0.5 0.75 1.0 --seeds 5

The search built what it stumbled on. Asked for one it usually returned two,
and which one you got depended on the seed, so the same request did not give
the same system twice. That is a calibration problem wearing a search's
clothes: if the winding geometry sets the count, the count should be
constructed, not hunted for.

This measures the map. One pair, one site, everything held fixed except how
far round the target the chain is taken, swept over several seeds. Counted
twice: as built, and again after minimisation, because a number that does not
survive is not a number you can order.
"""
from __future__ import annotations

import argparse
import collections
import os
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.paths import (  # noqa: E402
    Clearance,
    bridging_walk,
    loop_around,
    route_through,
    walk_through,
)
from tests.workflows.entangle_all import CASES  # noqa: E402
from tests.workflows.entangle_relaxed import (  # noqa: E402
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
    read_data,
    run_md,
    scale_for_design,
    write_system,
)


def ring_points(radius, span, bond=BOND, spacing=1.6):
    """How many waypoints an arc of this length needs.

    Waypoints are landed on exactly and the walk fills in between them, so
    spacing them much further than a bond or two hands the winding back to a
    random walk that does not know it is supposed to be winding.
    """
    arc = 2.0 * np.pi * radius * max(span, 0.05)
    return int(np.clip(round(arc / (spacing * bond)), 4, 24))


def build(routed_path, target_path, at, radius, span, phase, avoid, rng,
          dp=DP, bond=BOND, taut=False):
    """One routed chain taken `span` turns around the target at `at`.

    With ``taut`` the legs are zigzags rather than random walks, so the chain's
    leftover contour is spent in a fixed shape instead of wandering back over
    the target and adding crossings of its own.
    """
    i = int(np.clip(round(at * len(target_path)), 1, len(target_path) - 2))
    ring = loop_around(target_path, i, radius, ring_points(radius, span),
                       phase, avoid, span)
    if taut:
        return route_through(routed_path[0], routed_path[-1], list(ring),
                             dp + 1, bond, target_path[i], avoid)
    return walk_through(routed_path[0], routed_path[-1], list(ring),
                        dp + 1, bond, rng, avoid)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spans", type=float, nargs="*",
                    default=[0.35, 0.5, 0.65, 0.8, 1.0, 1.25, 1.5, 2.0])
    ap.add_argument("--radius", type=float, default=2.0)
    ap.add_argument("--at", type=float, default=0.5)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--dims", type=int, default=4)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--clearance", type=float, default=0.9)
    ap.add_argument("--taut", action="store_true",
                    help="spend the leftover contour as a zigzag rather than "
                         "a random walk, so the count belongs to the design")
    ap.add_argument("--scan-pairs", type=int, default=0,
                    help="instead of sweeping span, hold span fixed and scan "
                         "this many partners. Tests whether the count's "
                         "parity belongs to the pair rather than the winding.")
    ap.add_argument("--margin", type=float, default=2.0,
                    help="contour headroom. The chord-based estimate is made "
                         "before relaxation, and relaxation moves chains, so "
                         "the real detour is longer than the one the box was "
                         "sized for.")
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (args.dims,) * 3
    graph = build_network(spec)

    box_g = np.asarray(graph.graph["box"], float)
    raw = {n: np.asarray(d["pos"], float) for n, d in graph.nodes(data=True)}
    edges = sorted(graph.edges())

    def chord_pts(k, span=(0.3, 0.7)):
        u, v = edges[k]
        a = raw[u]
        mic = (raw[v] - a) - box_g * np.round((raw[v] - a) / box_g)
        return a + np.linspace(span[0], span[1], 24)[:, None] * mic

    routed = 0
    best, target = np.inf, None
    A = chord_pts(routed, (0.0, 1.0))
    for b in range(len(edges)):
        if b == routed or set(edges[routed]) & set(edges[b]):
            continue
        d = A[:, None, :] - chord_pts(b)[None, :, :]
        d -= box_g * np.round(d / box_g)
        gap = float(np.linalg.norm(d, axis=2).min())
        if gap < best:
            best, target = gap, b

    sc = scale_for_design(graph, [(routed, target)], dp=args.dp, bond=BOND,
                          radius=args.radius, site_span=(0.3, 0.7),
                          margin=args.margin)
    geo = geometry(graph, dp=args.dp, scale=sc)
    keys = sorted(geo["chords"])
    n_beads = graph.number_of_edges() * args.dp + graph.number_of_nodes()
    rho = n_beads / float(np.prod(geo["L"]))
    print(f"  chain {routed} around chain {target}, box {geo['L'][0]:.1f} "
          f"sigma, density {rho:.3f}")

    rng0 = np.random.default_rng(42)
    paths0 = {k: bridging_walk(c0, c1, args.dp + 1, BOND, rng0)
              for k, (c0, c1) in geo["chords"].items()}
    root = OUT / "calibrate"
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

    # Both chains are unwrapped, but not necessarily into the same periodic
    # image, and a route to the wrong copy of the target is a route across the
    # whole box: measured, 166 sigma of travel for a chain carrying 77.
    paths[target] = paths[target] + box * np.round(
        (paths[routed].mean(0) - paths[target].mean(0)) / box)

    mine = set(seq[routed])
    avoid = Clearance(np.array([xyz0[i] for i in sorted(xyz0)
                                if i not in mine]), box, args.clearance)

    pair = (min(routed, target), max(routed, target))
    work = OUT / f"calibrate_work_{os.getpid()}"
    cand = OUT / f"calibrate_cand_{os.getpid()}"
    (cand / "03_Conformation").mkdir(parents=True, exist_ok=True)
    for sub in ("02_Chemistry", "04_Simulation"):
        (cand / sub).mkdir(parents=True, exist_ok=True)
        for f in (root / sub).glob("*"):
            if f.is_file() and (sub == "02_Chemistry"
                                or f.suffix == ".in"):
                (cand / sub / f.name).write_text(f.read_text())

    def relax_and_count(data_file):
        (cand / "03_Conformation" / "system_relaxed.data").write_text(
            Path(data_file).read_text())
        for stale in ("system_after_soft.data", "system_ramped.data"):
            q = cand / "04_Simulation" / stale
            if q.exists():
                q.unlink()
        try:
            run_md(cand / "04_Simulation", args.stages, quiet=True)
        except TypeError:
            run_md(cand / "04_Simulation", args.stages)
        except Exception:
            return None
        out = cand / "04_Simulation" / after
        if not out.exists():
            return None
        return measure_pairs(out, seq, [pair], work).get(pair)

    if args.scan_pairs:
        # Is the count's parity a property of the pair?
        #
        # A chain that leaves one junction, detours out to a partner and comes
        # back to the other junction on the same side of it crosses that
        # partner an even number of times, whatever it does when it gets
        # there. An odd count needs the partner to lie between the two
        # junctions. If that is what is going on, parity should sort by pair
        # and ignore the span, which would explain why sweeping the span
        # never produced a one.
        A = paths[routed]
        cands = []
        for b in keys:
            if b == routed or set(geo["ends"][routed]) & set(geo["ends"][b]):
                continue
            q = paths[b] + box * np.round((A.mean(0) - paths[b].mean(0)) / box)
            d = A[:, None, :] - q[None, :, :]
            d -= box * np.round(d / box)
            cands.append((float(np.linalg.norm(d, axis=2).min()), b, q))
        cands.sort(key=lambda t: t[0])

        chord = A[-1] - A[0]
        span_len = max(float(np.linalg.norm(chord)), 1e-9)
        n_hat = chord / span_len
        print(f"\n  span held at {args.spans[0]}, "
              f"{args.scan_pairs} nearest partners")
        print(f"\n  {'pair':>10} {'gap':>6} {'straddle':>9} "
              f"{'counts':>18} {'parity':>7}")
        for gap, b, q in cands[:args.scan_pairs]:
            # Where along the routed chord the partner sits. Near 0 means it
            # runs beside the chain; near 1 means it cuts clean across between
            # the two junctions.
            t = (q - A[0]) @ n_hat / span_len
            straddle = float(np.clip(t.max(), 0, 1) - np.clip(t.min(), 0, 1))
            got = []
            for sd in range(args.seeds):
                rng = np.random.default_rng(1000 + sd)
                ph = float(rng.uniform(0.0, 2.0 * np.pi))
                try:
                    pth = build(A, q, args.at, args.radius, args.spans[0],
                                ph, avoid, rng, args.dp, taut=args.taut)
                except ValueError:
                    continue
                work.mkdir(parents=True, exist_ok=True)
                tmp = work / "cand.data"
                rewrite_coords(relaxed, tmp,
                               {a: x for a, x in zip(seq[routed], pth)})
                pr = (min(routed, b), max(routed, b))
                got.append(measure_pairs(tmp, seq, [pr], work).get(pr))
            g = [x for x in got if x is not None]
            odd = sum(1 for x in g if x % 2)
            tag = ("mixed" if g and 0 < odd < len(g)
                   else "odd" if odd else "even" if g else "-")
            print(f"  {f'{routed}-{b}':>10} {gap:6.2f} {straddle:9.2f} "
                  f"{str(g):>18} {tag:>7}")
        return 0

    print(f"\n  {'span':>6} {'pts':>4} {'clear':>6} {'as built':>22} "
          f"{'after MD':>22}")
    rows = []
    reported = False
    for span in args.spans:
        built, survived, clears = [], [], []
        for s in range(args.seeds):
            rng = np.random.default_rng(1000 + s)
            phase = float(rng.uniform(0.0, 2.0 * np.pi))
            try:
                p = build(paths[routed], paths[target], args.at, args.radius,
                          span, phase, avoid, rng, args.dp, taut=args.taut)
            except ValueError as e:
                if not reported:
                    print(f"    span {span}: {e}")
                    reported = True
                continue
            clears.append(avoid.worst(p[1:-1]))
            tmp = work / "cand.data"
            work.mkdir(parents=True, exist_ok=True)
            rewrite_coords(relaxed, tmp,
                           {a: x for a, x in zip(seq[routed], p)})
            built.append(measure_pairs(tmp, seq, [pair], work).get(pair))
            survived.append(relax_and_count(tmp))
        def show(v):
            v = [x for x in v if x is not None]
            if not v:
                return "not measured"
            c = collections.Counter(v)
            return " ".join(f"{k}x{n}" for k, n in sorted(c.items()))
        print(f"  {span:6.2f} {ring_points(args.radius, span):4d} "
              f"{np.median(clears) if clears else 0:6.2f} "
              f"{show(built):>22} {show(survived):>22}")
        rows.append((span, built, survived))

    print("\n  a span is usable when one count appears every time in both "
          "columns:")
    for span, built, surv in rows:
        b = set(x for x in built if x is not None)
        s = set(x for x in surv if x is not None)
        if len(b) == 1 and b == s:
            print(f"    span {span:.2f} -> {b.pop()} every time")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
