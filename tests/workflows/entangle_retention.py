"""Does clearance predict whether a designed entanglement is kept?

    python tests/workflows/entangle_retention.py
    python tests/workflows/entangle_retention.py --ranks 6 --phases 6

Ordering exact placements by clearance lost both designed pairs, while the
rule that ignores it kept both, one of them sitting at a 0.14 sigma contact.
That is an overlap by the same standard that explains why designs survive at
all, so it reads as a paradox.

It was also measured on two placements, which is not enough to tell a
mechanism from a coincidence. Choosing by clearance also chose a different
site and a different phase, so the site may be what mattered and the clearance
may have come along for the ride.

This enumerates placements for one pair with everything else held fixed and
records, for each, what it built, what survived, how much room it had, and how
much of the routed chain lies against the target. If clearance drives
retention the survival column sorts by it. If it does not, the earlier reading
was a confound and something else is doing the work.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.paths import Clearance, bridging_walk  # noqa: E402
from tests.workflows.entangle_all import CASES  # noqa: E402
from tests.workflows.entangle_relaxed import (  # noqa: E402
    construct,
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


def grip(path, target, box, cutoff=1.5):
    """How many beads of the routed chain lie against the target strand.

    A proxy for how committed a winding is, as opposed to how much room it
    has. The two are not the same thing and may not even point the same way.
    """
    d = path[:, None, :] - target[None, :, :]
    d -= box * np.round(d / box)
    return int((np.linalg.norm(d, axis=2).min(axis=1) < cutoff).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ranks", type=int, default=6)
    ap.add_argument("--phases", type=int, default=4)
    ap.add_argument("--want", type=int, default=2)
    ap.add_argument("--ring", type=float, default=2.0)
    ap.add_argument("--margin", type=float, default=2.0)
    ap.add_argument("--dims", type=int, default=4)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--site-lo", type=float, default=0.30)
    ap.add_argument("--site-hi", type=float, default=0.70)
    ap.add_argument("--clearance", type=float, default=0.9)
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (args.dims,) * 3
    graph = build_network(spec)

    box_g = np.asarray(graph.graph["box"], float)
    raw = {n: np.asarray(d["pos"], float) for n, d in graph.nodes(data=True)}
    edges = sorted(graph.edges())

    def chord_pts(k, span=(args.site_lo, args.site_hi)):
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
                          radius=args.ring, margin=args.margin,
                          site_span=(args.site_lo, args.site_hi))
    geo = geometry(graph, dp=args.dp, scale=sc)
    keys = sorted(geo["chords"])
    n_beads = graph.number_of_edges() * args.dp + graph.number_of_nodes()
    print(f"  chain {routed} around chain {target}, box {geo['L'][0]:.1f} "
          f"sigma, density {n_beads / float(np.prod(geo['L'])):.3f}")

    rng0 = np.random.default_rng(42)
    paths0 = {k: bridging_walk(c0, c1, args.dp + 1, BOND, rng0)
              for k, (c0, c1) in geo["chords"].items()}
    root = OUT / "retention"
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

    mine = set(seq[routed])
    avoid = Clearance(np.array([xyz0[i] for i in sorted(xyz0)
                                if i not in mine]), box, args.clearance)
    tgt = paths[target] + box * np.round(
        (paths[routed].mean(0) - paths[target].mean(0)) / box)

    pair = (min(routed, target), max(routed, target))
    work = OUT / f"retention_work_{os.getpid()}"
    cand = OUT / f"retention_cand_{os.getpid()}"
    (cand / "03_Conformation").mkdir(parents=True, exist_ok=True)
    for sub in ("02_Chemistry", "04_Simulation"):
        (cand / sub).mkdir(parents=True, exist_ok=True)
        for f in (root / sub).glob("*"):
            if f.is_file() and (sub == "02_Chemistry" or f.suffix == ".in"):
                (cand / sub / f.name).write_text(f.read_text())

    print(f"\n  {'rank':>4} {'phase':>6} {'clear':>6} {'grip':>5} "
          f"{'built':>6} {'after':>6} {'moved':>6}  kept")
    rows = []
    for rank in range(args.ranks):
        for k in range(args.phases):
            phase = 2.0 * np.pi * k / args.phases
            try:
                p = construct(paths, routed, target, args.want, box, avoid,
                              radius=args.ring, phase=phase, dp=args.dp,
                              span=(args.site_lo, args.site_hi), rank=rank)
            except ValueError:
                continue
            work.mkdir(parents=True, exist_ok=True)
            trial = dict(zip(seq[routed], p))
            rewrite_coords(relaxed, work / "try.data", trial)
            built = measure_pairs(work / "try.data", seq, [pair], work)[pair]
            if built is None:
                continue

            (cand / "03_Conformation" / "system_relaxed.data").write_text(
                (work / "try.data").read_text())
            try:
                run_md(cand / "04_Simulation", args.stages)
            except Exception:
                continue
            out = cand / "04_Simulation" / after
            if not out.exists():
                continue
            surv = measure_pairs(out, seq, [pair], work)[pair]

            _b, xyz2, _m = read_data(out)
            d = np.array([xyz2[a] - trial[a] for a in seq[routed]])
            d -= box * np.round(d / box)
            moved = float(np.median(np.linalg.norm(d, axis=1)))

            room = avoid.worst(p[1:-1])
            g = grip(p, tgt, box)
            kept = (surv == built and built == args.want)
            rows.append((room, g, built, surv, moved, kept))
            print(f"  {rank:4d} {phase:6.2f} {room:6.2f} {g:5d} "
                  f"{str(built):>6} {str(surv):>6} {moved:6.2f}"
                  f"  {'yes' if kept else 'no'}")

    if not rows:
        print("\n  nothing built")
        return 1

    keep = [r for r in rows if r[5]]
    lose = [r for r in rows if not r[5]]
    print(f"\n  {len(keep)} of {len(rows)} placements gave {args.want} "
          f"and kept it")
    if keep and lose:
        print(f"\n  {'':>18} {'kept':>8} {'lost':>8}")
        for i, name in ((0, "clearance"), (1, "grip"), (4, "moved")):
            print(f"  {name:>18} {np.median([r[i] for r in keep]):8.2f} "
                  f"{np.median([r[i] for r in lose]):8.2f}")
        cl = np.array([r[0] for r in rows])
        kp = np.array([1.0 if r[5] else 0.0 for r in rows])
        if cl.std() > 1e-9 and kp.std() > 1e-9:
            print(f"\n  clearance vs kept, correlation: "
                  f"{np.corrcoef(cl, kp)[0, 1]:+.2f}")
            gr = np.array([float(r[1]) for r in rows])
            if gr.std() > 1e-9:
                print(f"  grip      vs kept, correlation: "
                      f"{np.corrcoef(gr, kp)[0, 1]:+.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
