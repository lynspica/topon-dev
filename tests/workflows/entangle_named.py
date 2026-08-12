"""Entangle named chains with named partners, at named counts.

    python tests/workflows/entangle_named.py --pairs 0:29:2 20:51:2
    python tests/workflows/entangle_named.py --pairs 0:29:1 --verify

Requirement 4, and the one place survival of an individual pair is the whole
point: nothing else can stand in for the pair that was asked for.

Two things make it hold rather than hope.

The path is drawn around the beads already there. A path that lands on top of
them starts inside the WCA hard core, and the minimisation that follows
resolves that by pushing chains through each other, which is the only thing
that can change a topology once it is built. Measured, routing blind took the
closest pair in the system from 0.502 sigma to 0.195 and put 153 beads inside
0.5 sigma.

And the count is verified rather than assumed. `--verify` relaxes each
candidate placement and keeps the first whose count survives, which is the
difference between a design that usually holds and one that holds by
construction: enumerated over 24 placements, 17 built the number asked for and
only 15 of those kept it, with neither clearance nor contact length predicting
which.
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
    ap.add_argument("--pairs", nargs="+", default=["0:29:2", "20:51:2"],
                    help="chain:partner:count, e.g. 0:29:2")
    ap.add_argument("--dims", type=int, default=4)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--ring", type=float, default=2.0)
    ap.add_argument("--margin", type=float, default=2.0)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--ranks", type=int, default=5)
    ap.add_argument("--phases", type=int, default=4)
    ap.add_argument("--passes", type=int, default=3)
    ap.add_argument("--clearance", type=float, default=0.9)
    ap.add_argument("--site-lo", type=float, default=0.30)
    ap.add_argument("--site-hi", type=float, default=0.70)
    ap.add_argument("--verify", action="store_true", default=True,
                    help="keep the placement whose count survives relaxation")
    ap.add_argument("--no-verify", dest="verify", action="store_false")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tag", default="_named")
    args = ap.parse_args()

    plan = []
    for spec in args.pairs:
        a, b, n = spec.split(":")
        plan.append((int(a), int(b), int(n)))
    if not plan:
        print("  nothing asked for")
        return 1

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (args.dims,) * 3
    graph = build_network(spec)

    sc = scale_for_design(graph, [(a, b) for a, b, _n in plan], dp=args.dp,
                          bond=BOND, radius=args.ring, margin=args.margin,
                          site_span=(args.site_lo, args.site_hi))
    geo = geometry(graph, dp=args.dp, scale=sc)
    keys = sorted(geo["chords"])
    for a, b, _n in plan:
        if a not in keys or b not in keys:
            print(f"  chain {a} or {b} is not in this network "
                  f"(0 to {len(keys) - 1})")
            return 1

    n_beads = graph.number_of_edges() * args.dp + graph.number_of_nodes()
    print("  asked: " + ", ".join(f"{a}-{b} at {n}" for a, b, n in plan))
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
    print("\n  --- relaxing a plain melt ---")
    run_md(sim, args.stages)
    after = {1: "system_after_soft.data", 2: "system_ramped.data",
             3: "system_equilibrated.data"}[args.stages]
    relaxed = root / "04_Simulation" / after
    if not relaxed.exists():
        print("  relaxation produced no output")
        return 1
    box, paths, xyz0 = paths_from(relaxed, keys, seq)

    work = OUT / f"named_work_{os.getpid()}"
    pairs = sorted({(min(a, b), max(a, b)) for a, b, _n in plan})
    base = measure_pairs(relaxed, seq, pairs, work)

    cand = OUT / f"named_cand_{os.getpid()}"
    (cand / "03_Conformation").mkdir(parents=True, exist_ok=True)
    for sub in ("02_Chemistry", "04_Simulation"):
        (cand / sub).mkdir(parents=True, exist_ok=True)
        for f in (root / sub).glob("*"):
            if f.is_file() and (sub == "02_Chemistry" or f.suffix == ".in"):
                (cand / sub / f.name).write_text(f.read_text())

    def relax(data_file):
        (cand / "03_Conformation" / "system_relaxed.data").write_text(
            Path(data_file).read_text())
        for stale in ("system_after_soft.data", "system_ramped.data",
                      "system_equilibrated.data"):
            q = cand / "04_Simulation" / stale
            if q.exists():
                q.unlink()
        try:
            run_md(cand / "04_Simulation", args.stages)
        except Exception:
            return None
        out = cand / "04_Simulation" / after
        return out if out.exists() else None

    current = root / "04_Simulation" / "current.data"
    xyz_now = dict(xyz0)
    rewrite_coords(relaxed, current, xyz_now)

    def place(a, b, n):
        mine = set(seq[a])
        avoid = Clearance(np.array([xyz_now[i] for i in sorted(xyz_now)
                                    if i not in mine]), box, args.clearance)
        try:
            p, got = construct_exact(
                paths, a, b, n, box, avoid, seq, current, work,
                (args.site_lo, args.site_hi), args.ring, args.dp,
                ranks=args.ranks, phases=args.phases,
                relax=(relax if args.verify else None))
        except ValueError as e:
            print(f"    {a}-{b}: {e}")
            return
        if p is None:
            print(f"    {a}-{b}: nothing built")
            return
        paths[a] = p
        for aid, xyzp in zip(seq[a], p):
            xyz_now[aid] = xyzp
        rewrite_coords(relaxed, current, xyz_now)
        print(f"    {a}-{b}: {got}" + ("" if got == n else f", asked {n}"))

    print("\n  --- placing ---")
    for a, b, n in plan:
        place(a, b, n)

    # Chains are placed one at a time against what is already committed, so
    # the first never sees the last. Re-place only those that drifted, and end
    # on a verification rather than a placement.
    settled = False
    for sweep in range(1, max(1, args.passes)):
        probe = relax(current)
        if probe is None:
            break
        now = measure_pairs(probe, seq, pairs, work)
        off = [(a, b, n) for a, b, n in plan
               if now.get((min(a, b), max(a, b))) != n]
        if not off:
            settled = True
            print(f"  every count holds after pass {sweep}")
            break
        print(f"  pass {sweep + 1}: re-placing "
              + ", ".join(f"{a}-{b}" for a, b, _n in off))
        for a, b, n in off:
            place(a, b, n)
    if not settled and args.passes > 1:
        print(f"  did not settle in {args.passes} passes")

    designed = root / "04_Simulation" / "designed.data"
    rewrite_coords(relaxed, designed, xyz_now)
    built = measure_pairs(designed, seq, pairs, work)

    root2 = OUT / f"relaxed{args.tag}_again"
    (root2 / "03_Conformation").mkdir(parents=True, exist_ok=True)
    for sub in ("02_Chemistry", "04_Simulation"):
        (root2 / sub).mkdir(parents=True, exist_ok=True)
        for f in (root / sub).glob("*"):
            if f.is_file() and (sub == "02_Chemistry" or f.suffix == ".in"):
                (root2 / sub / f.name).write_text(f.read_text())
    (root2 / "03_Conformation" / "system_relaxed.data").write_text(
        Path(designed).read_text())
    print("\n  --- minimising the finished system ---")
    run_md(root2 / "04_Simulation", args.stages)
    final = root2 / "04_Simulation" / after
    if not final.exists():
        print("  the finished system did not relax")
        return 1
    surv = measure_pairs(final, seq, pairs, work)

    print(f"\n  {'pair':>10} {'asked':>6} {'before':>7} {'as built':>9} "
          f"{'after MD':>9}")
    ok = 0
    for a, b, n in plan:
        q = (min(a, b), max(a, b))
        ok += (surv[q] == n)
        print(f"  {f'{a}-{b}':>10} {n:>6} {str(base[q]):>7} "
              f"{str(built[q]):>9} {str(surv[q]):>9}"
              + ("   ok" if surv[q] == n else ""))
    print(f"\n  {ok} of {len(plan)} carry exactly the count they were asked "
          f"for, after minimisation")
    print(f"  system: {final}")
    return 0 if ok == len(plan) else 1


if __name__ == "__main__":
    raise SystemExit(main())
