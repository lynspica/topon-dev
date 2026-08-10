"""Build a named entanglement topology by search, several chains at once.

    python tests/workflows/entangle_design.py
    python tests/workflows/entangle_design.py --rounds 8 --per-round 24

Takes a plan -- "chain A entangled with B and C, chain D entangled with B
twice and E once" -- and delivers it, one routed chain at a time, each
verified by primitive-path analysis before the next is attempted.

Chains are routed in sequence rather than together, and each search runs on
the conformation left by the ones before it. That ordering matters: routing
a chain changes what its neighbours are threaded by, so a plan settled in
parallel would not survive being assembled. Doing it in sequence means every
measurement is of the system as it will actually be.

The final report is the whole plan measured at once: every requested pair,
plus everything that appeared without being asked for.
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
from tests.workflows.entangle_all import CASES  # noqa: E402
from tests.workflows.entangle_search import (  # noqa: E402
    Wish,
    _both,
    measure_batch,
    propose,
    cost,
)
from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    DP,
    LATTICE,
    OUT,
    build_network,
    conform_and_script,
    geometry,
    report_bonds,
    run_md,
    write_system,
    write_z1,
)


def measure_one(paths, keys, geo, work, tag="probe"):
    """Z1+ on a single configuration. Returns {chain: Counter(partner)}."""
    L = geo["L"]
    work.mkdir(parents=True, exist_ok=True)
    for old in list(work.glob("*.Z1")) + list(work.glob("SP_*.dat")):
        old.unlink()
    arr = [paths[k] for k in keys]
    ref = arr[0].mean(0)
    write_z1(work / f"{tag}.Z1",
             [p + L * np.round((ref - p.mean(0)) / L) for p in arr], L)
    return (measure_batch(work) or {}).get(tag, {})


def route_one(paths, keys, geo, routed, wish, rounds, per_round, rng, work):
    """Search for a path for ``routed`` that satisfies ``wish``.

    Returns the winning path and what it measured, or None.
    """
    L = geo["L"]
    a0, a1 = geo["chords"][routed]
    targets = sorted(wish.want)
    best = None
    around, spread = None, 1.0

    for _ in range(rounds):
        # Every requested partner in one path. A chain wanting two of them
        # needs to loop around both; taking the best of attempts aimed at one
        # at a time delivers whichever was aimed at and never the pair.
        tgts = []
        for t in targets:
            q = paths[keys[t - 1]]
            tgts.append(q + L * np.round((paths[routed].mean(0)
                                          - q.mean(0)) / L))
        cands = propose(a0, a1, tgts, L, per_round, rng, around, spread)
        if not cands:
            continue

        work.mkdir(parents=True, exist_ok=True)
        for old in list(work.glob("*.Z1")) + list(work.glob("SP_*.dat")):
            old.unlink()
        for i, (_knobs, path) in enumerate(cands):
            trial = dict(paths)
            trial[routed] = path
            arr = [trial[k] for k in keys]
            ref = arr[0].mean(0)
            write_z1(work / f"c{i:03d}.Z1",
                     [p + L * np.round((ref - p.mean(0)) / L) for p in arr], L)

        res = measure_batch(work)
        if not res:
            return None

        scored = []
        for i, (knobs, path) in enumerate(cands):
            p = res.get(f"c{i:03d}")
            if p is None:
                continue
            c, miss, added = cost(p, wish)
            scored.append((c, miss, added, knobs, path, p))
        if not scored:
            continue
        scored.sort(key=lambda t: t[0])
        if best is None or scored[0][0] < best[0]:
            best = scored[0]
            around, spread = best[3], max(0.35, spread * 0.6)
        else:
            spread = min(1.5, spread * 1.4)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rounds", type=int, default=6)
    ap.add_argument("--per-round", type=int, default=20)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--n-routed", type=int, default=2,
                    help="how many chains carry a designed pair; above 2 a "
                         "wider plan is generated")
    ap.add_argument("--run-md", action="store_true",
                    help="write the system and run minimize 1 and 2")
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--passes", type=int, default=3,
                    help="sweeps over the plan; a chain routed later can "
                         "knock an earlier pair off target")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (4, 4, 4)
    graph = build_network(spec)
    geo = geometry(graph, dp=args.dp, density=0.85)
    ch, ends, L = geo["chords"], geo["ends"], geo["L"]
    keys = sorted(ch)
    idx = {k: i + 1 for i, k in enumerate(keys)}

    rng = np.random.default_rng(args.seed)
    paths = {k: bridging_walk(c0, c1, args.dp + 1, BOND, rng)
             for k, (c0, c1) in ch.items()}

    # The plan, in the shape the original request described: one chain with
    # two partners, another with a repeat and a single.
    def far_from(a, rank):
        m = 0.5 * (ch[a][0] + ch[a][1])
        order = sorted(
            (float(np.linalg.norm(
                (lambda v: v - L * np.round(v / L))(
                    0.5 * (ch[b][0] + ch[b][1]) - m))), b)
            for b in keys if b != a and not set(ends[a]) & set(ends[b]))
        return order[rank][1]

    if args.n_routed <= 2:
        A = keys[0]
        B, C = far_from(A, len(keys) // 3), far_from(A, len(keys) // 5)
        D = keys[20]
        E = far_from(D, len(keys) // 4)
        plan = [(A, {B: 1, C: 1}), (D, {B: 2, E: 1})]
        print(f"  plan: chain {A} with {B} and {C}; "
              f"chain {D} with {B} twice and {E} once")
    else:
        # A wider plan, for finding where the collateral from one routed
        # chain starts undoing another. Routed chains are spread through the
        # network and each is given one partner a third of the way down its
        # own distance list.
        step = max(1, len(keys) // args.n_routed)
        plan = []
        for i in range(args.n_routed):
            a = keys[(i * step) % len(keys)]
            if any(a == r for r, _ in plan):
                continue
            plan.append((a, {far_from(a, len(keys) // 3): 1}))
        print(f"  plan: {len(plan)} chains, one partner each, "
              f"{sum(len(w) for _, w in plan)} requested pairs")
    print(f"  {len(keys)} chains, DP {args.dp}, melt density")
    print(f"  {args.rounds} rounds x {args.per_round} candidates per chain\n")

    work = OUT / "design_work"
    base = measure_one(paths, keys, geo, work, "base")

    # Passes, not one sweep. Routing a chain changes what its neighbours are
    # threaded by, so a pair settled early can be knocked off by a chain
    # routed later: measured, chain 0 was routed with both its targets exact
    # and then routing chain 20, which shares one of those targets, took the
    # 0-15 count from 1 to 2. Re-routing whatever has drifted, against the
    # conformation as it now stands, is what closes that.
    for sweep in range(args.passes):
        state = base if sweep == 0 else measure_one(paths, keys, geo, work,
                                                    f"s{sweep}")
        todo = [(r, w) for r, w in plan
                if any(_both(state, idx[r], idx[t]) != n for t, n in w.items())]
        if not todo:
            print(f"  pass {sweep + 1}: everything already on target")
            break
        print(f"  pass {sweep + 1}: {len(todo)} chain(s) to route")
        for routed, want in todo:
            w = {idx[t]: n for t, n in want.items()}
            got0 = state.get(idx[routed], collections.Counter())
            baseline = sum(v for q, v in got0.items() if q not in w)
            wish = Wish(chain=idx[routed], want=w, penalty=1.0,
                        baseline=baseline)
            best = route_one(paths, keys, geo, routed, wish,
                             args.rounds, args.per_round, rng, work)
            if best is None:
                print(f"    chain {routed}: nothing built")
                continue
            c, miss, added, knobs, path, p = best
            paths[routed] = path
            detail = ", ".join(f"{t}:{_both(p, idx[routed], idx[t])}/{n}"
                               for t, n in want.items())
            print(f"    chain {routed} -- {detail}, {added} added elsewhere")

    final = measure_one(paths, keys, geo, work, "final")
    print("\n  " + "-" * 58)
    print("  the finished system, every requested pair measured at once")
    print("  " + "-" * 58)
    hit = 0
    for routed, want in plan:
        for t, n in want.items():
            got = _both(final, idx[routed], idx[t])
            ok = got == n
            hit += ok
            print(f"  {routed}-{t}: asked {n}, got {got}"
                  + ("   ok" if ok else ""))
    total = sum(len(w) for _, w in plan)
    print("  " + "-" * 58)
    print(f"  {hit} of {total} requested pairs delivered exactly")

    routed_ids = {idx[r] for r, _ in plan}
    asked = {(idx[r], idx[t]) for r, w in plan for t in w}
    extra = 0
    for a in routed_ids:
        for b, v in final.get(a, {}).items():
            if (a, b) not in asked and (b, a) not in asked:
                extra += v
    base_extra = 0
    for a in routed_ids:
        for b, v in base.get(a, {}).items():
            if (a, b) not in asked and (b, a) not in asked:
                base_extra += v
    print(f"  entanglements on the routed chains that nobody asked for: "
          f"{extra}, against {base_extra} before routing")
    bonds = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                            for p in paths.values()])
    print(f"  bonds {bonds.min():.3f} to {bonds.max():.3f}")

    if args.run_md:
        root = OUT / "design_md"
        n_atoms, node_atom, chain_atoms = write_system(graph, geo, paths, root)
        sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                                 protocol="hardcore")
        print(f"\n  wrote {n_atoms} beads to {root.name}/")
        print(f"  --- LAMMPS, stages 1 to {args.stages} ---")
        run_md(sim, args.stages)
        print()
        report_bonds(root)

        # Does the designed topology survive the protocol? Measure the same
        # pairs on the minimised coordinates, not on what was built.
        from tests.workflows.entangle_steps import chain_ids, read_data, unwrap_chain
        after = {1: "system_after_soft.data", 2: "system_ramped.data",
                 3: "system_equilibrated.data"}[args.stages]
        f = root / "04_Simulation" / after
        if f.exists():
            box, xyz, _ = read_data(f)
            md = {}
            for k in keys:
                seq = chain_ids(k, node_atom, chain_atoms, geo["ends"])
                md[k] = unwrap_chain(seq, xyz, box)
            geo_md = dict(geo)
            geo_md["L"] = box
            post = measure_one(md, keys, geo_md, work, "post")
            print("\n  after minimisation:")
            ok = 0
            for routed, want in plan:
                for t, n in want.items():
                    got = _both(post, idx[routed], idx[t])
                    ok += (got == n)
                    print(f"    {routed}-{t}: asked {n}, got {got}"
                          + ("   ok" if got == n else ""))
            print(f"    {ok} of {total} still exact after stages 1-"
                  f"{args.stages}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
