"""Design entanglements into an already-relaxed melt, so they survive.

    python tests/workflows/entangle_relaxed.py
    python tests/workflows/entangle_relaxed.py --rounds 8 --per-round 24

Building a designed entanglement into a fresh random-walk conformation and
then minimising destroys about half of it. The reason is not subtle: the
beads move 3.34 sigma on average during minimisation, up to 12.45, while the
designed loop is only 1.2 to 3.4 sigma across. The feature is the same size
as the motion, so no amount of care in building it helps.

That motion is not the entanglement's fault. A random walk drawn without
excluded volume overlaps itself and its neighbours everywhere, and the first
minimisation stage is mostly relieving that. So relax first, and design into
the relaxed melt:

    1. build a plain melt, no designed entanglements
    2. minimise it -- this absorbs all the overlap relief
    3. route chains through the relaxed coordinates
    4. minimise again, which now has little left to do
    5. measure each requested pair on its own

Step 4 is the one that matters. A system already at equilibrium has no reason
to move far, so the designed feature is no longer competing with a global
rearrangement.
"""
from __future__ import annotations

import argparse
import collections
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.paths import bridging_walk  # noqa: E402
from tests.workflows.entangle_all import CASES  # noqa: E402
from tests.workflows.entangle_design import route_one  # noqa: E402
from tests.workflows.entangle_search import (  # noqa: E402
    Wish,
    _both,
    measure_batch,
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
    report_bonds,
    run_md,
    unwrap_chain,
    write_system,
    z1_export,
)


def rewrite_coords(src, dst, new_xyz):
    """Copy a LAMMPS data file, replacing the coordinates of some atoms.

    Everything else -- header, masses, bonds, angles -- is carried through
    untouched, so the result is the same system in a different conformation.
    Writing a fresh file through the chemistry and conformation stages would
    also work but would rebuild the topology, and the point here is to change
    only where the atoms are.
    """
    out, section = [], None
    for line in Path(src).read_text().splitlines(True):
        s = line.strip()
        if s and s.split()[0] in ("Atoms", "Bonds", "Masses", "Velocities",
                                  "Angles"):
            section = s.split()[0]
            out.append(line)
            continue
        if section == "Atoms" and s and not s.startswith("#"):
            p = s.split()
            if len(p) >= 7 and int(p[0]) in new_xyz:
                x, y, z = new_xyz[int(p[0])]
                p[4], p[5], p[6] = f"{x:.6f}", f"{y:.6f}", f"{z:.6f}"
                out.append(" ".join(p) + "\n")
                continue
        out.append(line)
    Path(dst).write_text("".join(out))
    return dst


def paths_from(data_file, keys, seq):
    box, xyz, _ = read_data(data_file)
    return box, {k: unwrap_chain(seq[k], xyz, box) for k in keys}, xyz


def measure_pairs(data_file, seq, pairs, work):
    """Each pair on its own -- the only export Z1+ takes reliably."""
    work.mkdir(parents=True, exist_ok=True)
    for old in list(work.glob("*.Z1")) + list(work.glob("SP_*.dat")):
        old.unlink()
    for a, b in pairs:
        z1_export(data_file, [seq[a], seq[b]], work / f"p{a}_{b}.Z1")
    res = measure_batch(work) or {}
    return {(a, b): (_both(res[f"p{a}_{b}"], 1, 2)
                     if f"p{a}_{b}" in res else None) for a, b in pairs}


def route_pairwise(paths, keys, seq, geo, routed, target, want, template,
                   rounds, per_round, rng, work, relax=None):
    """Search for a path, scoring each candidate on the pair alone.

    The whole-system score used elsewhere is unavailable here: Z1+ crashes on
    a relaxed melt read back from LAMMPS output, so every candidate comes back
    empty and the search has no signal. Measuring just the two chains works
    reliably and scores the thing being asked for.

    The cost is collateral blindness -- what else the routed chain picks up is
    not seen, so it cannot be traded against. The primary objective is served;
    the secondary one is not.
    """
    from tests.workflows.entangle_search import propose

    L = geo["L"]
    a0, a1 = paths[routed][0], paths[routed][-1]
    best = None
    around, spread = None, 1.0

    for _ in range(rounds):
        tgt = paths[target]
        tgt = tgt + L * np.round((paths[routed].mean(0) - tgt.mean(0)) / L)
        cands = propose(a0, a1, [tgt], L, per_round, rng, around, spread)
        if not cands:
            continue

        work.mkdir(parents=True, exist_ok=True)
        for old in list(work.glob("*.Z1")) + list(work.glob("SP_*.dat")):
            old.unlink()
        tmp = work / "cand.data"
        for i, (_knobs, path) in enumerate(cands):
            new = {aid: xyzp for aid, xyzp in zip(seq[routed], path)}
            rewrite_coords(template, tmp, new)
            # Score after a short minimisation, not as built.
            #
            # A loop placed into a relaxed melt is a local strain, and
            # relaxing it is what opens or tightens the winding. Measured on
            # candidates chosen as-built: a pair designed at 1 came back at 0
            # and one designed at 2 came back at 5, with beads moving only
            # 0.57 sigma. Nothing about the built geometry predicts that, so
            # the only way to select for it is to relax each candidate and
            # look. In an already-relaxed system stage 1 costs a fraction of
            # a second, which is what makes this affordable at all.
            scored_file = tmp
            if relax is not None:
                r = relax(tmp, i)
                if r is not None:
                    scored_file = r
            z1_export(scored_file, [seq[routed], seq[target]],
                      work / f"c{i:03d}.Z1")
        res = measure_batch(work) or {}

        scored = []
        for i, (knobs, path) in enumerate(cands):
            r = res.get(f"c{i:03d}")
            if r is None:
                continue
            got = _both(r, 1, 2)
            scored.append((abs(got - want), got, knobs, path))
        if not scored:
            continue
        scored.sort(key=lambda t: t[0])
        if best is None or scored[0][0] < best[0]:
            best = scored[0]
            around, spread = best[2], max(0.35, spread * 0.6)
        else:
            spread = min(1.5, spread * 1.4)
        if best[0] == 0:
            break
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rounds", type=int, default=6)
    ap.add_argument("--per-round", type=int, default=20)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (4, 4, 4)
    graph = build_network(spec)
    geo = geometry(graph, dp=args.dp, density=0.85)
    ch, ends = geo["chords"], geo["ends"]
    keys = sorted(ch)
    idx = {k: i + 1 for i, k in enumerate(keys)}

    rng = np.random.default_rng(args.seed)
    paths0 = {k: bridging_walk(c0, c1, args.dp + 1, BOND, rng)
              for k, (c0, c1) in ch.items()}

    # ---- 1 and 2: a plain melt, relaxed ---------------------------------
    root = OUT / "relaxed"
    _n, node_atom, chain_atoms = write_system(graph, geo, paths0, root)
    seq = {k: chain_ids(k, node_atom, chain_atoms, ends) for k in keys}
    sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                             protocol="hardcore")
    print("  --- relaxing a plain melt, no designed entanglements ---")
    run_md(sim, args.stages)
    after = {1: "system_after_soft.data", 2: "system_ramped.data",
             3: "system_equilibrated.data"}[args.stages]
    relaxed = root / "04_Simulation" / after
    if not relaxed.exists():
        print("  relaxation produced no output")
        return 1

    box, paths, xyz0 = paths_from(relaxed, keys, seq)
    geo_r = dict(geo)
    geo_r["L"] = box
    geo_r["chords"] = {k: (paths[k][0], paths[k][-1]) for k in keys}
    print(f"  relaxed: box {box[0]:.1f} sigma")

    # ---- 3: route into the relaxed coordinates --------------------------
    def far_from(a, rank):
        m = 0.5 * (ch[a][0] + ch[a][1])
        order = sorted(
            (float(np.linalg.norm(
                (lambda v: v - geo["L"] * np.round(v / geo["L"]))(
                    0.5 * (ch[b][0] + ch[b][1]) - m))), b)
            for b in keys if b != a and not set(ends[a]) & set(ends[b]))
        return order[rank][1]

    A, D = keys[0], keys[20]
    plan = [(A, {far_from(A, len(keys) // 3): 1}),
            (D, {far_from(D, len(keys) // 4): 1})]
    pairs = sorted({(min(r, t), max(r, t)) for r, w in plan for t in w})
    print(f"  plan: " + "; ".join(f"{r} with {list(w)[0]}" for r, w in plan))

    work = OUT / f"relaxed_work_{os.getpid()}"
    base = measure_pairs(relaxed, seq, pairs, work)
    print(f"  before routing: "
          + ", ".join(f"{a}-{b}={v}" for (a, b), v in base.items()))

    # A one-stage relaxation of a candidate, reusing the scripts already
    # generated for this system.
    cand_root = OUT / f"relaxed_cand_{os.getpid()}"
    (cand_root / "03_Conformation").mkdir(parents=True, exist_ok=True)
    (cand_root / "04_Simulation").mkdir(parents=True, exist_ok=True)
    (cand_root / "02_Chemistry").mkdir(parents=True, exist_ok=True)
    for f in (root / "02_Chemistry").glob("*"):
        (cand_root / "02_Chemistry" / f.name).write_text(f.read_text())
    for f in (root / "04_Simulation").glob("*.in"):
        (cand_root / "04_Simulation" / f.name).write_text(f.read_text())

    def relax_candidate(data_file, i):
        (cand_root / "03_Conformation" / "system_relaxed.data").write_text(
            Path(data_file).read_text())
        for stale in ("system_after_soft.data",):
            q = cand_root / "04_Simulation" / stale
            if q.exists():
                q.unlink()
        try:
            run_md(cand_root / "04_Simulation", 1)
        except Exception:
            return None
        out = cand_root / "04_Simulation" / "system_after_soft.data"
        return out if out.exists() else None

    for routed, want in plan:
        target, n = list(want.items())[0]
        best = route_pairwise(paths, keys, seq, geo_r, routed, target, n,
                              relaxed, args.rounds, args.per_round, rng, work,
                              relax=relax_candidate)
        if best is None:
            print(f"    chain {routed}: nothing built")
            continue
        paths[routed] = best[3]
        print(f"    chain {routed} -> {target}: got {best[1]}, wanted {n}")

    # ---- 4: write it back and minimise again ----------------------------
    new_xyz = {}
    for k in keys:
        for aid, xyzp in zip(seq[k], paths[k]):
            new_xyz[aid] = xyzp
    designed = root / "04_Simulation" / "designed.data"
    rewrite_coords(relaxed, designed, new_xyz)

    built = measure_pairs(designed, seq, pairs, work)
    print(f"\n  as designed:    "
          + ", ".join(f"{a}-{b}={v}" for (a, b), v in built.items()))

    root2 = OUT / "relaxed_again"
    (root2 / "04_Simulation").mkdir(parents=True, exist_ok=True)
    for f in (root / "04_Simulation").glob("*.in"):
        (root2 / "04_Simulation" / f.name).write_text(f.read_text())
    (root2 / "03_Conformation").mkdir(parents=True, exist_ok=True)
    (root2 / "03_Conformation" / "system_relaxed.data").write_text(
        Path(designed).read_text())
    for sub in ("02_Chemistry",):
        (root2 / sub).mkdir(parents=True, exist_ok=True)
        for f in (root / sub).glob("*"):
            (root2 / sub / f.name).write_text(f.read_text())

    print(f"\n  --- minimising again, from the designed relaxed system ---")
    run_md(root2 / "04_Simulation", args.stages)
    final = root2 / "04_Simulation" / after
    if not final.exists():
        print("  second minimisation produced no output")
        return 1
    print()
    report_bonds(root2)

    _b, xyz2, _ = read_data(final)
    ids = sorted(xyz0)
    d = np.array([xyz2[i] - new_xyz[i] for i in ids])
    d -= box * np.round(d / box)
    dist = np.linalg.norm(d, axis=1)
    print(f"\n  bead motion in this second minimisation: median "
          f"{np.median(dist):.2f}, max {dist.max():.2f} sigma")

    surv = measure_pairs(final, seq, pairs, work)
    print(f"\n  {'pair':>10} {'before':>7} {'designed':>9} {'after MD':>9}")
    ok = 0
    for a, b in pairs:
        want = next(n for r, w in plan for t, n in w.items()
                    if {r, t} == {a, b})
        ok += (surv[(a, b)] == want)
        print(f"  {f'{a}-{b}':>10} {str(base[(a, b)]):>7} "
              f"{str(built[(a, b)]):>9} {str(surv[(a, b)]):>9}"
              + ("   ok" if surv[(a, b)] == want else ""))
    print(f"\n  {ok} of {len(pairs)} survived exactly")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
