"""Does the new construction reproduce the old Gaussian kink, statistically?

    python tests/workflows/entangle_reproduce.py --per-chain 2.0

Requirement 3: asking for the same density with second neighbours only should
give, statistically, what the old kink gives on the same lattice.

Both systems are built here, from the same network, the same chains, the same
density and the same LAMMPS protocol, and both are measured the same way. That
matters more than it sounds: the old path's `avg_crosslinks_per_chain` is a
*request*, never checked against the built system, so "e = 2 in the old code"
and "2 delivered" are different quantities. Only the delivered ones can be
compared, and they are what this reports.

What is compared:

    entanglement points per chain     the density itself
    partners per chain                spread over neighbours, not doubled up
    shell of the delivered pairs      whether they land where the kink puts them
    points per entangled pair         how concentrated each contact is

Nothing here decides which method is better. It reports where the two agree and
where they do not.
"""
from __future__ import annotations

import argparse
import collections
import os
import random
import shutil
import sys
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.assignment.entanglements import (  # noqa: E402
    chain_distances,
    find_crossing_candidates,
    neighbour_shells,
    select_by_shells,
    select_entanglements,
)
from topon.config.schema import EntanglementsConfig  # noqa: E402
from topon.conformation.paths import Clearance, bridging_walk  # noqa: E402
from tests.workflows.entangle_all import CASES  # noqa: E402
from tests.workflows.entangle_density import (  # noqa: E402
    measure_system,
    watch_set,
)
from tests.workflows.entangle_relaxed import (  # noqa: E402
    construct,
    measure_pairs,
    paths_from,
    rewrite_coords,
)
from tests.workflows.entangle_spatial import kinked_paths  # noqa: E402
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
    write_system,
)


def relax_system(graph, geo, paths, root, stages, tag):
    """Build and minimise one system. Returns (final data file, seq)."""
    shutil.rmtree(root, ignore_errors=True)
    _n, node_atom, chain_atoms = write_system(graph, geo, paths, root)
    sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                             protocol="hardcore")
    print(f"  --- {tag} ---")
    run_md(sim, stages)
    out = root / "04_Simulation" / {
        1: "system_after_soft.data", 2: "system_ramped.data",
        3: "system_equilibrated.data"}[stages]
    keys = sorted(geo["chords"])
    seq = {k: chain_ids(k, node_atom, chain_atoms, geo["ends"]) for k in keys}
    return (out if out.exists() else None), seq


def stats(data_file, seq, keys, work, watch, shell_of):
    """The four numbers both methods are compared on."""
    from tests.workflows.entangle_relaxed import measure_pairs as mp

    per_chain, partners, _m, blind = measure_system(data_file, seq, keys,
                                                    work, watch)
    pairs, _scale = watch
    got = mp(data_file, seq, pairs, work)
    live = {q: v for q, v in got.items() if v}
    by_shell = collections.Counter()
    for q, v in live.items():
        s = shell_of.get(q)
        if s is not None:
            by_shell[s] += v
    conc = (sum(live.values()) / len(live)) if live else 0.0
    return {
        "per_chain": per_chain,
        "partners": partners,
        "per_pair": conc,
        "shells": by_shell,
        "blind": blind,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--per-chain", type=float, default=2.0,
                    help="the density both methods are asked for")
    ap.add_argument("--dims", type=int, default=4)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--density", type=float, default=0.85)
    ap.add_argument("--ring", type=float, default=2.0)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--clearance", type=float, default=0.9)
    ap.add_argument("--cutoff", type=float, default=3.0)
    ap.add_argument("--sample", type=int, default=250)
    ap.add_argument("--pair-yield", type=float, default=0.1,
                    help="entanglements per chain a designed pair delivers "
                         "here; entangle_density.py measures it, 0.095 to "
                         "0.147 on this system")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (args.dims,) * 3
    graph = build_network(spec)
    geo = geometry(graph, dp=args.dp, density=args.density)
    dims = np.asarray(graph.graph["box"], float)
    keys = sorted(geo["chords"])
    edges = sorted(graph.edges())
    idx_of = {frozenset(e): i for i, e in enumerate(edges)}

    G = nx.MultiGraph()
    G.add_nodes_from(graph.nodes(data=True))
    for u, v in graph.edges():
        G.add_edge(u, v)

    dist = chain_distances(G, dims)
    shells = neighbour_shells(G, dims, max_shell=4, distances=dist)
    shell_of = {}
    for chain, by in shells.items():
        a = idx_of.get(frozenset((chain[0], chain[1])))
        if a is None:
            continue
        for s, others in by.items():
            for o in others:
                b = idx_of.get(frozenset((o[0], o[1])))
                if b is not None and a != b:
                    shell_of.setdefault((min(a, b), max(a, b)), s)

    print(f"  SC {args.dims}^3, DP {args.dp}, density {args.density}, "
          f"{len(keys)} chains, asking {args.per_chain} per chain")

    # ---- the plain melt, which both are measured against -----------------
    rng0 = np.random.default_rng(args.seed)
    paths0 = {k: bridging_walk(c0, c1, args.dp + 1, BOND, rng0)
              for k, (c0, c1) in geo["chords"].items()}
    plain_f, seq = relax_system(graph, geo, paths0, OUT / "repro_plain",
                                args.stages, "plain melt, no entanglements")
    if plain_f is None:
        print("  the plain melt did not relax")
        return 1
    box, paths, xyz0 = paths_from(plain_f, keys, seq)

    work = OUT / f"repro_work_{os.getpid()}"
    rng = np.random.default_rng(args.seed)

    # ---- the old way ------------------------------------------------------
    cands = find_crossing_candidates(G, dims)
    random.seed(args.seed)
    cfg = EntanglementsConfig(enabled=True,
                              avg_crosslinks_per_chain=args.per_chain)
    sel = select_entanglements(G, cfg, dims, candidates=list(cands),
                               num_chains=G.number_of_edges())
    old_pairs = []
    for e1, e2, c in sel:
        a = idx_of.get(frozenset((e1[0], e1[1])))
        b = idx_of.get(frozenset((e2[0], e2[1])))
        if a is not None and b is not None and a != b:
            old_pairs.append((min(a, b), max(a, b)))
    kinked, _p, _s = kinked_paths(graph, geo, sel, dims, args.dp,
                                  {"overshoot": 0.2, "z_amp": 0.5,
                                   "sigma": 0.15})
    old_f, _sq = relax_system(graph, geo, kinked, OUT / "repro_old",
                              args.stages, "old Gaussian kink")

    # ---- the new way, second shell only -----------------------------------
    new_sel = select_by_shells(G, args.per_chain, {2: 1.0}, dims,
                               shells=shells, rng=rng)
    new_pairs = []
    for e1, e2, c in new_sel:
        a = idx_of.get(frozenset((e1[0], e1[1])))
        b = idx_of.get(frozenset((e2[0], e2[1])))
        if a is not None and b is not None and a != b:
            new_pairs.append((min(a, b), max(a, b), c))

    watch = watch_set(paths, box, old_pairs
                      + [(a, b) for a, b, _c in new_pairs],
                      args.cutoff, args.sample, rng)
    base = stats(plain_f, seq, keys, work, watch, shell_of)
    print(f"\n  the plain melt carries {base['per_chain']:.2f} per chain "
          f"over {base['partners']:.1f} partners")

    ids = sorted(xyz0)
    row = {a: i for i, a in enumerate(ids)}
    xyz_arr = np.array([xyz0[a] for a in ids])
    xyz_now = dict(xyz0)
    built = 0

    # Only as many pairs as the density needs.
    #
    # Routing every selected pair puts a hundred chains through a melt that
    # only has room for a few dozen, and they start landing on each other: the
    # pair energy came back at 1.8e11 and LAMMPS stopped. A designed pair is
    # worth about 0.1 entanglements per chain here, measured, so the number
    # that matches the request is what gets routed.
    n_route = max(1, int(round(args.per_chain / max(args.pair_yield, 1e-9))))
    if n_route < len(new_pairs):
        print(f"  routing {n_route} of {len(new_pairs)} selected pairs, at "
              f"{args.pair_yield:.3f} entanglements per chain per pair")
        new_pairs = new_pairs[:n_route]

    for a, b, count in new_pairs:
        keep = np.ones(len(ids), bool)
        keep[[row[i] for i in seq[a]]] = False
        avoid = Clearance(xyz_arr[keep], box, args.clearance)
        try:
            p = construct(paths, a, b, max(0.5, 0.5 * count), box, avoid,
                          radius=args.ring, dp=args.dp, span=(0.3, 0.7))
        except ValueError:
            continue
        paths[a] = p
        for aid, xyzp in zip(seq[a], p):
            xyz_now[aid] = xyzp
            xyz_arr[row[aid]] = xyzp
        built += 1

    new_root = OUT / "repro_new"
    (new_root / "04_Simulation").mkdir(parents=True, exist_ok=True)
    (new_root / "03_Conformation").mkdir(parents=True, exist_ok=True)
    for sub in ("02_Chemistry", "04_Simulation"):
        (new_root / sub).mkdir(parents=True, exist_ok=True)
        for f in (OUT / "repro_plain" / sub).glob("*"):
            if f.is_file() and (sub == "02_Chemistry" or f.suffix == ".in"):
                (new_root / sub / f.name).write_text(f.read_text())
    rewrite_coords(plain_f, new_root / "03_Conformation" /
                   "system_relaxed.data", xyz_now)
    print(f"  --- new construction, second shell, {built} pairs routed ---")
    try:
        run_md(new_root / "04_Simulation", args.stages)
    except Exception:
        pass
    new_f = new_root / "04_Simulation" / {
        1: "system_after_soft.data", 2: "system_ramped.data",
        3: "system_equilibrated.data"}[args.stages]

    # ---- compare ----------------------------------------------------------
    old_s = stats(old_f, seq, keys, work, watch, shell_of) if old_f else None
    new_s = (stats(new_f, seq, keys, work, watch, shell_of)
             if new_f.exists() else None)

    print(f"\n  {'':<26} {'plain melt':>12} {'old kink':>12} {'new':>12}")

    def row_of(label, key, fmt="{:.2f}"):
        def one(d):
            return fmt.format(d[key]) if d else "not built"
        print(f"  {label:<26} {fmt.format(base[key]):>12} "
              f"{one(old_s):>12} {one(new_s):>12}")

    row_of("points per chain", "per_chain")
    row_of("partners per chain", "partners")
    row_of("points per entangled pair", "per_pair")

    print(f"\n  added over the plain melt:")
    for name, d in (("old kink", old_s), ("new", new_s)):
        if d:
            print(f"    {name:<10} {d['per_chain'] - base['per_chain']:+.2f} "
                  f"per chain, asked {args.per_chain:.2f}")

    print(f"\n  where the entanglements sit, by shell:")
    print(f"  {'shell':>6} {'plain melt':>12} {'old kink':>12} {'new':>12}")
    allsh = sorted(set(base["shells"])
                   | set(old_s["shells"] if old_s else ())
                   | set(new_s["shells"] if new_s else ()))
    for s in allsh:
        def frac(d):
            if not d:
                return "-"
            t = sum(d["shells"].values())
            return f"{d['shells'].get(s, 0) / t:.2f}" if t else "0.00"
        print(f"  {s:>6} {frac(base):>12} {frac(old_s):>12} "
              f"{frac(new_s):>12}")
    print("\n  the same network, chains, density and protocol throughout, "
          "and both measured the same way")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
