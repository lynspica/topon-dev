"""One pair, both methods: do they give the same thing?

    python tests/workflows/entangle_onepair.py --count 3

The system-wide comparison in `entangle_reproduce.py` asked a different
question and its answer was easy to misread. Building a hundred kinks changes
the whole conformation, and most of the difference it reports is that
second-order effect on the melt, not a statement about what one designed
entanglement is.

This asks the narrow question instead. Take two parallel chains that are
neighbours in an SC lattice, ask each method for the same count on that one
pair, change nothing else, and measure that pair. If the two constructions mean
the same thing by "entangle these two chains n times", this is where it shows.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import networkx as nx
import numpy as np

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


def axis_of(graph, e):
    """Which lattice direction a chain runs along, or None if it is diagonal.

    Node identifiers here are integers, not coordinates, so the direction has
    to come from the stored positions.
    """
    u, v = e
    pu = np.asarray(graph.nodes[u]["pos"], float)
    pv = np.asarray(graph.nodes[v]["pos"], float)
    box = np.asarray(graph.graph["box"], float)
    d = pv - pu
    d -= box * np.round(d / box)
    nz = np.flatnonzero(np.abs(d) > 1e-9)
    return int(nz[0]) if len(nz) == 1 else None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--count", type=int, default=3,
                    help="entanglements asked for on the one pair")
    ap.add_argument("--shell", type=int, default=2,
                    help="which neighbour shell the pair is drawn from")
    ap.add_argument("--dims", type=int, default=4)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--density", type=float, default=0.85)
    ap.add_argument("--coil", type=float, default=None,
                    help="contour over chord, the scale entangle_all.py uses "
                         "for the legacy work (1.8). Overrides --density.")
    ap.add_argument("--ring", type=float, default=2.0)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--ranks", type=int, default=5)
    ap.add_argument("--phases", type=int, default=4)
    ap.add_argument("--clearance", type=float, default=0.9)
    ap.add_argument("--old-protocol", default="generated",
                    choices=("generated", "hardcore"),
                    help="which LAMMPS protocol the old kink is built under. "
                         "Its own is 'generated', whose soft push expands the "
                         "collapsed chords it starts from.")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (args.dims,) * 3
    graph = build_network(spec)
    # The scale the legacy path was built for, not one picked here.
    #
    # At density 0.85 a DP-80 chain laid on its chord has beads 0.067 sigma
    # apart against a 0.95 sigma bond, so the soft push has to expand every
    # chain fourteenfold and no designed bulge survives that. The coil route
    # gives 0.528 sigma instead. Measuring the old kink at 0.85 was measuring
    # it somewhere it was never meant to run.
    geo = (geometry(graph, dp=args.dp, bond=BOND, coil=args.coil)
           if args.coil else
           geometry(graph, dp=args.dp, density=args.density))
    dims = np.asarray(graph.graph["box"], float)
    keys = sorted(geo["chords"])
    edges = sorted(graph.edges())
    idx_of = {frozenset(e): i for i, e in enumerate(edges)}

    G = nx.MultiGraph()
    G.add_nodes_from(graph.nodes(data=True))
    for u, v in graph.edges():
        G.add_edge(u, v)
    dist = chain_distances(G, dims)
    shells = neighbour_shells(G, dims, max_shell=max(2, args.shell),
                              distances=dist)

    # A parallel pair in the requested shell: two chains running along the same
    # lattice direction, which is the geometry the old kink was built for.
    pick = None
    for chain, by in sorted(shells.items()):
        a = idx_of.get(frozenset((chain[0], chain[1])))
        if a is None or axis_of(graph, edges[a]) is None:
            continue
        for other in by.get(args.shell, ()):
            b = idx_of.get(frozenset((other[0], other[1])))
            if b is None or b == a:
                continue
            if axis_of(graph, edges[a]) == axis_of(graph, edges[b]):
                pick = (min(a, b), max(a, b))
                break
        if pick:
            break
    if pick is None:
        print(f"  no parallel pair found in shell {args.shell}")
        return 1
    a, b = pick
    # `dist` is keyed by the MultiGraph's (u, v, key) triples, not by the
    # 2-tuples `graph.edges()` yields, so the gap has to be looked up through
    # the index map rather than with the edge itself.
    gap_of = {}
    for (ca, cb), r in dist.items():
        ia = idx_of.get(frozenset((ca[0], ca[1])))
        ib = idx_of.get(frozenset((cb[0], cb[1])))
        if ia is not None and ib is not None:
            gap_of[(min(ia, ib), max(ia, ib))] = r
    gap = gap_of.get(pick, float("nan"))
    print(f"  chains {a} and {b}, parallel along axis "
          f"{'xyz'[axis_of(graph, edges[a])]}, shell {args.shell}, "
          f"{gap:.2f} lattice units apart")
    print(f"  asking both methods for {args.count} entanglements on this pair "
          f"and nothing else\n")

    # ---- the plain melt ---------------------------------------------------
    rng0 = np.random.default_rng(args.seed)
    paths0 = {k: bridging_walk(c0, c1, args.dp + 1, BOND, rng0)
              for k, (c0, c1) in geo["chords"].items()}
    root = OUT / "onepair_plain"
    shutil.rmtree(root, ignore_errors=True)
    _n, node_atom, chain_atoms = write_system(graph, geo, paths0, root)
    seq = {k: chain_ids(k, node_atom, chain_atoms, geo["ends"]) for k in keys}
    sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                             protocol="hardcore")
    print("  --- plain melt ---")
    run_md(sim, args.stages)
    after = {1: "system_after_soft.data", 2: "system_ramped.data",
             3: "system_equilibrated.data"}[args.stages]
    plain_f = root / "04_Simulation" / after
    if not plain_f.exists():
        print("  the plain melt did not relax")
        return 1
    box, paths, xyz0 = paths_from(plain_f, keys, seq)

    work = OUT / f"onepair_work_{os.getpid()}"
    base = measure_pairs(plain_f, seq, [pick], work)[pick]

    # ---- the old kink, this pair only -------------------------------------
    e1 = (edges[a][0], edges[a][1], 0)
    e2 = (edges[b][0], edges[b][1], 0)
    kinked, _p, _s = kinked_paths(graph, geo, [(e1, e2, args.count)], dims,
                                  args.dp,
                                  {"overshoot": 0.2, "z_amp": 0.5,
                                   "sigma": 0.15})
    old_root = OUT / "onepair_old"
    shutil.rmtree(old_root, ignore_errors=True)
    _n2, na2, ca2 = write_system(graph, geo, kinked, old_root)
    seq_old = {k: chain_ids(k, na2, ca2, geo["ends"]) for k in keys}
    # The old kink gets its own protocol by default.
    #
    # It lays chains along their chords, which at DP 80 means bonds near 0.07
    # sigma, and the soft push in the generated scripts is what expands them.
    # Judging it under the hard-core protocol instead would be judging it on a
    # conformation it was never built for.
    sim2 = conform_and_script(old_root, graph, geo, pair_style="repulsive",
                              protocol=args.old_protocol)
    print(f"  --- old kink, one pair, {args.old_protocol} protocol ---")
    run_md(sim2, args.stages)
    old_f = old_root / "04_Simulation" / after
    old_got = (measure_pairs(old_f, seq_old, [pick], work)[pick]
               if old_f.exists() else None)

    # ---- the new construction, this pair only -----------------------------
    mine = set(seq[a])
    avoid = Clearance(np.array([xyz0[i] for i in sorted(xyz0)
                                if i not in mine]), box, args.clearance)
    cand = OUT / f"onepair_cand_{os.getpid()}"
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

    print("  --- new construction, one pair ---")
    p, built = construct_exact(paths, a, b, args.count, box, avoid, seq,
                               plain_f, work, (0.3, 0.7), args.ring, args.dp,
                               ranks=args.ranks, phases=args.phases,
                               relax=relax)
    new_got = None
    if p is not None:
        xyz_now = dict(xyz0)
        for aid, xyzp in zip(seq[a], p):
            xyz_now[aid] = xyzp
        new_root = OUT / "onepair_new"
        (new_root / "03_Conformation").mkdir(parents=True, exist_ok=True)
        for sub in ("02_Chemistry", "04_Simulation"):
            (new_root / sub).mkdir(parents=True, exist_ok=True)
            for f in (root / sub).glob("*"):
                if f.is_file() and (sub == "02_Chemistry"
                                    or f.suffix == ".in"):
                    (new_root / sub / f.name).write_text(f.read_text())
        rewrite_coords(plain_f, new_root / "03_Conformation" /
                       "system_relaxed.data", xyz_now)
        run_md(new_root / "04_Simulation", args.stages)
        new_f = new_root / "04_Simulation" / after
        new_got = (measure_pairs(new_f, seq, [pick], work)[pick]
                   if new_f.exists() else None)

    def show(v):
        return "not measured" if v is None else str(v)

    print(f"\n  {'':<32} {'entanglements on this pair':>28}")
    print(f"  {'asked for':<32} {args.count:>28}")
    print(f"  {'the melt already had':<32} {show(base):>28}")
    print(f"  {f'old kink ({args.old_protocol})':<32} {show(old_got):>28}")
    print(f"  {'new construction delivered':<32} {show(new_got):>28}")
    if old_got is not None and new_got is not None:
        same = old_got == new_got
        print(f"\n  the two methods {'agree' if same else 'do not agree'} "
              f"on this pair")
    print("\n  one pair, one seed, one lattice. Nothing else was entangled, "
          "so the melt's own entanglement is the only other contributor.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
