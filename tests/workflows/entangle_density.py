"""Ask for e entanglements per chain, get e entanglements per chain.

    python tests/workflows/entangle_density.py --want 2.0
    python tests/workflows/entangle_density.py --want 1.0 --shells 2:1.0

A system-averaged target, not a named pair. That difference is what makes it
robust: for a density it does not matter *which* pairs carry the entanglements,
only how many there are, so a pair that fails to survive can be replaced by
another and the collateral a routed chain picks up on its way counts toward the
total rather than against it.

So the loop closes on the measured value:

    select -> route -> relax -> measure the whole system -> top up or stop

Each round adds only the shortfall, and the measurement is always of the
relaxed system, so what is reported is what survives rather than what was
built. Survival of an individual pair stops mattering; the delivered density is
what is being controlled.
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
    select_by_shells,
)
from topon.conformation.paths import Clearance, bridging_walk  # noqa: E402
from tests.workflows.entangle_all import CASES  # noqa: E402
from tests.workflows.entangle_relaxed import (  # noqa: E402
    construct,
    paths_from,
    rewrite_coords,
)
from tests.workflows.entangle_search import measure_batch  # noqa: E402
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
    z1_export,
)


def close_pairs(paths, box, cutoff=3.0):
    """Chain pairs that come within ``cutoff`` of each other.

    Two chains that never come close cannot be entangled, so restricting the
    measurement to these pairs is exact rather than a sample, and it is what
    makes the measurement affordable: a 648-chain system has 210k pairs and
    only a few thousand are ever in contact.

    Found with a KD-tree over the beads rather than chain by chain. The
    pairwise form is quadratic in the chain count and quadratic again in the
    beads per chain, which is minutes at 648 chains; this is seconds.
    """
    from scipy.spatial import cKDTree

    keys = sorted(paths)
    owner, pts = [], []
    for k in keys:
        q = paths[k]
        q = q - box * np.floor(q / box)
        q = np.clip(q, 0.0, np.nextafter(box, 0.0))
        pts.append(q)
        owner += [k] * len(q)
    tree = cKDTree(np.vstack(pts), boxsize=box)
    owner = np.asarray(owner)

    out = set()
    for i, j in tree.query_pairs(cutoff):
        a, b = owner[i], owner[j]
        if a != b:
            out.add((min(a, b), max(a, b)))
    return sorted(out)


def measure_system(data_file, seq, keys, work, paths, box, cutoff=3.0):
    """Entanglement points per chain, summed over every pair in contact.

    Built from the per-pair export rather than a whole-system one. The
    whole-system export is not usable here: Z1+ refuses a configuration once a
    chain is longer than the periodic cell, which at the low densities routing
    needs is routine, and it has a size ceiling besides -- 8544 beads measured,
    52056 did not. The per-pair export has worked throughout.

    Returns ``(points_per_chain, partners_per_chain, pairs_measured,
    pairs_unmeasured)``.
    """
    from tests.workflows.entangle_relaxed import measure_pairs

    pairs = close_pairs(paths, box, cutoff)
    if not pairs:
        return 0.0, 0.0, 0, 0
    got = measure_pairs(data_file, seq, pairs, work)
    live = [v for v in got.values() if v is not None]
    blind = sum(1 for v in got.values() if v is None)
    pts = sum(live)
    partners = sum(1 for v in live if v)
    n = len(keys)
    return 2.0 * pts / n, 2.0 * partners / n, len(pairs), blind


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--want", type=float, default=2.0,
                    help="target entanglements per chain")
    ap.add_argument("--shells", default="1:0.5,2:0.5",
                    help="shell mix, e.g. '1:0.2,2:0.5,3:0.25,4:0.05'")
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--tol", type=float, default=0.15,
                    help="stop when within this many entanglements per chain")
    ap.add_argument("--dims", type=int, default=4)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--ring", type=float, default=2.0)
    ap.add_argument("--density", type=float, default=None,
                    help="melt route: pack to this density instead of sizing "
                         "by coil ratio. 0.85 is the physical melt the LAMMPS "
                         "scripts are calibrated for.")
    ap.add_argument("--coil", type=float, default=6.0,
                    help="contour over chord, which sets the box. 6 puts "
                         "neighbouring chains 1.5 sigma apart with real free "
                         "volume to route through; 1.8 leaves them 35 sigma "
                         "apart and unreachable, 12 leaves no room at all.")
    ap.add_argument("--margin", type=float, default=2.0)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--clearance", type=float, default=0.9)
    ap.add_argument("--cutoff", type=float, default=3.0,
                    help="chains further apart than this everywhere cannot be "
                         "entangled, so they are not measured")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tag", default="_density")
    args = ap.parse_args()

    mix = {}
    for part in args.shells.split(","):
        s, f = part.split(":")
        mix[int(s)] = float(f)

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
    shells = neighbour_shells(G, dims, max_shell=max(mix), distances=dist)
    edges = sorted(graph.edges())
    idx_of = {frozenset(e): i for i, e in enumerate(edges)}

    def to_idx(sel):
        out = []
        for e1, e2, c in sel:
            a = idx_of.get(frozenset((e1[0], e1[1])))
            b = idx_of.get(frozenset((e2[0], e2[1])))
            if a is not None and b is not None and a != b:
                out.append((min(a, b), max(a, b), c))
        return out

    rng = np.random.default_rng(args.seed)
    plan = to_idx(select_by_shells(G, args.want, mix, dims, shells=shells,
                                   rng=rng))
    if not plan:
        print("  selection produced nothing")
        return 1
    print(f"  target {args.want:.2f} per chain, mix "
          f"{ {s: mix[s] for s in sorted(mix)} }")
    print(f"  selected {len(plan)} pairs, "
          f"{sum(c for _a, _b, c in plan)} entanglements")

    # Sized by coil ratio, not by the design.
    #
    # `scale_for_design` sizes the box so every requested route fits, which is
    # right for a handful of named pairs and wrong here: with a hundred pairs
    # some chain has several partners, its ring and detour budget dominates,
    # and the box collapses -- measured, 2.2 sigma at density 843. A density
    # target does not need every pair to build. It needs a physically sensible
    # box and enough pairs, and the loop replaces whatever does not fit.
    geo = (geometry(graph, dp=args.dp, density=args.density)
           if args.density else
           geometry(graph, dp=args.dp, bond=BOND, coil=args.coil))
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
    print("\n  --- relaxing a plain melt ---")
    run_md(sim, args.stages)
    after = {1: "system_after_soft.data", 2: "system_ramped.data",
             3: "system_equilibrated.data"}[args.stages]
    relaxed = root / "04_Simulation" / after
    if not relaxed.exists():
        print("  relaxation produced no output")
        return 1
    box, paths, xyz0 = paths_from(relaxed, keys, seq)

    work = OUT / f"density_work_{os.getpid()}"
    base = measure_system(relaxed, seq, keys, work, paths, box,
                          args.cutoff)
    print(f"  the plain melt already carries {base[0]:.2f} per chain over "
          f"{base[1]:.1f} partners ({base[2]} pairs in contact"
          + (f", {base[3]} unmeasured)" if base[3] else ")"))

    cand = OUT / f"density_cand_{os.getpid()}"
    (cand / "03_Conformation").mkdir(parents=True, exist_ok=True)
    for sub in ("02_Chemistry", "04_Simulation"):
        (cand / sub).mkdir(parents=True, exist_ok=True)
        for f in (root / sub).glob("*"):
            if f.is_file() and (sub == "02_Chemistry" or f.suffix == ".in"):
                (cand / sub / f.name).write_text(f.read_text())

    def relax(xyz):
        rewrite_coords(relaxed, cand / "03_Conformation" /
                       "system_relaxed.data", xyz)
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

    xyz_now = dict(xyz0)
    done, target = set(), args.want + base[0]
    print(f"\n  {'round':>6} {'routed':>7} {'per chain':>10} "
          f"{'partners':>9} {'designed':>9}")
    for rnd in range(1, args.rounds + 1):
        for a, b, count in plan:
            if (a, b) in done:
                continue
            mine = set(seq[a])
            avoid = Clearance(np.array([xyz_now[i] for i in sorted(xyz_now)
                                        if i not in mine]), box,
                              args.clearance)
            try:
                p = construct(paths, a, b, max(0.5, 0.5 * count), box, avoid,
                              radius=args.ring, dp=args.dp, span=(0.3, 0.7))
            except ValueError:
                continue
            paths[a] = p
            for aid, xyzp in zip(seq[a], p):
                xyz_now[aid] = xyzp
            done.add((a, b))

        out = relax(xyz_now)
        if out is None:
            print(f"  {rnd:>6}   relaxation failed")
            return 1
        got = measure_system(out, seq, keys, work, paths, box,
                             args.cutoff)
        print(f"  {rnd:>6} {len(done):>7} {got[0]:>10.2f} {got[1]:>9.1f} "
              f"{len(done):>9}"
              + (f"   ({got[3]} of {got[2]} pairs unmeasured)"
                 if got[3] else ""))

        added = got[0] - base[0]
        if abs(added - args.want) <= args.tol:
            print(f"\n  delivered {added:.2f} per chain against a target of "
                  f"{args.want:.2f}, within {args.tol}")
            break
        if added > args.want:
            print(f"\n  overshot: {added:.2f} against {args.want:.2f}")
            break
        # Top up only the shortfall, so the loop converges instead of
        # doubling.
        short = args.want - added
        more = to_idx(select_by_shells(G, short, mix, dims, shells=shells,
                                       rng=rng))
        plan = [q for q in more if (q[0], q[1]) not in done]
        if not plan:
            print(f"\n  no unused pair left; delivered {added:.2f} against "
                  f"{args.want:.2f}")
            break
    else:
        print(f"\n  did not reach the target in {args.rounds} rounds")

    Path(root / "04_Simulation" / "designed.data").write_text(
        Path(rewrite_coords(relaxed,
                            root / "04_Simulation" / "designed.data",
                            xyz_now)).read_text())
    print(f"  system: {root / '04_Simulation' / 'designed.data'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
