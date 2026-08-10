"""Select entanglements after conformation, ranking pairs by proximity.

    python tests/workflows/entangle_by_proximity.py
    python tests/workflows/entangle_by_proximity.py --bias region
    python tests/workflows/entangle_by_proximity.py --no-rank    # control

The pipeline places entanglements during assignment, before any coordinates
exist, so the only thing available to choose pairs by is the distance between
their crosslinks. That is a property of the network rather than of the
chains, and it is the wrong quantity: capacity for entanglement is set by how
much of two chains lies alongside, and two chains can be nearest neighbours
by crosslink and barely touch.

This runs the selection in passes instead:

    1. draw a provisional conformation with no entanglements in it
    2. score every candidate pair on that conformation
    3. select, with the score multiplying into whatever spatial or shell
       bias was configured
    4. redraw the chosen chains with their kinks

Measured on a melt of 106 chains, ranking every pair this way and checking
against a primitive-path analysis: the top 20 average 3.70 entanglements and
the bottom 50 carry none at all.
"""
from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.assignment.entanglements import (  # noqa: E402
    compute_proximity_weights,
    find_crossing_candidates,
    select_entanglements,
)
from topon.config.schema import EntanglementsConfig  # noqa: E402
from topon.conformation.paths import bridging_walk  # noqa: E402
from tests.workflows.entangle_spatial import BIAS, kinked_paths  # noqa: E402
from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    DP,
    OUT,
    build_network,
    conform_and_script,
    geometry,
    report_bonds,
    run_md,
    write_system,
)


def provisional(geo, dp, bond, seed):
    """Pass 1: a conformation with no entanglements, to rank on.

    Coiled rather than straight. Straight chains lie on their chords, and two
    chords being close says nothing about whether the chains meet -- which is
    the whole reason for doing this after conformation rather than before.
    """
    rng = np.random.default_rng(seed)
    return {k: bridging_walk(c0, c1, dp + 1, bond, rng)
            for k, (c0, c1) in geo["chords"].items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bias", default="uniform", choices=sorted(BIAS))
    ap.add_argument("--per-chain", type=float, default=0.20)
    ap.add_argument("--cutoff", type=float, default=2.0,
                    help="bead pairs closer than this count toward a score")
    ap.add_argument("--no-rank", action="store_true",
                    help="skip the proximity pass, for comparison")
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--density", type=float, default=0.85)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--no-md", action="store_true")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    graph = build_network()
    geo = geometry(graph, dp=args.dp, density=args.density)
    dims = np.asarray(graph.graph["box"], float)

    G = nx.MultiGraph()
    G.add_nodes_from(graph.nodes(data=True))
    for u, v in graph.edges():
        G.add_edge(u, v)
    cands = find_crossing_candidates(G, dims)

    # Pass 1 ---------------------------------------------------------------
    paths0 = provisional(geo, args.dp, BOND, args.seed)
    print(f"  provisional conformation: {len(paths0)} coiled chains, "
          f"box {geo['L'][0]:.1f} sigma, density {geo['density']:.3f}")

    cfg = EntanglementsConfig(enabled=True,
                              avg_crosslinks_per_chain=args.per_chain,
                              placement_bias_kind=args.bias,
                              placement_bias_params=BIAS[args.bias])

    # Assignment works in lattice units, the same as the node positions, so
    # the conformation is handed over in those units and the cutoff with it.
    scale = geo["scale"]
    by_pair = {frozenset(v): paths0[k] / scale for k, v in geo["ends"].items()}
    cut = args.cutoff / scale
    allw = compute_proximity_weights(cands, by_pair, box=dims, cutoff=cut)

    # Pass 2 ---------------------------------------------------------------
    if not args.no_rank:
        live = [x for x in allw if x > 0]
        print(f"  scored {len(cands)} candidates: {len(live)} have chains "
              f"that come within {args.cutoff} sigma")
        if live:
            print(f"    score range {min(live):.0f} to {max(live):.0f}, "
                  f"median {np.median(live):.0f}")


    # Pass 3 ---------------------------------------------------------------
    random.seed(args.seed)
    sel = select_entanglements(G, cfg, dims, candidates=list(cands),
                               num_chains=G.number_of_edges(),
                               chain_paths=None if args.no_rank else by_pair,
                               proximity_cutoff=cut)
    print(f"  selected {len(sel)} entanglements")

    if sel:
        chosen = compute_proximity_weights(
            [(a, b) for a, b, _ in sel], by_pair, box=dims, cutoff=cut)
        pool = [x for x in allw if x > 0]
        print(f"  chosen pairs score {np.median(chosen):.0f} median, "
              f"against {np.median(pool):.0f} for the candidate pool "
              f"({sum(1 for x in chosen if x == 0)} of {len(chosen)} chosen "
              f"never come within range)")

    # Pass 4 ---------------------------------------------------------------
    paths, partner, sites = kinked_paths(
        graph, geo, sel, dims, args.dp,
        {"overshoot": 0.2, "z_amp": 0.5, "sigma": 0.15})
    bonds = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                            for p in paths.values()])
    print(f"  {len(partner)} chains kinked, bonds "
          f"{bonds.min():.3f} to {bonds.max():.3f}")

    tag = "unranked" if args.no_rank else "ranked"
    root = OUT / f"proximity_{tag}_{args.bias}"
    write_system(graph, geo, paths, root)
    sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                             protocol="hardcore")
    if not args.no_md:
        print(f"\n  --- LAMMPS, stages 1 to {args.stages} ---")
        run_md(sim, args.stages)
        print()
        report_bonds(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
