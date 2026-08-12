"""Does the legacy kink deliver the pair counts it selects?

    python tests/workflows/entangle_legacy_check.py --per-chain 2.0

`select_entanglements` picks pairs and assigns each a count, and
`calculate_entangled_kink` builds a bulge sized for that count. Whether the
built system actually carries those entanglements has never been measured --
the count is an input to the geometry, not an output of it.

This measures it. Build the system the legacy way, run the minimisation it is
designed for, then ask Z1+ what is actually entangled with what.

The measurement has to be after minimisation, not as built. At melt density
the legacy path lays chains along their chords, which at DP 80 means bonds
near 0.07 sigma -- a collapsed chain that Z1+ refuses outright. That is not a
fault of the method; the soft push is what expands it, and the expanded state
is the one that means anything.

The whole-system export used to fail here and the reason is now known: Z1+
rejects a configuration read back from LAMMPS output whenever a chain is
longer than the periodic cell, which after minimisation is routine rather than
exceptional. Measuring one pair at a time avoids it, and that is what the
designed-entanglement workflows do. Pairs that still cannot be measured are
reported and excluded, never counted as zero, because a failed measurement and
an absent entanglement print identically.
"""
from __future__ import annotations

import argparse
import collections
import os
import random
import sys
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.assignment.entanglements import (  # noqa: E402
    find_crossing_candidates,
    select_entanglements,
)
from topon.config.schema import EntanglementsConfig  # noqa: E402
from tests.workflows.entangle_all import CASES  # noqa: E402
from tests.workflows.entangle_relaxed import measure_pairs  # noqa: E402
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
    report_bonds,
    run_md,
    write_system,
    z1_export,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--per-chain", type=float, default=2.0)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--protocol", default="generated",
                    choices=("generated", "hardcore"))
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (4, 4, 4)
    graph = build_network(spec)
    geo = geometry(graph, dp=args.dp, density=0.85)
    dims = np.asarray(graph.graph["box"], float)

    G = nx.MultiGraph()
    G.add_nodes_from(graph.nodes(data=True))
    for u, v in graph.edges():
        G.add_edge(u, v)

    cands = find_crossing_candidates(G, dims)
    random.seed(args.seed)
    cfg = EntanglementsConfig(enabled=True,
                              avg_crosslinks_per_chain=args.per_chain)
    sel = select_entanglements(G, cfg, dims, candidates=list(cands),
                               num_chains=G.number_of_edges())
    paths, partner, _sites = kinked_paths(
        graph, geo, sel, dims, args.dp,
        {"overshoot": 0.2, "z_amp": 0.5, "sigma": 0.15})

    asked = sum(c for _, _, c in sel)
    print(f"  {len(cands)} candidates, {len(sel)} pairs selected, "
          f"{asked} entanglements requested")

    root = OUT / f"legacy_check_{args.protocol}"
    _n, node_atom, chain_atoms = write_system(graph, geo, paths, root)
    sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                             protocol=args.protocol)
    print(f"\n  --- LAMMPS, stages 1 to {args.stages} ---")
    run_md(sim, args.stages)
    print()
    report_bonds(root)

    final = root / "04_Simulation" / {
        1: "system_after_soft.data", 2: "system_ramped.data",
        3: "system_equilibrated.data"}[args.stages]
    if not final.exists():
        print("\n  minimisation produced no output")
        return 1

    keys = sorted(geo["chords"])
    seq = {k: chain_ids(k, node_atom, chain_atoms, geo["ends"]) for k in keys}

    # One pair at a time, not the whole system.
    #
    # Exporting every chain at once is what stalled this measurement. Z1+
    # rejects a configuration read back from LAMMPS output whenever a chain is
    # longer than the periodic cell, and after minimisation that is routine.
    # Per pair it works, which is how the designed-entanglement workflows
    # measure survival, and it is the same question being asked here.
    key_of = {frozenset(v): k for k, v in geo["ends"].items()}
    wanted = []
    for e1, e2, cnt in sel:
        k1 = key_of.get(frozenset((e1[0], e1[1])))
        k2 = key_of.get(frozenset((e2[0], e2[1])))
        if k1 is not None and k2 is not None and k1 != k2:
            wanted.append((min(k1, k2), max(k1, k2), cnt))
    if not wanted:
        print("\n  no selected pair maps onto a chain, nothing to check")
        return 1

    work = OUT / f"legacy_z1_{os.getpid()}"
    pl = [(a, b) for a, b, _c in wanted]
    got = measure_pairs(final, seq, pl, work)
    base = measure_pairs(base_final, seq, pl, work) if base_final else {
        q: None for q in pl}

    blind = sum(1 for a, b, _c in wanted if got[(a, b)] is None)
    if blind:
        print(f"\n  {blind} of {len(wanted)} pairs could not be measured. "
              f"They are left out below rather than counted as zero, since a "
              f"failed measurement and an absent entanglement look identical.")

    rows = [(a, b, c, got[(a, b)]) for a, b, c in wanted
            if got[(a, b)] is not None]
    if not rows:
        print("  nothing measured at all; no conclusion can be drawn")
        return 1

    hit = sum(1 for _a, _b, c, g in rows if g == c)
    hist = collections.Counter((c, g) for _a, _b, c, g in rows)
    print(f"\n  {hit} of {len(rows)} requested pairs carry exactly the "
          f"count they were given ({100.0 * hit / len(rows):.0f}%)")
    print(f"  {sum(1 for _a, _b, _c, g in rows if g == 0)} carry none at all")
    print(f"  asked for {sum(c for _a, _b, c, _g in rows)} entanglements "
          f"over these pairs, measured "
          f"{sum(g for _a, _b, _c, g in rows)}")
    print("\n  asked -> got, most common:")
    for (a, b), n in hist.most_common(8):
        print(f"    asked {a}, got {b}: {n} pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
