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

STATUS: does not yet produce an answer. Z1+ rejects the minimised
configuration as well, and it is not clear why. The file is well formed --
106 chains, 82 beads each, box 21.58, bonds 0.373 to 1.636 with none over 2
sigma, coordinates spanning 35 sigma because chains are unwrapped across the
boundary. The same export path works on a two-chain pair, and the same
whole-system shape works when the paths are written directly rather than read
back from a data file, so the failing combination is whole-system *and* read
from LAMMPS output. Until that is resolved, whether the legacy kink delivers
the counts it selects is unmeasured -- which is worth stating plainly, since
the count has always been an input to the geometry and never an output
checked against it.
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
from tests.workflows.entangle_search import _both, measure_batch  # noqa: E402
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
    idx = {k: i + 1 for i, k in enumerate(keys)}
    seqs = [chain_ids(k, node_atom, chain_atoms, geo["ends"]) for k in keys]

    work = OUT / f"legacy_z1_{os.getpid()}"
    work.mkdir(parents=True, exist_ok=True)
    for old in list(work.glob("*.Z1")) + list(work.glob("SP_*.dat")):
        old.unlink()
    z1_export(final, seqs, work / "final.Z1")
    res = (measure_batch(work) or {}).get("final", {})
    total = sum(sum(v.values()) for v in res.values())
    print(f"\n  Z1+ after minimisation: {total} entanglement points "
          f"over {len(keys)} chains")
    if not total:
        print("  nothing measured; the numbers below mean nothing")
        return 1

    key_of = {frozenset(v): k for k, v in geo["ends"].items()}
    hit = 0
    rows = []
    hist = collections.Counter()
    for e1, e2, cnt in sel:
        k1 = key_of.get(frozenset((e1[0], e1[1])))
        k2 = key_of.get(frozenset((e2[0], e2[1])))
        if k1 is None or k2 is None:
            continue
        got = _both(res, idx[k1], idx[k2])
        rows.append((k1, k2, cnt, got))
        hist[(cnt, got)] += 1
        hit += (got == cnt)

    print(f"\n  {hit} of {len(rows)} requested pairs carry exactly the "
          f"count they were given ({100.0 * hit / max(len(rows), 1):.0f}%)")
    print(f"  {sum(1 for _, _, _, g in rows if g == 0)} carry none at all")
    print("\n  asked -> got, most common:")
    for (a, b), n in hist.most_common(8):
        print(f"    asked {a}, got {b}: {n} pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
