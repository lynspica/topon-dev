"""Is the relaxed state a fixed point of the minimiser?

    python tests/workflows/entangle_fixedpoint.py

A designed entanglement can only survive if the state it is designed into is
one the minimiser is already happy with. If the system is still relaxing when
the design goes in, the residual relaxation sweeps the design away, and no
amount of care in placing it matters.

The test is to relax twice and compare. Feed the relaxed coordinates straight
back into the same protocol with nothing changed: whatever moves the second
time is motion the design would have had to survive for no reason at all.
Anything the design does is on top of that.

Reported per coil ratio, because coil sets the box, and box sets both how much
free volume there is to route through and how much the neighbours constrain a
chain from rearranging. Those pull in opposite directions and the point of
this is to find where they balance.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.paths import Clearance, bridging_walk  # noqa: E402
from tests.workflows.entangle_all import CASES  # noqa: E402
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
    write_system,
)


def relax(root, graph, geo, paths, stages, tag):
    """One pass of the protocol. Returns the resulting coordinates."""
    shutil.rmtree(root, ignore_errors=True)
    _n, node_atom, chain_atoms = write_system(graph, geo, paths, root)
    sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                             protocol="hardcore")
    run_md(sim, stages)
    out = root / "04_Simulation" / {
        1: "system_after_soft.data", 2: "system_ramped.data",
        3: "system_equilibrated.data"}[stages]
    if not out.exists():
        return None, None, None
    box, xyz, _ = read_data(out)
    return box, xyz, (node_atom, chain_atoms)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coils", type=float, nargs="*",
                    default=[4.0, 6.0, 9.0, 12.0])
    ap.add_argument("--density", type=float, nargs="*", default=[0.85])
    ap.add_argument("--dims", type=int, default=4)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (args.dims,) * 3
    graph = build_network(spec)
    n_beads = graph.number_of_edges() * args.dp + graph.number_of_nodes()

    print(f"  {'scale':>14} {'box':>7} {'density':>8} {'crowding':>9} "
          f"{'pass 1 moves':>13} {'pass 2 moves':>13} {'settled':>8}")
    runs = ([("coil", c) for c in args.coils]
            + [("density", d) for d in args.density])
    for how, val in runs:
        geo = (geometry(graph, dp=args.dp, bond=BOND, coil=val)
               if how == "coil"
               else geometry(graph, dp=args.dp, density=val))
        L = geo["L"]
        rho = n_beads / float(np.prod(L))
        rng = np.random.default_rng(args.seed)
        paths = {k: bridging_walk(c0, c1, args.dp + 1, BOND, rng)
                 for k, (c0, c1) in geo["chords"].items()}

        built = {}
        root = OUT / f"fixedpoint_{how}_{val}"
        box, xyz1, meta = relax(root, graph, geo, paths, args.stages, "1")
        if xyz1 is None:
            print(f"  {how} {val:<9} relaxation produced no output")
            continue
        node_atom, chain_atoms = meta
        for k in sorted(geo["chords"]):
            built[k] = chain_ids(k, node_atom, chain_atoms, geo["ends"])
        start = {aid: p for k in built
                 for aid, p in zip(built[k], paths[k])}

        # Second pass: same protocol, starting from where the first ended.
        paths2 = {k: np.array([xyz1[a] for a in built[k]]) for k in built}
        _b, xyz2, _m = relax(OUT / f"fixedpoint_{how}_{val}_again", graph, geo,
                             paths2, args.stages, "2")
        if xyz2 is None:
            print(f"  {how} {val:<9} second pass produced no output")
            continue

        def moved(a, b):
            ids = sorted(set(a) & set(b))
            d = np.array([b[i] - a[i] for i in ids])
            d -= L * np.round(d / L)
            return float(np.median(np.linalg.norm(d, axis=1)))

        m1, m2 = moved(start, xyz1), moved(xyz1, xyz2)
        crowd = rho * 4.0 / 3.0 * np.pi * 0.9 ** 3
        print(f"  {how + ' ' + str(val):>14} {L[0]:7.1f} {rho:8.3f} "
              f"{crowd:9.2f} {m1:13.2f} {m2:13.2f} "
              f"{'yes' if m2 < 0.3 else 'no':>8}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
