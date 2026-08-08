"""How often does the delivered entanglement count match the asked one.

    python tests/workflows/entangle_validate.py
    python tests/workflows/entangle_validate.py --cases 12 --stages 3

Runs a set of designed pairs through the whole pipeline and measures each
one with Z1+, the network stripped so only that pair is counted. Reports the
hit rate rather than a single example, because a single example is what a
distribution looks like when you only draw once: the same configuration was
seen to give 3/3 one run and 2/4 the next before the conformation noise was
seeded.

Nothing here is tuned per case. Every case gets the same coil, the same
bond, the same protocol, and the sites are placed evenly along the chain.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    COIL,
    DP,
    OUT,
    Site,
    build_group,
    build_network,
    chain_ids,
    conform_and_script,
    geometry,
    report_bonds,
    run_md,
    run_z1,
    separation_bands,
    write_system,
    z1_export,
)


def cases(geo, want, bands=(1, 2, 3)):
    """One pair per band, several site counts each."""
    out = []
    found = separation_bands(geo)
    for b in bands:
        if b > len(found):
            continue
        gap_u, ka, kb = found[b - 1][0]
        for n in want:
            out.append((b, gap_u, ka, kb, n))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sites", type=int, nargs="*", default=[1, 2, 3])
    ap.add_argument("--bands", type=int, nargs="*", default=[1, 2, 3])
    ap.add_argument("--coil", type=float, default=COIL)
    ap.add_argument("--bond", type=float, default=BOND)
    ap.add_argument("--stages", type=int, default=3)
    ap.add_argument("--reach", type=float, default=0.45)
    args = ap.parse_args()

    graph = build_network()
    probe = geometry(graph, dp=DP, bond=args.bond, coil=args.coil)
    todo = cases(probe, args.sites, args.bands)
    print(f"  {len(todo)} cases: bands {args.bands} x sites {args.sites}, "
          f"coil {args.coil}, {args.stages} stage(s)\n")

    z1_dir = OUT / "validate_z1"
    if z1_dir.exists():
        for f in z1_dir.glob("*.Z1"):
            f.unlink()
    z1_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for i, (band, gap_u, ka, kb, n) in enumerate(todo):
        tag = f"c{i:02d}_b{band}_n{n}"
        sites = [Site(at=(j + 1) / (n + 1), turns=1) for j in range(n)]
        plan = [(ka, kb, sites)]
        try:
            geo, paths, info = build_group(graph, plan, args.bond, DP,
                                           args.reach, args.coil)
        except ValueError as exc:
            rows.append((tag, band, n, None, str(exc)[:44]))
            print(f"  {tag}: refused, {str(exc)[:60]}")
            continue

        root = OUT / f"validate_{tag}"
        for stale in ("04_Simulation", "03_Conformation"):
            pass
        _, node_atom, chain_atoms = write_system(graph, geo, paths, root)
        seqs = [chain_ids(k, node_atom, chain_atoms, geo["ends"])
                for k in (ka, kb)]
        sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                                 protocol="hardcore")
        run_md(sim, args.stages)

        final = {1: "system_after_soft.data", 2: "system_ramped.data",
                 3: "system_equilibrated.data"}[args.stages]
        out_file = root / "04_Simulation" / final
        if not out_file.exists():
            rows.append((tag, band, n, None, "MD produced no output"))
            continue
        z1_export(out_file, seqs, z1_dir / f"{tag}.Z1")
        bonds = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                                for p in paths.values()])
        rows.append((tag, band, n, float(bonds.max()), None))
        print(f"  {tag}: built and run, longest bond {bonds.max():.2f}")

    print("\n  measuring with Z1+ ...")
    z = run_z1(z1_dir)

    print()
    print("  " + "-" * 60)
    print(f"  {'case':>12} {'band':>5} {'asked':>6} {'Z':>8} {'max bond':>9} "
          f"{'verdict':>8}")
    print("  " + "-" * 60)
    hits = total = 0
    for tag, band, n, bond, err in rows:
        if err:
            print(f"  {tag:>12} {band:5d} {n:6d} {'-':>8} {'-':>9} "
                  f"{'refused':>8}   {err}")
            continue
        got = z.get(tag) if z else None
        total += 1
        if got:
            shown = "/".join(str(v) for v in got)
            ok = all(v == n for v in got)
            hits += ok
            verdict = "hit" if ok else "miss"
        else:
            shown, verdict = "-", "?"
        print(f"  {tag:>12} {band:5d} {n:6d} {shown:>8} {bond:9.2f} "
              f"{verdict:>8}")
    print("  " + "-" * 60)
    if total:
        print(f"  {hits} of {total} delivered exactly what was asked "
              f"({100.0 * hits / total:.0f}%)")
    print()
    print("  Z is per chain, both chains of the pair, measured on the data")
    print("  file after the last stage with the network removed. A hit means")
    print("  both chains read exactly the number of sites asked for.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
