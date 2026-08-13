"""Verification: does a requested shell distribution come back as itself?

    python tests/workflows/entangle_shellsuite.py
    python tests/workflows/entangle_shellsuite.py --dims 5 --want 2.0

Runs several distributions end to end and tabulates asked against delivered for
each. One run proves nothing about a controller -- it can hit one target by
luck, or by a bias that happens to point the right way for that particular
mix. A set of them, including the degenerate single-shell cases and a uniform
one, is what shows the mix is being controlled rather than approximated.

Each case is a full build: select, route, minimise with LAMMPS, and measure
with Z1+ on the minimised system. Nothing here reads the built coordinates.

The distributions are chosen to probe different failure modes:

    all in one shell        can it concentrate at all
    even over two          the simplest mix
    the four-shell request  the case actually asked for
    uniform over four      no shell favoured, so any bias shows up
    weighted to the far shell   outer shells deliver more per pair, so this
                                is where over-delivery would show
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent

CASES = [
    ("all shell 1", "1:1.0"),
    ("all shell 2", "2:1.0"),
    ("even, 1 and 2", "1:0.5,2:0.5"),
    ("the four-shell request", "1:0.2,2:0.5,3:0.25,4:0.05"),
    ("uniform over four", "1:0.25,2:0.25,3:0.25,4:0.25"),
    ("weighted outward", "1:0.05,2:0.15,3:0.3,4:0.5"),
]


def parse(out):
    """Pull the final asked/delivered table out of a run's output."""
    rows = {}
    in_table = False
    for line in out.splitlines():
        if re.search(r"shell\s+asked\s+delivered", line):
            in_table = True
            continue
        if in_table:
            m = re.match(r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)", line)
            if m:
                rows[int(m.group(1))] = (float(m.group(2)),
                                         float(m.group(3)))
            elif rows:
                break
    got = re.search(r"delivered ([\d.-]+) per chain against a target of "
                    r"([\d.]+)", out)
    every = "every shell within" in out
    return rows, (float(got.group(1)) if got else None), every


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--want", type=float, default=2.0)
    ap.add_argument("--dims", type=int, default=4)
    ap.add_argument("--density", type=float, default=0.85)
    ap.add_argument("--rounds", type=int, default=8)
    ap.add_argument("--tol-pct", type=float, default=5.0)
    ap.add_argument("--only", type=int, default=None,
                    help="run just one case, by index")
    args = ap.parse_args()

    cases = CASES if args.only is None else [CASES[args.only]]
    results = []
    for i, (name, mix) in enumerate(cases):
        print(f"\n{'=' * 70}\n  case {i}: {name}   --shells {mix}\n{'=' * 70}",
              flush=True)
        cmd = [sys.executable,
               str(ROOT / "tests/workflows/entangle_density.py"),
               "--want", str(args.want), "--shells", mix,
               "--dims", str(args.dims), "--density", str(args.density),
               "--rounds", str(args.rounds), "--tol-pct", str(args.tol_pct)]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        rows, delivered, every = parse(r.stdout)
        results.append((name, mix, rows, delivered, every))
        if rows:
            for sh in sorted(rows):
                a, d = rows[sh]
                print(f"    shell {sh}: asked {a:.2f}, delivered {d:.2f}")
        else:
            print("    no table produced; the run did not finish")
            tail = [l for l in r.stdout.splitlines()[-6:] if l.strip()]
            for l in tail:
                print(f"      {l[:90]}")

    print(f"\n\n{'=' * 70}\n  SUMMARY\n{'=' * 70}")
    print(f"\n  {'case':<24} {'shell':>6} {'asked':>7} {'delivered':>10} "
          f"{'error':>7}")
    worst_overall = 0.0
    for name, _mix, rows, delivered, every in results:
        if not rows:
            print(f"  {name:<24} {'-':>6} {'-':>7} {'did not finish':>10}")
            continue
        first = True
        worst = 0.0
        for sh in sorted(rows):
            a, d = rows[sh]
            worst = max(worst, abs(d - a))
            print(f"  {name if first else '':<24} {sh:>6} {a:>7.2f} "
                  f"{d:>10.2f} {d - a:>+7.2f}")
            first = False
        worst_overall = max(worst_overall, worst)
        print(f"  {'':<24} {'':>6} {'':>7} {'worst':>10} {worst:>7.2f}"
              + ("   all shells in tolerance" if every else ""))
    print(f"\n  largest single-shell error across every case: "
          f"{worst_overall:.2f}")
    print("  measured by Z1+ on the LAMMPS-minimised system, per case")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
