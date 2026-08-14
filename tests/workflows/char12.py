"""Twelve plain networks, equilibrated, and their characteristic ratios.

    python tests/workflows/char12.py build
    (LAMMPS runs driven by the shell between the two phases)
    python tests/workflows/char12.py analyze

Six lattices -- SC, BCC, FCC and three SC:BCC:FCC mixtures -- at DP 40 and
DP 80, dims 5, no entanglements, built at density 0.30 because that is the
density the reference protocol (`equil_demo.in`) equilibrates at.

The protocol file is copied verbatim and driven through its own `-var` knobs:
soft push-off, harmonic pre-min, FENE, then NVT at T=1. Two departures from
the file as shipped, both stated rather than silent:

  * the data files are post-processed before LAMMPS sees them: the Angles
    sections are stripped (the protocol defines no angle_style and read_data
    would refuse the file), the placeholder box is replaced by the real one,
    and coordinates are wrapped into it. Bond styles apply the minimum image,
    so wrapped coordinates with bonds across the boundary are correct.
  * RUN is shortened from 800000 to fit twelve systems in the session; the
    actual number of steps is recorded on every plot.

Analysis is `characteristic_ratio.py` imported from its own location, not a
copy, so the numbers are computed by exactly the user's code.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

REF = Path("E:/PhD/topology_datasets")
EQUIL = REF / "demo_downturn" / "equil_demo.in"
OUTDIR = ROOT / "tests/output/char12"

CASES = [
    ("SC",     "SC",  None),
    ("BCC",    "BCC", None),
    ("FCC",    "FCC", None),
    ("MIX801010", "MIX", (0.8, 0.1, 0.1)),
    ("MIXeven",   "MIX", (0.34, 0.33, 0.33)),
    ("MIX101080", "MIX", (0.1, 0.1, 0.8)),
]
DPS = (40, 80)


def fix_data(src, dst, L):
    """Strip Angles, set the real box, wrap the coordinates.

    The build writes a placeholder box of -50..50; the real cell is ``L``.
    Coordinates land wherever the walks wandered, so they are wrapped in.
    Bond styles apply the minimum image, so a bond across the boundary is
    computed correctly from wrapped coordinates with zero image flags.
    """
    lines = Path(src).read_text().splitlines()
    out, section, skip = [], None, False
    for ln in lines:
        s = ln.strip()
        if s.endswith("angles") and s.split()[0].isdigit():
            out.append("0 angles")
            continue
        if s.endswith("angle types") and s.split()[0].isdigit():
            out.append("0 angle types")
            continue
        if s.endswith("xlo xhi"):
            out.append(f"0.0 {L[0]:.6f} xlo xhi")
            continue
        if s.endswith("ylo yhi"):
            out.append(f"0.0 {L[1]:.6f} ylo yhi")
            continue
        if s.endswith("zlo zhi"):
            out.append(f"0.0 {L[2]:.6f} zlo zhi")
            continue
        if s and s[0].isalpha():
            section = s.split()[0] if not s.startswith("Angle") else "Angle"
            skip = section == "Angle"
            if skip:
                continue
            out.append(ln)
            continue
        if skip:
            continue
        if section == "Atoms" and s and not s.startswith("#"):
            p = s.split()
            if len(p) >= 7:
                xyz = np.array([float(p[4]), float(p[5]), float(p[6])])
                xyz -= np.asarray(L) * np.floor(xyz / np.asarray(L))
                p[4], p[5], p[6] = (f"{xyz[0]:.6f}", f"{xyz[1]:.6f}",
                                    f"{xyz[2]:.6f}")
                out.append(" ".join(p))
                continue
        out.append(ln)
    Path(dst).write_text("\n".join(out) + "\n")


def build_all(dims, density, seed):
    from tests.workflows.entangle_all import CASES as LATCASES
    from tests.workflows.entangle_relaxed import rewrite_coords
    from tests.workflows.entangle_steps import (BOND, LATTICE, build_network,
                                                chain_ids, geometry,
                                                write_system)
    from topon.conformation.paths import bridging_walk

    OUTDIR.mkdir(parents=True, exist_ok=True)
    manifest = []
    for tag, lat, mix in CASES:
        for dp in DPS:
            name = f"{tag}_dp{dp}"
            run = OUTDIR / name
            if (run / "net.data").exists():
                print(f"  {name:16s} already built, kept", flush=True)
                manifest.append(name)
                continue
            spec = dict(LATTICE)
            spec.update(LATCASES[lat])
            spec["dims"] = (dims,) * 3
            if mix is not None:
                spec["mix"] = {"SC": mix[0], "BCC": mix[1], "FCC": mix[2]}
                spec["lattice"] = "MIX"
                # Full site functionality for custom mixes.
                # The MIX case's degree list demands degree 7 and 8
                # nodes, and a mix with few BCC/FCC sites has almost
                # none that can carry them: the rejection sampler then
                # never converges, and both this build and the lattice
                # catalogue hung on it for half an hour or more.
                spec["degree_dist"] = "0:0,1:0"
            graph = build_network(spec)
            geo = geometry(graph, dp=dp, density=density)
            rng = np.random.default_rng(seed)
            paths = {k: bridging_walk(c0, c1, dp + 1, BOND, rng)
                     for k, (c0, c1) in geo["chords"].items()}
            shutil.rmtree(run, ignore_errors=True)
            _n, node_atom, chain_atoms = write_system(graph, geo, paths, run)
            src = run / "01_Topology" / "system.data"
            if not src.exists():
                src = next(run.rglob("system.data"))

            # write_system writes zero coordinates plus displacement files;
            # the conform stage this build skips is what applies them. Every
            # atom in the raw file sits at the origin -- 8805 atoms in one
            # neighbour bin, which is the neighbour-list overflow LAMMPS
            # refused with. So the real coordinates go in here, from the
            # paths, the way the shell gallery always did. The tiny noise is
            # the same degeneracy-breaker the conform stage applies.
            xyz = {}
            for k in sorted(geo["chords"]):
                for aid, q in zip(chain_ids(k, node_atom, chain_atoms,
                                            geo["ends"]), paths[k]):
                    xyz[aid] = q
            noise = np.random.default_rng(seed).normal(0.0, 1e-4,
                                                       (len(xyz), 3))
            xyz = {a: q + noise[i] for i, (a, q) in
                   enumerate(sorted(xyz.items()))}
            rewrite_coords(src, src, xyz)
            fix_data(src, run / "net.data", geo["L"])
            shutil.copy2(EQUIL, run / "equil_demo.in")
            n = graph.number_of_edges()
            beads = n * dp + graph.number_of_nodes()
            print(f"  {name:16s} {n:5d} chains {beads:7d} beads "
                  f"box {geo['L'][0]:.1f}")
            manifest.append(name)
    (OUTDIR / "manifest.txt").write_text("\n".join(manifest) + "\n")


def analyze(run_steps):
    sys.path.insert(0, str(REF))
    from characteristic_ratio import char_ratio

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = (OUTDIR / "manifest.txt").read_text().split()
    colors = {t: c for (t, _l, _m), c in zip(
        CASES, ["#2C3E50", "#B03A2E", "#3A76B0", "#0e7c6b", "#8E6C3A",
                "#7B4FA0"])}

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9))
    rows = []
    for name in names:
        eq = OUTDIR / name / "eq.data"
        if not eq.exists():
            rows.append((name, None))
            print(f"  {name}: no eq.data, skipped")
            continue
        r = char_ratio(str(eq), bond=0.965)
        rows.append((name, r))
        tag, dp = name.rsplit("_dp", 1)
        row = 0 if dp == "40" else 1
        axes[row][0].plot(r["s"], r["prof"], "-", lw=1.4, ms=2,
                          color=colors[tag], label=tag)
        axes[row][1].plot(r["s"], r["C"], "-", lw=1.4, ms=2,
                          color=colors[tag], label=tag)
        # Per-system plot as the reference script would make it.
        f2, a2 = plt.subplots(1, 2, figsize=(11, 4.2))
        a2[0].plot(r["s"], r["prof"], "o-", ms=3)
        a2[0].set_xlabel("contour separation s")
        a2[0].set_ylabel(r"$\langle R^2(s)\rangle/s$")
        a2[1].plot(r["s"], r["C"], "s-", ms=3, color="tab:green")
        a2[1].set_xlabel("contour separation s")
        a2[1].set_ylabel("C(s)")
        f2.suptitle(f"{name}  ({run_steps} NVT steps)")
        f2.tight_layout()
        f2.savefig(OUTDIR / f"char_{name}.png", dpi=120)
        plt.close(f2)

    for row, dp in ((0, 40), (1, 80)):
        axes[row][0].set_ylabel(r"$\langle R^2(s)\rangle/s$"
                                + f"   (DP {dp})")
        axes[row][1].set_ylabel(f"C(s)   (DP {dp})")
        for col in (0, 1):
            axes[row][col].set_xlabel("contour separation s")
            axes[row][col].legend(fontsize=8, ncol=2)
    axes[0][0].set_title("internal distances (flat = ideal)")
    axes[0][1].set_title("characteristic ratio")
    fig.suptitle(f"12 networks, dims 5, density 0.30, no entanglements, "
                 f"{run_steps} NVT steps of the reference protocol")
    fig.tight_layout()
    fig.savefig(OUTDIR / "char12_summary.png", dpi=130)

    print(f"\n  {'system':16s} {'strands':>8} {'N':>4} {'Ree2/Rg2':>9} "
          f"{'C_app':>6} {'downturn':>9}")
    for name, r in rows:
        if r is None:
            print(f"  {name:16s} {'skipped':>8}")
            continue
        print(f"  {name:16s} {r['n_strands']:8d} {r['N']:4d} "
              f"{r['ree2rg2']:9.3f} {r['Cinf']:6.2f} {r['downturn']:9.3f}"
              + ("   FLAT/OK" if r["downturn"] < 0.05 else "   COMPRESSED"))
    print(f"\n  {OUTDIR / 'char12_summary.png'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("phase", choices=("build", "analyze"))
    ap.add_argument("--dims", type=int, default=5)
    ap.add_argument("--density", type=float, default=0.30)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--run-steps", type=int, default=150000,
                    help="recorded on the plots; the shell passes the same "
                         "number to LAMMPS")
    args = ap.parse_args()
    if args.phase == "build":
        build_all(args.dims, args.density, args.seed)
    else:
        analyze(args.run_steps)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
