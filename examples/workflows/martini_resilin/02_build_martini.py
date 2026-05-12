"""Step 2 of the MARTINI resilin example — build a LAMMPS system from a
checkpoint chosen in step 1.

Reads the CSV summary that 01_topology_sweep.py produced, picks one row by
``--system`` + ``--snapshot``, then drives the full chain through
``topon.protein_network.workflow.run_protein_network``:
  BFM topology -> MARTINI 3 protein chemistry -> water packing -> LAMMPS files.

The BFM step is re-run with the same parameters + seed as step 1 (so the
topology is bit-identical), then the snapshot_label is used to pick the
desired checkpoint.

Output (under runs/<system>__<snapshot>__W<density>/):
    protein_network.data            LAMMPS data file (atom_style full)
    protein_network.in.settings     pair / bond / angle / dihedral coeffs
    protein_network.in.groups       protein / water / ions groups
    protein_network_topology.json   re-emitted topology JSON
    relaxation/
      protein_network_stage1.in     soft-push overlap removal
      protein_network_stage2.in     LJ epsilon ramp
      protein_network_stage3.in     tight CG min + NVT + NPT @ 310 K

After stage 3, `system_equilibrated.data` sits in the parent folder.

Usage:
    python 02_build_martini.py --system natpro_100chain --snapshot pre_gel_conv0250
    python 02_build_martini.py --system natpro_50chain  --snapshot pre_gel_conv0250
    python 02_build_martini.py --system natpro_100chain --snapshot gel_point
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from topon.protein_network.workflow import run_protein_network


# Must match the COMMON block in 01_topology_sweep.py — keep these in sync.
# BFM knobs that are SHARED across systems. Per-system overrides
# (segs_per_block, target_packing, seed) live in the CSV from step 1.
BFM_SHARED = dict(
    n_repeats=18,
    equil_steps=500_000,
    n_extra_snapshots=12,
    snapshot_delta_conv=0.025,
    min_intrachain_sep=2,
    crosslink_method="adjacent",
    pre_gel_conversions=[0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50],
)


# "50 % water content with regular W-bead size":
# Bulk water ~33 H2O/nm^3, regular W bead = 4 H2O/bead -> 8.25 W/nm^3 at bulk.
# 50 % of bulk -> ~4 W/nm^3 (matches the v42 "W04nm3" medium-water reference).
DEFAULT_WATER_DENSITY_W = 4.0
DEFAULT_WATER_BEAD = "W"


HERE = Path(__file__).resolve().parent
DEFAULT_CSV = HERE / "topologies/summary.csv"
DEFAULT_OUT = HERE / "runs"


def load_systems(csv_path: Path) -> list[dict]:
    if not csv_path.exists():
        raise SystemExit(
            f"Step-1 summary not found: {csv_path}\n"
            f"Run 01_topology_sweep.py first."
        )
    with csv_path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--system", required=True,
                   help="System label from step 1 (e.g. natpro_50chain, natpro_100chain)")
    p.add_argument("--snapshot", required=True,
                   help="Snapshot label (e.g. pre_gel_conv0250, gel_point, post_gel_8)")
    p.add_argument("--water-density", type=float, default=DEFAULT_WATER_DENSITY_W,
                   help=f"W beads per nm^3 (default {DEFAULT_WATER_DENSITY_W} "
                        f"= ~50%% bulk water with W bead)")
    p.add_argument("--water-bead", default=DEFAULT_WATER_BEAD,
                   choices=["W", "SW", "TW"],
                   help="MARTINI 3 water bead type")
    p.add_argument("--csv", type=Path, default=DEFAULT_CSV,
                   help="Path to summary.csv from step 1")
    p.add_argument("--output-root", type=Path, default=DEFAULT_OUT,
                   help="Where to land the per-build output folder")
    args = p.parse_args()

    rows = load_systems(args.csv)
    row = next(
        (r for r in rows
         if r["system"] == args.system and r["snapshot_label"] == args.snapshot),
        None,
    )
    if row is None:
        avail = sorted({(r["system"], r["snapshot_label"]) for r in rows})
        raise SystemExit(
            f"No row for system='{args.system}' snapshot='{args.snapshot}'.\n"
            f"Available combos:\n" + "\n".join(f"  {s}  {sl}" for s, sl in avail)
        )

    n_chains = int(row["n_chains"])
    segs_per_block = int(row["segs_per_block"])
    target_packing = float(row["target_packing"])
    seed = int(row["seed"])
    conv = float(row["conversion"])
    n_xlinks = int(row["n_crosslinks"])
    max_xlinks = int(row["max_crosslinks"])
    print(f"Selected:")
    print(f"  system            : {args.system}")
    print(f"  snapshot          : {args.snapshot}")
    print(f"  n_chains          : {n_chains}")
    print(f"  segs_per_block    : {segs_per_block}")
    print(f"  target_packing    : {target_packing}")
    print(f"  seed              : {seed}")
    print(f"  conversion        : {conv:.4f}")
    print(f"  crosslinks        : {n_xlinks}/{max_xlinks}")
    print(f"  water_density (W) : {args.water_density} per nm^3")
    print()

    out_dir = args.output_root / f"{args.system}__{args.snapshot}__W{args.water_density:g}"
    out_dir.mkdir(parents=True, exist_ok=True)

    paths = run_protein_network(
        block_seq="GGRPSDSYGAPGGGN",
        n_chains=n_chains,
        n_repeats=BFM_SHARED["n_repeats"],
        snapshot_label=args.snapshot,
        output_dir=str(out_dir),
        segs_per_block=segs_per_block,
        target_packing=target_packing,
        equil_steps=BFM_SHARED["equil_steps"],
        n_extra_snapshots=BFM_SHARED["n_extra_snapshots"],
        snapshot_delta_conv=BFM_SHARED["snapshot_delta_conv"],
        min_intrachain_sep=BFM_SHARED["min_intrachain_sep"],
        crosslink_method=BFM_SHARED["crosslink_method"],
        pre_gel_conversions=BFM_SHARED["pre_gel_conversions"],
        water_density_w_per_nm3=args.water_density,
        water_bead_type=args.water_bead,
        seed=seed,
        hierarchical_stage1=True,
        verbose=True,
    )

    print("\nArtifacts:")
    for key, path in paths.items():
        print(f"  {key:25s} {Path(path).resolve()}")

    relax_dir = Path(paths.get("stage1", "")).parent
    print(f"\nTo run LAMMPS:")
    print(f"  cd {relax_dir}")
    print(f"  lmp -in protein_network_stage1.in   # soft overlap removal")
    print(f"  lmp -in protein_network_stage2.in   # LJ epsilon ramp")
    print(f"  lmp -in protein_network_stage3.in   # CG min + brief NVT/NPT")


if __name__ == "__main__":
    main()
