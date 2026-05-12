"""Sweep BFM topology parameters and record the gel-point conversion
for each configuration.

For a BFM crosslinked network, the *gel point* is the smallest crosslink
conversion at which all chains end up in a single connected cluster
(determined via Union-Find inside `topon.protein_network.bfm`). Where
the gel point lands depends primarily on:
  - n_chains          how many chains in the cell
  - n_repeats         chain length (Y-positions per chain)
  - segs_per_block    inter-Y spacing (controls intra-chain reach)
  - min_intrachain_sep how far apart two Y nodes on the same chain must
                      be to be eligible for a crosslink

This script sweeps each of those one at a time around a sensible default
and saves both per-run topology JSONs and a single summary CSV showing
the gel-point conversion (or 'no_gel' if the sweep didn't reach percolation
within the allotted equilibration). Reading the CSV side-by-side with
`docs/...` is the fastest way to develop intuition before committing to
a multi-hour CHARMM/MARTINI build.

Usage:
    python examples/workflows/bfm_gel_point_sweep/run.py
"""
from __future__ import annotations

import csv
import time
from pathlib import Path

from topon.protein_network.bfm import generate_topology
from topon.protein_network.charmm.topology_io import save_topology


OUTPUT_ROOT = Path(__file__).parent / "output"
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

# Sensible default that we vary one parameter at a time around
DEFAULT = dict(
    n_chains=8,
    n_repeats=8,
    segs_per_block=2,
    target_packing=0.45,
    equil_steps=10_000,
    n_extra_snapshots=2,
    snapshot_delta_conv=0.05,
    min_intrachain_sep=2,
    seed=42,
    verbose=False,
)

# A sweep is a list of dicts overriding DEFAULT. Add or remove rows freely.
SWEEPS = [
    {"label": "baseline"},
    {"label": "more_chains_12",        "n_chains": 12},
    {"label": "more_chains_16",        "n_chains": 16},
    {"label": "longer_chains_12",      "n_repeats": 12},
    {"label": "longer_chains_16",      "n_repeats": 16},
    {"label": "segs_per_block_3",      "segs_per_block": 3},
    {"label": "tight_intrachain_sep_1","min_intrachain_sep": 1},
    {"label": "loose_intrachain_sep_3","min_intrachain_sep": 3},
    {"label": "more_equilibration",    "equil_steps": 50_000},
]


def gel_point_of(topology: dict) -> tuple[str, float]:
    """Pull the gel snapshot's conversion out of a topology dict.

    Returns (label, conv). If no gel snapshot exists, returns the highest
    conversion reached (label='no_gel').
    """
    for snap in topology["snapshots"]:
        if snap.get("label") == "gel_point":
            return "gel_point", snap["conv"]
    # No gel reached — return the max conv we saw
    if topology["snapshots"]:
        worst = max(topology["snapshots"], key=lambda s: s["conv"])
        return worst.get("label", "no_gel"), worst["conv"]
    return "empty", 0.0


def main() -> None:
    rows: list[dict] = []
    for sweep in SWEEPS:
        label = sweep.pop("label")
        params = {**DEFAULT, **sweep}
        t0 = time.perf_counter()
        topo = generate_topology(**params)
        dt = time.perf_counter() - t0

        snap_label, gel_conv = gel_point_of(topo)
        n_snapshots = len(topo["snapshots"])

        topo_path = OUTPUT_ROOT / f"{label}.json"
        save_topology(topo, str(topo_path))

        row = {
            "label": label,
            "n_chains": params["n_chains"],
            "n_repeats": params["n_repeats"],
            "segs_per_block": params["segs_per_block"],
            "min_intrachain_sep": params["min_intrachain_sep"],
            "equil_steps": params["equil_steps"],
            "snapshot_label": snap_label,
            "gel_conv": round(gel_conv, 4),
            "n_snapshots": n_snapshots,
            "wall_s": round(dt, 2),
        }
        rows.append(row)
        print(f"  {label:30s}  {snap_label:10s}  conv={gel_conv:.4f}  ({dt:.1f}s)")

    csv_path = OUTPUT_ROOT / "summary.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {csv_path.relative_to(Path.cwd())} + {len(rows)} topology JSONs")
    print("Open `output/*.json` with `topon.protein_network.charmm.topology_io.load_topology`")
    print("then feed any of them into the bfm_to_martini or bfm_to_charmm workflows.")


if __name__ == "__main__":
    main()
