"""Step 1 of the MARTINI resilin example — BFM topology sweep.

Generates two crosslinked network topologies under the **natpro** sequence
setting (GGRPSDSYGAPGGGN x 18 = 270 residues per chain, segs_per_block=2)
at two system sizes:

    50  chains   (900 TYR sites,  max 450 dityrosine crosslinks)
    100 chains   (1800 TYR sites, max 900 dityrosine crosslinks)

For each system the BFM Monte Carlo loop runs until conversion sweeps
past the user's target of 25 % (= a quarter of all TYR sites reacted).
Snapshots are emitted at gel point + ~10 intermediate conversion values
(every 0.025 in conversion), so you can pick the right checkpoint for
the MARTINI build in step 2.

Outputs (under topologies/):
    {label}_gel_point.json
    {label}_post_gel_1.json, _2.json, ...
    summary.csv     one row per snapshot
    summary.png     crosslink-vs-conversion plot with the 25% target line

The same `seed=42` is used for both system sizes so the runs are
reproducible. Bump `equil_steps` if you want better-equilibrated chains
before crosslinking; the default 50 000 is enough at this scale.

Usage:
    python examples/workflows/martini_resilin/01_topology_sweep.py
"""
from __future__ import annotations

import csv
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no GUI required
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from topon.protein_network.bfm import generate_topology
from topon.protein_network.charmm.topology_io import save_topology


# ---- inputs ------------------------------------------------------------

OUTPUT_ROOT = Path(__file__).parent / "topologies"

# natpro sequence: 15-aa resilin repeat with Y (TYR) at position 8 (0-indexed: 7)
BLOCK_SEQ = "GGRPSDSYGAPGGGN"
N_REPEATS = 18                      # 18 x 15 = 270 residues per chain
SEGS_PER_BLOCK = 2                  # "two for BFM" — segments between Y nodes
TARGET_CONVERSION = 0.25            # 25 % of all TYR sites reacted

# Pre-gel snapshot targets — captured no matter whether gel was reached.
# Spans below + at + above the user's 25 % target so they can compare.
PRE_GEL_TARGETS = [0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]

SYSTEMS = [
    # Each system has its own BFM tuning that produces a fully-percolated
    # network whose first few post-gel snapshots land near the 25% target.
    #
    # 50-chain: segs_per_block=2, packing=0.45, equil=500k, seed=500
    #   -> gel at conv=0.202, post_gel_2 at 0.256 (fully percolated)
    # 100-chain: segs_per_block=3, packing=0.60, equil=500k, seed=1000
    #   -> gel at conv=0.198, post_gel_2 at ~0.25 (fully percolated)
    dict(label="natpro_50chain",  n_chains=50,  seed=500,
         segs_per_block=2, target_packing=0.45),
    dict(label="natpro_100chain", n_chains=100, seed=1000,
         segs_per_block=3, target_packing=0.60),
]

# BFM knobs (shared)
#
# Notes:
# - segs_per_block=2 puts a TYR-eligible "Y node" every other backbone bead,
#   matching the topro convention the user used previously.
# - target_packing=0.30 keeps the lattice tight; with 50–100 chains and
#   segs_per_block=2 (37 nodes/chain), 0.45 makes the box too sparse for
#   gel via lattice-adjacent crosslinks only.
# - equil_steps=200_000 gives the chains enough MC time to bring TYR
#   pairs into adjacent lattice sites.
# - crosslink_method="distance" with max_crosslink_distance_ang=6.0 lets
#   the algorithm pick up near-but-not-strictly-adjacent Y-Y pairs;
#   chemically still well within a dityrosine bond.
COMMON = dict(
    n_repeats=N_REPEATS,
    equil_steps=500_000,
    n_extra_snapshots=12,
    snapshot_delta_conv=0.025,
    min_intrachain_sep=2,
    crosslink_method="adjacent",
    pre_gel_conversions=PRE_GEL_TARGETS,
)


def lcc_stats(snap: dict, n_chains: int, n_y_per_chain: int) -> dict:
    """Compute largest-connected-component metrics for one snapshot.

    Builds a graph at TYR resolution:
      - one node per (chain_idx, y_idx) TYR site (whether reacted or not)
      - one edge per dityrosine reaction (TYR-TYR)
      - all TYR sites on the same chain are linked by the chain backbone
        (so a chain is one connected sub-block regardless of crosslink status)

    The "load-carrying" / active crosslinker count is the number of *reacted*
    TYR sites that fall inside the largest connected component. Below the
    gel point this can be much smaller than the raw reaction count, because
    crosslinks formed inside isolated sub-networks don't contribute to a
    percolating, load-bearing path.

    Returns:
        n_lcc_tyr           — total TYR sites in the LCC
        n_lcc_tyr_reacted   — *reacted* TYR sites in the LCC (active set)
        n_lcc_chains        — chains with >= 1 TYR in the LCC
        active_frac         — n_lcc_tyr_reacted / total TYR (the user's "active %")
        lcc_chain_frac      — n_lcc_chains / n_chains
        lcc_tyr_frac        — n_lcc_tyr / total TYR
    """
    total_tyr = n_chains * n_y_per_chain
    G = nx.Graph()
    # Backbone connectivity per chain (link all Y sites on the chain).
    for ci in range(n_chains):
        prev = None
        for yi in range(n_y_per_chain):
            node = (ci, yi)
            G.add_node(node)
            if prev is not None:
                G.add_edge(prev, node)
            prev = node
    # Crosslink edges + remember which (chain, y_idx_local) tuples are reacted.
    reacted: set = set()
    for r in snap.get("reactions") or []:
        (ci1, ni1), (ci2, ni2) = r[0], r[1]
        # ni1, ni2 are BFM node indices on the chain (not y-index). Map them
        # to a y-index by finding the position in the chain's y_positions.
        # snap["crosslinker_positions"] is sorted list of y positions on a chain.
        y_positions = snap.get("crosslinker_positions", [])
        try:
            yi1 = y_positions.index(ni1)
            yi2 = y_positions.index(ni2)
        except ValueError:
            continue
        G.add_edge((ci1, yi1), (ci2, yi2))
        reacted.add((ci1, yi1))
        reacted.add((ci2, yi2))

    if G.number_of_nodes() == 0:
        return dict(n_lcc_tyr=0, n_lcc_tyr_reacted=0, n_lcc_chains=0,
                    active_frac=0.0, lcc_chain_frac=0.0, lcc_tyr_frac=0.0)

    lcc = max(nx.connected_components(G), key=len)
    n_lcc_tyr = len(lcc)
    n_lcc_tyr_reacted = sum(1 for node in lcc if node in reacted)
    n_lcc_chains = len({ci for (ci, _yi) in lcc})

    return dict(
        n_lcc_tyr=n_lcc_tyr,
        n_lcc_tyr_reacted=n_lcc_tyr_reacted,
        n_lcc_chains=n_lcc_chains,
        active_frac=round(n_lcc_tyr_reacted / total_tyr, 4) if total_tyr else 0.0,
        lcc_chain_frac=round(n_lcc_chains / n_chains, 4),
        lcc_tyr_frac=round(n_lcc_tyr / total_tyr, 4) if total_tyr else 0.0,
    )


def gel_conv_of(topology: dict) -> float | None:
    """Return the conversion at which the system first formed a single
    connected cluster (the gel point), or None if it didn't gel."""
    for snap in topology["snapshots"]:
        if snap.get("label") == "gel_point":
            return snap["conv"]
    return None


def total_y_sites(n_chains: int) -> int:
    """For one TYR per repeat: n_chains * n_repeats."""
    return n_chains * N_REPEATS


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    series: dict[str, list[tuple[float, int]]] = {}  # label -> [(conv, n_crosslinks)]

    for sys_cfg in SYSTEMS:
        label = sys_cfg["label"]
        n_chains = sys_cfg["n_chains"]
        seed = sys_cfg["seed"]
        segs = sys_cfg["segs_per_block"]
        pack = sys_cfg["target_packing"]

        print(f"--- {label}: n_chains={n_chains}, n_repeats={N_REPEATS}, "
              f"segs_per_block={segs}, target_packing={pack}, seed={seed} ---")
        t0 = time.perf_counter()
        topo = generate_topology(
            n_chains=n_chains, seed=seed, segs_per_block=segs,
            target_packing=pack, verbose=False, **COMMON,
        )
        dt = time.perf_counter() - t0
        print(f"  generated in {dt:.1f} s; "
              f"{len(topo['snapshots'])} snapshot(s)")

        gel_conv = gel_conv_of(topo)
        if gel_conv is None:
            print(f"  WARN: {label} did not reach gel within "
                  f"{COMMON['equil_steps']} equil_steps; max conv = "
                  f"{topo['snapshots'][-1]['conv']:.3f}")
        else:
            print(f"  gel point at conv = {gel_conv:.3f}")

        # Save each snapshot as its own JSON for easy step-2 selection.
        max_y = total_y_sites(n_chains)
        n_y_per_chain = N_REPEATS
        series[label] = []
        for snap in topo["snapshots"]:
            snap_path = OUTPUT_ROOT / f"{label}_{snap['label']}.json"
            single = {"config": topo["config"], "snapshots": [snap]}
            save_topology(single, str(snap_path))

            n_xlinks = len(snap.get("reactions") or [])
            stats = lcc_stats(snap, n_chains, n_y_per_chain)

            row = {
                "system": label,
                "n_chains": n_chains,
                "segs_per_block": segs,
                "target_packing": pack,
                "seed": seed,
                "snapshot_label": snap["label"],
                "conversion": round(snap["conv"], 4),
                "n_crosslinks": n_xlinks,
                "max_crosslinks": max_y // 2,
                # active = reacted TYR in the largest connected component / total TYR
                "active_frac": stats["active_frac"],
                "lcc_chain_frac": stats["lcc_chain_frac"],
                "n_lcc_chains": stats["n_lcc_chains"],
                "n_lcc_tyr_reacted": stats["n_lcc_tyr_reacted"],
                "Nx": snap.get("Nx"),
                "gel_conv": (round(gel_conv, 4) if gel_conv is not None else "no_gel"),
                "topology_json": str(snap_path.relative_to(Path.cwd())),
            }
            rows.append(row)
            series[label].append((snap["conv"], n_xlinks, stats["active_frac"]))

        print()

    # ---- CSV ----
    csv_path = OUTPUT_ROOT / "summary.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {csv_path}")

    # ---- Plot: active fraction (load-carrying) vs. conversion ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Left: active fraction vs. conversion (the curve that matters for the user)
    ax = axes[0]
    for label, data in series.items():
        if not data:
            continue
        data_sorted = sorted(data, key=lambda t: t[0])
        c = [d[0] for d in data_sorted]
        a = [d[2] for d in data_sorted]
        ax.plot(c, a, "o-", label=label)
    ax.axhline(TARGET_CONVERSION, color="grey", ls="--",
               label=f"target active = {TARGET_CONVERSION:.0%}")
    ax.set_xlabel("conversion (fraction of TYR sites reacted)")
    ax.set_ylabel("active fraction (reacted TYR in LCC / total TYR)")
    ax.set_title("Active (load-carrying) fraction vs. conversion")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)

    # Right: total crosslinks vs. conversion (reference curve)
    ax = axes[1]
    for label, data in series.items():
        if not data:
            continue
        data_sorted = sorted(data, key=lambda t: t[0])
        c = [d[0] for d in data_sorted]
        x = [d[1] for d in data_sorted]
        ax.plot(c, x, "o-", label=label)
    ax.set_xlabel("conversion (fraction of TYR sites reacted)")
    ax.set_ylabel("number of dityrosine crosslinks")
    ax.set_title("Total crosslink count vs. conversion")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    plot_path = OUTPUT_ROOT / "summary.png"
    fig.savefig(plot_path, dpi=140)
    plt.close(fig)
    print(f"Wrote {plot_path}")

    # ---- Friendly summary at the end ----
    print()
    print("Snapshots (target: active_frac ~ 0.25 = 25% of TYR are load-carrying):")
    print(f"  {'system':<16s} {'snapshot':<22s} {'conv':>6s} {'xlinks':>7s} "
          f"{'active':>7s} {'chains_in_LCC':>14s}")
    for r in rows:
        active = float(r["active_frac"])
        marker = " <- near 25% active" if abs(active - TARGET_CONVERSION) < 0.025 else ""
        chains_lcc = f"{r['n_lcc_chains']}/{r['n_chains']}"
        print(f"  {r['system']:<16s} {r['snapshot_label']:<22s} "
              f"{r['conversion']:>6.3f} {r['n_crosslinks']:>4d}/{r['max_crosslinks']:<3d} "
              f"{active:>7.3f} {chains_lcc:>14s}{marker}")


if __name__ == "__main__":
    main()
