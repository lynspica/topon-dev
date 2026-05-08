"""JSON I/O for topology snapshots.

Forked from `legacy/subprojects/protein_network/topro/topro/utils/io.py`. The
on-disk envelope is byte-compatible with topro's topo_*.json files, so a topro
JSON loads here unchanged and vice versa.
"""
from __future__ import annotations

import json
import os
import warnings
from datetime import datetime

TOPOLOGY_VERSION = "1.0"


def save_topology(topology: dict, path: str) -> None:
    """Write a topology dict (output of `bfm.generate_topology`) to JSON."""
    data = {
        "version": TOPOLOGY_VERSION,
        "created": datetime.now().isoformat(),
        **topology,
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Topology saved: {path}  ({len(topology['snapshots'])} snapshots)")


def load_topology(path: str) -> dict:
    """Read a topology JSON file. Warns on version mismatch."""
    with open(path, "r") as f:
        data = json.load(f)
    ver = data.get("version", "0.0")
    if ver != TOPOLOGY_VERSION:
        warnings.warn(
            f"Topology version mismatch: file={ver}, expected={TOPOLOGY_VERSION}"
        )
    return {
        "config": data["config"],
        "snapshots": data["snapshots"],
        "created": data.get("created", "unknown"),
    }


def get_snapshot(topology: dict, label: str | int = "gel_point") -> dict:
    """Return one snapshot by label string or integer index."""
    snapshots = topology["snapshots"]
    if isinstance(label, int):
        return snapshots[label]
    for snap in snapshots:
        if snap["label"] == label:
            return snap
    available = [s["label"] for s in snapshots]
    raise KeyError(f"Snapshot {label!r} not found. Available: {available}")


def list_snapshots(topology: dict) -> None:
    """Print a summary table of available snapshots."""
    cfg = topology["config"]
    print(f"\nTopology created: {topology.get('created', 'unknown')}")
    print(
        f"Config : {cfg['n_chains']} chains | {cfg['n_repeats']} repeats | "
        f"segs_per_block={cfg['segs_per_block']} | "
        f"lattice {cfg['Nx']}x{cfg['Ny']}x{cfg['Nz']}"
    )
    print(f"\n{'#':<4} {'label':<22} {'conv':>7}  {'reactions':>10}")
    print("-" * 48)
    for i, snap in enumerate(topology["snapshots"]):
        print(
            f"[{i}]  {snap['label']:<22} {snap['conv']:7.4f}  "
            f"{len(snap['reactions']):>10}"
        )
