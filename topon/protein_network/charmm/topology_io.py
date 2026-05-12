"""topro.utils.io — Topology file I/O."""

import json
import os
from datetime import datetime

TOPOLOGY_VERSION = "1.0"


def save_topology(topology, path):
    """
    Save a topology dict (output of generate_topology) to a JSON file.

    Parameters
    ----------
    topology : dict
        Must contain 'config' and 'snapshots' keys.
    path : str
        Output file path (.json).
    """
    data = {
        "version": TOPOLOGY_VERSION,
        "created": datetime.now().isoformat(),
        **topology,
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Topology saved: {path}  ({len(topology['snapshots'])} snapshots)")


def load_topology(path):
    """
    Load a topology JSON file.

    Returns
    -------
    dict with 'config', 'snapshots', and 'created' keys.
    """
    with open(path, "r") as f:
        data = json.load(f)

    ver = data.get("version", "0.0")
    if ver != TOPOLOGY_VERSION:
        import warnings
        warnings.warn(
            f"Topology version mismatch: file={ver}, expected={TOPOLOGY_VERSION}"
        )

    return {
        "config": data["config"],
        "snapshots": data["snapshots"],
        "created": data.get("created", "unknown"),
    }


def get_snapshot(topology, label="gel_point"):
    """
    Retrieve a specific snapshot from a loaded topology.

    Parameters
    ----------
    label : str or int
        Snapshot label (e.g. 'gel_point', 'post_gel_1') or integer index.
    """
    snapshots = topology["snapshots"]

    if isinstance(label, int):
        return snapshots[label]

    for snap in snapshots:
        if snap["label"] == label:
            return snap

    available = [s["label"] for s in snapshots]
    raise KeyError(f"Snapshot '{label}' not found. Available: {available}")


def list_snapshots(topology):
    """Print a summary table of available snapshots."""
    cfg = topology["config"]
    print(f"\nTopology created: {topology.get('created', 'unknown')}")
    print(
        f"Config : {cfg['n_chains']} chains | {cfg['n_repeats']} repeats | "
        f"segs_per_block={cfg['segs_per_block']} | "
        f"lattice {cfg['Nx']}×{cfg['Ny']}×{cfg['Nz']}"
    )
    print(f"\n{'#':<4} {'label':<22} {'conv':>7}  {'reactions':>10}")
    print("-" * 48)
    for i, snap in enumerate(topology["snapshots"]):
        print(
            f"[{i}]  {snap['label']:<22} {snap['conv']:7.4f}  "
            f"{len(snap['reactions']):>10}"
        )
