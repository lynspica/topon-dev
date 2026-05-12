"""BFM (Bond-Fluctuation Model) topology generation — standalone demo.

Generates a small crosslinked-network topology JSON with the same recipe
that drives `topon.protein_network` MARTINI builds. Useful for inspecting
the gel-point structure or feeding a snapshot into a different chemistry
builder without going through the full LAMMPS-emitting workflow.

Defaults reproduce the topro V41/V42 reference (8 chains x 18 repeats x
segs_per_block=2, packing=0.45, equil=20k, adjacent crosslinks, seed=42)
which gels at conv=0.125 — visible in the snapshot list printed at the end.

Outputs (under output/):
    bfm_topology.json   all snapshots in topro JSON format
    *.nodes / *.edges   per-snapshot dual graph (one file per chain-end +
                        crosslinker node, one file per backbone + crosslink edge)
    *.png               3D visualisation of each snapshot (chains as lines,
                        TYR sites as scatter, crosslinks bolded)

Usage:
    python examples/demos/topology/bfm/run.py
"""
from __future__ import annotations

import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from topon.protein_network.bfm import generate_topology
from topon.protein_network.charmm.topology_io import save_topology


OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


def _idx_to_xyz(idx: int, Nx: int, Ny: int, Nz: int) -> tuple[int, int, int]:
    """Convert flat lattice index to (x, y, z) (matches BFM's row-major)."""
    z = idx // (Nx * Ny)
    rem = idx - z * Nx * Ny
    y = rem // Nx
    x = rem - y * Nx
    return x, y, z


def export_snapshot(snap: dict, cfg: dict, basename: Path) -> None:
    """Emit .nodes, .edges, and a 3D PNG visualisation for one snapshot."""
    Nx, Ny, Nz = snap["Nx"], snap["Ny"], snap["Nz"]
    chains = snap["chains"]
    n_chains = len(chains)
    y_positions = snap["crosslinker_positions"]

    # ---- .nodes / .edges ----
    # Each TYR site = one node. Backbone path = chain of backbone edges
    # between consecutive Y positions on a chain.
    # ID convention: (chain_idx * n_y_per_chain + y_local_idx)
    n_y_per_chain = len(y_positions)
    nodes_lines = ["# NodeID X Y Z Degree IsCrosslinker"]
    edges_lines = []

    for ci in range(n_chains):
        chain = chains[ci]
        for yi, y_node_idx in enumerate(y_positions):
            nid = ci * n_y_per_chain + yi
            lattice = chain[y_node_idx]
            x, y, z = _idx_to_xyz(lattice, Nx, Ny, Nz)
            deg = 1 if (yi == 0 or yi == n_y_per_chain - 1) else 2  # along-chain
            nodes_lines.append(f"{nid} {x} {y} {z} {deg} 1")
            if yi > 0:
                prev_nid = ci * n_y_per_chain + (yi - 1)
                edges_lines.append(f"{prev_nid} {nid}  # backbone")

    # Crosslink edges (between two different chain's TYR sites or same-chain ones)
    for (a, b) in snap.get("reactions", []):
        # 'reactions' format: [[ci1, y_node_idx1], [ci2, y_node_idx2]] where
        # y_node_idx is the BFM node index (not the local y-index). Map to local.
        ni1 = y_positions.index(a[1])
        ni2 = y_positions.index(b[1])
        u = a[0] * n_y_per_chain + ni1
        v = b[0] * n_y_per_chain + ni2
        edges_lines.append(f"{u} {v}  # crosslink")

    (basename.with_suffix(".nodes")).write_text("\n".join(nodes_lines) + "\n")
    (basename.with_suffix(".edges")).write_text(
        "# u v  # type\n" + "\n".join(edges_lines) + "\n"
    )

    # ---- 3D PNG ----
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")
    colours = plt.get_cmap("tab10", n_chains)

    # Chain backbones (line through Y sites)
    for ci in range(n_chains):
        chain = chains[ci]
        xs, ys, zs = [], [], []
        # Walk the full BFM chain to keep it visually continuous.
        for node_idx in chain:
            x, y, z = _idx_to_xyz(node_idx, Nx, Ny, Nz)
            xs.append(x); ys.append(y); zs.append(z)
        ax.plot(xs, ys, zs, color=colours(ci), alpha=0.6, linewidth=1.0)
        # Mark TYR sites
        for yi in y_positions:
            x, y, z = _idx_to_xyz(chain[yi], Nx, Ny, Nz)
            ax.scatter([x], [y], [z], color=colours(ci), s=18, edgecolor="k", linewidth=0.4)

    # Crosslinks (red lines)
    for (a, b) in snap.get("reactions", []):
        ci1, ni1 = a; ci2, ni2 = b
        x1, y1, z1 = _idx_to_xyz(chains[ci1][ni1], Nx, Ny, Nz)
        x2, y2, z2 = _idx_to_xyz(chains[ci2][ni2], Nx, Ny, Nz)
        ax.plot([x1, x2], [y1, y2], [z1, z2],
                color="red", linewidth=1.5, alpha=0.8)

    n_xlinks = len(snap.get("reactions") or [])
    ax.set_title(f"{snap['label']}  (conv={snap['conv']:.3f}, "
                 f"{n_xlinks} crosslinks, {n_chains} chains on {Nx}^3 lattice)")
    ax.set_xlim(0, Nx); ax.set_ylim(0, Ny); ax.set_zlim(0, Nz)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    fig.tight_layout()
    fig.savefig(basename.with_suffix(".png"), dpi=130)
    plt.close(fig)


def main() -> None:
    print("--- BFM topology generation (topro V41/V42 reference recipe) ---")
    t0 = time.perf_counter()
    topology = generate_topology(
        n_chains=8,
        n_repeats=18,
        segs_per_block=2,
        target_packing=0.45,
        equil_steps=20_000,
        n_extra_snapshots=4,
        snapshot_delta_conv=0.05,
        min_intrachain_sep=2,
        seed=42,
        verbose=False,
    )
    dt = time.perf_counter() - t0
    print(f"  generated in {dt:.2f} s")

    out_path = OUTPUT_DIR / "bfm_topology.json"
    save_topology(topology, str(out_path))
    print(f"  saved: {out_path}")

    # Per-snapshot dual-graph + visualization
    print("\nEmitting per-snapshot .nodes / .edges / .png ...")
    cfg = topology["config"]
    for snap in topology["snapshots"]:
        base = OUTPUT_DIR / f"snapshot_{snap['label']}"
        export_snapshot(snap, cfg, base)
        print(f"  {snap['label']:<20s} -> {base.name}.{{nodes,edges,png}}")

    cfg = topology["config"]
    print(f"\nConfig:")
    print(f"  n_chains={cfg['n_chains']}, n_repeats={cfg['n_repeats']}, "
          f"segs_per_block={cfg['segs_per_block']}")
    print(f"  Nx={cfg.get('Nx', '?')}^3, "
          f"packing={cfg.get('actual_packing', cfg.get('target_packing')):.3f}, "
          f"equil_steps={cfg['equil_steps']}")
    print(f"  crosslink_method={cfg.get('crosslink_method', 'adjacent')}")

    print(f"\nSnapshots ({len(topology['snapshots'])}):")
    total_y = cfg["n_chains"] * cfg["n_repeats"]
    print(f"  total TYR sites: {total_y}; max crosslinks: {total_y // 2}")
    print(f"  {'label':<20s} {'conv':>6s} {'n_xlinks':>10s}")
    for snap in topology["snapshots"]:
        n_xlinks = len(snap.get("reactions") or [])
        print(f"  {snap['label']:<20s} {snap['conv']:>6.3f} {n_xlinks:>10d}")

    print()
    print("Feed this JSON into a chemistry build, e.g.:")
    print("  topon> topro generate --block-seq GGRPSDSYGAPGGGN \\")
    print("          --n-chains 8 --n-repeats 18 --water-density 4 \\")
    print("          --snapshot-label gel_point --output runs/resilin_demo")
    print()
    print("Or load the JSON directly in Python:")
    print("  from topon.protein_network.charmm.topology_io import load_topology")
    print(f'  topo = load_topology("{out_path.relative_to(Path.cwd())}")')


if __name__ == "__main__":
    main()
