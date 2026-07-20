"""How topon prunes a lattice to a target degree distribution.

Drives the REAL strict-sculpting generator, replays its `move_history` (the exact
sequence of edge removals it recorded), and captures the graph at intervals.
Nodes are coloured by current degree so the distribution is visible in the network
itself; a companion degree histogram is rendered per frame.

This module writes per-frame LAMMPS data files (graph -> atoms+bonds) so the
network panel renders in the same OVITO Tachyon style as the rest of the gallery.
"""
import collections
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

# degree -> colour (deep blue = fully connected 6, warm = pruned/dangling)
DEG_COLOUR = {
    6: (0.16, 0.28, 0.52), 5: (0.20, 0.46, 0.70), 4: (0.10, 0.62, 0.55),
    3: (0.65, 0.68, 0.20), 2: (0.86, 0.52, 0.18), 1: (0.78, 0.25, 0.20),
    0: (0.80, 0.80, 0.83),
}


def sculpt_history(size="6x6x6", lattice="SC", target="e:420", max_func=6, trials=80):
    """Run the generator, return (base_graph, move_history)."""
    from topon.topology.generator_python import PythonTopologyGenerator
    dims = tuple(int(x) for x in size.split("x"))
    cfg = SimpleNamespace(lattice_size=size, lattice_type=lattice,
                          max_functionality=max_func, degree_distribution=target)
    gen = PythonTopologyGenerator(cfg)
    base = gen._create_lattice(dims, lattice)
    gs = gen.generate(trials=trials, max_saves=1, time_limit=30)
    if not gs:
        raise RuntimeError(f"generator failed for {size} {lattice} {target}")
    return base, dims, gs[0].graph.get("move_history", [])


def capture(base, move_history, n_frames=32):
    """Replay removals; return a list of edge-sets sampled evenly, plus the
    removal index of each captured frame (for the histogram/label)."""
    total = len(move_history)
    idxs = sorted(set(int(round(i)) for i in np.linspace(0, total, n_frames)))
    g = base.copy()
    frames, at = [], []
    ptr = 0
    for stop in idxs:
        while ptr < stop:
            u, v = move_history[ptr]["edge"]
            if g.has_edge(u, v):
                g.remove_edge(u, v)
            ptr += 1
        frames.append(list(g.edges()))
        at.append(stop)
    return frames, at


def write_data(path, dims, pos, edges, degrees):
    """Minimal LAMMPS data: nodes as atoms (type = degree+1), edges as bonds, in a
    periodic box sized to the lattice so wrapping edges draw as short stubs."""
    nx_, ny_, nz_ = dims
    box = [nx_, ny_, nz_]
    N = len(pos)
    lines = ["graph frame\n", f"{N} atoms", f"{len(edges)} bonds",
             "7 atom types", "1 bond types",
             f"0 {box[0]} xlo xhi", f"0 {box[1]} ylo yhi", f"0 {box[2]} zlo zhi",
             "", "Masses", ""]
    lines += [f"{t} 1.0" for t in range(1, 8)]
    lines += ["", "Atoms # full", ""]
    for i in range(N):
        x, y, z = pos[i]
        t = min(max(degrees[i], 0), 6) + 1        # type 1..7 = degree 0..6
        lines.append(f"{i+1} 1 {t} 0 {x:.4f} {y:.4f} {z:.4f}")
    lines += ["", "Bonds", ""]
    for k, (u, v) in enumerate(edges):
        lines.append(f"{k+1} 1 {u+1} {v+1}")
    Path(path).write_text("\n".join(lines) + "\n")


def build_all(outdir, size="6x6x6", lattice="SC", target="e:420", max_func=6,
              n_frames=30):
    """TOPON-python stage: generate, capture, write f{k}.data + meta.json.

    The render stage (OVITO python) needs no networkx -- each node's degree is
    encoded as its atom type (type = degree + 1), and base/target degree counts
    are saved to meta.json."""
    import json
    import networkx as nx
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    base, dims, mh = sculpt_history(size, lattice, target, max_func)
    total = len(mh)
    frames, at = capture(base, mh, n_frames=n_frames)
    pos = [base.nodes[i]["pos"] for i in range(base.number_of_nodes())]
    N = len(pos)
    base_counts = collections.Counter(dict(base.degree()).values())
    gf = nx.Graph(); gf.add_nodes_from(range(N)); gf.add_edges_from(frames[-1])
    target_counts = collections.Counter(dict(gf.degree()).values())
    for k, edges in enumerate(frames):
        g = nx.Graph(); g.add_nodes_from(range(N)); g.add_edges_from(edges)
        deg = [g.degree(i) for i in range(N)]
        write_data(outdir / f"f{k:03d}.data", dims, pos, edges, deg)
    (outdir / "meta.json").write_text(json.dumps({
        "dims": list(dims), "n_nodes": N, "n_frames": len(frames), "total": total,
        "at": at, "base_edges": base.number_of_edges(),
        "final_edges": len(frames[-1]),
        "base_counts": {int(k): int(v) for k, v in base_counts.items()},
        "target_counts": {int(k): int(v) for k, v in target_counts.items()},
        "size": size, "lattice": lattice, "target": target,
    }, indent=1))
    print(f"[sculpt_frames] {N} nodes, {base.number_of_edges()} -> "
          f"{len(frames[-1])} edges, {len(frames)} frames -> {outdir}")


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "C:/tmp/sculpt_frames"
    kw = {}
    for a in sys.argv[2:]:
        k, v = a.split("=")
        kw[k] = int(v) if k == "n_frames" else v
    build_all(out, **kw)
