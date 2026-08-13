"""One system per shell, so the shells can be looked at rather than read.

    python tests/workflows/entangle_shellshow.py
    python tests/workflows/entangle_shellshow.py --shells 1 2 3 4 --count 2

Writes a LAMMPS data file per shell, each showing the same routed chain wound
around a partner from that shell and nothing else designed. Every other chain
is drawn straight, so the picture is the design and not the melt.

Atom types carry the meaning, so OVITO colours them without any setup:

    type 1   every other chain, straight
    type 2   the routed chain
    type 3   its partner, the one it is wound around
    type 4   the crosslink junctions

A shell is every chain at the same closest approach to the routed one. On a
simple cubic lattice those distances are discrete, so shell 1 is the chains
sharing a crosslink at distance zero, shell 2 the chains one lattice unit away,
shell 3 sqrt(2), shell 4 sqrt(3). What the files show is how much further the
routed chain has to travel as the shell number goes up.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.assignment.entanglements import (  # noqa: E402
    chain_distances,
    neighbour_shells,
)
from topon.conformation.paths import Clearance, straight  # noqa: E402
from tests.workflows.entangle_all import CASES  # noqa: E402

from tests.workflows.entangle_relaxed import construct  # noqa: E402
from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    DP,
    LATTICE,
    OUT,
    build_network,
    chain_ids,
    geometry,
    write_system,
)


def chain_dir(graph, e):
    """Unit vector along a chain, minimum image.

    Not "which lattice axis": on BCC and FCC the chains run along body and
    face diagonals, so an axis test returns nothing for them and a request for
    parallel partners silently finds none. Two chains are parallel when their
    directions are aligned to within a tolerance, whatever direction that is.
    """
    u, v = e
    box = np.asarray(graph.graph["box"], float)
    d = (np.asarray(graph.nodes[v]["pos"], float)
         - np.asarray(graph.nodes[u]["pos"], float))
    d -= box * np.round(d / box)
    n = float(np.linalg.norm(d))
    return d / n if n > 1e-12 else np.array([1.0, 0.0, 0.0])


def retype(src, dst, types):
    """Copy a data file, setting each atom's type from ``types``.

    The type column is what OVITO colours by, so this is the whole labelling
    mechanism: no modifiers, no selection expressions, just open the file.

    Adding types means every per-type section has to grow with them. The
    header count, Masses and Pair Coeffs all have one entry per type, and a
    file whose header promises four types while Masses lists two is rejected:
    OVITO reads the "Pair Coeffs" header as the third mass and reports an
    invalid mass specification. The mass lines also carry a trailing comment,
    so the entry is not simply two fields.
    """
    lines = Path(src).read_text().splitlines()
    n_types = max(types.values())

    out, section, i = [], None, 0
    while i < len(lines):
        line = lines[i]
        s = line.strip()

        if s.endswith("atom types"):
            out.append(f"{n_types} atom types")
            i += 1
            continue

        if s and s[0].isalpha():
            section = s.split()[0] if not s.startswith("Pair") else "Pair"
            out.append(line)
            # Replay a per-type section with one entry per type, copying the
            # first entry's values so every type is physically identical.
            if section in ("Masses", "Pair"):
                i += 1
                while i < len(lines) and not lines[i].strip():
                    out.append(lines[i])
                    i += 1
                body = []
                while i < len(lines) and lines[i].strip() and                         not lines[i].strip()[0].isalpha():
                    body.append(lines[i].strip())
                    i += 1
                if body:
                    first = body[0].split("#")[0].split()
                    for t in range(1, n_types + 1):
                        out.append(" ".join([str(t)] + first[1:]))
                continue
            i += 1
            continue

        if section == "Atoms" and s and not s.startswith("#"):
            p = s.split()
            if len(p) >= 7:
                p[2] = str(types.get(int(p[0]), 1))
                out.append(" ".join(p))
                i += 1
                continue

        out.append(line)
        i += 1

    Path(dst).write_text("\n".join(out) + "\n")
    return dst


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shells", type=int, nargs="*", default=[1, 2, 3, 4])
    ap.add_argument("--count", type=int, default=2,
                    help="entanglements to build on the one designed pair")
    ap.add_argument("--dims", type=int, default=4)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--ring", type=float, default=2.0)
    ap.add_argument("--coil", type=float, default=4.0,
                    help="sets the box. Lower is roomier, which makes the "
                         "designed winding easier to see.")
    ap.add_argument("--lattice", default="SC",
                    choices=("SC", "BCC", "FCC", "MIX"))
    ap.add_argument("--orient", default="any",
                    choices=("any", "parallel", "perpendicular"),
                    help="a shell is a distance, and at a given distance "
                         "there are both parallel and perpendicular "
                         "neighbours -- shell 2 on SC is 740 parallel pairs "
                         "and 1898 perpendicular. This picks between them.")
    ap.add_argument("--mix", default=None,
                    help="SC,BCC,FCC fractions for a mixed lattice, "
                         "e.g. 0.6,0.2,0.2. Implies --lattice MIX.")
    ap.add_argument("--clearance", type=float, default=0.9)
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES[args.lattice])
    if args.mix:
        # SC:BCC:FCC as fractions, e.g. "0.6,0.2,0.2". The lattice a strand
        # belongs to is what sets its neighbour spacing, so the mix is
        # effectively a knob on how tightly the shells are packed.
        f = [float(x) for x in args.mix.split(",")]
        spec["mix"] = {"SC": f[0], "BCC": f[1], "FCC": f[2]}
        spec["lattice"] = "MIX"
    spec["dims"] = (args.dims,) * 3
    graph = build_network(spec)
    geo = geometry(graph, dp=args.dp, bond=BOND, coil=args.coil)
    dims = np.asarray(graph.graph["box"], float)
    keys = sorted(geo["chords"])
    edges = sorted(graph.edges())
    idx_of = {frozenset(e): i for i, e in enumerate(edges)}

    G = nx.MultiGraph()
    G.add_nodes_from(graph.nodes(data=True))
    for u, v in graph.edges():
        G.add_edge(u, v)
    dist = chain_distances(G, dims)
    shells = neighbour_shells(G, dims, max_shell=max(args.shells),
                              distances=dist)

    gap_of = {}
    for (ca, cb), r in dist.items():
        ia = idx_of.get(frozenset((ca[0], ca[1])))
        ib = idx_of.get(frozenset((cb[0], cb[1])))
        if ia is not None and ib is not None:
            gap_of[(min(ia, ib), max(ia, ib))] = r

    out_root = OUT / f"shellshow_{args.lattice}_{args.orient}"
    shutil.rmtree(out_root, ignore_errors=True)
    out_root.mkdir(parents=True, exist_ok=True)

    print(f"  {args.lattice} {args.dims}^3, box {geo['L'][0]:.1f} sigma, "
          f"{len(keys)} chains, {args.orient} partners, "
          f"all straight except the designed one")
    print(f"\n  {'shell':>6} {'partner':>8} {'gap':>7} {'built':>7}  file")

    # Which chain to route, per shell.
    #
    # Fixing on chain 0 finds nothing when that chain happens to have no
    # partner of the requested orientation: on MIX, asking for parallel
    # partners produced no file at all for any shell. Any chain with a
    # qualifying partner will do, so look for one.
    def find_pair(sh):
        for chain, by in sorted(shells.items()):
            ia = idx_of.get(frozenset((chain[0], chain[1])))
            if ia is None:
                continue
            mine_dir = chain_dir(graph, edges[ia])
            for o in by.get(sh, ()):
                ib = idx_of.get(frozenset((o[0], o[1])))
                if ib is None or ib == ia:
                    continue
                if args.orient != "any":
                    par = abs(float(mine_dir
                                    @ chain_dir(graph, edges[ib]))) > 0.9
                    if par != (args.orient == "parallel"):
                        continue
                return ia, ib
        return None, None

    for sh in args.shells:
        routed, partner = find_pair(sh)
        if partner is None:
            print(f"  {sh:>6} {'-':>8} {'-':>7} {'-':>7}  "
                  f"no partner in this shell")
            continue

        # Everything straight, which is what makes the design visible.
        paths = {k: straight(c0, c1, args.dp + 1)
                 for k, (c0, c1) in geo["chords"].items()}

        root = out_root / f"shell{sh}"
        _n, node_atom, chain_atoms = write_system(graph, geo, paths, root)
        seq = {k: chain_ids(k, node_atom, chain_atoms, geo["ends"])
               for k in keys}

        xyz = {}
        for k in keys:
            for aid, q in zip(seq[k], paths[k]):
                xyz[aid] = q
        ids = sorted(xyz)
        mine = set(seq[routed])
        avoid = Clearance(np.array([xyz[i] for i in ids if i not in mine]),
                          geo["L"], args.clearance)
        try:
            p = construct(paths, routed, partner,
                          max(0.5, 0.5 * args.count), geo["L"], avoid,
                          radius=args.ring, dp=args.dp, span=(0.3, 0.7))
        except ValueError as e:
            print(f"  {sh:>6} {partner:>8} "
                  f"{gap_of.get((min(routed, partner), max(routed, partner)), float('nan')):>7.2f} "
                  f"{'-':>7}  {str(e)[:40]}")
            continue

        for aid, q in zip(seq[routed], p):
            xyz[aid] = q

        # Rewrite with the routed chain in place, then label by type.
        src = root / "01_Topology" / "system.data"
        if not src.exists():
            src = next(root.rglob("system.data"))
        from tests.workflows.entangle_relaxed import rewrite_coords
        placed = root / "placed.data"
        rewrite_coords(src, placed, xyz)

        junctions = set(node_atom.values())
        types = {}
        for i in ids:
            if i in junctions:
                types[i] = 4
            elif i in mine:
                types[i] = 2
            elif i in set(seq[partner]):
                types[i] = 3
            else:
                types[i] = 1
        final = (out_root / f"{args.lattice}_{args.orient}_shell{sh}"
                 f"_chain{routed}_with{partner}.data")
        retype(placed, final, types)

        gap = gap_of.get((min(routed, partner), max(routed, partner)),
                         float("nan"))
        print(f"  {sh:>6} {partner:>8} {gap:>7.2f} {args.count:>7}  "
              f"{final.name}")

    print(f"\n  written to {out_root}")
    print("  open in OVITO and colour by Particle Type:")
    print("    1  the rest of the network, straight")
    print("    2  the routed chain")
    print("    3  the partner it winds around")
    print("    4  crosslink junctions")
    print("\n  the gap column is the closest approach between the two chains "
          "in lattice units,")
    print("  which is what defines the shell: 0 means they share a junction, "
          "then 1, sqrt(2), sqrt(3)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
