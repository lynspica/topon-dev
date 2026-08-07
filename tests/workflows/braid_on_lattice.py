"""Build a real lattice network with waypoint braids, for viewing in OVITO.

Run directly::

    python tests/workflows/braid_on_lattice.py --lattice SC --dims 4x4x4
    python tests/workflows/braid_on_lattice.py --lattice Diamond --dims 3x3x3
    python tests/workflows/braid_on_lattice.py --lattice MIX --mix 0.2,0.4,0.4

This is the honest test of the braid: a real lattice from the real
generator, sculpted to a real functionality, with real chains between real
junctions under periodic boundaries. Synthetic parallel chords in empty
space cannot show whether a braid survives contact with a lattice -- whether
the chains it braids are actually near each other, whether the braid volume
is clear, or whether the winding it produces is the one that was asked for.

Writes a LAMMPS data file per case. Open it in OVITO: entanglement is a
three-dimensional question and a flat projection cannot answer it, since
which strand passes in front is exactly what distinguishes a hook from a
near miss.

Bead type 1 is a junction, 2 is an ordinary chain bead, 3 is a bead inside a
braid. Colouring by type in OVITO shows immediately where the entanglements
are; the per-atom `molecule` id is the chain, so colouring by that instead
separates the partners.
"""
from __future__ import annotations

import argparse
import random
import sys
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.entanglement import (  # noqa: E402
    BraidShape,
    ContactRequest,
    allocate_contacts,
    compose_chain_path,
    far_closed_linking,
    min_separation,
)
from topon.topology.generator_python import PythonTopologyGenerator  # noqa: E402

OUT = ROOT / "tests/output/braid_lattice"


class _Cfg:
    def __init__(self, lattice, dims, max_func, degree_dist, mix, cutoff):
        self.lattice_type = lattice
        self.lattice_size = dims
        self.max_functionality = max_func
        self.degree_distribution = degree_dist
        self.periodicity = "111"
        self.mix_fractions = mix
        self.mix_cutoff = cutoff


def build_lattice(lattice, dims, max_func, degree_dist, mix, cutoff, seed):
    random.seed(seed)
    gen = PythonTopologyGenerator(
        _Cfg(lattice, dims, max_func, degree_dist, mix, cutoff))
    with redirect_stdout(StringIO()):
        graphs = gen.generate(trials=4000, max_saves=1, time_limit=120)
    if not graphs:
        raise SystemExit(f"sculpting produced no {lattice} graph for {dims}")
    return graphs[0]


def density_scale(graph, dp, density):
    """Lattice units -> sigma, from the bead count and the target density.

    This is what the pipeline does, and getting it from physics rather than
    picking a number matters here: it sets the junction spacing, which sets
    how far apart neighbouring strands are, which decides whether any pair is
    close enough to braid at all. An arbitrary scale of 12 sigma per lattice
    unit put every neighbour 12 sigma apart and produced zero candidates on a
    lattice whose strands are in fact adjacent.
    """
    n_beads = graph.number_of_edges() * dp + graph.number_of_nodes()
    volume = n_beads / density
    cells = float(np.prod(np.asarray(graph.graph["box"], float)))
    return (volume / cells) ** (1.0 / 3.0)


def chain_chords(graph, box, scale):
    """One chord per edge, in real units, unwrapped across the boundary.

    A chain between two junctions that sit on opposite sides of the box is
    physically short -- it crosses the boundary. Its chord has to be built
    from the minimum-image vector, or every wrapping chain is drawn as a
    line straight across the system.
    """
    pos = {n: np.asarray(d["pos"], float) * scale
           for n, d in graph.nodes(data=True)}
    L = np.asarray(box, float) * scale

    chords, wraps = {}, set()
    for k, (u, v) in enumerate(sorted(graph.edges())):
        a = pos[u]
        raw = pos[v] - a
        mic = raw - L * np.round(raw / L)
        if not np.allclose(mic, raw):
            wraps.add(k)
        chords[k] = (a, a + mic)
    return chords, pos, wraps


def candidate_requests(chords, max_gap, max_requests):
    """Pairs of chains close enough to be worth braiding.

    Ranked by how close they come, so genuinely adjacent strands win the
    room over distant ones. This stands in for the shell-weighted selection
    that belongs in the assignment stage; here it only has to produce a
    realistic candidate set.
    """
    from topon.conformation.entanglement import closest_approach, gap_at

    ids = sorted(chords)
    out = []
    for i, ka in enumerate(ids):
        a0, a1 = chords[ka]
        for kb in ids[i + 1:]:
            b0, b1 = chords[kb]
            # Cheap reject on chord midpoints before the exact test.
            if np.linalg.norm((a0 + a1) / 2 - (b0 + b1) / 2) > 3.0 * max_gap:
                continue
            s, _ = closest_approach(a0, a1, b0, b1)
            gap, _ = gap_at(a0, a1, b0, b1, s)
            if 1e-6 < gap <= max_gap:
                out.append((gap, ka, kb))

    out.sort()
    return [ContactRequest(ka, kb, windings=1, priority=-gap)
            for gap, ka, kb in out[:max_requests]]


def beads_per_chain(chords, dp):
    return {k: dp + 2 for k in chords}          # +2 for the two junctions


def write_lammps(path, chords, paths, alloc, box, scale):
    """Write a bead-spring data file: one molecule per chain."""
    braided = {}
    for a in alloc.accepted:
        for chain in (a.request.chain_a, a.request.chain_b):
            braided.setdefault(chain, []).append(a)

    atoms, bonds = [], []
    for chain in sorted(paths):
        p = paths[chain]
        first = len(atoms) + 1
        for i, xyz in enumerate(p):
            if i == 0 or i == len(p) - 1:
                t = 1                                   # junction
            else:
                t = 2
                for a in braided.get(chain, []):
                    u = float((xyz - a.contact.origin) @ a.contact.axis)
                    if abs(u) < a.half_span:
                        t = 3                           # inside a braid
                        break
            atoms.append((len(atoms) + 1, chain + 1, t, xyz))
        for i in range(len(p) - 1):
            bonds.append((len(bonds) + 1, 1, first + i, first + i + 1))

    L = np.asarray(box, float) * scale
    lo = -0.25 * L
    hi = 1.25 * L
    with open(path, "w") as f:
        f.write("LAMMPS data file: lattice network with waypoint braids\n\n")
        f.write(f"{len(atoms)} atoms\n{len(bonds)} bonds\n\n")
        f.write("3 atom types\n1 bond types\n\n")
        for k, ax in enumerate("xyz"):
            f.write(f"{lo[k]:.6f} {hi[k]:.6f} {ax}lo {ax}hi\n")
        f.write("\nMasses\n\n1 1.0\n2 1.0\n3 1.0\n\nAtoms # molecular\n\n")
        for aid, mol, t, xyz in atoms:
            f.write(f"{aid} {mol} {t} {xyz[0]:.5f} {xyz[1]:.5f} {xyz[2]:.5f}\n")
        f.write("\nBonds\n\n")
        for bid, bt, i, j in bonds:
            f.write(f"{bid} {bt} {i} {j}\n")
    return len(atoms), len(bonds)


def run_case(lattice, dims, max_func, degree_dist, mix, cutoff, dp, density,
             gap_factor, max_requests, seed, shape):
    graph = build_lattice(lattice, dims, max_func, degree_dist, mix, cutoff, seed)
    box = graph.graph["box"]
    scale = density_scale(graph, dp, density)
    chords, _, wraps = chain_chords(graph, box, scale)

    # Candidates are pairs closer than a fraction of the junction spacing,
    # so the threshold follows the lattice instead of being a fixed number.
    max_gap = gap_factor * scale
    reqs = candidate_requests(chords, max_gap, max_requests)
    alloc = allocate_contacts(reqs, chords, shape)

    nb = beads_per_chain(chords, dp)
    paths = {k: compose_chain_path(k, alloc, chords, nb[k]) for k in chords}

    # Verify against the built paths, not the plan.
    correct = mismatched = 0
    clearances = []
    for a in alloc.accepted:
        r = a.request
        lk = far_closed_linking(paths[r.chain_a], paths[r.chain_b], a.contact)
        clearances.append(min_separation(paths[r.chain_a], paths[r.chain_b]))
        # Compare magnitudes. The braid's handedness follows the frame,
        # which follows the chord orientations, so on a real lattice both
        # signs appear. A left-handed single winding is as much an
        # entanglement as a right-handed one; only the count is prescribed.
        if round(abs(lk)) == a.windings:
            correct += 1
        else:
            mismatched += 1

    reasons = {}
    for rj in alloc.rejected:
        reasons[rj.reason] = reasons.get(rj.reason, 0) + 1

    OUT.mkdir(parents=True, exist_ok=True)
    tag = f"{lattice}_{dims[0]}x{dims[1]}x{dims[2]}"
    n_atoms, n_bonds = write_lammps(OUT / f"{tag}.data", chords, paths,
                                    alloc, box, scale)

    print(f"{tag:18s} {graph.number_of_nodes():5d} junctions "
          f"{len(chords):5d} chains ({len(wraps)} wrap the boundary), "
          f"spacing {scale:.1f} sigma")
    print(f"{'':18s} {len(reqs):5d} candidates -> {len(alloc.accepted):4d} braids, "
          f"{correct} realised exactly, {mismatched} wrong")
    if clearances:
        print(f"{'':18s} clearance {min(clearances):.2f} to {max(clearances):.2f}, "
              f"multi-partner chains: "
              f"{sum(1 for v in alloc.partners.values() if len(v) > 1)}")
    for why, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"{'':18s}   {n:4d} refused: {why}")
    print(f"{'':18s} wrote {tag}.data  ({n_atoms} beads, {n_bonds} bonds)")
    print()
    return mismatched == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lattice", default="all",
                    help="SC, BCC, FCC, Diamond, MIX, or all")
    ap.add_argument("--dims", default="4x4x4")
    ap.add_argument("--max-func", type=int, default=4)
    ap.add_argument("--degree-dist", default="0:0,1:0")
    ap.add_argument("--mix", default="0.2,0.4,0.4")
    ap.add_argument("--mix-cutoff", type=float, default=1.0)
    ap.add_argument("--dp", type=int, default=60,
                    help="beads per chain between the junctions")
    ap.add_argument("--density", type=float, default=0.85,
                    help="bead density; sets the junction spacing")
    ap.add_argument("--gap-factor", type=float, default=0.6,
                    help="candidate cutoff, in units of junction spacing")
    ap.add_argument("--max-requests", type=int, default=60)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    dims = tuple(int(v) for v in args.dims.lower().split("x"))
    mix = dict(zip(("SC", "BCC", "FCC"),
                   (float(v) for v in args.mix.split(","))))
    shape = BraidShape()

    cases = (["SC", "BCC", "FCC", "Diamond", "MIX"]
             if args.lattice == "all" else [args.lattice])

    print("=" * 78)
    print("Waypoint braids on real lattices (open the .data files in OVITO)")
    print("=" * 78)
    ok = True
    for lat in cases:
        d = (3, 3, 3) if lat == "Diamond" and dims == (4, 4, 4) else dims
        mf = 4 if lat != "MIX" else 12
        try:
            ok &= run_case(lat, d, mf, args.degree_dist, mix, args.mix_cutoff,
                           args.dp, args.density, args.gap_factor,
                           args.max_requests, args.seed, shape)
        except SystemExit as exc:
            print(f"{lat}: {exc}\n")
            ok = False

    print("=" * 78)
    print(f"VERDICT: {'PASS' if ok else 'FAIL'}    files in {OUT}")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
