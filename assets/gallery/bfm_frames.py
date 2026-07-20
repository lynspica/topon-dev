"""How the BFM protein-network path reaches its gel point (topro).

Places self-avoiding-walk chains on a 6-neighbour cubic lattice (the REAL
`topon.protein_network.bfm` routines), then applies inter-chain crosslinks one by
one, capturing the state at intervals. Each captured frame records the chain
positions (fixed), the crosslinks so far, and each chain's connected component, so
the render can colour the growing giant cluster and mark the gel point (first
percolation to a single component).

TOPON-python stage: writes f{k}.data (chain nodes as atoms, backbone + crosslink
bonds, mol = component rank so the render needs no bfm/networkx) + meta.json with
the full conversion / largest-fraction / n-components series.
"""
import json
import sys
from pathlib import Path

import numpy as np


def run(n_chains=30, n_repeats=12, segs=2, pack=0.5, cut_ang=22.0, seed=7,
        equil=120000, n_equil_frames=16):
    """Place chains, EQUILIBRATE IN SLICES (capturing the Monte-Carlo chain
    motion -- end moves, kink/crankshaft, reptation), then crosslink to the gel.

    The chain motion is the defining BFM step, so it is animated rather than run
    once up front: `equilibrate` mutates `chains` in place, so calling it in
    slices and deep-copying between them yields real MC trajectory frames."""
    import copy
    from topon.protein_network import bfm
    rng = np.random.default_rng(seed)
    n_nodes = n_repeats * segs * 2 + 1
    Nx = bfm.compute_lattice_size(n_chains, n_nodes, pack)
    chains = bfm.place_chains(n_chains, n_nodes, Nx, Nx, Nx, rng)
    equil_frames = [copy.deepcopy(chains)]                 # as-placed
    per = max(1, equil // n_equil_frames)
    for _ in range(n_equil_frames):
        bfm.equilibrate(chains, Nx, Nx, Nx, per, rng)
        equil_frames.append(copy.deepcopy(chains))
    yps = bfm.get_y_positions(n_repeats, segs)
    cands = sorted(bfm.find_crosslink_candidates_distance(
        chains, yps, Nx, Nx, Nx, lattice_scale_ang=6.0, max_distance_ang=cut_ang,
        min_intrachain_sep=999), key=lambda c: c[2])   # inter-chain only, near-first
    # apply, recording the full series and the crosslink after which each forms
    uf = bfm.UnionFind(n_chains)
    reacted, reactions, series = set(), [], []
    total_y = len(yps) * n_chains
    gel_conv = None
    for (a, b, dist) in cands:
        (c1, i1), (c2, i2) = a, b
        if (c1, i1) in reacted or (c2, i2) in reacted:
            continue
        reacted.add((c1, i1)); reacted.add((c2, i2))
        reactions.append((c1, chains[c1][i1], c2, chains[c2][i2]))  # chain,flat,chain,flat
        if c1 != c2:
            uf.union(c1, c2)
        comps = uf.components()
        big = max(len(v) for v in comps.values())
        conv = len(reacted) / total_y
        series.append((conv, len(comps), big / n_chains))
        if gel_conv is None and len(comps) == 1:
            gel_conv = conv
    return bfm, Nx, chains, reactions, series, gel_conv, n_chains, equil_frames


def component_rank(bfm, n_chains, reactions_upto):
    """Return {chain_index: rank} where rank 0 = largest component, 1.. others by
    descending size, so the render can colour the giant cluster boldly."""
    uf = bfm.UnionFind(n_chains)
    for (c1, _, c2, _) in reactions_upto:
        if c1 != c2:
            uf.union(c1, c2)
    comps = uf.components()                         # root -> [chain ids]
    order = sorted(comps.values(), key=len, reverse=True)
    # Rank 0 is drawn as the bold "giant cluster". Before any chains have joined
    # up, every component is a single chain and there is no giant -- so shift the
    # ranks and let every chain take an ordinary cluster colour instead.
    offset = 1 if len(order[0]) == 1 else 0
    rank = {}
    for r, members in enumerate(order):
        for ci in members:
            rank[ci] = r + offset
    return rank


def write_data(path, Nx, chains, reactions_upto, rank):
    """chain nodes as atoms (mol = component rank + 1), backbone + crosslink bonds."""
    from topon.protein_network.bfm import flat_to_xyz
    # global atom id per (chain, node)
    aid = {}
    lines_atoms, lines_bonds = [], []
    nid = 0
    for ci, chain in enumerate(chains):
        molr = rank.get(ci, 0) + 1
        for ni, flat in enumerate(chain):
            nid += 1
            aid[(ci, ni)] = nid
            x, y, z = flat_to_xyz(flat, Nx, Nx)
            lines_atoms.append(f"{nid} {molr} 1 0 {x}.0 {y}.0 {z}.0")
    N = nid
    bid = 0
    # backbone bonds (type 1)
    for ci, chain in enumerate(chains):
        for ni in range(len(chain) - 1):
            bid += 1
            lines_bonds.append(f"{bid} 1 {aid[(ci, ni)]} {aid[(ci, ni+1)]}")
    # crosslink bonds (type 2) -- need atom ids from flat; find node index by flat
    flat_index = {}
    for ci, chain in enumerate(chains):
        for ni, flat in enumerate(chain):
            flat_index[(ci, flat)] = ni
    for (c1, f1, c2, f2) in reactions_upto:
        n1 = flat_index[(c1, f1)]; n2 = flat_index[(c2, f2)]
        bid += 1
        lines_bonds.append(f"{bid} 2 {aid[(c1, n1)]} {aid[(c2, n2)]}")
    hdr = ["bfm frame\n", f"{N} atoms", f"{bid} bonds",
           "1 atom types", "2 bond types",
           f"0 {Nx} xlo xhi", f"0 {Nx} ylo yhi", f"0 {Nx} zlo zhi",
           "", "Masses", "", "1 1.0", "", "Atoms # full", ""]
    Path(path).write_text("\n".join(hdr + lines_atoms + ["", "Bonds", ""] +
                                    lines_bonds) + "\n")


def build_all(outdir, n_frames=30, **kw):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    (bfm, Nx, chains, reactions, series, gel_conv, n_chains,
     equil_frames) = run(**kw)
    nrx = len(reactions)
    k = 0
    # --- phase A: Monte-Carlo chain motion (no crosslinks yet) ---
    rank0 = component_rank(bfm, n_chains, [])
    for ch in equil_frames:
        write_data(outdir / f"f{k:03d}.data", Nx, ch, [], rank0)
        k += 1
    n_equil = k
    # how far the chains actually moved during equilibration (sanity + caption)
    a = np.array([f for ch in equil_frames[0] for f in ch], dtype=float)
    b = np.array([f for ch in equil_frames[-1] for f in ch], dtype=float)
    moved = float((a != b).mean())
    # --- phase B: crosslinking to the gel (chains now fixed) ---
    idxs = sorted(set(int(round(i)) for i in np.linspace(0, nrx, n_frames)))
    for upto in idxs:
        rr = reactions[:upto]
        rank = component_rank(bfm, n_chains, rr)
        write_data(outdir / f"f{k:03d}.data", Nx, chains, rr, rank)
        k += 1
    conv_series = [s[0] for s in series]
    frac_series = [s[2] for s in series]
    ncomp_series = [s[1] for s in series]
    # conversion at each captured frame (0 crosslinks -> conv 0)
    frame_conv = [0.0] + [series[i - 1][0] for i in idxs[1:]]
    frame_frac = [1.0 / n_chains] + [series[i - 1][2] for i in idxs[1:]]
    frame_ncomp = [n_chains] + [series[i - 1][1] for i in idxs[1:]]
    (outdir / "meta.json").write_text(json.dumps({
        "Nx": Nx, "n_chains": n_chains, "n_crosslinks": nrx,
        "n_equil_frames": n_equil, "n_xlink_frames": len(idxs),
        "n_frames": k, "gel_conv": gel_conv, "equil_moved_frac": moved,
        "conv_series": conv_series, "frac_series": frac_series,
        "ncomp_series": ncomp_series,
        "frame_conv": frame_conv, "frame_frac": frame_frac,
        "frame_ncomp": frame_ncomp,
    }, indent=1))
    print(f"[bfm_frames] {n_chains} chains, {Nx}^3 lattice; "
          f"{n_equil} equilibration frames ({moved*100:.0f}% of nodes moved), "
          f"{len(idxs)} crosslink frames, {nrx} crosslinks, gel@{gel_conv} "
          f"-> {outdir}")


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "C:/tmp/bfm_frames"
    build_all(out)
