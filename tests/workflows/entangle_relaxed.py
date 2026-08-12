"""Design entanglements into an already-relaxed melt, so they survive.

    python tests/workflows/entangle_relaxed.py
    python tests/workflows/entangle_relaxed.py --rounds 8 --per-round 24

Building a designed entanglement into a fresh random-walk conformation and
then minimising destroys about half of it. The reason is not subtle: the
beads move 3.34 sigma on average during minimisation, up to 12.45, while the
designed loop is only 1.2 to 3.4 sigma across. The feature is the same size
as the motion, so no amount of care in building it helps.

That motion is not the entanglement's fault. A random walk drawn without
excluded volume overlaps itself and its neighbours everywhere, and the first
minimisation stage is mostly relieving that. So relax first, and design into
the relaxed melt:

    1. build a plain melt, no designed entanglements
    2. minimise it -- this absorbs all the overlap relief
    3. route chains through the relaxed coordinates
    4. minimise again, which now has little left to do
    5. measure each requested pair on its own

Step 4 is the one that matters. A system already at equilibrium has no reason
to move far, so the designed feature is no longer competing with a global
rearrangement.
"""
from __future__ import annotations

import argparse
import collections
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.paths import (  # noqa: E402
    Clearance,
    bridging_walk,
    loop_around,
    route_through,
)
from tests.workflows.entangle_all import CASES  # noqa: E402
from tests.workflows.entangle_design import route_one  # noqa: E402
from tests.workflows.entangle_search import (  # noqa: E402
    Wish,
    _both,
    measure_batch,
)
from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    COIL,
    DP,
    LATTICE,
    OUT,
    build_network,
    chain_ids,
    conform_and_script,
    geometry,
    scale_for_design,
    read_data,
    report_bonds,
    run_md,
    unwrap_chain,
    write_system,
    z1_export,
)


def rewrite_coords(src, dst, new_xyz):
    """Copy a LAMMPS data file, replacing the coordinates of some atoms.

    Everything else -- header, masses, bonds, angles -- is carried through
    untouched, so the result is the same system in a different conformation.
    Writing a fresh file through the chemistry and conformation stages would
    also work but would rebuild the topology, and the point here is to change
    only where the atoms are.
    """
    out, section = [], None
    for line in Path(src).read_text().splitlines(True):
        s = line.strip()
        if s and s.split()[0] in ("Atoms", "Bonds", "Masses", "Velocities",
                                  "Angles"):
            section = s.split()[0]
            out.append(line)
            continue
        if section == "Atoms" and s and not s.startswith("#"):
            p = s.split()
            if len(p) >= 7 and int(p[0]) in new_xyz:
                x, y, z = new_xyz[int(p[0])]
                p[4], p[5], p[6] = f"{x:.6f}", f"{y:.6f}", f"{z:.6f}"
                out.append(" ".join(p) + "\n")
                continue
        out.append(line)
    Path(dst).write_text("".join(out))
    return dst


def paths_from(data_file, keys, seq):
    box, xyz, _ = read_data(data_file)
    return box, {k: unwrap_chain(seq[k], xyz, box) for k in keys}, xyz


def measure_pairs(data_file, seq, pairs, work):
    """Each pair on its own -- the only export Z1+ takes reliably."""
    work.mkdir(parents=True, exist_ok=True)
    for old in list(work.glob("*.Z1")) + list(work.glob("SP_*.dat")):
        old.unlink()
    for a, b in pairs:
        z1_export(data_file, [seq[a], seq[b]], work / f"p{a}_{b}.Z1")
    res = measure_batch(work) or {}
    return {(a, b): (_both(res[f"p{a}_{b}"], 1, 2)
                     if f"p{a}_{b}" in res else None) for a, b in pairs}


def route_pairwise(paths, keys, seq, geo, routed, target, want, template,
                   rounds, per_round, rng, work, relax=None, avoid=None,
                   site_span=(0.05, 0.95)):
    """Search for a path, scoring each candidate on the pair alone.

    The whole-system score used elsewhere is unavailable here: Z1+ crashes on
    a relaxed melt read back from LAMMPS output, so every candidate comes back
    empty and the search has no signal. Measuring just the two chains works
    reliably and scores the thing being asked for.

    The cost is collateral blindness -- what else the routed chain picks up is
    not seen, so it cannot be traded against. The primary objective is served;
    the secondary one is not.

    ``avoid`` holds every bead except the routed chain's own, so candidates
    are drawn around what is already there rather than through it.
    """
    from tests.workflows.entangle_search import propose

    L = geo["L"]
    a0, a1 = paths[routed][0], paths[routed][-1]
    best = None
    around, spread = None, 1.0

    for _ in range(rounds):
        tgt = paths[target]
        tgt = tgt + L * np.round((paths[routed].mean(0) - tgt.mean(0)) / L)
        cands = propose(a0, a1, [tgt], L, per_round, rng, around, spread,
                        avoid=avoid, site_span=site_span)
        if not cands:
            continue

        work.mkdir(parents=True, exist_ok=True)
        for old in list(work.glob("*.Z1")) + list(work.glob("SP_*.dat")):
            old.unlink()
        tmp = work / "cand.data"
        for i, (_knobs, path) in enumerate(cands):
            new = {aid: xyzp for aid, xyzp in zip(seq[routed], path)}
            rewrite_coords(template, tmp, new)
            # Score after a short minimisation, not as built.
            #
            # A loop placed into a relaxed melt is a local strain, and
            # relaxing it is what opens or tightens the winding. Measured on
            # candidates chosen as-built: a pair designed at 1 came back at 0
            # and one designed at 2 came back at 5, with beads moving only
            # 0.57 sigma. Nothing about the built geometry predicts that, so
            # the only way to select for it is to relax each candidate and
            # look. In an already-relaxed system stage 1 costs a fraction of
            # a second, which is what makes this affordable at all.
            scored_file = tmp
            if relax is not None:
                r = relax(tmp, i)
                if r is not None:
                    scored_file = r
            z1_export(scored_file, [seq[routed], seq[target]],
                      work / f"c{i:03d}.Z1")
        res = measure_batch(work) or {}

        scored = []
        for i, (knobs, path) in enumerate(cands):
            r = res.get(f"c{i:03d}")
            if r is None:
                continue
            got = _both(r, 1, 2)
            scored.append((abs(got - want), got, knobs, path))
        if not scored:
            continue
        scored.sort(key=lambda t: t[0])
        if best is None or scored[0][0] < best[0]:
            best = scored[0]
            around, spread = best[2], max(0.35, spread * 0.6)
        else:
            spread = min(1.5, spread * 1.4)
        if best[0] == 0:
            break
    return best


def construct(paths, routed, target, want, box, avoid, at=None, radius=2.0,
              phase=0.0, dp=DP, bond=BOND, span=(0.3, 0.7), rank=0):
    """Build the requested count instead of searching for it.

    One full turn around a strand is two entanglement points, not one: pulled
    taut, a loop presses on the strand from both sides. Measured over four
    phases, one turn gave 2 as built and 2 after minimisation every time, so
    the count follows from the winding and does not have to be hunted for.

    Two things had to be true before that worked. The leftover contour is spent
    as a zigzag rather than a random walk -- a chain carrying 77 sigma for a
    21 sigma route will otherwise wander back over the target and add crossings
    of its own, which is what made the same design come back 4, 7 and 0 on
    three seeds. And the path is drawn around the beads already there, so the
    minimisation that follows has no overlap to resolve by pushing chains
    through one another.

    ``at`` names the site along the target, as a fraction of its length. Left
    unset it is chosen as the place inside ``span`` where the two chains
    actually come closest, which matters more than it sounds: a ring put at
    mid-strand regardless built with a clean 0.92 sigma of clearance and still
    measured zero, because the routed chain only reaches it by a long
    excursion and the primitive path retracts straight back out. A winding has
    to sit where the chain already passes. ``rank`` takes the second-closest
    such place, and so on, so two chains sent to one target do not stack.

    Odd counts are not available this way. Raises rather than quietly building
    the nearest even number, since a count you did not ask for is worse than a
    refusal.
    """
    if want <= 0 or want % 2:
        raise ValueError(
            f"{want} is not buildable: one turn is two entanglement points, "
            f"so this construction makes even counts only")
    turns = want // 2
    tgt = paths[target]
    tgt = tgt + box * np.round((paths[routed].mean(0) - tgt.mean(0)) / box)
    if at is not None:
        i = int(np.clip(round(at * len(tgt)), 1, len(tgt) - 2))
    else:
        lo = max(1, int(span[0] * len(tgt)))
        hi = min(len(tgt) - 2, int(span[1] * len(tgt)))
        d = tgt[lo:hi + 1, None, :] - paths[routed][None, :, :]
        d -= box * np.round(d / box)
        near = np.linalg.norm(d, axis=2).min(axis=1)
        order = np.argsort(near)
        i = int(lo + order[min(rank, len(order) - 1)])
    n_pts = int(np.clip(round(2.0 * np.pi * radius * turns / (1.6 * bond)),
                        4, 24))
    ring = loop_around(tgt, i, radius, n_pts, phase, avoid, float(turns))
    return route_through(paths[routed][0], paths[routed][-1], list(ring),
                         dp + 1, bond, tgt[i], avoid)


def construct_exact(paths, routed, target, want, box, avoid, seq, template,
                    work, span, radius, dp, ranks=8, phases=6, relax=None):
    """Build the exact count, by enumerating placements rather than drawing.

    One turn is two entanglement points where the routed chain passes its
    target once, and four where it passes twice, so the winding alone does not
    fix the count -- what else the route does with that target counts as well.
    Measured on two pairs built identically: one came back 2 and the other 4.

    So build, measure, and move the winding if the number is wrong. The
    placements are tried in a fixed order and the first exact match wins, so
    this is an enumeration and not a search: the same request gives the same
    system every time.

    With ``relax`` the count is measured *after* a minimisation rather than as
    built, and that is the difference between a design that usually holds and
    one that holds by construction. Enumerated over 24 placements of one pair:
    17 built the number asked for, and 15 of those kept it, so building it
    correctly still left one in eight to be lost later. Neither clearance
    (correlation -0.27) nor how much of the chain lay against the target
    (-0.18) predicted which. Since what survives cannot be predicted from the
    built geometry, the only sound thing to do is relax each placement and
    look, and 17 of the 24 give the right count afterwards, so a surviving
    placement is not hard to find.

    Returns ``(path, count)``, with the closest miss if no placement is exact,
    so the caller can say what it got rather than imply it got what was asked.
    """
    best = None
    for rank in range(ranks):
        for k in range(phases):
            phase = 2.0 * np.pi * k / phases
            try:
                p = construct(paths, routed, target, want, box, avoid,
                              radius=radius, phase=phase, dp=dp, span=span,
                              rank=rank)
            except ValueError:
                continue
            trial = dict(zip(seq[routed], p))
            rewrite_coords(template, work / "try.data", trial)
            scored = work / "try.data"
            if relax is not None:
                r = relax(scored)
                if r is None:
                    continue
                scored = r
            got = measure_pairs(scored, seq,
                                [(min(routed, target), max(routed, target))],
                                work)[(min(routed, target),
                                       max(routed, target))]
            if got is None:
                continue
            if got == want:
                return p, got
            if best is None or abs(got - want) < abs(best[1] - want):
                best = (p, got)
    return best if best is not None else (None, None)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rounds", type=int, default=6)
    ap.add_argument("--per-round", type=int, default=20)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--verify", action="store_true",
                    help="pick the placement whose count survives a "
                         "minimisation, rather than the one that builds it. "
                         "Building it right still loses one in eight, and "
                         "nothing about the built geometry says which.")
    ap.add_argument("--construct", action="store_true",
                    help="build the requested count from the winding instead "
                         "of searching for it. Deterministic: the same "
                         "request gives the same system.")
    ap.add_argument("--want", type=int, default=2,
                    help="entanglements per designed pair")
    ap.add_argument("--site-lo", type=float, default=0.30)
    ap.add_argument("--site-hi", type=float, default=0.70,
                    help="the stretch of the target a loop may sit on. Near a "
                         "tip it slides off during relaxation without any "
                         "chain crossing anything, and the count drops.")
    ap.add_argument("--ring", type=float, default=2.0,
                    help="ring radius the box is sized to hold")
    ap.add_argument("--margin", type=float, default=2.0,
                    help="how much chain contour to leave spare")
    ap.add_argument("--coil", type=float, default=None,
                    help="contour over chord. Sets the lattice scale, and so "
                         "the box: density is reported, not chosen.")
    ap.add_argument("--density", type=float, default=None,
                    help="use the melt route instead and pack to this "
                         "density. Leaves no free volume to route through: "
                         "at 0.85 a point is within 0.9 sigma of 2.6 beads on "
                         "average, so a routed path lands on beads whatever "
                         "it does, and the overlap relief that follows is "
                         "what destroys the design.")
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dims", type=int, default=6,
                    help="lattice cells per side. 4 puts the box at 21.6 "
                         "sigma, which is shorter than a routed chain's own "
                         "unwrapped extent, and Z1+ rejects those outright -- "
                         "so the search scores on refusals and is blind. The "
                         "spacing is set by density, so more cells is a bigger "
                         "box at the same physics.")
    ap.add_argument("--clearance", type=float, default=0.9,
                    help="how close a routed bead may come to an existing one")
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (args.dims,) * 3
    graph = build_network(spec)
    # The design comes first and the box is sized to it. Which chains are
    # partnered has to be decided in graph units, before there is a scale at
    # all, because the scale is what the decision determines.
    box_g = np.asarray(graph.graph["box"], float)
    raw = {n: np.asarray(d["pos"], float) for n, d in graph.nodes(data=True)}
    edges = sorted(graph.edges())

    def chord_pts(k, span=(0.0, 1.0)):
        u, v = edges[k]
        a = raw[u]
        mic = (raw[v] - a) - box_g * np.round((raw[v] - a) / box_g)
        return a + np.linspace(span[0], span[1], 24)[:, None] * mic

    def nearest_chord(a, taken=()):
        # Judged against the stretch of the partner a loop may actually sit
        # on, not its nearest approach, for the same reason the box is sized
        # that way.
        #
        # `taken` keeps two routed chains off the same partner. Sending both
        # to the nearest one puts the second ring against the first: measured,
        # the second pair built at a tightest contact of 0.34 sigma and read
        # zero, while the first built and survived cleanly.
        span = (args.site_lo, args.site_hi)
        best, out = np.inf, None
        A = chord_pts(a)
        for b in range(len(edges)):
            if b == a or b in taken or set(edges[a]) & set(edges[b]):
                continue
            d = A[:, None, :] - chord_pts(b, span)[None, :, :]
            d -= box_g * np.round(d / box_g)
            gap = float(np.linalg.norm(d, axis=2).min())
            if gap < best:
                best, out = gap, b
        return out

    plan, taken = [], set()
    for a in (0, 20):
        t = nearest_chord(a, taken)
        if t is None:
            continue
        taken.add(t)
        plan.append((a, {t: args.want}))
    pairs_g = [(r, t) for r, w in plan for t in w]

    if args.density:
        geo = geometry(graph, dp=args.dp, density=args.density)
    elif args.coil:
        geo = geometry(graph, dp=args.dp, bond=BOND, coil=args.coil)
    else:
        sc = scale_for_design(graph, pairs_g, dp=args.dp, bond=BOND,
                              radius=args.ring, margin=args.margin,
                              site_span=(args.site_lo, args.site_hi))
        geo = geometry(graph, dp=args.dp, scale=sc)
    n_beads = graph.number_of_edges() * args.dp + graph.number_of_nodes()
    rho = n_beads / float(np.prod(geo["L"]))
    print(f"  box {geo['L'][0]:.1f} sigma, density {rho:.3f}, "
          f"{n_beads} beads, crowding "
          f"{rho * 4 / 3 * np.pi * 0.9 ** 3:.2f} beads within 0.9 sigma")
    ch, ends = geo["chords"], geo["ends"]
    keys = sorted(ch)
    idx = {k: i + 1 for i, k in enumerate(keys)}

    rng = np.random.default_rng(args.seed)
    paths0 = {k: bridging_walk(c0, c1, args.dp + 1, BOND, rng)
              for k, (c0, c1) in ch.items()}

    # ---- 1 and 2: a plain melt, relaxed ---------------------------------
    root = OUT / "relaxed"
    _n, node_atom, chain_atoms = write_system(graph, geo, paths0, root)
    seq = {k: chain_ids(k, node_atom, chain_atoms, ends) for k in keys}
    sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                             protocol="hardcore")
    print("  --- relaxing a plain melt, no designed entanglements ---")
    run_md(sim, args.stages)
    after = {1: "system_after_soft.data", 2: "system_ramped.data",
             3: "system_equilibrated.data"}[args.stages]
    relaxed = root / "04_Simulation" / after
    if not relaxed.exists():
        print("  relaxation produced no output")
        return 1

    box, paths, xyz0 = paths_from(relaxed, keys, seq)
    geo_r = dict(geo)
    geo_r["L"] = box
    geo_r["chords"] = {k: (paths[k][0], paths[k][-1]) for k in keys}
    print(f"  relaxed: box {box[0]:.1f} sigma")

    # ---- 3: route into the relaxed coordinates --------------------------
    pairs = sorted({(min(r, t), max(r, t)) for r, w in plan for t in w})
    print(f"  plan: " + "; ".join(f"{r} with {list(w)[0]}" for r, w in plan))

    work = OUT / f"relaxed_work_{os.getpid()}"
    base = measure_pairs(relaxed, seq, pairs, work)
    print(f"  before routing: "
          + ", ".join(f"{a}-{b}={v}" for (a, b), v in base.items()))

    # A one-stage relaxation of a candidate, reusing the scripts already
    # generated for this system.
    cand_root = OUT / f"relaxed_cand_{os.getpid()}"
    (cand_root / "03_Conformation").mkdir(parents=True, exist_ok=True)
    (cand_root / "04_Simulation").mkdir(parents=True, exist_ok=True)
    (cand_root / "02_Chemistry").mkdir(parents=True, exist_ok=True)
    for f in (root / "02_Chemistry").glob("*"):
        (cand_root / "02_Chemistry" / f.name).write_text(f.read_text())
    for f in (root / "04_Simulation").glob("*.in"):
        (cand_root / "04_Simulation" / f.name).write_text(f.read_text())

    def relax_candidate(data_file, i=0):
        """Relax one candidate under the protocol it will be judged by.

        The same number of stages as the final check, not a cheaper screen.
        Selecting a placement for surviving one minimisation and then
        reporting it under a different one would be choosing against the wrong
        test, and the residual loss this is meant to remove is exactly the
        kind of thing that would slip through.
        """
        (cand_root / "03_Conformation" / "system_relaxed.data").write_text(
            Path(data_file).read_text())
        for stale in ("system_after_soft.data", "system_ramped.data",
                      "system_equilibrated.data"):
            q = cand_root / "04_Simulation" / stale
            if q.exists():
                q.unlink()
        try:
            run_md(cand_root / "04_Simulation", args.stages)
        except Exception:
            return None
        out = cand_root / "04_Simulation" / after
        return out if out.exists() else None

    # A template that carries the chains routed so far.
    #
    # Placements used to be judged against the pristine melt, so the screen for
    # the second pair never saw the first pair's new path. Measured: both
    # pairs screened at 2 on their own and the finished system read 3 and 1.
    # Designed windings interact, so each has to be chosen against what is
    # already committed.
    current = root / "04_Simulation" / "current.data"
    xyz_now = dict(xyz0)
    rewrite_coords(relaxed, current, xyz_now)
    # Two chains may be sent to the same partner, and winding both at the same
    # place puts the second ring inside the first. Measured: the second pair
    # built nothing at all. Spread the sites across the permitted stretch so
    # each has its own piece of the target.
    used = collections.Counter()
    for routed, want in plan:
        target, n = list(want.items())[0]
        # Everything except the chain being redrawn, by atom id rather than
        # by chain: a crosslink junction belongs to every chain meeting there,
        # so dropping whole chains still leaves the routed chain's own two
        # endpoints in the obstacle set. It then spends its accept test trying
        # to get away from the junctions it is pinned to, and the tightest
        # contact reads 0.00 no matter what the path does.
        mine = set(seq[routed])
        others = np.array([xyz_now[i] for i in sorted(xyz_now)
                           if i not in mine])
        avoid = Clearance(others, box, args.clearance)
        before = avoid.worst(paths[routed])
        if args.construct:
            k = used[target]
            used[target] += 1
            try:
                p, got = construct_exact(
                    paths, routed, target, n, box, avoid, seq, current, work,
                    (args.site_lo, args.site_hi), args.ring, args.dp,
                    relax=(relax_candidate if args.verify else None))
            except ValueError as e:
                print(f"    chain {routed} -> {target}: {e}")
                continue
            if p is None:
                print(f"    chain {routed} -> {target}: nothing built")
                continue
            if got != n:
                print(f"    chain {routed} -> {target}: asked {n}, "
                      f"closest placement gives {got}")
            paths[routed] = p
            for aid, xyzp in zip(seq[routed], p):
                xyz_now[aid] = xyzp
            rewrite_coords(relaxed, current, xyz_now)
            print(f"    chain {routed} -> {target}: built {got} by "
                  f"construction"
                  f"   (tightest contact {before:.2f} -> "
                  f"{avoid.worst(p[1:-1]):.2f} sigma)")
            continue
        best = route_pairwise(paths, keys, seq, geo_r, routed, target, n,
                              relaxed, args.rounds, args.per_round, rng, work,
                              relax=relax_candidate, avoid=avoid,
                              site_span=(args.site_lo, args.site_hi))
        if best is None:
            print(f"    chain {routed}: nothing built")
            continue
        paths[routed] = best[3]
        for aid, xyzp in zip(seq[routed], best[3]):
            xyz_now[aid] = xyzp
        print(f"    chain {routed} -> {target}: got {best[1]}, wanted {n}"
              f"   (tightest contact {before:.2f} -> "
              f"{avoid.worst(best[3]):.2f} sigma)")

    # ---- 4: write it back and minimise again ----------------------------
    new_xyz = {}
    for k in keys:
        for aid, xyzp in zip(seq[k], paths[k]):
            new_xyz[aid] = xyzp
    designed = root / "04_Simulation" / "designed.data"
    rewrite_coords(relaxed, designed, new_xyz)

    built = measure_pairs(designed, seq, pairs, work)
    print(f"\n  as designed:    "
          + ", ".join(f"{a}-{b}={v}" for (a, b), v in built.items()))

    root2 = OUT / "relaxed_again"
    (root2 / "04_Simulation").mkdir(parents=True, exist_ok=True)
    for f in (root / "04_Simulation").glob("*.in"):
        (root2 / "04_Simulation" / f.name).write_text(f.read_text())
    (root2 / "03_Conformation").mkdir(parents=True, exist_ok=True)
    (root2 / "03_Conformation" / "system_relaxed.data").write_text(
        Path(designed).read_text())
    for sub in ("02_Chemistry",):
        (root2 / sub).mkdir(parents=True, exist_ok=True)
        for f in (root / sub).glob("*"):
            (root2 / sub / f.name).write_text(f.read_text())

    print(f"\n  --- minimising again, from the designed relaxed system ---")
    run_md(root2 / "04_Simulation", args.stages)
    final = root2 / "04_Simulation" / after
    if not final.exists():
        print("  second minimisation produced no output")
        return 1
    print()
    report_bonds(root2)

    _b, xyz2, _ = read_data(final)
    ids = sorted(xyz0)
    d = np.array([xyz2[i] - new_xyz[i] for i in ids])
    d -= box * np.round(d / box)
    dist = np.linalg.norm(d, axis=1)
    print(f"\n  bead motion in this second minimisation: median "
          f"{np.median(dist):.2f}, max {dist.max():.2f} sigma")

    surv = measure_pairs(final, seq, pairs, work)
    print(f"\n  {'pair':>10} {'before':>7} {'designed':>9} {'after MD':>9}")
    ok = 0
    for a, b in pairs:
        want = next(n for r, w in plan for t, n in w.items()
                    if {r, t} == {a, b})
        ok += (surv[(a, b)] == want)
        print(f"  {f'{a}-{b}':>10} {str(base[(a, b)]):>7} "
              f"{str(built[(a, b)]):>9} {str(surv[(a, b)]):>9}"
              + ("   ok" if surv[(a, b)] == want else ""))
    print(f"\n  {ok} of {len(pairs)} survived exactly")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
