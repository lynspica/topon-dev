"""Ask for e entanglements per chain, get e entanglements per chain.

    python tests/workflows/entangle_density.py --want 2.0
    python tests/workflows/entangle_density.py --want 1.0 --shells 2:1.0

A system-averaged target, not a named pair. That difference is what makes it
robust: for a density it does not matter *which* pairs carry the entanglements,
only how many there are, so a pair that fails to survive can be replaced by
another and the collateral a routed chain picks up on its way counts toward the
total rather than against it.

So the loop closes on the measured value:

    select -> route -> relax -> measure the whole system -> top up or stop

Each round adds only the shortfall, and the measurement is always of the
relaxed system, so what is reported is what survives rather than what was
built. Survival of an individual pair stops mattering; the delivered density is
what is being controlled.

Under ``--mix-only`` the currency narrows: only the designed pairs are
watched, so collateral on undesigned pairs is not counted, the reported
total undercounts the system, and the flag exists for verifying the *mix*
cheaply rather than the density exactly.
"""
from __future__ import annotations

import argparse
import json
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
    select_by_shells,
)
from topon.conformation.paths import Clearance, bridging_walk  # noqa: E402
from tests.workflows.entangle_all import CASES  # noqa: E402
from tests.workflows.entangle_relaxed import (  # noqa: E402
    construct,
    paths_from,
    rewrite_coords,
)
from tests.workflows.entangle_search import measure_batch  # noqa: E402
from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    DP,
    LATTICE,
    OUT,
    build_network,
    chain_ids,
    conform_and_script,
    geometry,
    run_md,
    scale_for_design,
    write_system,
    z1_export,
)


def keep_pairs_with_counts(plan, wanted):
    """The plan entries for a set of pairs, in the plan's own order."""
    want = set(wanted)
    return [q for q in plan if (q[0], q[1]) in want]


def close_pairs(paths, box, cutoff=3.0):
    """Chain pairs that come within ``cutoff`` of each other.

    Two chains that never come close cannot be entangled, so restricting the
    measurement to these pairs is exact rather than a sample, and it is what
    makes the measurement affordable: a 648-chain system has 210k pairs and
    only a few thousand are ever in contact.

    Found with a KD-tree over the beads rather than chain by chain. The
    pairwise form is quadratic in the chain count and quadratic again in the
    beads per chain, which is minutes at 648 chains; this is seconds.
    """
    from scipy.spatial import cKDTree

    keys = sorted(paths)
    owner, pts = [], []
    for k in keys:
        q = paths[k]
        q = q - box * np.floor(q / box)
        q = np.clip(q, 0.0, np.nextafter(box, 0.0))
        pts.append(q)
        owner += [k] * len(q)
    tree = cKDTree(np.vstack(pts), boxsize=box)
    owner = np.asarray(owner)

    out = set()
    for i, j in tree.query_pairs(cutoff):
        a, b = owner[i], owner[j]
        if a != b:
            out.add((min(a, b), max(a, b)))
    return sorted(out)


def measure_system(data_file, seq, keys, work, watch):
    """Entanglement points per chain, over a fixed set of pairs.

    ``watch`` is decided once and reused for every round, so the before and
    after numbers are a paired comparison of the same pairs and the difference
    between them is the increment the design is responsible for.

    Measuring every pair in contact would be exact and is unaffordable: a
    melt at density 0.85 has 14754 of them, each needing its own export. The
    set is instead every designed pair, which is where the increment is meant
    to appear, plus a fixed random sample of the rest to carry the background
    and any collateral.

    Returns ``(points_per_chain, partners_per_chain, measured, unmeasured)``,
    scaled from the sample back to the whole system.
    """
    from tests.workflows.entangle_relaxed import measure_pairs

    pairs, scale = watch
    if not pairs:
        return 0.0, 0.0, 0, 0
    got = measure_pairs(data_file, seq, pairs, work)
    blind = sum(1 for v in got.values() if v is None)

    pts = partners = 0.0
    for q in pairs:
        v = got.get(q)
        if v is None:
            continue
        w = scale.get(q, 1.0)
        pts += w * v
        partners += w * (1 if v else 0)
    n = len(keys)
    return 2.0 * pts / n, 2.0 * partners / n, len(pairs) - blind, blind


def routed_watch(paths, box, chains, cutoff=3.0):
    """Every contact pair involving a chain that will be routed.

    Only routed chains change, so this is the whole signal and it is measured
    exactly. Sampling the rest and scaling it up is what made the estimate
    unusable: at melt density there are 14754 contact pairs, a 250-pair sample
    scales each one by 59, and over 106 chains a single fluctuation in the
    sample moves the reported average by 1.1 per chain. The first round of one
    run came back at -0.197 per pair, implying routing *removed* entanglement,
    when measuring a routed chain directly shows it going from 11 points over
    4 partners to 28 over 12.
    """
    want = set(chains)
    return [q for q in close_pairs(paths, box, cutoff)
            if q[0] in want or q[1] in want]


def watch_set(paths, box, designed, cutoff=3.0, sample=400, rng=None):
    """The pairs to measure, and what each one stands for.

    Every designed pair counts for itself. The rest of the contact pairs are
    represented by a random sample, each standing for ``len(rest)/len(sample)``
    of them, so the total is an unbiased estimate of the background rather
    than a count of the part that happened to be looked at.
    """
    rng = np.random.default_rng() if rng is None else rng
    contact = set(close_pairs(paths, box, cutoff))
    named = {(min(a, b), max(a, b)) for a, b in designed}
    contact |= named
    rest = sorted(contact - named)
    if len(rest) > sample:
        pick = rng.choice(len(rest), size=sample, replace=False)
        chosen = [rest[i] for i in sorted(pick)]
        w = len(rest) / float(sample)
    else:
        chosen, w = rest, 1.0
    pairs = sorted(named) + chosen
    scale = {q: 1.0 for q in named}
    scale.update({q: w for q in chosen})
    return pairs, scale


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--want", type=float, default=2.0,
                    help="target entanglements per chain")
    ap.add_argument("--shells", default="1:0.5,2:0.5",
                    help="shell mix, e.g. '1:0.2,2:0.5,3:0.25,4:0.05'")
    ap.add_argument("--bias", default="uniform",
                    choices=("uniform", "region", "anti_region", "gradient",
                             "clusters"),
                    help="spatial placement. The shell mix says which "
                         "neighbour a chain winds; this says where in the "
                         "box the wound pairs sit.")
    ap.add_argument("--bias-params", default="{}",
                    help="JSON for the bias, e.g. "
                         "'{\"axis\":\"z\",\"strength\":3}' for gradient, "
                         "'{\"center\":[0.5,0.5,0.5],\"radius\":0.3,"
                         "\"strength\":8}' for region")
    ap.add_argument("--rounds", type=int, default=6)
    ap.add_argument("--reject-below", type=float, default=0.25,
                    help="drop a routed path whose tightest contact is below "
                         "this. 0.55 was too strict: at melt density no path "
                         "clears it and every round dropped 25 of 25, while "
                         "paths built at 0.14 to 0.25 sigma have repeatedly "
                         "minimised fine. What actually blows LAMMPS up is "
                         "near-zero separation, and 0.25 keeps that out "
                         "without vetoing the whole density regime.")
    ap.add_argument("--headroom", type=float, default=3.0,
                    help="how many times the expected pair count to watch. "
                         "Covers the back-off restarting with a different "
                         "number without measuring the whole selection.")
    ap.add_argument("--min-pairs", type=int, default=25,
                    help="pairs that must be routed before the measured yield "
                         "is trusted. Below this the background sampling "
                         "error swamps the signal and the estimate can even "
                         "come out negative.")
    ap.add_argument("--pair-yield", type=float, default=0.10,
                    help="yield assumed before enough pairs are routed to "
                         "measure it")
    ap.add_argument("--calibrate", type=float, default=0.15,
                    help="fraction of the selected pairs to route in the "
                         "first batch, before what a pair is worth in this "
                         "system has been measured")
    ap.add_argument("--gain", type=float, default=0.7,
                    help="fraction of each shell's estimated need routed per "
                         "round. Below 1 approaches the target from below, "
                         "which matters because the per-pair yield is "
                         "ask-dependent: routing a shell heavily inflates "
                         "its own yield through same-shell collateral, so "
                         "meeting the full shortfall at the current estimate "
                         "overshoots and has to be walked back.")
    ap.add_argument("--tol-pct", type=float, default=5.0,
                    help="stop when within this percent of the target. Scales "
                         "with what was asked for, where an absolute "
                         "tolerance does not: 0.15 is 7.5 percent at e=2 and "
                         "3.75 at e=4. The floor is one entanglement point on "
                         "one pair, 2/num_chains per chain, so a larger "
                         "network can be asked for a tighter percentage.")
    ap.add_argument("--tol", type=float, default=None,
                    help="absolute tolerance in entanglements per chain, "
                         "overriding --tol-pct")
    ap.add_argument("--dims", type=int, default=4)
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--ring", type=float, default=2.0)
    ap.add_argument("--density", type=float, default=None,
                    help="melt route: pack to this density instead of sizing "
                         "by coil ratio. 0.85 is the physical melt the LAMMPS "
                         "scripts are calibrated for.")
    ap.add_argument("--max-density", type=float, default=0.85,
                    help="ceiling for auto-sizing. The box is sized from the "
                         "design, which gives the roomiest box it fits in, "
                         "and this only bites when a design would otherwise "
                         "be crushed into a smaller cell than the system is "
                         "meant to be built at.")
    ap.add_argument("--coil", type=float, default=None,
                    help="contour over chord, which sets the box. 6 puts "
                         "neighbouring chains 1.5 sigma apart with real free "
                         "volume to route through; 1.8 leaves them 35 sigma "
                         "apart and unreachable, 12 leaves no room at all.")
    ap.add_argument("--margin", type=float, default=2.0)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--clearance", type=float, default=0.9)
    ap.add_argument("--cutoff", type=float, default=3.0,
                    help="chains further apart than this everywhere cannot be "
                         "entangled, so they are not measured")
    ap.add_argument("--mix-only", action="store_true",
                    help="watch only the designed pairs. The e total then "
                         "counts designed increments rather than the whole "
                         "watched neighbourhood, so collateral on undesigned "
                         "pairs is not counted and the delivered *system* "
                         "density will exceed --want. The run is Z1-cheap: "
                         "the full watch at melt density is 4000+ pairs and "
                         "hours of Z1, the plan is a couple of hundred. The "
                         "shell table is unchanged -- it was always built "
                         "from the designed pairs alone.")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tag", default=None,
                    help="output directory suffix. Defaults to the lattice "
                         "size and target, so two runs at different sizes do "
                         "not delete each other's files -- which is exactly "
                         "what happened when a dims 5 run was started while a "
                         "dims 6 run was still going.")
    args = ap.parse_args()

    if args.tag is None:
        # The process id, not just the parameters.
        #
        # Keying on lattice size and target stopped two *different* runs
        # colliding and did nothing for two of the *same* kind, which promptly
        # happened: a second dims 5 run at e=2 deleted the first one's files
        # mid-flight and it died on FileNotFoundError. Twice in one session.
        args.tag = f"_d{args.dims}_e{args.want:g}_{os.getpid()}"

    mix = {}
    for part in args.shells.split(","):
        s, f = part.split(":")
        mix[int(s)] = float(f)

    # Vet the bias request before the expensive work, not after: a typo'd
    # axis used to bias along x (compute_bias_weights falls back silently)
    # and then crash in the report, after the whole build.
    bias_params = json.loads(args.bias_params)
    if args.bias == "gradient" and bias_params.get("axis", "x") not in (
            "x", "y", "z"):
        print(f"  unknown gradient axis "
              f"{bias_params.get('axis')!r}; use x, y or z")
        return 1
    if (args.bias in ("region", "anti_region")
            and float(bias_params.get("radius", 0.3)) > 0.5):
        # Beyond half the box the minimum-image "inside" test is no longer
        # a sphere and the volume formula exceeds the cell, which turns the
        # reported density ratio negative.
        print("  region radius is a fraction of min(box); above 0.5 the "
              "minimum-image sphere stops being a sphere. Ask for 0.5 or "
              "less.")
        return 1

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (args.dims,) * 3
    graph = build_network(spec)
    dims = np.asarray(graph.graph["box"], float)
    G = nx.MultiGraph()
    G.add_nodes_from(graph.nodes(data=True))
    for u, v in graph.edges():
        G.add_edge(u, v)

    dist = chain_distances(G, dims)
    shells = neighbour_shells(G, dims, max_shell=max(mix), distances=dist)
    edges = sorted(graph.edges())
    idx_of = {frozenset(e): i for i, e in enumerate(edges)}

    def to_idx(sel):
        out = []
        for e1, e2, c in sel:
            a = idx_of.get(frozenset((e1[0], e1[1])))
            b = idx_of.get(frozenset((e2[0], e2[1])))
            if a is not None and b is not None and a != b:
                out.append((min(a, b), max(a, b), c))
        return out

    rng = np.random.default_rng(args.seed)
    plan = to_idx(select_by_shells(G, args.want, mix, dims, shells=shells,
                                   rng=rng, bias_kind=args.bias,
                                   bias_params=bias_params))
    if not plan:
        print("  selection produced nothing")
        return 1
    print(f"  target {args.want:.2f} per chain, mix "
          f"{ {s: mix[s] for s in sorted(mix)} }"
          + (f", bias {args.bias} {bias_params}"
             if args.bias != "uniform" else ""))
    print(f"  selected {len(plan)} pairs, "
          f"{sum(c for _a, _b, c in plan)} entanglements")

    # Sized by coil ratio, not by the design.
    #
    # `scale_for_design` sizes the box so every requested route fits, which is
    # right for a handful of named pairs and wrong here: with a hundred pairs
    # some chain has several partners, its ring and detour budget dominates,
    # and the box collapses -- measured, 2.2 sigma at density 843. A density
    # target does not need every pair to build. It needs a physically sensible
    # box and enough pairs, and the loop replaces whatever does not fit.
    # Sized from the design, with the build density as a ceiling.
    #
    # This is the point of auto-sizing and it took three reminders to get
    # right: the largest box the design fits in is the *lowest* density, and
    # 0.85 is a cap that bites only when a design would crush the system.
    # Pinning every run to 0.85 instead put the routing in the worst case it
    # could be in -- a point at that density lies within 0.9 sigma of 2.6
    # beads, so there is no free volume to route through, and the result was
    # paths landing on beads, a pair energy of 4e9, and chains longer than the
    # periodic cell.
    #
    # Sized on the pairs that will actually be routed rather than the whole
    # selection, since the rest are spares the loop may never reach.
    n_likely = int(args.headroom * args.want / max(args.pair_yield, 1e-9))
    likely = plan[:max(n_likely, args.min_pairs)]
    if args.density:
        geo = geometry(graph, dp=args.dp, density=args.density)
    elif args.coil:
        geo = geometry(graph, dp=args.dp, bond=BOND, coil=args.coil)
    else:
        sc = scale_for_design(graph, [(a, b) for a, b, _c in likely],
                              dp=args.dp, bond=BOND, radius=args.ring,
                              margin=args.margin, site_span=(0.3, 0.7),
                              max_density=args.max_density)
        geo = geometry(graph, dp=args.dp, scale=sc)
    keys = sorted(geo["chords"])
    n_beads = graph.number_of_edges() * args.dp + graph.number_of_nodes()
    tol = (args.tol if args.tol is not None
           else args.want * args.tol_pct / 100.0)
    floor = 2.0 / len(sorted(geo["chords"]))
    print(f"  tolerance {tol:.3f} per chain"
          + (f" ({args.tol_pct:.1f} percent of {args.want})"
             if args.tol is None else "")
          + f"; one entanglement on one pair is {floor:.3f}, the finest step "
            f"this network can take")
    print(f"  box {geo['L'][0]:.1f} sigma, density "
          f"{n_beads / float(np.prod(geo['L'])):.3f}")

    rng0 = np.random.default_rng(args.seed)
    paths0 = {k: bridging_walk(c0, c1, args.dp + 1, BOND, rng0)
              for k, (c0, c1) in geo["chords"].items()}
    root = OUT / f"relaxed{args.tag}"
    shutil.rmtree(root, ignore_errors=True)
    _n, node_atom, chain_atoms = write_system(graph, geo, paths0, root)
    seq = {k: chain_ids(k, node_atom, chain_atoms, geo["ends"]) for k in keys}
    sim = conform_and_script(root, graph, geo, pair_style="repulsive",
                             protocol="hardcore")
    print("\n  --- relaxing a plain melt ---")
    run_md(sim, args.stages)
    after = {1: "system_after_soft.data", 2: "system_ramped.data",
             3: "system_equilibrated.data"}[args.stages]
    relaxed = root / "04_Simulation" / after
    if not relaxed.exists():
        print("  relaxation produced no output")
        return 1
    box, paths, xyz0 = paths_from(relaxed, keys, seq)

    work = OUT / f"density_work_{os.getpid()}"
    # Watch every pair that touches a chain the plan might route, and measure
    # them exactly. The rest of the system does not move, so it contributes the
    # same amount before and after and cancels out of the difference.
    # Only the chains that will plausibly be routed, not every selected pair.
    #
    # The plan selects far more pairs than the target needs -- 380 for a
    # request of 2.0 per chain that takes about 20 -- and watching both members
    # of all of them measured 12134 contact pairs to observe twenty chains. A
    # round then took longer than the whole run should. The pool is ordered, so
    # the pairs that will be used are at the front of it, and the headroom
    # covers the back-off restarting with a different count.
    may_route = {a for a, _b, _c in likely} | {b for _a, b, _c in likely}
    if args.mix_only:
        wp = sorted({(min(a, b), max(a, b)) for a, b, _c in plan})
        print(f"  mix-only: watching the {len(wp)} designed pairs and "
              f"nothing else")
    else:
        wp = routed_watch(paths, box, may_route, args.cutoff)
        print(f"  watching {len(wp)} contact pairs that touch a routable "
              f"chain, measured exactly")
    watch = (wp, {q: 1.0 for q in wp})
    base = measure_system(relaxed, seq, keys, work, watch)
    if args.mix_only:
        # Not a statement about the melt: the watch is the designed plan,
        # chosen by graph shell rather than by distance, so this is only the
        # baseline those particular pairs will be differenced against.
        print(f"  the designed pairs carry {base[0]:.2f} per chain in the "
              f"plain melt ({base[2]} measured"
              + (f", {base[3]} unmeasured)" if base[3] else ")"))
    else:
        print(f"  the plain melt already carries {base[0]:.2f} per chain "
              f"over {base[1]:.1f} partners ({base[2]} pairs in contact"
              + (f", {base[3]} unmeasured)" if base[3] else ")"))

    cand = OUT / f"density_cand_{os.getpid()}"
    (cand / "03_Conformation").mkdir(parents=True, exist_ok=True)
    for sub in ("02_Chemistry", "04_Simulation"):
        (cand / sub).mkdir(parents=True, exist_ok=True)
        for f in (root / sub).glob("*"):
            if f.is_file() and (sub == "02_Chemistry" or f.suffix == ".in"):
                (cand / sub / f.name).write_text(f.read_text())

    def relax(xyz):
        rewrite_coords(relaxed, cand / "03_Conformation" /
                       "system_relaxed.data", xyz)
        for stale in ("system_after_soft.data", "system_ramped.data",
                      "system_equilibrated.data"):
            q = cand / "04_Simulation" / stale
            if q.exists():
                q.unlink()
        try:
            run_md(cand / "04_Simulation", args.stages)
        except Exception:
            return None
        out = cand / "04_Simulation" / after
        return out if out.exists() else None

    # Coordinates as an array with a row index, not a dict.
    #
    # The obstacle set is rebuilt for every routed chain, and doing it by dict
    # comprehension over every atom costs more than the routing: 52k lookups
    # times 380 chains is 20 million, and the first run stalled in round one.
    # With rows, excluding a chain is a boolean mask.
    ids = sorted(xyz0)
    row = {a: i for i, a in enumerate(ids)}
    xyz_arr = np.array([xyz0[a] for a in ids])

    xyz_now = dict(xyz0)
    done, routed_chains = set(), set()

    # A calibrated controller, not a single shot.
    #
    # Routing every selected pair at once overshoots badly: measured, 102 pairs
    # asked for 2.00 per chain and delivered 5.48. A designed pair is worth
    # more than the entanglements it was asked for, because the routed chain
    # picks others up on its way and for a density those count.
    #
    # So route a small batch first, measure what a pair is actually worth in
    # this system, and size every later batch from that. The yield is
    # re-estimated each round, so it corrects itself rather than trusting the
    # first estimate.
    # Which shell each pair sits in, built before the loop so the mix can be
    # steered rather than only reported afterwards.
    shell_of = {}
    for chain, by in shells.items():
        ia = idx_of.get(frozenset((chain[0], chain[1])))
        if ia is None:
            continue
        for sh, others in by.items():
            for o in others:
                ib = idx_of.get(frozenset((o[0], o[1])))
                if ib is not None and ia != ib:
                    shell_of.setdefault((min(ia, ib), max(ia, ib)), sh)

    pool = list(likely) + [q for q in plan if q not in likely]

    # The full population of each requested shell, as a reserve.
    #
    # The plan's pool holds pairs in proportion to the ask, but the
    # controller consumes them in proportion to shortfall over yield, which
    # for a low-yield shell runs several times deeper: every first-generation
    # controller run ended on "pool dry" with the mix unconverged. The
    # reserve is every pair the geometry offers at that shell distance, so a
    # shell's pool now runs out only when the network itself does.
    rng_res = np.random.default_rng(args.seed + 1)
    in_plan = {(a, b) for a, b, _c in plan}
    reserve_by_shell = {}
    for q, sh_r in shell_of.items():
        if sh_r in mix and q not in in_plan:
            reserve_by_shell.setdefault(sh_r, []).append((q[0], q[1], 1))
    for sh_r, lst in reserve_by_shell.items():
        reserve_by_shell[sh_r] = [lst[i]
                                  for i in rng_res.permutation(len(lst))]

    added, yield_per, raw, reliable = 0.0, None, 0.0, False
    print(f"\n  {'round':>6} {'routed':>7} {'added':>7} {'per pair':>9} "
          f"{'target':>7}")
    # The mix is a target per shell, not just a selection rule.
    #
    # Selecting the right fraction of *pairs* does not deliver the right
    # fraction of *entanglements*: a pair in an outer shell is worth more,
    # because the routed chain travels further and picks up more on the way.
    # Asking 0.20 / 0.50 / 0.25 delivered 0.20 / 0.60 / 0.20. So each shell
    # gets its own target and its own shortfall, and a round routes from the
    # shells that are behind rather than from the pool in order.
    nrm = sum(mix.values())
    want_by_shell = {sh: args.want * f / nrm for sh, f in mix.items()}
    got_by_shell = {sh: 0.0 for sh in mix}

    # In-loop mix controller state. The ask-to-delivered map is nonlinear --
    # measured on the four-shell request, asking 20/50/25/5 delivers
    # 6/71/14/10, and dividing the ask by that yield over-corrects to
    # 52/24/21/2 -- so no single ask hits the target. What converged by hand
    # was iterating with damping (two calibration passes and a secant step
    # landed within 0.06); this loop is that iteration, run per round inside
    # one build: per-shell yields learned with damping, each round routing a
    # gain-scaled shortfall as a *mixed* batch, and overshooting shells
    # giving pairs back.
    per_shell_tol = tol / max(len(mix), 1)
    yield_pair = {sh: None for sh in mix}    # per routed pair, per-chain units
    inc_of = {}                              # pair -> measured increment
    base_cache = {}                          # pair -> melt baseline count
    dead = set()                             # shells measured undeliverable
    n_ch = len(keys)

    for rnd in range(1, args.rounds + 1):
        routed_by_shell = {sh: 0 for sh in mix}
        for q in done:
            s = shell_of.get(q)
            if s in routed_by_shell:
                routed_by_shell[s] += 1

        # Give back from shells that overshot, before asking for more.
        #
        # Un-routing restores the chain's original melt path, which was
        # clear when everything else was in place; other routed chains have
        # moved since, so the restore is checked against the current
        # obstacle set and skipped if it would land closer than the same
        # threshold a new route must clear.
        unrouted = 0
        if rnd > 1 and done:
            # Overshoot judged on fractions of the delivered total in
            # mix-only mode, absolutes otherwise. The deliverable total in
            # the roomy regime is far below the absolute want, so no shell
            # ever exceeded its absolute target and the give-back never
            # fired -- while shell 1 sat at 0.33 of the total against an
            # asked 0.20.
            frac_tol = max(0.03, per_shell_tol / max(args.want, 1e-9))

            def over(sh):
                if args.mix_only:
                    t = sum(got_by_shell.values())
                    if t <= 0:
                        return False
                    return (got_by_shell[sh] / t
                            > want_by_shell[sh] / args.want + frac_tol)
                return got_by_shell[sh] - want_by_shell[sh] > per_shell_tol

            for sh in sorted(mix):
                if not over(sh):
                    continue
                # Damped: remove gain-worth of the excess, not all of it.
                # Walking the shell all the way down chased its target
                # fraction of a shrinking total -- one round stripped shell
                # 2 from 0.34 to 0.06 per chain, the allocation flooded the
                # others, and the mix oscillated instead of settling.
                t0 = sum(got_by_shell.values())
                excess = (got_by_shell[sh]
                          - want_by_shell[sh] / args.want * t0
                          if args.mix_only and t0 > 0
                          else got_by_shell[sh] - want_by_shell[sh])
                budget = args.gain * max(excess, 0.0)
                # Not named `cand`: that is the candidate working directory
                # relax() closes over, and shadowing it here turned the
                # first give-back into a crash on the next relax call.
                gb = sorted(((inc_of.get(q, 0.0), q) for q in list(done)
                             if shell_of.get(q) == sh), reverse=True)
                for v, q in gb:
                    if budget <= 0 or not over(sh):
                        break
                    if v <= 0:
                        continue
                    a = q[0]
                    keep = np.ones(len(ids), bool)
                    keep[[row[i] for i in seq[a]]] = False
                    avoid = Clearance(xyz_arr[keep], box, args.clearance)
                    if avoid.worst(paths0[a][1:-1]) < args.reject_below:
                        continue
                    paths[a] = paths0[a].copy()
                    for aid, xyzp in zip(seq[a], paths0[a]):
                        xyz_now[aid] = xyzp
                        xyz_arr[row[aid]] = xyzp
                    done.discard(q)
                    routed_chains.discard(a)
                    routed_by_shell[sh] -= 1
                    pool += keep_pairs_with_counts(plan, [q])
                    got_by_shell[sh] -= v
                    budget -= v
                    unrouted += 1
            if unrouted:
                print(f"  {'':>6} gave back {unrouted} pairs from "
                      f"overshooting shells")

        # Each shell's ask: its shortfall over what one of its pairs is
        # worth, gain-scaled, routed as one mixed batch.
        floor_y = 2.0 / n_ch
        short_amt = {sh: want_by_shell[sh] - got_by_shell[sh] for sh in mix}
        tot_del = sum(got_by_shell.values())
        need = {}
        for sh in mix:
            if sh in dead or short_amt[sh] <= 0.5 * per_shell_tol:
                continue
            # A shell at or above its fraction of the delivered total does
            # not get more routing in mix-only mode, even though it is
            # below its absolute want: routing it fights the give-back.
            if (args.mix_only and tot_del > 0
                    and got_by_shell[sh] / tot_del
                    >= want_by_shell[sh] / args.want):
                continue
            y = (yield_pair[sh] if yield_pair[sh] and yield_pair[sh] > 0
                 else max(args.pair_yield, floor_y))
            need[sh] = max(args.gain * short_amt[sh] / max(y, 1e-4), 1.0)
        cap = (max(args.min_pairs, int(args.calibrate * len(pool)))
               if rnd == 1 else 50)
        total_need = sum(need.values())
        scale = min(1.0, cap / total_need) if total_need > 0 else 0.0
        quota = {sh: int(need[sh] * scale) for sh in need}
        for sh in sorted(need, key=lambda s: -(need[s] * scale
                                               - quota[s])):
            if sum(quota.values()) >= min(cap, int(total_need)):
                break
            quota[sh] += 1

        # Both chains of every delivered pair are off limits, not just the
        # routed one. Routing the *partner* of a delivered pair for some
        # other pair moves it and disturbs the winding it already carries:
        # the outward run's shell 4 fell from 0.49 of the total to 0.38
        # across the late rounds with no give-back ever touching it.
        locked = {c for q in done for c in q} | routed_chains
        take, seen_q = [], set()
        pool_dry = []
        for sh, nq in sorted(quota.items()):
            got_n = 0
            for q in pool:
                if got_n >= nq:
                    break
                if q in seen_q or shell_of.get((q[0], q[1])) != sh:
                    continue
                if q[0] in locked or q[1] in locked:
                    continue
                take.append(q)
                seen_q.add(q)
                got_n += 1
            # Top up from the shell's reserve when the plan's pool is dry.
            extra = reserve_by_shell.get(sh, [])
            while extra and got_n < nq:
                q = extra.pop()
                if (q in seen_q or (q[0], q[1]) in done
                        or q[0] in locked or q[1] in locked):
                    continue
                take.append(q)
                seen_q.add(q)
                got_n += 1
            if got_n < nq:
                pool_dry.append((sh, nq - got_n))
        if pool_dry:
            print(f"  {'':>6} pool dry: "
                  + ", ".join(f"shell {s} short {m} pairs"
                              for s, m in pool_dry))
        pool = [q for q in pool if q not in seen_q]
        skipped = 0
        if not take and not unrouted:
            print(f"  no pairs left to route; delivered {added:.2f} against "
                  f"{args.want:.2f}")
            break
        if quota:
            print(f"  {'':>6} routing "
                  + ", ".join(f"shell {s}: {quota[s]}"
                              for s in sorted(quota)))

        for a, b, count in take:
            # One routing per chain, ever.
            #
            # `done` tracks pairs, so a chain appearing in two of them was
            # routed twice and the second path replaced the first, destroying
            # the entanglement it carried and moving coordinates that later
            # chains had been routed around. Including junction-sharing pairs
            # made this common, since a chain then has five more partners:
            # the yield fell 0.247, 0.101, 0.038 over three rounds and the
            # relaxation then failed with beads at zero separation.
            if a in locked or b in locked:
                continue
            keep = np.ones(len(ids), bool)
            keep[[row[i] for i in seq[a]]] = False
            avoid = Clearance(xyz_arr[keep], box, args.clearance)
            try:
                p = construct(paths, a, b, max(0.5, 0.5 * count), box, avoid,
                              radius=args.ring, dp=args.dp, span=(0.3, 0.7))
            except ValueError:
                continue
            # Refuse a path that would make the minimisation unstable.
            #
            # At melt density there is very little free volume -- a point lies
            # within 0.9 sigma of 2.6 beads on average -- so `taut_leg` often
            # cannot find a clearing step and falls back to the roomiest one
            # available, which may still be an overlap. One such path is
            # enough: LAMMPS came back with a pair energy of 4e9 and stopped,
            # losing the whole round. Committing only paths that clear costs a
            # few candidates and keeps the system runnable.
            room = avoid.worst(p[1:-1])
            if room < args.reject_below:
                skipped += 1
                continue

            paths[a] = p
            for aid, xyzp in zip(seq[a], p):
                xyz_now[aid] = xyzp
                xyz_arr[row[aid]] = xyzp
            done.add((a, b))
            routed_chains.add(a)
            locked.add(a)
            locked.add(b)

        if skipped:
            print(f"  {'':>6} {skipped} of {len(take)} paths dropped for "
                  f"landing closer than {args.reject_below} sigma")
        out = relax(xyz_now)
        if out is None:
            print(f"  {rnd:>6}   relaxation failed")
            return 1
        got = measure_system(out, seq, keys, work, watch)
        added = got[0] - base[0]

        # Per-shell delivery, so the next round can correct the mix.
        #
        # Baselines are measured once per pair and cached: the melt does not
        # change between rounds, and re-measuring it doubled the Z1 cost of
        # every round.
        from tests.workflows.entangle_relaxed import measure_pairs as _mp
        if done:
            fresh = [q for q in sorted(done) if q not in base_cache]
            if fresh:
                base_cache.update(_mp(relaxed, seq, fresh, work))
            now = _mp(out, seq, sorted(done), work)
            got_by_shell = {sh: 0.0 for sh in mix}
            inc_of = {}
            for q in sorted(done):
                sh = shell_of.get(q)
                v, v0 = now.get(q), base_cache.get(q)
                if sh in got_by_shell and v is not None and v0 is not None:
                    inc_of[q] = 2.0 * max(0, v - v0) / n_ch
                    got_by_shell[sh] += inc_of[q]

            # What a pair in each shell is worth here, damped: half the new
            # measurement, half the running estimate, so one noisy round
            # does not slew the allocation.
            routed_by_shell = {sh: 0 for sh in mix}
            for q in done:
                s = shell_of.get(q)
                if s in routed_by_shell:
                    routed_by_shell[s] += 1
            for sh in mix:
                if routed_by_shell[sh] > 0:
                    y_new = got_by_shell[sh] / routed_by_shell[sh]
                    yield_pair[sh] = (y_new if yield_pair[sh] is None
                                      else 0.5 * y_new + 0.5 * yield_pair[sh])
                    if (sh not in dead and routed_by_shell[sh] >= 8
                            and got_by_shell[sh] <= 0):
                        dead.add(sh)
                        print(f"  {'':>6} shell {sh}: {routed_by_shell[sh]} "
                              f"pairs routed, nothing delivered; no longer "
                              f"asking")

        # A yield estimate has to be a measurement, not one noisy difference.
        #
        # Three things went wrong when it was not. The first round measured
        # -0.197 per pair -- routing eight chains *removed* entanglement,
        # because rearranging a chain destroys melt entanglement as well as
        # adding designed entanglement, and over a handful of pairs that can
        # dominate. The loop then divided the shortfall by that negative
        # number, and the back-off shrank the pool until nothing was left:
        # -0.197, 0.186, 0.480, 0.623, 0.464 over five rounds, ending on three
        # pairs and a mix table describing those three.
        #
        # So: no estimate until enough pairs have been routed for the signal
        # to exceed the sampling error of the background, and no dividing by a
        # yield that is not positive.
        raw = added / max(len(done), 1)
        reliable = len(done) >= args.min_pairs and added > 0
        if reliable:
            yield_per = raw
        elif yield_per is None:
            yield_per = args.pair_yield
        print(f"  {rnd:>6} {len(done):>7} {added:>7.2f} {raw:>9.3f} "
              f"{args.want:>7.2f}"
              + ("" if reliable else "   (too few pairs to trust)")
              + (f"   ({got[3]} unmeasured)" if got[3] else ""))

        miss = {sh: want_by_shell[sh] - got_by_shell[sh] for sh in mix}
        tot_del = sum(got_by_shell.values())
        print(f"  {'':>6} by shell: "
              + ", ".join(
                  f"{sh}: {got_by_shell[sh]:.2f}/{want_by_shell[sh]:.2f}"
                  + (f" ({got_by_shell[sh] / tot_del:.2f})"
                     if tot_del > 0 else "")
                  for sh in sorted(mix)))

        # In mix-only mode the run is done when the *fractions* are right:
        # the absolute total is regime-bound and no absolute test can pass.
        if args.mix_only and tot_del > 0 and len(done) >= args.min_pairs:
            frac_tol = max(0.03, per_shell_tol / max(args.want, 1e-9))
            if all(abs(got_by_shell[sh] / tot_del
                       - want_by_shell[sh] / args.want) <= frac_tol
                   for sh in mix if sh not in dead):
                print(f"\n  mix converged: every live shell within "
                      f"{frac_tol:.2f} of its fraction over {len(done)} "
                      f"pairs, {tot_del:.2f} per chain delivered")
                break

        # Done when every shell is within tolerance, not just the total.
        if all(abs(m) <= per_shell_tol for m in miss.values()):
            print(f"\n  every shell within {per_shell_tol:.3f} per chain")
            break
        if abs(added - args.want) <= tol and not any(
                m > per_shell_tol for m in miss.values()):
            print(f"\n  delivered {added:.2f} per chain against a target of "
                  f"{args.want:.2f}, within {tol:.3f}")
            break
        if (added > args.want + tol and reliable
                and not any(m < -per_shell_tol for m in miss.values())):
            # Total over but no shell over: the excess is collateral outside
            # the designed pairs (full-watch mode only; in mix-only the
            # total is the sum of the shells). Per-shell give-back cannot
            # reach it, so start again from the plain melt with the right
            # number of pairs, which the measured yield now gives exactly.
            n_want = max(1, int(args.want / max(yield_per, 1e-9)))
            if n_want >= len(done) or rnd >= args.rounds:
                print(f"\n  overshot: {added:.2f} against {args.want:.2f}; a "
                      f"pair is worth {yield_per:.3f} here")
                break
            print(f"  {'':>6} backing off to {n_want} pairs, since a pair is "
                  f"worth {yield_per:.3f}")
            keep_pairs = sorted(done)[:n_want]
            paths.update({k: v.copy() for k, v in paths0.items()})
            xyz_now = dict(xyz0)
            xyz_arr[:] = np.array([xyz0[a] for a in ids])
            done, routed_chains = set(), set()
            inc_of, base_cache = {}, {}
            got_by_shell = {sh: 0.0 for sh in mix}
            pool = keep_pairs_with_counts(plan, keep_pairs) + pool
            continue
    else:
        print(f"\n  did not reach the target in {args.rounds} rounds; "
              f"delivered {added:.2f} against {args.want:.2f}")

    # The delivered shell mix, against the requested one.
    #
    # Selection produces the requested mix by construction; whether the built
    # and relaxed system carries it is a separate question, and the only one
    # that matters. Counted over the designed pairs that survived, since a
    # background entanglement was not placed in any shell by anybody.
    from tests.workflows.entangle_relaxed import measure_pairs
    shell_of = {}
    for chain, by in shells.items():
        a = idx_of.get(frozenset((chain[0], chain[1])))
        if a is None:
            continue
        for s, others in by.items():
            for o in others:
                b = idx_of.get(frozenset((o[0], o[1])))
                if b is not None and a != b:
                    shell_of.setdefault((min(a, b), max(a, b)), s)

    if done:
        out = relax(xyz_now)
        final = measure_pairs(out, seq, sorted(done), work) if out else {}
        base_d = measure_pairs(relaxed, seq, sorted(done), work)
        got_mix, blind = {}, 0
        for q in sorted(done):
            s = shell_of.get(q)
            v, b0 = final.get(q), base_d.get(q)
            if s is None or v is None or b0 is None:
                blind += 1
                continue
            got_mix[s] = got_mix.get(s, 0) + max(0, v - b0)
        tot = sum(got_mix.values())
        print(f"\n  shell mix over {len(done)} designed pairs"
              + (f", {blind} unmeasured" if blind else ""))
        print(f"\n  {'shell':>6} {'asked':>8} {'delivered':>10} "
              f"{'yield':>8} {'ask next':>9}")

        # What to ask for next time, and why asking for the delivered mix
        # would not work.
        #
        # A pair in an outer shell delivers about twice what a first-shell pair
        # does, so a mix of pairs is not a mix of entanglements. That makes the
        # map from asked to delivered not a fixed point: feeding the delivered
        # mix back in drifts further out rather than settling. Measured,
        # asking 0.11/0.57/0.32 returns 0.05/0.58/0.37.
        #
        # Dividing the target by the measured yield is the correction, and it
        # is what `select_by_shells(yield_by_shell=...)` consumes. A shell with
        # zero yield gets no request, because no weighting reaches a shell that
        # delivers nothing.
        yields, next_ask = {}, {}
        for s in sorted(set(mix) | set(got_mix)):
            want_f = mix.get(s, 0.0) / sum(mix.values())
            have_f = got_mix.get(s, 0) / tot if tot else 0.0
            y = have_f / want_f if want_f > 0 else 0.0
            yields[s] = y
            next_ask[s] = want_f / y if y > 0 else 0.0
        nrm = sum(next_ask.values())
        for s in sorted(set(mix) | set(got_mix)):
            want_f = mix.get(s, 0.0) / sum(mix.values())
            have_f = got_mix.get(s, 0) / tot if tot else 0.0
            nxt = next_ask[s] / nrm if nrm else 0.0
            print(f"  {s:>6} {want_f:>8.2f} {have_f:>10.2f} "
                  f"{yields[s]:>8.2f} {nxt:>9.2f}"
                  + ("   unreachable" if yields[s] == 0 and want_f > 0
                     else ""))
        if nrm:
            print(f"\n  to deliver the mix asked for, run again with "
                  f"--shells "
                  + ",".join(f"{s}:{next_ask[s] / nrm:.2f}"
                             for s in sorted(next_ask)
                             if next_ask[s] > 0))

        # Per-pair dump for rendering: which chains, which shell, what
        # survived, and where the pair actually touches in the final
        # coordinates. A renderer then needs no topon imports and no
        # re-measurement -- the site is the closest approach of the two
        # chains in the coordinates designed.data carries.
        pairs_out = []
        for q in sorted(done):
            v, b0 = final.get(q), base_d.get(q)
            d = (max(0, v - b0) if v is not None and b0 is not None
                 else None)
            pa = np.array([xyz_now[i] for i in seq[q[0]]])
            pb = np.array([xyz_now[i] for i in seq[q[1]]])
            diff = pa[:, None, :] - pb[None, :, :]
            diff -= box * np.round(diff / box)
            d2 = np.einsum("ijk,ijk->ij", diff, diff)
            i0, j0 = np.unravel_index(int(np.argmin(d2)), d2.shape)
            site = pa[i0] - 0.5 * diff[i0, j0]
            site -= box * np.floor(site / box)
            pairs_out.append({
                "a": int(q[0]), "b": int(q[1]),
                "shell": int(shell_of.get(q, 0)),
                "delivered": (None if d is None else int(d)),
                "site": [round(float(x), 3) for x in site],
                "atoms_a": [int(i) for i in seq[q[0]]],
                "atoms_b": [int(i) for i in seq[q[1]]]})
        pj = root / "04_Simulation" / "designed_pairs.json"
        pj.write_text(json.dumps({
            "tag": args.tag, "want": args.want,
            "asked": {str(s): mix[s] / sum(mix.values()) for s in mix},
            "delivered": {str(s): (got_mix.get(s, 0) / tot if tot else 0.0)
                          for s in sorted(set(mix) | set(got_mix))},
            "box": [float(x) for x in box],
            "pairs": pairs_out}))
        print(f"  pair dump: {pj}")

        # Where the delivered entanglements ended up, against where the bias
        # asked for them. Selection placing the *pairs* right is checked in
        # unit tests; this is the built system, so it is the claim itself.
        # Binned on the bias's own geometry, because that is what was
        # promised: a gradient promises a rising profile, a region promises
        # an in-against-out density ratio.
        if args.bias != "uniform":
            def centre_of(q):
                mids = []
                for ci in q:
                    u, v = edges[ci]
                    pu = np.asarray(graph.nodes[u]["pos"], float)
                    pv = np.asarray(graph.nodes[v]["pos"], float)
                    d = pv - pu
                    d -= dims * np.round(d / dims)
                    mids.append(pu + 0.5 * d)
                d = mids[1] - mids[0]
                d -= dims * np.round(d / dims)
                c = mids[0] + 0.5 * d
                return c - dims * np.floor(c / dims)

            spots, vals = [], []
            for q in sorted(done):
                v, b0 = final.get(q), base_d.get(q)
                if v is None or b0 is None:
                    continue
                spots.append(centre_of(q))
                vals.append(max(0, v - b0))
            spots, vals = np.asarray(spots), np.asarray(vals, float)
            if not len(vals) or vals.sum() <= 0:
                print(f"\n  spatial check ({args.bias}): no delivered "
                      f"designed entanglements to place, nothing to check")
            if len(vals) and vals.sum() > 0:
                print(f"\n  spatial check ({args.bias}), over the delivered "
                      f"designed entanglements:")
                if args.bias == "gradient":
                    axn = bias_params.get("axis", "x")
                    ax = {"x": 0, "y": 1, "z": 2}[axn]
                    frac = spots[:, ax] / dims[ax]
                    cut = np.linspace(0.0, 1.0, 5)
                    for i in range(4):
                        m = (frac >= cut[i]) & (frac < cut[i + 1] + (i == 3))
                        print(f"    {axn} in {cut[i]:.2f}..{cut[i + 1]:.2f}: "
                              f"{vals[m].sum() / vals.sum():.2f} of the "
                              f"total ({int(m.sum())} pairs)")
                elif args.bias in ("region", "anti_region"):
                    c = (np.asarray(bias_params.get("center", [0.5] * 3),
                                    float) * dims)
                    r = (float(bias_params.get("radius", 0.3))
                         * float(dims.min()))
                    d = spots - c
                    d -= dims * np.round(d / dims)
                    inside = np.linalg.norm(d, axis=1) <= r
                    vol = (4.0 / 3.0) * np.pi * r ** 3 / float(np.prod(dims))
                    fin = vals[inside].sum() / vals.sum()
                    # fin == 0 is total depletion (the anti_region success
                    # case), not infinite concentration; only fin == 1 is.
                    if fin <= 0:
                        ratio = 0.0
                    elif fin >= 1:
                        ratio = float("inf")
                    else:
                        ratio = (fin / vol) / ((1 - fin) / (1 - vol))
                    print(f"    sphere holds {vol:.0%} of the box and "
                          f"{fin:.0%} of the entanglements: density ratio "
                          f"{ratio:.1f} in against out "
                          f"(asked strength "
                          f"{bias_params.get('strength', 10.0):g})")
                else:                                          # clusters
                    ccs = [np.asarray(x, float) * dims for x in
                           bias_params.get("centers",
                                           [[0.25] * 3, [0.75] * 3])]
                    sg = (float(bias_params.get("sigma", 0.15))
                          * float(dims.min()))
                    dmin = np.min([np.linalg.norm(
                        (spots - cc) - dims * np.round((spots - cc) / dims),
                        axis=1) for cc in ccs], axis=0)
                    near = dmin <= 2.0 * sg
                    # The null expectation, or the fraction reads as an
                    # absolute when it is only meaningful against the volume
                    # the balls cover. Monte-Carlo because the balls overlap
                    # the cell and each other in general.
                    mc = np.random.default_rng(0).uniform(0, 1, (20000, 3))
                    mc = mc * dims
                    mmin = np.min([np.linalg.norm(
                        (mc - cc) - dims * np.round((mc - cc) / dims),
                        axis=1) for cc in ccs], axis=0)
                    cover = float((mmin <= 2.0 * sg).mean())
                    print(f"    within 2 sigma of a cluster centre: "
                          f"{vals[near].sum() / vals.sum():.0%} of the "
                          f"entanglements ({int(near.sum())} pairs), in "
                          f"{cover:.0%} of the box")

    Path(root / "04_Simulation" / "designed.data").write_text(
        Path(rewrite_coords(relaxed,
                            root / "04_Simulation" / "designed.data",
                            xyz_now)).read_text())
    print(f"  system: {root / '04_Simulation' / 'designed.data'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
