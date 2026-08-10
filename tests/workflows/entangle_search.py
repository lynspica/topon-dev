"""Search for a chain path that delivers a named entanglement topology.

    python tests/workflows/entangle_search.py                    # one pair
    python tests/workflows/entangle_search.py --rounds 3
    python tests/workflows/entangle_search.py --ambitious        # a whole plan

Every geometric attempt in this work to predict whether a path will entangle
a particular partner has failed, and always the same way: proximity is not
entanglement, and no distance measure separates "passed alongside" from
"threaded through". So this does not predict. It proposes paths, measures
each one with a primitive-path analysis, and keeps what actually worked.

The loop:

    1. draw a base conformation
    2. propose N paths for the chain being routed, each encircling its
       target strand at a different radius, phase and bead
    3. write all N as separate configurations and measure them in one batch
    4. score: did the named partner appear the requested number of times,
       and how many partners appeared that nobody asked for
    5. keep the best, then propose again around it

Batching step 3 is what makes this affordable. Z1+ takes a few seconds on a
hundred-chain system, so a few dozen candidates per round is a minute or two,
and the loop is measuring the objective rather than a proxy for it.
"""
from __future__ import annotations

import argparse
import collections
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.conformation.paths import (  # noqa: E402
    bridging_walk,
    loop_around,
    walk_through,
)
from tests.workflows.entangle_all import CASES  # noqa: E402
from tests.workflows.entangle_steps import (  # noqa: E402
    BOND,
    DP,
    LATTICE,
    OUT,
    build_network,
    geometry,
    write_z1,
)

RUNNER = Path(__file__).resolve().parent / "run_z1_batch.sh"


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

def _wsl(p):
    p = Path(p).resolve()
    return f"/mnt/{p.drive.rstrip(':').lower()}{str(p)[len(p.drive):]}".replace(
        "\\", "/")


def measure_batch(directory):
    """Z1+ on every .Z1 in ``directory``. Returns {name: {chain: Counter}}.

    One subprocess for the whole batch. Z1+ writes its output into the
    working directory, so each configuration needs its own, which the runner
    handles.
    """
    try:
        subprocess.run(["wsl.exe", "-d", "Ubuntu", "--", "bash",
                        _wsl(RUNNER), _wsl(directory)],
                       capture_output=True, text=True, timeout=3600)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None

    out = {}
    for f in sorted(Path(directory).glob("SP_*.dat")):
        out[f.stem[3:]] = _partners(f)
    return out


def _partners(path):
    """Parse Z1+'s shortest-path file into {chain: Counter(partner: count)}.

    Rows carry ``x y z s E`` and, with the ``+`` option, ``chain node`` as
    well, so an entanglement row is seven fields with E of 1. Chain and
    partner numbering is 1-based.
    """
    rows = [ln.split() for ln in Path(path).read_text().splitlines()
            if ln.strip()]
    i = 0
    n_chains = int(float(rows[i][0]))
    i += 2                                   # skip the box line
    out = collections.defaultdict(collections.Counter)
    for c in range(1, n_chains + 1):
        n = int(float(rows[i][0]))
        i += 1
        for _ in range(n):
            r = rows[i]
            i += 1
            if len(r) >= 7 and int(float(r[4])) == 1:
                p = int(float(r[5]))
                if p > 0:
                    out[c][p] += 1
    return out


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

@dataclass
class Wish:
    """What the search is trying to achieve for one routed chain."""

    chain: int                  # 1-based index in the Z1+ file
    want: dict                  # {partner index: how many entanglements}
    penalty: float = 1.0        # cost of one entanglement nobody asked for
    baseline: int = 0           # entanglements the chain had before routing


def cost(partners, wish):
    """Lower is better: distance from the requested count, plus whatever the
    routing *added* elsewhere.

    Collateral is counted against the baseline, not against zero. A chain in
    a melt already carries several entanglements before anything is routed --
    three, on the chain used here -- so charging for all of them makes the
    cheapest path the one that entangles nothing, and the search converges on
    doing nothing. Measured before this was fixed: a candidate reading 0 of 1
    on the target beat one reading 1 of 1, because the second had picked up
    three more elsewhere than the first.
    """
    got = partners.get(wish.chain, collections.Counter())
    miss = sum(abs(wish.want.get(p, 0) - _both(partners, wish.chain, p))
               for p in wish.want)
    others = sum(v for p, v in got.items() if p not in wish.want)
    added = max(0, others - wish.baseline)

    # Lexicographic: hit the target first, then minimise collateral. A single
    # weighted sum cannot express "deliver what I asked for", because any
    # weight small enough to tolerate collateral also makes missing the
    # target outright the cheapest move -- measured twice while tuning this,
    # the search converged on entangling nothing, which scores well and is
    # useless.
    return (miss, wish.penalty * added), miss, added


def _both(partners, a, b):
    """Entanglements between a and b, counted from either side.

    Z1+ attributes a kink to the chain responsible for it, so a link between
    two chains is reported on one side or the other and not reliably on both.
    """
    return partners.get(a, collections.Counter())[b] + \
        partners.get(b, collections.Counter())[a]


# ---------------------------------------------------------------------------
# Proposals
# ---------------------------------------------------------------------------

def propose(a0, a1, target_path, box, n, rng, around=None, spread=1.0,
            dp=DP, bond=BOND):
    """``n`` candidate paths, each encircling the target strand somewhere.

    ``around`` is a previous winner's ``(bead, radius, points, phase)``; when
    given, proposals cluster near it instead of covering the whole range.
    That is the difference between the first round and the ones after it.
    """
    out = []
    for _ in range(n):
        if around is None:
            bead = int(rng.integers(4, len(target_path) - 4))
            radius = float(rng.uniform(1.2, 3.2))
            pts = int(rng.choice([4, 5, 6, 8]))
            phase = float(rng.uniform(0.0, 2.0 * np.pi))
        else:
            b0, r0, p0, ph0 = around
            bead = int(np.clip(b0 + rng.integers(-6, 7) * spread,
                               4, len(target_path) - 5))
            radius = float(np.clip(r0 + rng.normal(0.0, 0.35 * spread),
                                   1.0, 3.5))
            pts = int(np.clip(p0 + rng.integers(-1, 2), 4, 8))
            phase = ph0 + float(rng.normal(0.0, 0.8 * spread))

        ring = loop_around(target_path, bead, radius, pts, phase)
        try:
            path = walk_through(a0, a1, list(ring), dp + 1, bond, rng)
        except ValueError:
            continue
        out.append(((bead, radius, pts, phase), path))
    return out


# ---------------------------------------------------------------------------
# The loop
# ---------------------------------------------------------------------------

def search(base, keys, routed, wish, geo, rounds, per_round, rng, work):
    """Route ``routed`` to satisfy ``wish``, keeping the best each round."""
    L = geo["L"]
    a0, a1 = geo["chords"][routed]
    target_key = keys[list(wish.want)[0] - 1]

    best = None
    around, spread = None, 1.0
    for rnd in range(rounds):
        # The target in the image nearest the chain being routed.
        tgt = base[target_key]
        tgt = tgt + L * np.round((base[routed].mean(0) - tgt.mean(0)) / L)
        cands = propose(a0, a1, tgt, L, per_round, rng, around, spread)
        if not cands:
            print(f"  round {rnd + 1}: no candidate path could be built")
            continue

        # Clear the files rather than the directory: OneDrive holds a lock
        # on the directory itself and rmtree fails with a permission error.
        work.mkdir(parents=True, exist_ok=True)
        for old in list(work.glob("*.Z1")) + list(work.glob("SP_*.dat")):
            old.unlink()
        for i, (_, path) in enumerate(cands):
            trial = dict(base)
            trial[routed] = path
            arr = [trial[k] for k in keys]
            ref = arr[0].mean(0)
            write_z1(work / f"c{i:03d}.Z1",
                     [p + L * np.round((ref - p.mean(0)) / L) for p in arr], L)

        res = measure_batch(work)
        if not res:
            print("  Z1+ unavailable")
            return best

        want_of = list(wish.want)[0]
        scored = []
        for i, (knobs, path) in enumerate(cands):
            p = res.get(f"c{i:03d}")
            if p is None:
                continue
            c, miss, unw = cost(p, wish)
            scored.append((c, miss, unw, _both(p, wish.chain, want_of),
                           knobs, path))
        if not scored:
            continue
        scored.sort(key=lambda t: t[0])
        c, miss, unw, hits, knobs, path = scored[0]
        print(f"  round {rnd + 1}: best ({c[0]:.0f}, {c[1]:.0f}) "
              f"(target {hits}/{wish.want[want_of]}, {unw} added) "
              f"from {len(scored)} candidates")

        if best is None or c < best[0]:
            best = (c, miss, unw, hits, knobs, path)
            around, spread = knobs, max(0.35, spread * 0.6)
        else:
            spread = min(1.5, spread * 1.4)      # widen when a round fails
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--per-round", type=int, default=16)
    ap.add_argument("--want", type=int, default=1,
                    help="entanglements wanted with the named partner")
    ap.add_argument("--penalty", type=float, default=1.0)
    ap.add_argument("--rank", type=int, default=None,
                    help="which partner by distance; default is a third of "
                         "the way down the list")
    ap.add_argument("--dp", type=int, default=DP)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    spec = dict(LATTICE)
    spec.update(CASES["SC"])
    spec["dims"] = (4, 4, 4)
    graph = build_network(spec)
    geo = geometry(graph, dp=args.dp, density=0.85)
    ch, ends, L = geo["chords"], geo["ends"], geo["L"]
    keys = sorted(ch)

    rng = np.random.default_rng(args.seed)
    base = {k: bridging_walk(c0, c1, args.dp + 1, BOND, rng)
            for k, (c0, c1) in ch.items()}

    A = keys[0]
    mid_a = 0.5 * (ch[A][0] + ch[A][1])
    far = sorted(
        (float(np.linalg.norm(
            (lambda v: v - L * np.round(v / L))(
                0.5 * (ch[b][0] + ch[b][1]) - mid_a))), b)
        for b in keys[1:] if not set(ends[A]) & set(ends[b]))
    gap, B = far[args.rank if args.rank is not None else len(far) // 3]
    chord = float(np.linalg.norm(ch[A][1] - ch[A][0]))

    # What chain A carries before anything is routed, so collateral is
    # measured as what the routing added.
    probe = OUT / "search_probe"
    probe.mkdir(parents=True, exist_ok=True)
    for old in list(probe.glob("*.Z1")) + list(probe.glob("SP_*.dat")):
        old.unlink()
    arr = [base[k] for k in keys]
    ref = arr[0].mean(0)
    write_z1(probe / "base.Z1",
             [p + L * np.round((ref - p.mean(0)) / L) for p in arr], L)
    res0 = measure_batch(probe) or {}
    p0 = res0.get("base", {})
    a_idx, b_idx = keys.index(A) + 1, keys.index(B) + 1
    baseline = sum(v for q, v in p0.get(a_idx, {}).items() if q != b_idx)

    wish = Wish(chain=a_idx, want={b_idx: args.want},
                penalty=args.penalty, baseline=baseline)
    print(f"  routing chain {A} to entangle chain {B} exactly {args.want} "
          f"time(s)")
    print(f"  they are {gap / chord:.2f} chord-lengths apart, "
          f"{len(keys)} chains total")
    print(f"  {args.rounds} rounds x {args.per_round} candidates\n")

    work = OUT / "search_work"
    best = search(base, keys, A, wish, geo, args.rounds, args.per_round,
                  rng, work)
    if best is None:
        print("\n  nothing built")
        return 1
    c, miss, unw, hits, knobs, path = best
    bead, radius, pts, phase = knobs
    bonds = np.linalg.norm(np.diff(path, axis=0), axis=1)
    verdict = "exactly as asked" if miss == 0 else f"{miss} off target"
    print(f"\n  best: {hits} entanglement(s) with chain {B} "
          f"(wanted {args.want}) -- {verdict}, {unw} added elsewhere")
    print(f"  loop at bead {bead} of the target, radius {radius:.2f}, "
          f"{pts} points")
    print(f"  bonds {bonds.min():.3f} to {bonds.max():.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
