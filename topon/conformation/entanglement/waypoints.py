"""Draw a chain through points you choose, winding it round a partner.

You say where along the chain an entanglement goes and how many turns it
makes. This returns the two paths. It does not refuse: geometry that is
tight or ugly is still built, on the view that a minimisation will tidy it
and that design control matters more than a pretty starting structure.

    sites = [Site(at=0.2, turns=1), Site(at=0.5, turns=2), Site(at=0.8, turns=1)]
    path_a, path_b = entangled_pair(a0, a1, b0, b1, sites)

Both chains are splines through waypoints. At each site the waypoints spiral
about the line midway between the two chords, the partners in antiphase, so
a site with ``turns`` turns contributes exactly that many crossings. Away
from a site the waypoints are just the chord, so the chain runs straight.

The sites are given as fractions of chain A's own length, which is the
coordinate you actually want to place them in. Nothing here works in a
shared frame, so nothing here compresses the range you can ask for.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["Site", "catmull_rom", "site_frame", "winding_waypoints",
           "chain_through", "entangled_pair", "entangled_group",
           "meander_to_length", "resample_path", "MIN_REACH"]

# Below this the site has no radius worth the name: both chains sit on the
# midline together, which is not a winding, though it still reads as a
# plausible linking number because they are on top of each other.
MIN_REACH = 0.12


@dataclass(frozen=True)
class Site:
    """One entanglement, placed by hand.

    at      fraction along chain A, 0 at its first junction, 1 at its second
    turns   how many times the two chains wind about each other here
    radius  how far each chain swings from the midline; None takes a share
            of the local gap, so partners that are far apart reach further
    span    how much of the chain the site occupies, as a fraction; None
            scales it with the turn count
    """

    at: float
    turns: int = 1
    radius: float | None = None
    span: float | None = None


def catmull_rom(points, n_out: int, closed: bool = False) -> np.ndarray:
    """Smooth curve through every one of ``points``.

    Interpolating rather than approximating: a waypoint is a place the chain
    must go, so a spline that merely passes near it is the wrong tool. The
    Catmull-Rom tangent at each point is set by its neighbours, which keeps
    the curve from overshooting between widely spaced waypoints.
    """
    p = np.asarray(points, float)
    if len(p) < 2:
        return np.repeat(p[:1], n_out, axis=0)
    if len(p) == 2:
        t = np.linspace(0.0, 1.0, n_out)[:, None]
        return p[0] + t * (p[1] - p[0])

    # Duplicate the ends so the first and last segments have tangents too.
    q = np.vstack([p[0] + (p[0] - p[1]), p, p[-1] + (p[-1] - p[-2])])
    n_seg = len(p) - 1
    per = max(2, int(np.ceil(n_out / n_seg)) + 1)

    out = []
    for i in range(n_seg):
        p0, p1, p2, p3 = q[i], q[i + 1], q[i + 2], q[i + 3]
        t = np.linspace(0.0, 1.0, per)[:, None]
        t2, t3 = t * t, t * t * t
        seg = 0.5 * ((2 * p1)
                     + (-p0 + p2) * t
                     + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2
                     + (-p0 + 3 * p1 - 3 * p2 + p3) * t3)
        out.append(seg if i == n_seg - 1 else seg[:-1])
    curve = np.vstack(out)

    # Re-space at equal arc length and hit the requested count exactly.
    seg = np.linalg.norm(np.diff(curve, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] < 1e-12:
        return np.repeat(curve[:1], n_out, axis=0)
    want = np.linspace(0.0, s[-1], n_out)
    return np.column_stack([np.interp(want, s, curve[:, d]) for d in range(3)])


def site_frame(a0, a1, b0, b1, at: float, bias: float = 0.5):
    """Local frame where chain A is at fraction ``at`` of its own chord.

    Returns ``(mid, axis, toward, across, gap)``. ``mid`` sits a fraction
    ``bias`` of the way from A to B, so at the default of 0.5 both chains
    travel the same distance to meet.

    Moving the meeting point matters once a chain has more than one partner.
    Its excursion toward one of them sweeps through space, and if it reaches
    halfway every time, that sweep is large enough to catch the others in
    passing: measured on a chain with two partners, both asked for one
    winding, the isolated pairs read two. A small bias keeps the busy chain
    near its own chord and sends the partners to it instead.

    Everything is computed at the position asked for. The alternative --
    finding where the pair comes closest and working outward from there --
    is what made placement unusable: the site's projection back onto the
    chain covered only about a fifth of it, so two sites could not be told
    apart however far apart they were asked to be.
    """
    a0, a1 = np.asarray(a0, float), np.asarray(b0, float) * 0 + np.asarray(a1, float)
    b0, b1 = np.asarray(b0, float), np.asarray(b1, float)

    pa = a0 + float(np.clip(at, 0.0, 1.0)) * (a1 - a0)

    d = b1 - b0
    L2 = float(d @ d)
    u = float(np.clip(((pa - b0) @ d) / L2, 0.0, 1.0)) if L2 > 1e-12 else 0.0
    pb = b0 + u * d
    bias = float(np.clip(bias, 0.05, 0.95))

    toward = pb - pa
    gap = float(np.linalg.norm(toward))
    toward = toward / gap if gap > 1e-9 else np.array([1.0, 0.0, 0.0])

    axis = (a1 - a0)
    axis = axis / (np.linalg.norm(axis) + 1e-12)
    # Orthogonalise so the three directions are a proper frame.
    toward = toward - (toward @ axis) * axis
    n = np.linalg.norm(toward)
    if n < 1e-9:
        trial = np.array([0.0, 0.0, 1.0])
        if abs(trial @ axis) > 0.9:
            trial = np.array([1.0, 0.0, 0.0])
        toward = trial - (trial @ axis) * axis
        n = np.linalg.norm(toward)
    toward /= n
    across = np.cross(axis, toward)

    return pa + bias * (pb - pa), axis, toward, across, gap


def winding_waypoints(a0, a1, b0, b1, site: Site, per_turn: int = 8,
                      reach: float = 0.45, bias: float = 0.5):
    """Points the two chains must pass through to wind at one site.

    Antiphase spiral about the midline: where A is on one side, B is on the
    other. ``turns`` full revolutions therefore give ``turns`` crossings,
    prescribed rather than hoped for.
    """
    mid, axis, toward, across, gap = site_frame(a0, a1, b0, b1, site.at, bias)
    chord = float(np.linalg.norm(np.asarray(a1, float) - np.asarray(a0, float)))

    radius = site.radius if site.radius is not None else reach * gap

    # Span follows the radius, not the chain. A turn of radius r needs about
    # 2r of axial length to be a turn at all; give it less and the spiral is
    # wider than it is long, which is a flat loop that does not wind. That is
    # what a fixed fraction produces on a short chord: at density 0.30 the
    # chords come out near 10 sigma, so a span of 0.12 is 1.2 sigma against a
    # radius of 1.3, and the site read as no entanglement at all.
    #
    # An explicit span still wins, since the caller may be packing several
    # sites onto one chain and know what it is trading away.
    if site.span is not None:
        half = 0.5 * site.span * chord
    else:
        # A turn of radius r needs about 2.2r of axial length to be a turn
        # rather than a flat loop, with a floor in the chain's own length so
        # a site on a very close pair does not vanish. Widening this to cover
        # the sideways travel as well was tried and is worse: the span costs
        # contour of its own, and 39 pairs went from wanting 91.9 sigma to
        # 120.8.
        half = 0.5 * max(2.2 * radius * site.turns, 0.06 * chord)
    half = min(half, 0.45 * chord)

    n = max(3, per_turn * site.turns + 1)
    phase = np.linspace(0.0, 2.0 * np.pi * site.turns, n)
    along = np.linspace(-half, half, n)

    # Phase 0 puts each chain on its own side, A at -toward and B at
    # +toward, since toward points from A to B. Starting them the other way
    # round has both diving for the midline at the same axial position on
    # the way in, and they collide there before the winding begins: measured
    # 0.22 apart on a pair whose steady state held them at 2*radius.
    wa, wb = [], []
    for ph, u in zip(phase, along):
        radial = radius * (np.cos(ph) * toward + np.sin(ph) * across)
        wa.append(mid + u * axis - radial)
        wb.append(mid + u * axis + radial)
    return np.array(wa), np.array(wb), mid, axis, half


def chain_through(start, end, waypoints, n_beads: int) -> np.ndarray:
    """One chain: its two junctions, everything it must pass through between."""
    pts = [np.asarray(start, float)]
    for w in waypoints:
        pts.extend(np.asarray(w, float))
    pts.append(np.asarray(end, float))
    return catmull_rom(np.array(pts), n_beads)


def _frames(path):
    """A unit normal at every point, varying smoothly along the path."""
    d = np.gradient(path, axis=0)
    n = np.linalg.norm(d, axis=1, keepdims=True)
    tan = d / np.where(n < 1e-12, 1.0, n)

    ref = np.array([0.0, 0.0, 1.0])
    if abs(float(tan[0] @ ref)) > 0.9:
        ref = np.array([1.0, 0.0, 0.0])

    # Parallel transport, so the normal does not spin where the path turns.
    out = np.empty_like(path)
    v = ref - float(ref @ tan[0]) * tan[0]
    v /= np.linalg.norm(v) + 1e-12
    out[0] = v
    for i in range(1, len(path)):
        v = v - float(v @ tan[i]) * tan[i]
        m = np.linalg.norm(v)
        v = v / m if m > 1e-9 else out[i - 1]
        out[i] = v
    return tan, out


def _away_from(path, others, tan):
    """At each point, a unit vector pointing away from the nearest neighbour.

    The wave has to go somewhere, and an arbitrary direction is as likely to
    carry the chain into a neighbour as away from it. Waving *away* is what
    keeps the slack from turning into entanglements nobody asked for:
    measured on a pair asked for one winding, a blind wave delivered two to
    three, all of the excess made by the wave rather than by the site.

    Falls back to the transported normal where there is no neighbour near
    enough to matter.
    """
    _, normal = _frames(path)
    if not others:
        return normal

    pts = np.vstack([np.asarray(o, float) for o in others])
    out = np.array(normal, copy=True)
    for i, q in enumerate(path):
        d = pts - q
        j = int(np.argmin((d * d).sum(axis=1)))
        v = q - pts[j]
        v = v - float(v @ tan[i]) * tan[i]      # keep it off the tangent
        n = np.linalg.norm(v)
        if n > 1e-9:
            out[i] = v / n
    # Smooth, or the direction flips point to point and the wave becomes a
    # zigzag that is longer than it looks and does not clear anything.
    k = 9
    ker = np.ones(k) / k
    for d in range(3):
        out[:, d] = np.convolve(out[:, d], ker, mode="same")
    n = np.linalg.norm(out, axis=1, keepdims=True)
    return np.where(n < 1e-9, normal, out / np.where(n < 1e-9, 1.0, n))


def meander_to_length(path, target: float, protect=(), waves: float = 6.0,
                      tol: float = 1e-3, iters: int = 60,
                      avoid=()) -> np.ndarray:
    """Wave the path sideways until it is ``target`` long.

    A chain has a fixed number of beads and wants a fixed spacing, so its
    path has a length it must be. Drawing it straight and hoping the chord
    obliges is what fails on a uniform lattice, where every chord is
    identical and there is no slack anywhere: bonds came out at 4.039 sigma
    against a limit of 1.5, or the whole box had to be shrunk until the
    network was crushed. Neither is necessary. The path is ours to draw, so
    it is drawn at the length the beads need.

    The wave vanishes at both junctions and across every protected stretch,
    so the entanglements keep the geometry they were given and only the free
    run between them takes up the slack.

    ``protect`` is a list of (lo, hi) fractions of the path to leave alone.
    """
    p = np.asarray(path, float)
    have = float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum())
    if target <= have or have <= 0.0:
        return p

    t = np.linspace(0.0, 1.0, len(p))
    # Zero at the ends, zero across anything protected, smooth in between.
    w = np.sin(np.pi * t) ** 2
    for lo, hi in protect:
        mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
        if half <= 0:
            continue
        d = np.clip(np.abs(t - mid) / (1.6 * half), 0.0, 1.0)
        w = w * (0.5 - 0.5 * np.cos(np.pi * d))

    tan, _ = _frames(p)
    direction = _away_from(p, avoid, tan)
    phase = np.sin(2.0 * np.pi * waves * t)[:, None]
    shape = (w[:, None] * phase) * direction

    def length(a):
        q = p + a * shape
        return float(np.linalg.norm(np.diff(q, axis=0), axis=1).sum())

    lo_a, hi_a = 0.0, max(1e-6, 0.05 * have / max(waves, 1.0))
    for _ in range(40):
        if length(hi_a) >= target:
            break
        hi_a *= 1.8
    for _ in range(iters):
        mid = 0.5 * (lo_a + hi_a)
        if length(mid) < target:
            lo_a = mid
        else:
            hi_a = mid
        if abs(length(mid) - target) <= tol * target:
            break
    return p + 0.5 * (lo_a + hi_a) * shape


def entangled_pair(a0, a1, b0, b1, sites, n_beads: int = 200,
                   per_turn: int = 8, reach: float = 0.45,
                   bond: float | None = None):
    """Both chains, wound at every site asked for.

    With ``bond`` given, each path is drawn at exactly the length its beads
    need -- ``(n_beads - 1) * bond`` -- by waving the free stretches between
    the entanglements. Every bond then comes out at ``bond`` whatever the
    chord happens to be, so a lattice whose chords are all identical is no
    harder than one with a spread of them.

    Returns ``(path_a, path_b, info)`` where ``info`` carries each site's
    midpoint, axis and half span, for drawing and for measurement.
    """
    sites = sorted(sites, key=lambda s: s.at)
    wa_all, wb_all, info = [], [], []
    for s in sites:
        wa, wb, mid, axis, half = winding_waypoints(
            a0, a1, b0, b1, s, per_turn, reach)
        wa_all.append(wa)
        wb_all.append(wb)
        info.append(dict(at=s.at, turns=s.turns, mid=mid, axis=axis,
                         half=half))

    # Draw dense first: a wave has to be resolved before the beads are placed
    # on it, or the resampling cuts corners and loses the length it added.
    dense = max(n_beads * 6, 600)
    pa = chain_through(a0, a1, wa_all, dense)

    # Chain B's waypoints are ordered along A, which is the order B meets
    # them only when the two run the same way. Reverse if they do not, or
    # B doubles back between sites.
    same_way = float((np.asarray(a1, float) - np.asarray(a0, float))
                     @ (np.asarray(b1, float) - np.asarray(b0, float))) >= 0.0
    wb_ordered = wb_all if same_way else [w[::-1] for w in wb_all[::-1]]
    pb = chain_through(b0, b1, wb_ordered, dense)

    if bond is None:
        return resample_path(pa, n_beads), resample_path(pb, n_beads), info

    # Solve the site size so the path is the length the beads need, then wave
    # the free stretches to make up whatever is still missing.
    #
    # Both directions have to be handled and only one of them is meandering.
    # A short chord leaves the path shorter than the beads need, and the wave
    # takes up the slack. A long chord leaves it *longer*, since three
    # detours on a nearly-extended chain overshoot, and no amount of waving
    # shortens anything -- measured on an 76 sigma chord with 80 beads at
    # 0.95, bonds came out at 1.20. Shrinking the sites is what fixes that,
    # so reach is solved for rather than given.
    target = (n_beads - 1) * float(bond)
    a0v, a1v = np.asarray(a0, float), np.asarray(a1, float)
    chord_a = float(np.linalg.norm(a1v - a0v))

    def draw(r):
        wa, wb, nfo = [], [], []
        for s in sites:
            xa, xb, mid, axis, half = winding_waypoints(
                a0, a1, b0, b1, s, per_turn, r)
            wa.append(xa)
            wb.append(xb)
            nfo.append(dict(at=s.at, turns=s.turns, mid=mid, axis=axis,
                            half=half))
        ordered = wb if same_way else [w[::-1] for w in wb[::-1]]
        return (chain_through(a0, a1, wa, dense),
                chain_through(b0, b1, ordered, dense), nfo)

    def plen(p):
        return float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum())

    def worst(drawn):
        # Both chains, not just A. A site is placed by its position along A,
        # but B has to reach it too, and where it lands on B is whatever the
        # geometry says. Two sites spread evenly along A were measured
        # landing on the same point of B, so B doubled back and its path ran
        # to 163 sigma against A's 138 -- and checking only A never saw it.
        return max(plen(drawn[0]), plen(drawn[1]))

    if worst((pa, pb)) > target:
        lo_r, hi_r = 0.0, reach
        for _ in range(40):
            mid_r = 0.5 * (lo_r + hi_r)
            if worst(draw(mid_r)) > target:
                hi_r = mid_r
            else:
                lo_r = mid_r
            if hi_r - lo_r < 1e-4:
                break
        reach = 0.5 * (lo_r + hi_r)
        if reach < MIN_REACH:
            # Same guard entangled_group has. Without it this built sites of
            # almost no radius and said nothing: a three-site request came
            # out at reach 0.086, which is two chains lying against each
            # other rather than winding, and still measured as three
            # entanglements because they were on top of one another.
            raise ValueError(
                f"not enough contour: solving for the sites gives a reach of "
                f"{reach:.3f}, below {MIN_REACH}, which is not a winding. "
                f"Give the chains more slack (a larger coil), fewer sites, "
                f"or a closer partner.")
        pa, pb, info = draw(reach)

    protect = []
    for s, d in zip(sites, info):
        half_frac = min(0.45, d["half"] / max(chord_a, 1e-9))
        protect.append((s.at - half_frac, s.at + half_frac))
    pa = meander_to_length(pa, target, protect)
    pb = meander_to_length(pb, target, protect)

    for d in info:
        d["reach"] = reach
    return resample_path(pa, n_beads), resample_path(pb, n_beads), info


def entangled_group(chords, plan, n_beads: int = 200, per_turn: int = 8,
                    reach: float = 0.45, bond: float | None = None,
                    dropped: list | None = None):
    """Many chains at once, any of them carrying several partners.

    ``plan`` is a list of ``(chain_a, chain_b, sites)``. A chain appearing in
    more than one entry collects the waypoints from all of them, ordered
    along its own length, and is drawn as one path through the lot. That is
    the whole reason this exists: a chain entangled with two different
    partners has one path, not two, and building the pairs separately gives
    two paths that disagree about where the chain goes.

    Chains not named in ``plan`` are drawn straight, at the same length.

    Returns ``(paths, info)``; ``info`` lists every site with the pair it
    belongs to and where it landed on each partner.
    """
    dense = max(n_beads * 6, 600)
    target = None if bond is None else (n_beads - 1) * float(bond)
    plan = list(plan)
    if dropped is None:
        dropped = []

    load = {}
    for ka, kb, sites in plan:
        n = len(sites)
        load[ka] = load.get(ka, 0) + n
        load[kb] = load.get(kb, 0) + n

    def lay(r):
        # r is either one number for every pair, or a dict keyed by pair.
        marks = {k: [] for k in chords}
        nfo = []
        for ka, kb, sites in plan:
            a0, a1 = chords[ka]
            b0, b1 = chords[kb]
            b0v = np.asarray(b0, float)
            db = np.asarray(b1, float) - b0v
            Lb2 = float(db @ db)
            # Send the meeting point toward whichever chain has fewer
            # partners, so the busy one stays near its own chord.
            na, nb = load.get(ka, 1), load.get(kb, 1)
            bias = float(np.clip(nb / float(na + nb), 0.15, 0.85))
            rp = r[(ka, kb)] if isinstance(r, dict) else r
            for s in sorted(sites, key=lambda x: x.at):
                wa, wb, mid, axis, half = winding_waypoints(
                    a0, a1, b0, b1, s, per_turn, rp, bias)
                # Where this site falls along B, so B meets its waypoints in
                # the order B travels rather than the order A does.
                at_b = (float((mid - b0v) @ db) / Lb2) if Lb2 > 1e-12 else 0.5
                # A site is placed along A, and where it falls on B is
                # whatever the geometry gives. Off B's end there is no chain
                # to wind around, so B is sent past its own junction to meet
                # a partner that is not there: measured at -0.05, which read
                # as no entanglement. Clamp into B and record the shortfall
                # so the caller can see the site is not where it was asked
                # for.
                at_b_raw = at_b
                at_b = float(np.clip(at_b, 0.02, 0.98))
                marks[ka].append((s.at, wa))
                marks[kb].append((at_b, wb))
                nfo.append(dict(pair=(ka, kb), at_a=s.at, at_b=at_b,
                                at_b_raw=at_b_raw,
                                off_end=abs(at_b_raw - at_b) > 1e-9,
                                turns=s.turns, mid=mid, axis=axis, half=half,
                                reach=rp))

        out = {}
        for k, (c0, c1) in chords.items():
            items = sorted(marks[k], key=lambda x: x[0])
            if items:
                out[k] = chain_through(c0, c1, [p for _, p in items], dense)
            else:
                out[k] = resample_path(
                    np.stack([np.asarray(c0, float), np.asarray(c1, float)]),
                    dense)
        return out, nfo, marks

    def longest(paths):
        return max(float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum())
                   for p in paths.values())

    def plen_of(paths, k):
        return float(np.linalg.norm(np.diff(paths[k], axis=0), axis=1).sum())

    paths, info, marks = lay(reach)

    if target is not None:
        # Solve the site size so no path overshoots. A chain carrying two
        # partners detours twice and overshoots where a chain carrying one
        # would not, and every chain has the same bead count, so the
        # busiest chain sets the size. Without this, a chain with two
        # partners came out at 2.30 to 2.50 sigma bonds against a limit of
        # 1.5.
        if longest(paths) > target:
            # Solve each pair's reach on its own rather than one number for
            # all of them. A plan drawn from several separation bands mixes
            # gaps that differ several-fold, and a single reach is set by the
            # widest pair, so the close ones are shrunk to nothing for a
            # problem they do not have: 39 pairs across three bands wanted
            # 119.3 sigma against the 77 available, on a plan where every
            # chain carried exactly one site.
            keys = [(ka, kb) for ka, kb, _ in plan]
            solved_per = {}
            for key in keys:
                lo_r, hi_r = 0.0, reach
                for _ in range(30):
                    mid_r = 0.5 * (lo_r + hi_r)
                    trial = dict(solved_per)
                    trial[key] = mid_r
                    for other in keys:
                        trial.setdefault(other, reach)
                    drawn = lay(trial)[0]
                    a_len = plen_of(drawn, key[0])
                    b_len = plen_of(drawn, key[1])
                    if max(a_len, b_len) > target:
                        hi_r = mid_r
                    else:
                        lo_r = mid_r
                    if hi_r - lo_r < 1e-4:
                        break
                solved_per[key] = 0.5 * (lo_r + hi_r)
            # Drop what cannot be afforded rather than refusing the plan.
            # At scale a plan mixes bands whose gaps differ several-fold, and
            # one pair too far apart for its chains' contour should cost that
            # pair, not the other thirty-eight.
            dropped[:] = [k for k, v in solved_per.items() if v < MIN_REACH]
            if dropped and len(dropped) < len(keys):
                keep = [(a, b, st) for a, b, st in plan
                        if (a, b) not in set(dropped)]
                plan[:] = keep
                keys = [(a, b) for a, b, _ in plan]
                solved_per = {k: v for k, v in solved_per.items()
                              if k not in set(dropped)}
            solved = min(solved_per.values()) if solved_per else 0.0
            if solved >= MIN_REACH:
                paths, info, marks = lay(solved_per)
            if solved < MIN_REACH:
                # Even a site of no size overshoots, so there is not enough
                # contour for what was asked. Building it anyway produces a
                # site of zero radius, which is not an entanglement at all
                # but still reads as a plausible linking number because the
                # two chains end up on top of each other.
                bare = longest(lay(0.0)[0])
                at_min = longest(lay(MIN_REACH)[0])
                raise ValueError(
                    f"not enough contour: the busiest chain has "
                    f"{target:.1f} sigma, but a site big enough to be a "
                    f"winding (reach {MIN_REACH}) would need {at_min:.1f}. "
                    f"Sites of no size at all would need {bare:.1f}. "
                    f"Give the chains more slack (a larger coil), fewer "
                    f"sites, or partners that are closer.")
            paths, info, marks = lay(solved)

        # Each chain waves away from the ones it is entangled with, so the
        # slack does not become entanglements nobody asked for.
        partners = {k: set() for k in chords}
        for ka, kb, _ in plan:
            partners[ka].add(kb)
            partners[kb].add(ka)
        settled = dict(paths)
        for k, p in paths.items():
            protect = [(at - 0.08, at + 0.08)
                       for at, _ in sorted(marks[k], key=lambda x: x[0])]
            near = [settled[o] for o in partners.get(k, ()) if o in settled]
            settled[k] = meander_to_length(p, target, protect, avoid=near)
        paths = settled

    return {k: resample_path(p, n_beads) for k, p in paths.items()}, info


def resample_path(p, n: int) -> np.ndarray:
    """Place ``n`` points at equal arc length along ``p``."""
    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] < 1e-12:
        return np.repeat(np.asarray(p, float)[:1], n, axis=0)
    want = np.linspace(0.0, s[-1], n)
    return np.column_stack([np.interp(want, s, np.asarray(p, float)[:, d])
                            for d in range(3)])
