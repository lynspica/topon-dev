"""Chain paths between two fixed endpoints.

Coordinate generation, so it belongs to the conformation stage. Kept apart
from ``manager.py`` because that stage reads and rewrites a LAMMPS data file,
while these are plain geometry: give them two junctions and a bead count and
they hand back a path.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

__all__ = ["Clearance", "bridging_walk", "straight", "walk_via",
           "walk_through", "loop_around", "zigzag", "taut_leg",
           "route_through"]


class Clearance:
    """The beads already in the box, so a new path can be drawn around them.

    A path drawn without regard to what is already there lands beads on top of
    beads. Measured on a relaxed melt at density 0.85, routing a single chain
    took the closest pair in the system from 0.502 sigma to 0.195 and put 153
    beads inside 0.5 sigma where there had been none.

    That is not a small mistake to leave to the minimiser. At 0.195 sigma the
    WCA energy is of order 1e5 kT, so the next minimisation does not relax the
    contact, it shoves -- and the shove is large enough to drag chains through
    each other, which rewrites the topology that was just designed. It shows
    up as designed entanglements being lost and undesigned ones appearing, and
    no amount of care in choosing the winding survives it. The cure is to not
    make the overlap in the first place.

    Distances are minimum-image when a box is given.
    """

    def __init__(self, points, box=None, radius: float = 0.9):
        self.radius = float(radius)
        self.box = None if box is None else np.asarray(box, float).reshape(3)
        pts = np.asarray(points, float).reshape(-1, 3)
        self.tree = None if len(pts) == 0 else cKDTree(
            self._wrap(pts), boxsize=self.box)

    def _wrap(self, p):
        p = np.asarray(p, float).reshape(-1, 3)
        if self.box is None:
            return p
        p = p - self.box * np.floor(p / self.box)
        # cKDTree's periodic box is half-open, and it rejects the whole tree
        # over a single point sitting exactly on the upper face.
        return np.clip(p, 0.0, np.nextafter(self.box, 0.0))

    def near(self, pts) -> np.ndarray:
        """Distance from each point to the nearest bead already there."""
        q = self._wrap(pts)
        if self.tree is None:
            return np.full(len(q), np.inf)
        return np.atleast_1d(self.tree.query(q, k=1)[0])

    def worst(self, pts) -> float:
        """The tightest contact a path makes. Larger is better."""
        return float(self.near(pts).min())

    def ok(self, pts) -> bool:
        return self.worst(pts) >= self.radius


def straight(start, end, n_beads: int) -> np.ndarray:
    """A straight line, ``n_beads`` points including both ends."""
    t = np.linspace(0.0, 1.0, n_beads)[:, None]
    a = np.asarray(start, float)
    return a + t * (np.asarray(end, float) - a)


def bridging_walk(start, end, n_bonds: int, bond: float = 0.97,
                  rng=None, avoid: "Clearance | None" = None,
                  tries: int = 64) -> np.ndarray:
    """A random walk of fixed bond length that lands exactly on ``end``.

    Every step is drawn from the cone of directions that still leaves the far
    junction reachable in the bonds remaining, so the walk closes on it
    without any bond being rescaled afterwards. Returns ``n_bonds + 1``
    points.

    This is the shape a chain in a melt actually has, and the difference from
    a straight line is not cosmetic. A straight chain lies on its chord, so
    whether two chains meet is decided by where their crosslinks sit; a
    coiled one wanders, and chains whose crosslinks are far apart routinely
    run alongside each other while nearest neighbours by crosslink may never
    touch. Anything that needs to know which chains are near which -- the
    entanglement candidate ranking, for one -- needs the coiled shape.

    The walk is unbiased only in the limit of plenty of slack. As the chain
    approaches full extension the reachable cone narrows to nothing and the
    path straightens, which is correct: there is only one way to span a chord
    of ``n_bonds * bond``.

    ``avoid`` keeps the walk off beads that are already there. Each step is
    drawn from the reachable cone as before and kept if it clears; ``tries``
    sets how many draws before settling for the roomiest of them. Every draw
    comes from the same cone as before and the accept test is all that is
    added, so a walk with room everywhere has the same distribution as one
    with no obstacles at all -- though not the same sequence, since drawing
    ``tries`` at a time uses the generator differently. The last bead is a
    junction and is placed wherever the junction is, clear or not.

    Measured drawing one chain through a relaxed melt at density 0.85, the
    tightest contact the path makes: 0.081 to 0.227 sigma with no ``avoid``,
    0.822 with it at 16 draws or fewer only intermittently, 0.822 every time
    at 64. That 0.822 is the melt's own nearest-neighbour distance at the
    junction the walk is pinned to, so it is the best a path between those two
    junctions can do, and the default is the smallest number of draws that
    reaches it reliably.
    """
    rng = np.random.default_rng() if rng is None else rng
    start = np.asarray(start, float)
    end = np.asarray(end, float)
    n_try = 1 if avoid is None else max(1, int(tries))

    pts = [start]
    for k in range(n_bonds - 1):
        d = end - pts[-1]
        r = float(np.linalg.norm(d))
        remaining = n_bonds - k - 1

        # Cosine of the widest angle from which the end is still reachable.
        cos_min = ((r * r + bond * bond - (remaining * bond) ** 2)
                   / (2.0 * r * bond) if r > 1e-12 else -1.0)
        cos_min = min(1.0, max(-1.0, cos_min))

        c = rng.uniform(cos_min, 1.0, n_try)
        s = np.sqrt(np.maximum(0.0, 1.0 - c * c))
        head = d / r if r > 1e-12 else np.array([0.0, 0.0, 1.0])

        ref = [0.0, 0.0, 1.0] if abs(head[2]) < 0.9 else [1.0, 0.0, 0.0]
        t = np.cross(head, ref)
        t /= np.linalg.norm(t)
        u = np.cross(head, t)
        phi = rng.uniform(0.0, 2.0 * np.pi, n_try)

        cand = pts[-1] + bond * (c[:, None] * head
                                 + s[:, None] * (np.cos(phi)[:, None] * t
                                                 + np.sin(phi)[:, None] * u))
        if avoid is None:
            pts.append(cand[0])
            continue
        # First acceptable draw, not the roomiest of them: always taking the
        # roomiest would walk the chain down the middle of whatever void it
        # can find and stop looking like a melt chain. Falling back to the
        # roomiest only matters where nothing clears, which is where the cone
        # has closed and there is no choice left to make anyway.
        gap = avoid.near(cand)
        clear = np.flatnonzero(gap >= avoid.radius)
        pts.append(cand[clear[0]] if len(clear) else cand[int(gap.argmax())])
    pts.append(end)
    return np.array(pts)


def walk_via(start, end, via, n_bonds: int, bond: float = 0.97,
             rng=None, at: float = 0.5,
             avoid: "Clearance | None" = None) -> np.ndarray:
    """A bridging walk from ``start`` to ``end`` that passes through ``via``.

    Two bridging walks joined at the waypoint: ``at`` sets what fraction of
    the bonds are spent getting there. The chain still lands exactly on both
    junctions and every bond is still ``bond`` long.

    This is how a chain reaches a partner its crosslinks are nowhere near. A
    chain carries far more contour than its chord needs -- 77 sigma on a 5.4
    sigma chord at melt density -- and that slack is enough to visit a
    neighbour one or two chord-lengths away and come back. Measured with
    blind draws, partners at 1.0 to 1.7 chord-lengths were reached by 2 to 16
    of 50 attempts; routing through a point on the partner reaches them by
    construction.

    Raises when the detour does not fit: the two legs together cannot be
    shorter than the distance they have to cover.
    """
    rng = np.random.default_rng() if rng is None else rng
    start = np.asarray(start, float)
    end = np.asarray(end, float)
    via = np.asarray(via, float)

    n1 = max(1, int(round(at * n_bonds)))
    n2 = n_bonds - n1
    if n2 < 1:
        raise ValueError("no bonds left for the second leg")

    need1 = float(np.linalg.norm(via - start))
    need2 = float(np.linalg.norm(end - via))
    if need1 > n1 * bond or need2 > n2 * bond:
        raise ValueError(
            f"waypoint out of reach: legs need {need1:.1f} and {need2:.1f} "
            f"sigma but carry {n1 * bond:.1f} and {n2 * bond:.1f}. Move the "
            f"waypoint, shift `at`, or give the chain more contour.")

    first = bridging_walk(start, via, n1, bond, rng, avoid)
    second = bridging_walk(via, end, n2, bond, rng, avoid)
    return np.vstack([first, second[1:]])


def walk_through(start, end, waypoints, n_bonds: int, bond: float = 0.97,
                 rng=None, avoid: "Clearance | None" = None) -> np.ndarray:
    """A bridging walk visiting each of ``waypoints`` in order.

    Bonds are shared between the legs in proportion to how far each has to
    travel, so no leg is left short of what it needs. Raises when the whole
    route is longer than the chain.
    """
    rng = np.random.default_rng() if rng is None else rng
    pts = [np.asarray(start, float)]
    pts += [np.asarray(w, float) for w in waypoints]
    pts.append(np.asarray(end, float))

    legs = [float(np.linalg.norm(pts[i + 1] - pts[i]))
            for i in range(len(pts) - 1)]
    total = sum(legs)
    if total > n_bonds * bond:
        raise ValueError(
            f"route is {total:.1f} sigma but the chain carries "
            f"{n_bonds * bond:.1f}")

    # One bond minimum per leg, the rest shared by distance.
    share = [max(1, int(round(n_bonds * L / total))) for L in legs]
    while sum(share) > n_bonds:
        share[int(np.argmax(share))] -= 1
    while sum(share) < n_bonds:
        share[int(np.argmin([s / max(L, 1e-9) for s, L in zip(share, legs)]))] += 1

    out = [bridging_walk(pts[0], pts[1], share[0], bond, rng, avoid)]
    for i in range(1, len(legs)):
        out.append(
            bridging_walk(pts[i], pts[i + 1], share[i], bond, rng, avoid)[1:])
    return np.vstack(out)


def loop_around(target, i: int, radius: float, n_pts: int = 6,
                phase: float = 0.0,
                avoid: "Clearance | None" = None,
                span: float = 1.0) -> np.ndarray:
    """Waypoints that encircle ``target``'s strand at bead ``i``.

    Returns ``n_pts`` points on a circle of ``radius`` about the target's
    local tangent, so a chain routed through them in order passes once around
    that strand.

    Going *around* is the thing. Routing a chain to a point *on* its intended
    partner puts the two side by side and creates no link between them:
    measured, twelve such attempts raised the routed chain's own entanglement
    count from 3 to 8 while leaving the count with the named partner at zero,
    because at melt density the arriving chain is caught by whichever
    neighbour is topologically in the way. Encircling is what cannot be
    undone by pulling the two taut.

    ``span`` is how much of a turn to make, in units of a full circle. One
    full turn is not the smallest thing that links: it puts two crossings into
    the primitive path, not one, which is why asking for a count of one and
    only ever generating full turns comes back with two every time. Values
    under one make a hook rather than a loop, and somewhere above a half turn
    is where it starts to catch.

    With ``avoid``, the ring keeps its winding but not its exact shape: each
    waypoint is placed inside its own slice of the circle, wherever in that
    slice is clear of the beads already there. The waypoints are landed on
    exactly, so a ring laid across occupied sites puts beads on top of beads
    however well the walk between them behaves.
    """
    target = np.asarray(target, float)
    i = int(np.clip(i, 1, len(target) - 2))
    tan = target[i + 1] - target[i - 1]
    n = float(np.linalg.norm(tan))
    tan = tan / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])

    ref = np.array([0.0, 0.0, 1.0])
    if abs(float(tan @ ref)) > 0.9:
        ref = np.array([1.0, 0.0, 0.0])
    u = np.cross(tan, ref)
    u /= np.linalg.norm(u)
    v = np.cross(tan, u)

    full = 2.0 * np.pi * float(span)

    def ring(ph):
        th = np.linspace(0.0, full, n_pts, endpoint=(span != 1.0)) + ph
        return target[i] + radius * (np.cos(th)[:, None] * u
                                     + np.sin(th)[:, None] * v)

    if avoid is None:
        return ring(phase)

    # Place each waypoint on its own, inside its own slice of the circle.
    #
    # Rotating the ring rigidly does not work. At melt density a point lies
    # within 0.9 sigma of 2.6 beads on average, so there is almost no free
    # volume to rotate into, and one shared angle is not enough freedom to
    # clear every point at once: measured, it left the tightest contact at
    # 0.10 sigma, no better than not trying at all. Nor can the walk between
    # the points make up for it, because at radius 1.2 sigma with eight points
    # they sit 0.92 sigma apart, one bond, so there is nothing in between to
    # move.
    #
    # A path still goes once around as long as each waypoint stays in its own
    # sector, which leaves a slice of angle and a range of radius free per
    # point. Candidates are tried nearest-to-nominal first, so the ring keeps
    # the shape it was asked for wherever it can.
    th0 = np.linspace(0.0, full, n_pts, endpoint=(span != 1.0)) + phase
    d_th = np.linspace(-1.0, 1.0, 7) * (0.6 * abs(full) / max(n_pts, 2))
    d_r = radius * np.array([1.0, 1.15, 0.87, 1.35, 0.75, 1.6])
    rr, tt = (x.ravel() for x in np.meshgrid(d_r, d_th, indexing="ij"))

    out = []
    for a in th0:
        cand = target[i] + rr[:, None] * (np.cos(a + tt)[:, None] * u
                                          + np.sin(a + tt)[:, None] * v)
        nominal = target[i] + radius * (np.cos(a) * u + np.sin(a) * v)
        cand = cand[np.argsort(np.linalg.norm(cand - nominal, axis=1))]
        gap = avoid.near(cand)
        clear = np.flatnonzero(gap >= avoid.radius)
        out.append(cand[clear[0]] if len(clear) else cand[int(gap.argmax())])
    return np.array(out)


def _lateral(direction, hint):
    """A unit vector perpendicular to ``direction``, as near ``hint`` as it
    can be."""
    d = np.asarray(direction, float)
    n = float(np.linalg.norm(d))
    d = d / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])
    h = np.asarray(hint, float)
    h = h - d * float(h @ d)
    n = float(np.linalg.norm(h))
    if n < 1e-9:
        ref = np.array([0.0, 0.0, 1.0])
        if abs(float(d @ ref)) > 0.9:
            ref = np.array([1.0, 0.0, 0.0])
        h = np.cross(d, ref)
        n = float(np.linalg.norm(h))
    return h / n


def zigzag(start, end, n_bonds: int, bond: float = 0.97, hint=None):
    """Exact bonds from ``start`` to ``end``, using up every one of them.

    A chain carries far more contour than its route needs. Measured on one
    designed pair: 77 sigma of contour for a route of about 21. A random walk
    disposes of that slack by wandering, and the wandering crosses the target
    again on its own account, so the count stops being a property of the design
    -- the same design, same site, same winding, drawn with three different
    seeds, came back 4, 7 and 0.

    A zigzag spends the same slack in a fixed local shape. It folds against one
    lateral direction and stays near the straight line between its ends, so the
    contour goes somewhere known instead of somewhere random, and the count
    becomes reproducible. Nothing here draws from a generator.

    ``hint`` is the lateral direction to fold against, and only its component
    perpendicular to the leg is used. Point it away from whatever the path is
    meant not to touch.
    """
    start = np.asarray(start, float)
    end = np.asarray(end, float)
    d = end - start
    dist = float(np.linalg.norm(d))
    if n_bonds < 1:
        raise ValueError("a leg needs at least one bond")
    if dist > n_bonds * bond + 1e-9:
        raise ValueError(
            f"leg is {dist:.1f} sigma but carries {n_bonds * bond:.1f}")
    if n_bonds == 1:
        return np.vstack([start, end])

    axis = d / dist if dist > 1e-12 else _lateral([0, 0, 1], [1, 0, 0])
    lat = _lateral(axis, hint if hint is not None else axis + 1.0)

    # Offsets alternate about the straight line and vanish at both ends, so
    # the leg starts and finishes exactly where it was told to.
    sign = np.array([0.0] + [1.0 if i % 2 else -1.0
                             for i in range(1, n_bonds)] + [0.0])

    def reach(h):
        off = sign * h
        step = np.sqrt(np.maximum(bond ** 2 - np.diff(off) ** 2, 0.0))
        return float(step.sum())

    # One bisection on the fold height: taller folds eat more contour, so
    # there is exactly one height that lands the leg on its far end.
    if reach(0.0) < dist:
        h = 0.0
    else:
        lo, hi = 0.0, bond
        while reach(hi) > dist and hi < 50.0 * bond:
            hi *= 2.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if reach(mid) > dist:
                lo = mid
            else:
                hi = mid
        h = 0.5 * (lo + hi)

    off = sign * h
    step = np.sqrt(np.maximum(bond ** 2 - np.diff(off) ** 2, 0.0))
    total = step.sum()
    step = step * (dist / total) if total > 1e-12 else step
    along = np.concatenate([[0.0], np.cumsum(step)])
    return start + along[:, None] * axis + off[:, None] * lat


def route_through(start, end, waypoints, n_bonds: int, bond: float = 0.97,
                  avoid: "Clearance | None" = None):
    """Visit every waypoint in order, spending all the contour deterministically.

    The same job as ``walk_through`` and the same guarantees -- exact bonds,
    lands on both junctions, hits each waypoint -- but the legs come from
    ``taut_leg`` rather than a random walk, so the same design gives the same
    path every time. That is what makes a requested entanglement count
    reproducible: with random legs the leftover contour wanders back across the
    target and adds crossings of its own.

    ``avoid`` keeps the legs off the beads already in the box, and each leg
    also avoids the ones laid down before it.
    """
    pts = [np.asarray(start, float)]
    pts += [np.asarray(w, float) for w in waypoints]
    pts.append(np.asarray(end, float))

    legs = [float(np.linalg.norm(pts[i + 1] - pts[i]))
            for i in range(len(pts) - 1)]
    total = sum(legs)
    if total > n_bonds * bond:
        raise ValueError(
            f"route is {total:.1f} sigma but the chain carries "
            f"{n_bonds * bond:.1f}")

    # Each leg needs enough bonds to span itself; the rest is shared by
    # length. Reserving a spare bond per leg looks harmless and is not: a ring
    # of seventeen waypoints is eighteen legs, so it quietly costs eighteen
    # bonds of contour and turns routes that fit into routes that do not.
    need = [max(1, int(np.ceil(L / bond - 1e-9))) for L in legs]
    share = list(need)
    spare = n_bonds - sum(share)
    if spare < 0:
        raise ValueError(
            f"route needs {sum(need)} bonds but the chain has {n_bonds}")
    for i in np.argsort([-L for L in legs]):
        take = int(round(spare * legs[i] / total)) if total > 1e-12 else 0
        share[i] += take
    share[int(np.argmax(legs))] += n_bonds - sum(share)

    out, placed = [], []
    for i in range(len(legs)):
        leg = taut_leg(pts[i], pts[i + 1], share[i], bond, avoid, placed)
        placed.extend(leg[:-1])
        out.append(leg if not out else leg[1:])
    return np.vstack(out)


def _sphere(n=72):
    """A fixed, evenly spread set of unit directions.

    Deterministic by construction, which is the point: it stands in for the
    random draw so the same design gives the same path every time.
    """
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi),
                            np.cos(phi)])


_DIRS = _sphere()


def taut_leg(start, end, n_bonds: int, bond: float = 0.97,
             avoid: "Clearance | None" = None, placed=None):
    """A deterministic leg of exact bonds that avoids itself and its
    surroundings.

    Same guarantees as ``bridging_walk`` -- every step from the cone that keeps
    the far end reachable, so it lands there exactly with no bond rescaled --
    but the direction is chosen from a fixed set rather than drawn, so the path
    is reproducible.

    It also avoids the beads it has already laid down. A planar fold cannot:
    given far more bonds than the leg needs, its axial step goes to zero and
    every second bead lands exactly on the one two before it. That is not a
    near miss, it is a zero separation, and LAMMPS reports an infinite pair
    energy and stops. Chains carry three or four times the contour their route
    needs, so this is the normal case rather than the corner one.
    """
    start = np.asarray(start, float)
    end = np.asarray(end, float)
    if n_bonds < 1:
        raise ValueError("a leg needs at least one bond")
    dist = float(np.linalg.norm(end - start))
    if dist > n_bonds * bond + 1e-9:
        raise ValueError(
            f"leg is {dist:.1f} sigma but carries {n_bonds * bond:.1f}")

    mine = [] if placed is None else [np.asarray(p, float) for p in placed]
    pts = [start]
    mine.append(start)
    for k in range(n_bonds - 1):
        d = end - pts[-1]
        r = float(np.linalg.norm(d))
        remaining = n_bonds - k - 1
        cos_min = ((r * r + bond * bond - (remaining * bond) ** 2)
                   / (2.0 * r * bond) if r > 1e-12 else -1.0)
        cos_min = min(1.0, max(-1.0, cos_min))
        head = d / r if r > 1e-12 else np.array([0.0, 0.0, 1.0])

        ok = _DIRS[(_DIRS @ head) >= cos_min - 1e-12]
        if not len(ok):
            ok = head[None, :]
        cand = pts[-1] + bond * ok

        gap = (avoid.near(cand) if avoid is not None
               else np.full(len(cand), np.inf))
        if len(mine) > 1:
            # Its own beads, minus the one it is stepping off.
            prev = np.array(mine[:-1])
            own = np.linalg.norm(cand[:, None, :] - prev[None, :, :], axis=2)
            gap = np.minimum(gap, own.min(axis=1))

        # Spend the slack evenly, and prefer a step that clears.
        #
        # The two obvious rules both fail. Roomiest walks the chain down the
        # middle of whatever void it can find, and a chain carrying three
        # times the contour its route needs wanders far enough to cross the
        # target again on its own account: a pair built for 2 came back 3, and
        # 0 after minimisation. Straightest is worse in a quieter way -- it
        # arrives early and has to burn every leftover bond in the small space
        # left, folding on itself to 0.13 sigma, which is the overlap that
        # makes minimisation push chains through each other.
        #
        # Closing the gap at a constant rate spreads the slack along the whole
        # leg instead, which is what an unbiased walk does by accident and
        # this does on purpose.
        want = float(np.clip(r / (bond * max(remaining, 1)), cos_min, 1.0))
        room = np.flatnonzero(gap >= (avoid.radius if avoid is not None
                                      else 0.9))
        pool = room if len(room) else np.arange(len(cand))
        pick = int(pool[np.abs((ok[pool] @ head) - want).argmin()])
        nxt = cand[pick]

        if remaining == 1:
            # The last free bead has to sit exactly one bond from the junction
            # as well as one bond from here, and a fixed set of directions
            # cannot land on that by chance: it closed at 0.841 sigma instead
            # of 0.95. The two conditions describe a circle, so put the bead on
            # the point of it nearest the direction already chosen.
            ax = end - pts[-1]
            dd = float(np.linalg.norm(ax))
            if 1e-9 < dd <= 2.0 * bond:
                ax = ax / dd
                mid = pts[-1] + 0.5 * dd * ax
                rho = np.sqrt(max(bond ** 2 - 0.25 * dd * dd, 0.0))
                v = nxt - mid
                v = v - ax * float(v @ ax)
                n = float(np.linalg.norm(v))
                v = v / n if n > 1e-9 else _lateral(ax, ax + 1.0)
                nxt = mid + rho * v

        pts.append(nxt)
        mine.append(pts[-1])
    pts.append(end)
    return np.array(pts)
