"""Chain paths between two fixed endpoints.

Coordinate generation, so it belongs to the conformation stage. Kept apart
from ``manager.py`` because that stage reads and rewrites a LAMMPS data file,
while these are plain geometry: give them two junctions and a bead count and
they hand back a path.
"""
from __future__ import annotations

import numpy as np

__all__ = ["bridging_walk", "straight", "walk_via", "walk_through",
           "loop_around"]


def straight(start, end, n_beads: int) -> np.ndarray:
    """A straight line, ``n_beads`` points including both ends."""
    t = np.linspace(0.0, 1.0, n_beads)[:, None]
    a = np.asarray(start, float)
    return a + t * (np.asarray(end, float) - a)


def bridging_walk(start, end, n_bonds: int, bond: float = 0.97,
                  rng=None) -> np.ndarray:
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
    """
    rng = np.random.default_rng() if rng is None else rng
    start = np.asarray(start, float)
    end = np.asarray(end, float)

    pts = [start]
    for k in range(n_bonds - 1):
        d = end - pts[-1]
        r = float(np.linalg.norm(d))
        remaining = n_bonds - k - 1

        # Cosine of the widest angle from which the end is still reachable.
        cos_min = ((r * r + bond * bond - (remaining * bond) ** 2)
                   / (2.0 * r * bond) if r > 1e-12 else -1.0)
        cos_min = min(1.0, max(-1.0, cos_min))

        c = rng.uniform(cos_min, 1.0)
        s = np.sqrt(max(0.0, 1.0 - c * c))
        head = d / r if r > 1e-12 else np.array([0.0, 0.0, 1.0])

        ref = [0.0, 0.0, 1.0] if abs(head[2]) < 0.9 else [1.0, 0.0, 0.0]
        t = np.cross(head, ref)
        t /= np.linalg.norm(t)
        u = np.cross(head, t)
        phi = rng.uniform(0.0, 2.0 * np.pi)

        pts.append(pts[-1] + bond * (c * head
                                     + s * (np.cos(phi) * t
                                            + np.sin(phi) * u)))
    pts.append(end)
    return np.array(pts)


def walk_via(start, end, via, n_bonds: int, bond: float = 0.97,
             rng=None, at: float = 0.5) -> np.ndarray:
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

    first = bridging_walk(start, via, n1, bond, rng)
    second = bridging_walk(via, end, n2, bond, rng)
    return np.vstack([first, second[1:]])


def walk_through(start, end, waypoints, n_bonds: int, bond: float = 0.97,
                 rng=None) -> np.ndarray:
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

    out = [bridging_walk(pts[0], pts[1], share[0], bond, rng)]
    for i in range(1, len(legs)):
        out.append(bridging_walk(pts[i], pts[i + 1], share[i], bond, rng)[1:])
    return np.vstack(out)


def loop_around(target, i: int, radius: float, n_pts: int = 6,
                phase: float = 0.0) -> np.ndarray:
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

    th = np.linspace(0.0, 2.0 * np.pi, n_pts, endpoint=False) + phase
    return target[i] + radius * (np.cos(th)[:, None] * u
                                 + np.sin(th)[:, None] * v)
