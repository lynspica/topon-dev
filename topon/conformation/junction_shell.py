"""Spread the chains leaving a junction so their first beads do not overlap.

A junction of functionality ``f`` has ``f`` chains leaving one point. Built
by walking each chain along its own chord, every one of those chains puts
its first bead about one bond length from the junction, in whatever
direction its chord happens to point. Lattice chords are not evenly spread
-- they cluster along a few crystallographic directions -- so the first
beads land on top of each other.

Measured on a 4x4x4 MIX network with functionality up to 12: 14861
non-bonded pairs closer than 1.0 sigma, the closest at 0.060, and 64% of
them within 1 sigma of a junction.

That matters beyond a bad starting energy. The soft potential in the first
minimisation stage exists to resolve exactly such overlaps, and the way it
resolves them is by letting the two chains pass through one another. For a
network with prescribed entanglements that is the one move that undoes the
work: a braid measured at linking +0.77 as built came back from stage 1 at
-1.00 with the two strands 0.26 sigma apart.

The fix is geometric. Seat the first bead of each chain on a shell of
points spread evenly over the sphere, sized so neighbouring points are at
least ``spacing`` apart, and blend back onto the chain's own path over the
next few beads. The chain still goes where it was going; it just leaves in
a direction chosen with its siblings in mind.

The radius has to grow with functionality. Points spread on a sphere of
radius ``r`` sit roughly ``r * sqrt(8 pi / (sqrt(3) f))`` apart, so seating
12 chains at 1 sigma needs about 1.0 sigma of radius where 4 chains need
0.6. Holding the radius fixed and raising ``f`` just puts the overlap back.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "spread_points",
    "shell_radius",
    "junction_directions",
    "apply_junction_shells",
]


def spread_points(n: int, iters: int = 400, step: float = 0.05,
                  seed: int = 0) -> np.ndarray:
    """``n`` unit vectors spread as evenly as possible over the sphere.

    A Fibonacci lattice for the starting arrangement, then Coulomb-style
    repulsion. The relaxation matters most at the small counts this is used
    for, where the Fibonacci spiral is visibly uneven.
    """
    if n <= 0:
        return np.zeros((0, 3))
    if n == 1:
        return np.array([[0.0, 0.0, 1.0]])
    if n == 2:
        return np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]])

    i = np.arange(n) + 0.5
    phi = np.arccos(np.clip(1.0 - 2.0 * i / n, -1.0, 1.0))
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    p = np.column_stack([np.cos(theta) * np.sin(phi),
                         np.sin(theta) * np.sin(phi),
                         np.cos(phi)])

    for _ in range(iters):
        d = p[:, None, :] - p[None, :, :]
        r2 = (d * d).sum(-1) + np.eye(n)
        force = (d / r2[..., None] ** 1.5).sum(1)
        p = p + step * force
        p /= np.linalg.norm(p, axis=1, keepdims=True)
    return p


def _min_chord(p: np.ndarray) -> float:
    if len(p) < 2:
        return 2.0
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=-1)
    d[np.diag_indices_from(d)] = np.inf
    return float(d.min())


def shell_radius(n: int, spacing: float = 1.0,
                 points: np.ndarray | None = None) -> float:
    """Radius that seats ``n`` chains at least ``spacing`` apart.

    Solved from the actual spread rather than the asymptotic formula, since
    the counts here are small enough that the two disagree.
    """
    if n <= 1:
        return 0.0
    p = spread_points(n) if points is None else points
    return spacing / _min_chord(p)


def junction_directions(chord_dirs: np.ndarray) -> np.ndarray:
    """Match each chain to a shell point, keeping chains near their chord.

    A chain forced to leave in an unrelated direction has to turn back, and
    the turn costs bond length near the junction. Assignment is by total
    angular cost, so the arrangement stays as close to the chords as the
    even spread allows.
    """
    n = len(chord_dirs)
    shell = spread_points(n)
    if n < 2:
        return shell

    cost = -(chord_dirs @ shell.T)              # maximise alignment
    try:
        from scipy.optimize import linear_sum_assignment
        _, col = linear_sum_assignment(cost)
        return shell[col]
    except ImportError:
        # Greedy fallback: same intent, occasionally a worse pairing.
        out = np.empty_like(chord_dirs)
        free = list(range(n))
        for i in np.argsort(cost.min(axis=1)):
            j = min(free, key=lambda c: cost[i, c])
            out[i] = shell[j]
            free.remove(j)
        return out


def apply_junction_shells(paths, ends, spacing: float = 1.0,
                          blend: int = 4, max_radius: float | None = None,
                          carry: bool = True):
    """Seat the beads next to each junction on a spread shell.

    ``paths`` maps chain id to its bead array, first and last bead sitting
    on the two junctions. ``ends`` maps chain id to ``(u, v)``.

    With ``carry`` false the shell offset decays to nothing by bead
    ``blend`` and the chain rejoins its own path. That is enough when the
    chains leaving a junction diverge, and not enough when they do not:
    measured on a 4x4x4 MIX network, the smallest angle between two chains
    at one junction is 0.0 degrees, because the mix puts two lattice sites
    on the same ray from a third. Chains along collinear chords re-converge
    as soon as the blend releases them, and 3826 of the remaining overlaps
    were between chains sharing a junction.

    With ``carry`` true, the default, the offset is instead carried the
    whole length of the chain: it rises over one bond, holds while
    interpolating between the two ends' seats, and falls over one bond at
    the far junction. Collinear siblings then stay a shell diameter apart
    for their whole length rather than only at the ends.

    Junctions themselves never move. They carry the network's topology and
    the box, and a junction that drifts changes the network rather than
    its conformation.

    Returns a new dict; the input is not modified.
    """
    out = {k: np.array(p, dtype=float, copy=True) for k, p in paths.items()}
    seat_at = {}                       # (chain, side) -> offset from junction

    # Which chains meet at each junction, and at which end.
    at = {}
    for k, (u, v) in ends.items():
        at.setdefault(u, []).append((k, 0))
        at.setdefault(v, []).append((k, -1))

    for node, members in at.items():
        f = len(members)
        if f < 2:
            continue

        # Each chain carries its own copy of this junction. Paths are built
        # unwrapped from their own chord, so a chain that crosses the
        # boundary holds the junction in a different periodic image than its
        # neighbour does. Directions are local to each chain and safe to
        # compare; absolute positions are not, and mixing them up displaces
        # a chain clear across the box.
        dirs = []
        for k, side in members:
            step = 1 if side == 0 else -2
            d = out[k][step] - out[k][side]
            n = np.linalg.norm(d)
            dirs.append(d / n if n > 1e-9 else np.array([0.0, 0.0, 1.0]))
        dirs = np.asarray(dirs)

        pts = spread_points(f)
        r = spacing / _min_chord(pts)
        if max_radius is not None:
            r = min(r, max_radius)
        seats = junction_directions(dirs) * r

        for (k, side), seat in zip(members, seats):
            p = out[k]
            first = 1 if side == 0 else -2
            shift = (p[side] + seat) - p[first]
            if carry:
                seat_at[(k, side)] = shift
                continue
            n_blend = min(blend, len(p) - 2)
            if n_blend < 1:
                continue
            idx = (range(1, n_blend + 1) if side == 0
                   else range(len(p) - 2, len(p) - n_blend - 2, -1))
            for step, i in enumerate(idx):
                p[i] = p[i] + shift * (1.0 - step / float(n_blend))

    if not carry:
        return out

    # Carry each end's offset along the chain. The weight is 0 at the
    # junctions, full one bond in, and flat between, so the whole chain is
    # translated sideways rather than pinched near its ends. The two ends
    # generally get different seats, so the offset interpolates between
    # them across the chain.
    for k, p in out.items():
        n = len(p)
        if n < 4:
            continue
        head = seat_at.get((k, 0))
        tail = seat_at.get((k, -1))
        if head is None and tail is None:
            continue
        if head is None:
            head = np.zeros(3)
        if tail is None:
            tail = np.zeros(3)

        i = np.arange(n)
        rise = np.clip(i, 0, 1).astype(float)
        fall = np.clip(n - 1 - i, 0, 1).astype(float)
        w = (rise * fall)[:, None]
        t = (i / (n - 1.0))[:, None]
        p += w * ((1.0 - t) * head + t * tail)

    return out
