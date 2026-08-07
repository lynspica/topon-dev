"""Waypoint braid: prescribe an entanglement's winding count exactly.

The legacy kink (:func:`topon.utils.network_helpers.calculate_entangled_kink`)
gives each chain a Gaussian bulge peaked at the midpoint of its *own* span and
aimed at where it believes the partner to be. Three consequences follow, all
measured:

* the bulge is pinned to ``t = 0.5``, so a pair whose true closest approach
  sits elsewhere simply does not hook (100% hook rate at zero axial offset,
  18% past 0.25);
* each chain scales its reach by its own length, so unequal partners aim at
  different places and miss;
* the winding count is whatever the shape happens to produce.

This module inverts the construction. A *contact* is defined once from the
pair -- a shared axis midway between the two chords, and a shared frame -- and
both chains follow anti-phase ellipses about that one axis, each taking its
own side. Neither chain ever consults the other's path, only the shared frame.

Two anti-phase helices with ``e`` turns have linking number ``e``, so the
winding count is *prescribed* rather than emergent. Verified by closing each
path with its own chord and taking the Gauss linking number: ``e = 1..5``
measure 1..5 exactly, and a pair with no braid measures 0.

Nothing here knows about lattices, periodic images or graphs. It takes two
chords and returns two paths, which keeps it testable in isolation; the
callers in :mod:`topon.assignment` and :mod:`topon.conformation` supply the
minimum-image geometry.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

__all__ = [
    "BraidShape",
    "Contact",
    "make_contact",
    "plan_braid",
    "braid_path",
    "braid_pair",
    "closest_approach",
    "linking_number",
    "chord_closed_linking",
    "far_closed_linking",
]


# ---------------------------------------------------------------------------
# Shape
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BraidShape:
    """Geometry of one braid. All lengths in the same units as the chords.

    ``n_radius`` (toward the partner) is what sets clearance; ``m_radius``
    (across) is deliberately small. A circular braid swings a full half-gap
    sideways, which on a lattice reaches into neighbouring strand corridors
    and creates entanglements nobody prescribed, while costing several times
    the contour for no gain in clearance -- m-radius was measured to have no
    effect on partner clearance at fixed pitch.

    ``ramp`` is an absolute axial length, not a fraction of the span. As a
    fraction it scales with the turn count, so a 3-turn braid would spend
    three times the axial room on blending as a 1-turn braid -- and axial
    room is the budget that binds on a lattice chord. The blend is a radial
    move of order ``n_radius``, so its cost should track that instead.
    """

    n_radius: float = 0.9
    m_radius: float = 1.5
    pitch: float = 3.2
    ramp: float = 2.0

    def span(self, e: int) -> float:
        """Axial length a braid of ``e`` turns needs, blends included."""
        return e * self.pitch + 2.0 * self.ramp

    def turns_within(self, half_span: float) -> int:
        """Turns that fit in ``±half_span`` of axial room."""
        usable = 2.0 * half_span - 2.0 * self.ramp
        return int(max(0.0, usable) // self.pitch) if self.pitch > 0 else 0

    def fit_to_gap(self, gap: float, reach: float = 0.4) -> "BraidShape":
        """Shrink the braid so neither chain reaches past its partner's chord.

        A chain swings ``±n_radius`` about the shared axis, which sits at the
        midpoint of the gap, so its excursion toward the partner reaches
        ``gap/2 + n_radius`` from its own chord. Once ``n_radius`` exceeds
        ``gap/2`` that excursion passes *beyond* the partner's chord line.

        Two things go wrong there, and the first is the serious one:

        * the partners have effectively swapped sides, so the braid no
          longer separates them and nothing prevents the chains sliding
          through one another during minimisation;
        * the chord-closed loops used to verify the braid genuinely
          intersect, and the linking number of intersecting curves is
          undefined -- measured readings of -5.9 and +9.9 for a braid asked
          for one winding.

        ``reach`` is the fraction of the gap the radius may take. The default
        0.4 keeps a visible margin below the 0.5 where the chords touch.
        """
        if gap <= 0.0:
            return self
        n_max = reach * gap
        if self.n_radius <= n_max:
            return self
        scale = n_max / self.n_radius
        return BraidShape(n_radius=n_max,
                          m_radius=self.m_radius * scale,
                          pitch=self.pitch,
                          ramp=self.ramp)


# ---------------------------------------------------------------------------
# Contact
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Contact:
    """A braid site, defined once from a chord pair.

    Attributes
    ----------
    origin : the point the braid axis passes through, midway between the two
        chords at their closest approach.
    axis, toward, across : orthonormal frame. ``toward`` points from chain A's
        contact point to chain B's, projected perpendicular to the axis.
    s_a, s_b : where the contact falls along each chord, as a fraction.
    gap : distance between the chords at closest approach.
    """

    origin: np.ndarray
    axis: np.ndarray
    toward: np.ndarray
    across: np.ndarray
    s_a: float
    s_b: float
    gap: float


def _unit(v: np.ndarray, fallback: Optional[np.ndarray] = None) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n > 1e-12:
        return v / n
    if fallback is not None:
        return fallback
    return np.array([1.0, 0.0, 0.0])


def _perp_to(v: np.ndarray) -> np.ndarray:
    """Any unit vector perpendicular to ``v``."""
    seed = np.array([0.0, 0.0, 1.0])
    if abs(float(v @ seed)) > 0.9:
        seed = np.array([1.0, 0.0, 0.0])
    return _unit(np.cross(v, seed))


def closest_approach(a0, a1, b0, b1) -> tuple[float, float]:
    """Parameters ``(s, t)`` of closest approach between two segments.

    Both are clamped to ``[0, 1]``, so a contact that would fall beyond a
    chord's end is reported at the end -- which is the honest answer: those
    pairs meet at a junction, where any hook is held by the crosslink rather
    than by topology.
    """
    a0, a1, b0, b1 = (np.asarray(v, float) for v in (a0, a1, b0, b1))
    d1, d2 = a1 - a0, b1 - b0
    r = a0 - b0
    A, B, C = float(d1 @ d1), float(d1 @ d2), float(d2 @ d2)
    D, E = float(d1 @ r), float(d2 @ r)
    denom = A * C - B * B

    if abs(denom) > 1e-12:
        s = (B * E - C * D) / denom
        t = (A * E - B * D) / denom
        return float(np.clip(s, 0.0, 1.0)), float(np.clip(t, 0.0, 1.0))

    # Parallel chords: every point of the overlap is equally close, so the
    # useful answer is the middle of the overlap. Taking s = 0.5 regardless
    # (the obvious shortcut) puts the contact at A's midpoint even when the
    # partner only reaches A's far end, which then sizes the braid against
    # room that is not shared with the partner at all.
    if A <= 1e-12:
        return 0.0, float(np.clip(-E / C, 0.0, 1.0)) if C > 1e-12 else 0.0

    # B's endpoints projected onto A's parameter line.
    t0 = float(-D / A)                       # where b0 lands on A
    t1 = float((float(d1 @ (b1 - a0))) / A)  # where b1 lands on A
    lo, hi = (t0, t1) if t0 <= t1 else (t1, t0)
    lo, hi = max(lo, 0.0), min(hi, 1.0)

    if lo <= hi:                             # genuine overlap
        s = 0.5 * (lo + hi)
    else:                                    # disjoint: nearest ends
        s = 0.0 if hi < 0.0 else 1.0

    t = (float(d2 @ (a0 + s * d1 - b0)) / C) if C > 1e-12 else 0.0
    return float(np.clip(s, 0.0, 1.0)), float(np.clip(t, 0.0, 1.0))


def gap_at(a0, a1, b0, b1, s_a: float) -> tuple[float, float]:
    """Gap between the chords when the contact is pinned at ``s_a`` on A.

    Returns ``(gap, s_b)`` where ``s_b`` is where the partner comes closest
    to that point.
    """
    a0, a1, b0, b1 = (np.asarray(v, float) for v in (a0, a1, b0, b1))
    p = a0 + s_a * (a1 - a0)
    d = b1 - b0
    L2 = float(d @ d)
    s_b = float(np.clip(float((p - b0) @ d) / L2, 0.0, 1.0)) if L2 > 1e-12 else 0.0
    return float(np.linalg.norm(b0 + s_b * d - p)), s_b


def feasible_window(a0, a1, b0, b1, tolerance: float = 0.35,
                    samples: int = 65) -> tuple[float, float]:
    """Range of ``s_a`` where the pair is close enough to braid.

    A braid can only be placed where the chains genuinely approach each
    other, but that is rarely a single point. Two parallel chords are
    equidistant along their whole overlap, so the closest approach is
    degenerate and every pair on a lattice row would otherwise be assigned
    the same midpoint -- which is exactly what stops several partners from
    sharing a chain.

    ``tolerance`` is how much further than the minimum gap is still
    acceptable, as a fraction of that minimum. The window is the contiguous
    run around the closest approach that stays inside it, so a skew pair
    still yields a narrow window and only genuinely close stretches qualify.
    """
    s_star, _ = closest_approach(a0, a1, b0, b1)
    g_min, _ = gap_at(a0, a1, b0, b1, s_star)
    limit = g_min * (1.0 + tolerance) + 1e-9

    ss = np.linspace(0.0, 1.0, samples)
    ok = np.array([gap_at(a0, a1, b0, b1, s)[0] <= limit for s in ss])

    k = int(np.argmin(np.abs(ss - s_star)))
    lo = k
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    hi = k
    while hi < samples - 1 and ok[hi + 1]:
        hi += 1
    return float(ss[lo]), float(ss[hi])


def make_contact(a0, a1, b0, b1,
                 s_a: Optional[float] = None,
                 s_b: Optional[float] = None) -> Contact:
    """Build the shared frame for a chord pair.

    Pass ``s_a`` / ``s_b`` to place the contact deliberately (for several
    contacts along one chain); omit them to use the closest approach.
    """
    a0, a1, b0, b1 = (np.asarray(v, float) for v in (a0, a1, b0, b1))
    if s_a is None or s_b is None:
        s_a, s_b = closest_approach(a0, a1, b0, b1)

    ca = a0 + s_a * (a1 - a0)
    cb = b0 + s_b * (b1 - b0)
    origin = 0.5 * (ca + cb)

    ua, ub = _unit(a1 - a0), _unit(b1 - b0)
    # Flip B if the chords run opposite ways, so the mean is a real mean and
    # not a near-cancellation.
    if float(ua @ ub) < 0.0:
        ub = -ub
    axis = _unit(ua + ub, fallback=ua)

    sep = cb - ca
    toward = sep - float(sep @ axis) * axis
    toward = _unit(toward, fallback=_perp_to(axis))
    across = _unit(np.cross(axis, toward))

    return Contact(origin=origin, axis=axis, toward=toward, across=across,
                   s_a=float(s_a), s_b=float(s_b),
                   gap=float(np.linalg.norm(sep)))


# ---------------------------------------------------------------------------
# Planning
# ---------------------------------------------------------------------------

def axial_room(a0, a1, contact: Contact, margin: float = 0.06
               ) -> tuple[float, float]:
    """Axial extent available to a braid before it reaches a junction.

    Reported in shared-frame axial units, not chord fractions. The two differ
    for a skew chord, and differ in *scale* between chords of unequal length
    -- reconciling that is exactly what lets unequal partners braid about one
    axis.
    """
    a0, a1 = np.asarray(a0, float), np.asarray(a1, float)
    lo_end = a0 + margin * (a1 - a0)
    hi_end = a1 - margin * (a1 - a0)
    u_lo = float((lo_end - contact.origin) @ contact.axis)
    u_hi = float((hi_end - contact.origin) @ contact.axis)
    return (u_lo, u_hi) if u_lo <= u_hi else (u_hi, u_lo)


def plan_braid(a0, a1, b0, b1, contact: Contact, e: int,
               shape: Optional[BraidShape] = None,
               margin: float = 0.06) -> tuple[float, int]:
    """Size a braid to what *both* chains can give it.

    Returns ``(half_span, e_max)``. The half-span is shared: one braid, one
    axis, one phase map. Sizing per chain is what makes unequal partners
    collide -- they reach opposite phase at different axial positions, which
    is not anti-phase at all.

    ``e_max`` is what the room actually supports. Asking beyond it is not an
    error, it compresses the braid and costs clearance, so the caller is told
    rather than silently served.
    """
    shape = shape or BraidShape()
    a_lo, a_hi = axial_room(a0, a1, contact, margin)
    b_lo, b_hi = axial_room(b0, b1, contact, margin)
    lo, hi = max(a_lo, b_lo), min(a_hi, b_hi)
    have = min(abs(lo), abs(hi))            # symmetric about the contact

    half = min(0.5 * shape.span(e), have)
    return float(half), shape.turns_within(have)


# ---------------------------------------------------------------------------
# Path
# ---------------------------------------------------------------------------

def braid_path(a0, a1, contact: Contact, e: int, side: int, n_beads: int,
               half_span: float, shape: Optional[BraidShape] = None
               ) -> np.ndarray:
    """Route one chain through a braid about the shared axis.

    ``side`` is ``+1`` or ``-1``; partners take opposite sides, which is what
    makes them anti-phase and yields linking number ``e``.

    Phase and blend are keyed to the *axial* coordinate in the shared frame,
    never to the chain's own parameter. Keying them to each chain's own ``t``
    is subtly wrong: two chords of different length map the same ``t`` to
    different axial positions, so the partners reach opposite phase in
    different places and close to contact instead of braiding around one
    another. Measured before this was fixed, a 24-vs-8 pair fell to 0.05
    clearance; after, 1.38.

    Phase advances only across the plateau. Letting it advance over the whole
    span lets the ramps absorb the rotation and the realised linking number
    collapses -- symptom: ``e = 1, 2, 3`` all measure 1.
    """
    shape = shape or BraidShape()
    a0, a1 = np.asarray(a0, float), np.asarray(a1, float)
    chord = a1 - a0

    if n_beads <= 0:
        return np.empty((0, 3))
    if half_span <= 0.0:
        ts = np.linspace(0.0, 1.0, n_beads) if n_beads > 1 else np.array([0.5])
        return a0 + ts[:, None] * chord

    ramp_u = min(shape.ramp, 0.45 * 2.0 * half_span)
    p_lo, p_hi = -half_span + ramp_u, half_span - ramp_u

    pts = np.empty((n_beads, 3))
    for k in range(n_beads):
        t = k / (n_beads - 1) if n_beads > 1 else 0.5
        on_chord = a0 + t * chord
        u = float((on_chord - contact.origin) @ contact.axis)

        if u <= -half_span or u >= half_span:
            pts[k] = on_chord
            continue

        if u < p_lo:
            frac, w = 0.0, (u + half_span) / max(p_lo + half_span, 1e-12)
        elif u > p_hi:
            frac, w = 1.0, (half_span - u) / max(half_span - p_hi, 1e-12)
        else:
            frac = (u - p_lo) / max(p_hi - p_lo, 1e-12)
            w = 1.0

        phi = 2.0 * np.pi * e * frac
        w = float(np.clip(w, 0.0, 1.0))
        w = w * w * (3.0 - 2.0 * w)                       # smoothstep

        target = (contact.origin + u * contact.axis
                  + side * shape.n_radius * np.cos(phi) * contact.toward
                  + side * shape.m_radius * np.sin(phi) * contact.across)
        pts[k] = (1.0 - w) * on_chord + w * target

    return pts


def braid_pair(a0, a1, b0, b1, e: int = 1,
               shape: Optional[BraidShape] = None,
               n_beads: int = 128,
               contact: Optional[Contact] = None,
               fit_gap: bool = True
               ) -> tuple[np.ndarray, np.ndarray, int]:
    """Braid two chords. Returns ``(path_a, path_b, e_max)``.

    ``fit_gap`` shrinks an over-wide braid to the partner gap; see
    :meth:`BraidShape.fit_to_gap`. Pass False only to study the unclamped
    geometry, never to build a system.
    """
    shape = shape or BraidShape()
    contact = contact or make_contact(a0, a1, b0, b1)
    if fit_gap:
        shape = shape.fit_to_gap(contact.gap)
    half, e_max = plan_braid(a0, a1, b0, b1, contact, e, shape)
    pa = braid_path(a0, a1, contact, e, -1, n_beads, half, shape)
    pb = braid_path(b0, b1, contact, e, +1, n_beads, half, shape)
    return pa, pb, e_max


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def linking_number(loop_a: np.ndarray, loop_b: np.ndarray) -> float:
    """Gauss linking integral for two closed polylines."""
    a, b = np.asarray(loop_a, float), np.asarray(loop_b, float)
    s1, s2 = a[1:] - a[:-1], b[1:] - b[:-1]
    m1, m2 = 0.5 * (a[1:] + a[:-1]), 0.5 * (b[1:] + b[:-1])
    r = m1[:, None, :] - m2[None, :, :]
    d = np.linalg.norm(r, axis=-1)
    d = np.where(d < 1e-9, 1e-9, d)
    cr = np.cross(s1[:, None, :], s2[None, :, :])
    return float(((r * cr).sum(-1) / d ** 3).sum() / (4.0 * np.pi))


def _close_with_chord(path: np.ndarray, n_fill: int = 40) -> np.ndarray:
    back = np.linspace(path[-1], path[0], n_fill)[1:-1]
    return np.vstack([path, back, path[:1]])


def _close_far(path: np.ndarray, direction: np.ndarray, distance: float,
               n_fill: int = 24) -> np.ndarray:
    """Close a path via a detour far away along ``direction``."""
    d = _unit(np.asarray(direction, float)) * distance
    out = np.linspace(path[-1], path[-1] + d, n_fill)
    across = np.linspace(path[-1] + d, path[0] + d, n_fill)
    back = np.linspace(path[0] + d, path[0], n_fill)
    return np.vstack([path, out[1:], across[1:], back[1:]])


def chord_closed_linking(path_a: np.ndarray, path_b: np.ndarray) -> float:
    """Linking number of two open paths, each closed by its own chord.

    Closing with the chord means the result reads as "does this path wind
    around its partner, relative to the un-entangled straight reference".
    A pair that pulls apart freely measures 0.

    **Validity.** The closure chords are part of the closed loops, so the
    measure requires the braided paths to stay clear of *both* chords. A
    braid whose radius approaches half the partner gap sweeps close to the
    partner's chord and the integral diverges: measured on a single-winding
    braid, gap 3 reads 1.38, gap 2 reads 3.47, gap 1 reads 6.21, while the
    actual chain-chain clearance is unchanged and healthy throughout. Use
    :func:`far_closed_linking` when the gap is tight; it is accurate across
    the whole range and is what the tests use.

    Not valid as a post-MD entanglement count either: once chains crumple,
    |Lk| is not a topological invariant for open curves. Use primitive-path
    analysis there.
    """
    return linking_number(_close_with_chord(np.asarray(path_a, float)),
                          _close_with_chord(np.asarray(path_b, float)))


def far_closed_linking(path_a: np.ndarray, path_b: np.ndarray,
                       contact: Optional[Contact] = None) -> float:
    """Linking number of two open paths, closed through distant detours.

    Each path is closed by a loop that runs far away, the two loops in
    opposite directions, so neither closure passes near either curve. That
    removes the divergence that limits :func:`chord_closed_linking` at tight
    gaps and leaves the winding of the braid itself, which is the quantity
    being verified.

    The detour direction defaults to the contact's ``across`` axis -- the
    thin direction of the braid ellipse, where there is most room -- or to a
    perpendicular of the paths' overall extent when no contact is given.
    """
    a = np.asarray(path_a, float)
    b = np.asarray(path_b, float)

    if contact is not None:
        d = contact.across
    else:
        span = np.vstack([a, b])
        principal = _unit(span[-1] - span[0])
        d = _perp_to(principal)

    reach = float(np.linalg.norm(np.vstack([a, b]).ptp(axis=0))) + 1.0
    far = 25.0 * reach
    return linking_number(_close_far(a, d, far), _close_far(b, -d, far))


def min_separation(path_a: np.ndarray, path_b: np.ndarray) -> float:
    """Closest approach between two paths."""
    a, b = np.asarray(path_a, float), np.asarray(path_b, float)
    return float(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1).min())
