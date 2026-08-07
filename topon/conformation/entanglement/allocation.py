"""Allocate braid contacts so several can share one chain.

A braid is not a property of one chain, it is a property of a pair, and it
occupies a stretch of *both* partners' chords. That makes multi-partner
entanglement an interval-packing problem rather than a per-chain choice:

* chain A may want its own contacts, while also being someone else's partner;
* two contacts landing on the same stretch of A produce blends that fight,
  and the realised winding of both collapses;
* the axial room a contact needs is set by the pair, so a request that fits
  one partner may not fit the other.

Deciding this per chain is the documented failure: a chain that appears in
another chain's partner list *and* carries its own contact ends up with two
braids on the same piece of chord. This module allocates over the whole set
at once.

The formulation: each candidate contact occupies an interval on each of its
two chords. Accept a set of candidates such that no chord carries two
overlapping intervals, and no chain exceeds its partner budget. Greedy by
priority, which is the usual choice for interval packing and is stable and
explainable -- exactness is not worth much here because the candidate set is
itself a physical approximation.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Hashable, Iterable, Optional, Sequence

import numpy as np

from topon.conformation.entanglement.braid import (
    BraidShape,
    Contact,
    axial_room,
    braid_path,
    closest_approach,
    far_closed_linking,
    feasible_window,
    gap_at,
    make_contact,
    plan_braid,
)

__all__ = [
    "ContactRequest",
    "AllocatedContact",
    "Rejection",
    "Allocation",
    "allocate_contacts",
    "compose_chain_path",
]

ChainId = Hashable


@dataclass(frozen=True)
class ContactRequest:
    """A wish for one braid between two chains.

    ``priority`` orders the greedy pass: higher goes first. The natural
    choice is a decreasing function of the partner gap, so genuinely close
    pairs win the room over distant ones, but the caller decides.
    """

    chain_a: ChainId
    chain_b: ChainId
    windings: int = 1
    priority: float = 0.0


@dataclass(frozen=True)
class AllocatedContact:
    """A request that was granted room, with the geometry to build it."""

    request: ContactRequest
    contact: Contact
    half_span: float
    windings: int
    """Windings actually allocated. May be below ``request.windings`` when the
    room allowed fewer -- see ``Rejection`` for requests that got none."""
    shape: BraidShape = field(default_factory=BraidShape)
    """The shape this contact was planned with, already fitted to the pair's
    gap. Build with this, not the caller's original -- a braid wider than the
    gap reaches past its partner's chord."""


@dataclass(frozen=True)
class Rejection:
    """A request that could not be granted, and why."""

    request: ContactRequest
    reason: str


@dataclass
class Allocation:
    accepted: list[AllocatedContact] = field(default_factory=list)
    rejected: list[Rejection] = field(default_factory=list)

    @property
    def partners(self) -> dict[ChainId, list[ChainId]]:
        """Which chains each chain ended up braided with."""
        out: dict[ChainId, list[ChainId]] = defaultdict(list)
        for a in self.accepted:
            out[a.request.chain_a].append(a.request.chain_b)
            out[a.request.chain_b].append(a.request.chain_a)
        return dict(out)

    def summary(self) -> str:
        n_ch = len(self.partners)
        tot = sum(a.windings for a in self.accepted)
        return (f"{len(self.accepted)} contacts over {n_ch} chains, "
                f"{tot} windings, {len(self.rejected)} rejected")


def _interval_on(chord_start, chord_end, contact: Contact,
                 half_span: float) -> tuple[float, float]:
    """Axial interval a braid occupies, expressed along this chord.

    Returned as fractions of the chord so intervals from chords of different
    length and direction can be compared on their own chain.
    """
    a0 = np.asarray(chord_start, float)
    a1 = np.asarray(chord_end, float)
    chord = a1 - a0
    L2 = float(chord @ chord)
    if L2 <= 1e-12:
        return 0.0, 1.0

    # Project the two ends of the axial span onto the chord.
    lo_pt = contact.origin - half_span * contact.axis
    hi_pt = contact.origin + half_span * contact.axis
    t_lo = float((lo_pt - a0) @ chord) / L2
    t_hi = float((hi_pt - a0) @ chord) / L2
    if t_lo > t_hi:
        t_lo, t_hi = t_hi, t_lo
    return max(t_lo, 0.0), min(t_hi, 1.0)


def _span_fraction(iv: tuple[float, float]) -> float:
    """How much of a chord an interval covers, after clamping to it."""
    lo, hi = iv
    return max(0.0, hi - lo)


def _overlaps(iv: tuple[float, float],
              taken: Sequence[tuple[float, float]],
              clearance: float) -> bool:
    lo, hi = iv
    return any(not (hi + clearance <= t_lo or lo - clearance >= t_hi)
               for t_lo, t_hi in taken)


def _min_separation(pa, pb) -> float:
    """Closest approach between two bead paths."""
    d = np.linalg.norm(np.asarray(pa)[:, None, :] - np.asarray(pb)[None, :, :],
                       axis=-1)
    return float(d.min())


def _obstructing_chain(contact: Contact, half_span: float, radius: float,
                       chords, exclude, samples: int = 21):
    """A third chain lying inside the braid volume, if there is one.

    The braid sweeps a tube of roughly ``radius`` about its axis over
    ``±half_span``. Any other chord passing through that tube is in the way:
    the two partners would wind straight through it.

    This is the hazard that makes long-range entanglement different from
    nearest-neighbour entanglement. The axis of a contact sits midway between
    the partners, so on a lattice a second-neighbour contact puts its axis
    *exactly where the first neighbour lies*. Measured on a row of chains 3
    apart: a 1-3 contact placed its axis at the position of chain 2, chain 1
    swung to within 0.90 of chain 2, and the pair picked up an unwanted
    linking number of 1.00 with a chain nobody had asked it to entangle --
    while the prescribed 1-2 braid over-counted from 2 to 3.
    """
    us = np.linspace(-half_span, half_span, samples)
    tube = contact.origin[None, :] + us[:, None] * contact.axis[None, :]

    for cid, (c0, c1) in chords.items():
        if cid in exclude:
            continue
        c0 = np.asarray(c0, float)
        c1 = np.asarray(c1, float)
        d = c1 - c0
        L2 = float(d @ d)
        if L2 <= 1e-12:
            continue
        t = np.clip(((tube - c0) @ d) / L2, 0.0, 1.0)
        closest = c0[None, :] + t[:, None] * d[None, :]
        if float(np.linalg.norm(tube - closest, axis=1).min()) < radius:
            return cid
    return None


def _candidate_positions(a0, a1, b0, b1, samples: int,
                         taken: Sequence[tuple[float, float]] = (),
                         tolerance: float = 0.35) -> list[float]:
    """Positions to try for a contact, best first.

    With nothing on the chord yet, the closest approach comes first and the
    rest spread outward from it, so a lone braid sits at the tightest part
    of the approach.

    Once the chord carries braids, the order flips to farthest-first from
    what is already there. That is what lets several entanglements between
    the same pair sit at, say, a fifth, a half and four fifths of the chain
    instead of piling up against each other in the middle. Spreading them is
    also the only way to raise the count on a pair at all: clearance is set
    by how tight the helix is, so more turns at one site means a tighter
    spiral and a smaller margin, while more sites costs nothing.
    """
    lo, hi = feasible_window(a0, a1, b0, b1, tolerance=tolerance)
    s_star, _ = closest_approach(a0, a1, b0, b1)
    if hi - lo < 1e-6:
        return [s_star]

    grid = np.linspace(lo, hi, max(samples, 3))
    if not taken:
        return [float(s) for s in sorted(grid, key=lambda s: abs(s - s_star))]

    centres = [0.5 * (t_lo + t_hi) for t_lo, t_hi in taken]

    def spread(s):
        # Farthest from any existing braid; the closest-approach distance
        # breaks ties, so among equally free positions the tightest wins.
        return (-min(abs(s - c) for c in centres), abs(s - s_star))

    return [float(s) for s in sorted(grid, key=spread)]


def allocate_contacts(
    requests: Iterable[ContactRequest],
    chords: dict[ChainId, tuple[Sequence[float], Sequence[float]]],
    shape: Optional[BraidShape] = None,
    max_partners: Optional[int] = None,
    min_clearance: float = 1.0,
    verify_windings: bool = True,
    window_tolerance: float = 0.35,
    separation: float = 0.02,
    window_samples: int = 25,
    check_obstruction: bool = True,
    clearance: float = 0.5,
    min_overlap: float = 0.05,
) -> Allocation:
    """Grant room to as many requests as the chords can carry.

    Parameters
    ----------
    requests : the wished-for contacts.
    chords : ``{chain_id: (start, end)}`` for every chain named in requests.
    shape : braid geometry; sets how much axial room a winding needs.
    max_partners : cap on distinct partners per chain. None means no cap
        beyond what the geometry allows.
    verify_windings : reject a position whose built braid does not measure
        the winding count it was planned for. The measure is the Gauss
        integral over chords closed far away, which is exact for a pair
        running parallel and unreliable for one whose chords diverge, since
        the closure then contributes. Worth keeping on for parallel work and
        worth turning off when the arbiter is a primitive-path analysis.
    window_tolerance : how far the pair's gap may grow, as a fraction of its
        minimum, for a position still to be considered. The default keeps
        every braid near the tightest part of the approach, which is right
        for one site and is exactly what stops a second: measured on a pair
        whose chords diverge, a tolerance of 0.35 left a window covering 42%
        of the chord while one site occupies 24%, so two could not fit
        however much chord remained. Raising it to 1.0 opens 73% and room
        for three, at the cost of siting the outer braids where the chains
        are further apart and the detour to reach across is longer.
    min_clearance : how close the two partners of a braid may come, in the
        units of the chords. A braid squeezed into short chord is a braid
        whose partners are closer together, so this is the floor on how much
        squeezing is allowed; past it the winding count is cut instead.
    separation : minimum gap between two braids on one chord, as a fraction
        of that chord. Prevents adjacent blends from overlapping.
    check_obstruction : refuse a contact whose braid volume contains a third
        chain. Matters most for long-range pairs: the braid axis sits midway
        between the partners, which on a lattice is exactly where the chain
        between them lies.
    clearance : margin added to the braid radius when testing for an
        obstructing chain.
    min_overlap : least fraction of each chord the braid must occupy. Guards
        against a contact whose origin falls in empty space between two
        chains that never come near one another.

    Notes
    -----
    A request whose room allows fewer windings than asked is *granted at the
    lower count* rather than refused, and the shortfall is visible in
    ``AllocatedContact.windings``. Refusing outright would silently drop
    entanglements the user asked for; serving the full count would silently
    compress the braid and cost clearance. Reporting the reduction is the
    only option that keeps both facts.
    """
    shape = shape or BraidShape()
    taken: dict[ChainId, list[tuple[float, float]]] = defaultdict(list)
    partner_of: dict[ChainId, set[ChainId]] = defaultdict(set)
    out = Allocation()

    ordered = sorted(requests, key=lambda r: -r.priority)

    for req in ordered:
        if req.chain_a not in chords or req.chain_b not in chords:
            out.rejected.append(Rejection(req, "chain has no chord"))
            continue
        if req.chain_a == req.chain_b:
            out.rejected.append(Rejection(req, "self-pair"))
            continue

        a0, a1 = chords[req.chain_a]
        b0, b1 = chords[req.chain_b]

        if max_partners is not None:
            # A repeat pairing costs no new partner slot on either side.
            over_budget = any(
                other not in partner_of[me] and len(partner_of[me]) >= max_partners
                for me, other in ((req.chain_a, req.chain_b),
                                  (req.chain_b, req.chain_a))
            )
            if over_budget:
                out.rejected.append(Rejection(req, "partner budget full"))
                continue

        # Try the natural contact first, then slide along the stretch where
        # the pair stays close. Parallel chords are equidistant along their
        # whole overlap, so without sliding every pair on a lattice row lands
        # on the same midpoint and only the first is ever placed.
        placed = None
        blocked_by = None
        too_far = False
        too_close = False
        unrealised = False
        for s_a in _candidate_positions(a0, a1, b0, b1, window_samples,
                                        taken[req.chain_a], window_tolerance):
            gap, s_b = gap_at(a0, a1, b0, b1, s_a)
            contact = make_contact(a0, a1, b0, b1, s_a=s_a, s_b=s_b)
            # Two independent fits. Narrow the braid to this pair's gap so
            # it cannot reach past the partner's chord, then shorten it to
            # the axial room this chord actually has. A braid carrying one
            # absolute size is wrong at both ends of a lattice's range: too
            # long for a melt chord of 4.4 sigma, too narrow to reach across
            # a dilute gap of 26.
            fitted = shape.fit_to_gap(contact.gap)
            # Room is what BOTH chords can spare: a braid is one shared
            # axis, so the shorter partner's budget is the binding one.
            la, ha = axial_room(a0, a1, contact)
            lb, hb = axial_room(b0, b1, contact)
            room = 2.0 * min(ha, hb, -la, -lb)
            if room > 0.0:
                fitted = fitted.fit_to_room(room, req.windings,
                                            min_clearance)
            half, e_max = plan_braid(a0, a1, b0, b1, contact,
                                     req.windings, fitted)
            if e_max < 1:
                continue

            granted = min(req.windings, e_max)
            half, _ = plan_braid(a0, a1, b0, b1, contact, granted, fitted)

            # Verify on the built paths rather than trusting the plan. Two
            # things are checked, and both were caught by measurement rather
            # than predicted.
            #
            # Clearance, because pitch_for_clearance inverts a relation
            # measured on long chords and reads optimistic when the braid
            # spans most of a short one, where the two ramps sit close
            # enough to interact: planned for 1.0, one such braid measured
            # 0.81. Dropping a winding lengthens the pitch, which is what
            # clearance depends on, so that converges.
            #
            # The winding count, because a position can satisfy every budget
            # and still not link. Away from the closest approach the pair is
            # further apart, and a braid whose radius is capped for the
            # tightest part of the window no longer spans the gap out at its
            # edge: granted two, realised one. That is not recoverable by
            # dropping a turn, so the position is abandoned for the next.
            while granted >= 1:
                pa = braid_path(a0, a1, contact, granted, -1, 400, half, fitted)
                pb = braid_path(b0, b1, contact, granted, +1, 400, half, fitted)
                if _min_separation(pa, pb) >= min_clearance:
                    break
                granted -= 1
                if granted >= 1:
                    half, _ = plan_braid(a0, a1, b0, b1, contact,
                                         granted, fitted)
            if granted < 1:
                too_close = True
                continue
            if verify_windings and round(
                    abs(far_closed_linking(pa, pb, contact))) != granted:
                unrealised = True
                continue

            iv_a = _interval_on(a0, a1, contact, half)
            iv_b = _interval_on(b0, b1, contact, half)

            # The braid must actually land on both chords. When two chains
            # pass nowhere near each other -- perpendicular strands a full
            # lattice spacing apart, say -- the closest approach is at their
            # far ends and the contact origin falls in empty space between
            # them. Planning still succeeds and the braid is built, but not
            # one bead of either chain lies inside its span, so the chains
            # are left straight and the pair reads as unlinked. Measured on a
            # 4x4x4 SC lattice: ten of fourteen accepted contacts realised
            # nothing at all, every one of them this case.
            if (_span_fraction(iv_a) < min_overlap
                    or _span_fraction(iv_b) < min_overlap):
                too_far = True
                continue

            if (_overlaps(iv_a, taken[req.chain_a], separation)
                    or _overlaps(iv_b, taken[req.chain_b], separation)):
                continue

            if check_obstruction:
                blocker = _obstructing_chain(
                    contact, half, fitted.n_radius + clearance,
                    chords, {req.chain_a, req.chain_b})
                if blocker is not None:
                    blocked_by = blocker
                    continue

            placed = (contact, half, granted, iv_a, iv_b, fitted)
            break

        if placed is None:
            if blocked_by is not None:
                reason = "a third chain lies in the braid volume"
            elif too_far:
                reason = "chains do not approach along a shared stretch"
            elif too_close:
                reason = "no winding fits without breaking the clearance floor"
            elif unrealised:
                reason = "no position where the braid realises its windings"
            else:
                reason = "no free stretch with room for a winding"
            out.rejected.append(Rejection(req, reason))
            continue

        contact, half, granted, iv_a, iv_b, fitted = placed
        taken[req.chain_a].append(iv_a)
        taken[req.chain_b].append(iv_b)
        partner_of[req.chain_a].add(req.chain_b)
        partner_of[req.chain_b].add(req.chain_a)
        out.accepted.append(
            AllocatedContact(request=req, contact=contact,
                             half_span=half, windings=granted, shape=fitted))

    return out


def compose_chain_path(chain: ChainId, alloc: Allocation,
                       chords: dict[ChainId, tuple[Sequence[float], Sequence[float]]],
                       n_beads: int,
                       shape: Optional[BraidShape] = None) -> np.ndarray:
    """Build one chain's full path, carrying every braid allocated to it.

    A chain with three partners needs all three braids in a *single* path.
    Building one path per pair and checking those in isolation is a trap: it
    verifies geometry the simulation will never see, because the chain that
    reaches the writer has only one set of coordinates.

    The allocator guarantees the braid intervals on a chord do not overlap,
    so composition is a per-bead choice: a bead inside a braid's axial span
    follows that braid, and every other bead stays on the chord.
    """
    shape = shape or BraidShape()
    a0, a1 = (np.asarray(v, float) for v in chords[chain])
    chord = a1 - a0

    mine = [a for a in alloc.accepted
            if a.request.chain_a == chain or a.request.chain_b == chain]

    ts = np.linspace(0.0, 1.0, n_beads) if n_beads > 1 else np.array([0.5])
    path = a0 + ts[:, None] * chord
    if not mine:
        return path

    for a in mine:
        # `side` distinguishes the partners: the request's first chain takes
        # -1, the second +1, which is what makes them anti-phase.
        side = -1 if a.request.chain_a == chain else +1
        braided = braid_path(a0, a1, a.contact, a.windings, side,
                             n_beads, a.half_span, a.shape)
        u = (path - a.contact.origin) @ a.contact.axis
        inside = np.abs(u) < a.half_span
        path[inside] = braided[inside]

    return path
