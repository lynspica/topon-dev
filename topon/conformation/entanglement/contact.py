"""Wind two chains around each other where their paths already meet.

STATUS: works on clean geometry, does NOT yet work on real coiled chains.
Two parallel chords 1.2 sigma apart, asked for 1, 2 and 3 turns, measure
1.00, 2.04 and 3.13, at a contour cost of 0.41, 1.58 and 3.32 sigma and with
their separation unchanged. On random-walk chains from a real lattice the
same code produces correct bond lengths and no windings at all: 0 of 14 by
the Gauss integral, 1 of 14 by Z1+. Four fixes were tried -- putting B back
in its own periodic image, an adaptive window that follows the contact, a
minimum window derived from the bond limit, and one axis for the window
instead of a per-bead tangent -- and each fixed what it targeted without
moving the winding count. Do not use it for production paths yet.


The other construction in this package (:mod:`waypoints`) sites an
entanglement between two *chords* and sends both chains to the midpoint. That
costs contour in proportion to how far apart the chords are, and it is the
reason that construction refuses on uniform lattices and on any pair that is
not already a near neighbour.

This one does not move the chains anywhere. It finds a place where the two
paths, as built, come close, and rotates their separation about the local
axis. The midpoint of the pair is left exactly where it was; only the
*relative* position of the two chains turns.

Three things follow from rotating the offset rather than displacing the
chains:

* the cost is the circumference of a circle of radius ``sep/2``, which for
  chains already 1 sigma apart is about 3 sigma of contour, against the tens
  of sigma the midpoint construction needs;
* the rotation runs from 0 to ``2*pi*turns`` across the window, and both
  ends are whole turns, so the offset returns to itself and the path outside
  the window is untouched -- no blend, no ramp, nothing to tune;
* the winding count is exactly ``turns``, because that is how many times the
  two chains swap sides.

The one requirement is that the chains meet at all, which is a question about
how much slack they carry. A chain at 1.8 times its chord is nearly extended
and meets nothing; from about 6 times, chains interpenetrate and contacts are
everywhere. That ratio is the caller's to choose and is the whole ballgame.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["Contact", "find_contacts", "wind_at", "wind_all"]


@dataclass(frozen=True)
class Contact:
    """A place where two built paths come close enough to wind."""

    chain_a: object
    chain_b: object
    i_a: int
    i_b: int
    sep: float

    @property
    def pair(self):
        return (self.chain_a, self.chain_b)


def _min_image(d, box):
    return d if box is None else d - np.asarray(box, float) * np.round(
        d / np.asarray(box, float))


def find_contacts(paths, box=None, max_sep=2.5, margin=6, exclude=(),
                  min_sep=0.3):
    """Where every pair of paths comes closest, if that is close enough.

    ``margin`` beads are ignored at each end: chains that share a crosslink
    touch there by construction, and a winding placed on top of a junction is
    not an entanglement between the two chains.

    ``exclude`` is a set of frozenset pairs to skip, normally the pairs that
    share a crosslink.

    ``min_sep`` rejects pairs that are effectively on top of one another. The
    frame is built from their separation, so at a few hundredths of a sigma
    it is numerical noise and the rotation is meaningless.
    """
    keys = sorted(paths, key=str)
    out = []
    for i, ka in enumerate(keys):
        pa = paths[ka]
        if len(pa) <= 2 * margin + 2:
            continue
        a = pa[margin:-margin]
        ca = pa.mean(axis=0)
        for kb in keys[i + 1:]:
            if frozenset((ka, kb)) in exclude:
                continue
            pb = paths[kb]
            if len(pb) <= 2 * margin + 2:
                continue
            if np.linalg.norm(_min_image(ca - pb.mean(axis=0), box)) > 60.0:
                # Cheap reject; the real test is per bead below.
                pass
            b = pb[margin:-margin]
            d = _min_image(b[None, :, :] - a[:, None, :], box)
            n = np.linalg.norm(d, axis=-1)
            j = int(n.argmin())
            sep = float(n.flat[j])
            if not (min_sep < sep <= max_sep):
                continue
            ia, ib = np.unravel_index(j, n.shape)
            out.append(Contact(ka, kb, int(ia) + margin, int(ib) + margin, sep))
    out.sort(key=lambda c: c.sep)
    return out


def _pair_window(pa, pb, i_a, i_b, half, box, widen=3.0, turns=1):
    """Bead indices on each chain either side of the contact, and the sense.

    The two chains may run the same way or opposite ways through the contact.
    Walking B backwards when they are antiparallel is what keeps the two
    windows alongside each other rather than crossing.

    The window stops growing once the pair has drifted ``widen`` times its
    contact separation apart. Two chains that cross at an angle rather than
    running alongside separate quickly, and a fixed window then reaches out
    to where they are far apart; rotating a large offset flings beads around
    and produced 5.4 sigma bonds. Following the contact keeps the offset
    small, which is what makes the rotation cheap in the first place.
    """
    ta = pa[min(i_a + 1, len(pa) - 1)] - pa[max(i_a - 1, 0)]
    tb = pb[min(i_b + 1, len(pb) - 1)] - pb[max(i_b - 1, 0)]
    sense = 1 if float(ta @ tb) >= 0.0 else -1

    sep0 = float(np.linalg.norm(_min_image(pb[i_b] - pa[i_a], box)))
    cap = max(widen * sep0, 1.5)

    def reach(step):
        n = 0
        while n < half:
            ja = i_a + step * (n + 1)
            jb = i_b + sense * step * (n + 1)
            if not (0 <= ja < len(pa) and 0 <= jb < len(pb)):
                break
            if float(np.linalg.norm(_min_image(pb[jb] - pa[ja], box))) > cap:
                break
            n += 1
        return n

    lo, hi = reach(-1), reach(+1)
    # A full turn of an offset r spread over N beads moves neighbouring beads
    # about 2*pi*r/N apart, so N sets the bond length the winding costs. With
    # r near 1 sigma, keeping that under a few tenths needs roughly 20 beads
    # per turn. Accepting whatever the adaptive walk returned, sometimes 4,
    # is what left 3.6 sigma bonds.
    if lo + hi < _min_window(turns, max(sep0, 0.5)):
        return None
    return np.arange(-lo, hi + 1), sense


def _min_window(turns=1, offset=1.0, max_step=0.30):
    """Beads a winding needs so it does not stretch the bonds it rides on."""
    return int(np.ceil(2.0 * np.pi * turns * offset / max_step))


def wind_at(pa, pb, contact, turns=1, half=10, box=None):
    """Rotate the two chains' separation about the local axis ``turns`` times.

    Returns new copies of both paths. Beads outside the window are untouched,
    and the pair's midpoint is unchanged everywhere, so neither chain moves
    away from where the conformation put it.
    """
    pa = np.array(pa, dtype=float, copy=True)
    pb = np.array(pb, dtype=float, copy=True)

    got = _pair_window(pa, pb, contact.i_a, contact.i_b, half, box,
                       turns=turns)
    if got is None:
        return pa, pb
    ks, sense = got

    ia = contact.i_a + ks
    ib = contact.i_b + sense * ks

    A = pa[ia]
    B = pb[ib]
    # Work in the image of B nearest A, then put B back in its own image
    # afterwards. Rotating in A's frame and writing the result straight into
    # B moves every bead of B a box length: a contact found across the
    # boundary produced a 52 sigma bond that way.
    off = 0.5 * _min_image(B - A, box)
    mid = A + off
    b_shift = B - (A + 2.0 * off)

    # One axis for the whole window, from end to end of the shared stretch.
    # Taking the local tangent per bead sounds better and is not: on a
    # random-walk chain that direction is noise, so the offset turns about a
    # different axis every step and never accumulates a winding. Measured
    # with the per-bead axis: bonds perfect, 0 of 14 pairs linked.
    span = mid[-1] - mid[0]
    if float(np.linalg.norm(span)) < 1e-9:
        span = np.gradient(mid, axis=0).sum(axis=0)
    axis = span / (np.linalg.norm(span) + 1e-12)
    # Keep only the part of the offset perpendicular to it, or the component
    # along the axis rides around unchanged and dilutes the turn.
    off = off - np.outer(off @ axis, axis)
    axis = np.broadcast_to(axis, off.shape)

    # 0 to 2*pi*turns across the window. Both ends are whole turns, so the
    # offset comes back to itself and the join needs no blending.
    t = (ks - ks[0]) / float(ks[-1] - ks[0])
    ang = 2.0 * np.pi * turns * t

    # Rodrigues, per bead, about that bead's own axis.
    c = np.cos(ang)[:, None]
    s = np.sin(ang)[:, None]
    dot = np.sum(axis * off, axis=1)[:, None]
    rot = off * c + np.cross(axis, off) * s + axis * dot * (1.0 - c)

    pa[ia] = mid - rot
    pb[ib] = mid + rot + b_shift
    return pa, pb


def wind_all(paths, contacts, turns=1, half=10, box=None, min_gain=0.0):
    """Apply :func:`wind_at` to a list of contacts, one chain at a time.

    Contacts are taken in order. A chain may appear in several, and each is
    applied to the current paths, so later windings see the earlier ones.
    Windows that would overlap on the same chain are skipped, since two
    rotations over one stretch do not compose into a prescribed count.
    """
    out = {k: np.array(p, dtype=float, copy=True) for k, p in paths.items()}
    taken = {}
    applied = []
    for c in contacts:
        span_a = (c.i_a - half, c.i_a + half)
        span_b = (c.i_b - half, c.i_b + half)
        busy = False
        for chain, span in ((c.chain_a, span_a), (c.chain_b, span_b)):
            for lo, hi in taken.get(chain, ()):
                if not (span[1] < lo or span[0] > hi):
                    busy = True
                    break
            if busy:
                break
        if busy:
            continue
        na, nb = wind_at(out[c.chain_a], out[c.chain_b], c, turns, half, box)
        out[c.chain_a] = na
        out[c.chain_b] = nb
        taken.setdefault(c.chain_a, []).append(span_a)
        taken.setdefault(c.chain_b, []).append(span_b)
        applied.append(c)
    return out, applied
