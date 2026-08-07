"""Global contact allocation: several braids on one chain, without collisions.

The documented failure this prevents: a chain that is a partner in someone
else's list *and* carries its own contact gets two braids on the same stretch
of chord, and the blends fight. Deciding per chain cannot see that, because
each decision looks locally fine.
"""

import numpy as np
import pytest

from topon.conformation.entanglement import (
    compose_chain_path,
    BraidShape,
    ContactRequest,
    allocate_contacts,
    braid_path,
    far_closed_linking,
    min_separation,
)

SHAPE = BraidShape()


def _parallel_chords(n, length=60.0, spacing=3.0):
    """n parallel chords in a row, all the same length."""
    return {i: (np.array([0.0, i * spacing, 0.0]),
                np.array([length, i * spacing, 0.0]))
            for i in range(n)}


def _build(alloc, chords, n_beads=400):
    """Realise every accepted contact. Returns {(pair): (path_a, path_b)}."""
    out = {}
    for a in alloc.accepted:
        r = a.request
        a0, a1 = chords[r.chain_a]
        b0, b1 = chords[r.chain_b]
        pa = braid_path(a0, a1, a.contact, a.windings, -1, n_beads, a.half_span)
        pb = braid_path(b0, b1, a.contact, a.windings, +1, n_beads, a.half_span)
        out[(r.chain_a, r.chain_b)] = (pa, pb)
    return out


# ---------------------------------------------------------------------------
# The property the module exists for
# ---------------------------------------------------------------------------

def test_two_contacts_on_one_chain_do_not_share_chord_stretch():
    """Chain 1 partners with both 0 and 2. The braids must sit apart."""
    chords = _parallel_chords(3)
    alloc = allocate_contacts(
        [ContactRequest(0, 1, windings=1, priority=2.0),
         ContactRequest(1, 2, windings=1, priority=1.0)],
        chords, SHAPE)

    assert len(alloc.accepted) == 2, alloc.summary()
    # Both braids live on chain 1; their axial centres must differ.
    centres = [a.contact.origin[0] for a in alloc.accepted]
    spans = [a.half_span for a in alloc.accepted]
    assert abs(centres[0] - centres[1]) > (spans[0] + spans[1]) * 0.5


def test_a_chain_can_carry_several_partners():
    chords = _parallel_chords(4, length=90.0)
    alloc = allocate_contacts(
        [ContactRequest(1, 0, priority=3.0),
         ContactRequest(1, 2, priority=2.0),
         ContactRequest(1, 3, priority=1.0)],
        chords, SHAPE)
    assert len(alloc.partners.get(1, [])) >= 2, alloc.summary()


def test_every_allocated_contact_realises_its_winding_count():
    """Allocation must not quietly break the geometry it is scheduling."""
    chords = _parallel_chords(3, length=90.0)
    alloc = allocate_contacts(
        [ContactRequest(0, 1, windings=2, priority=2.0),
         ContactRequest(1, 2, windings=1, priority=1.0)],
        chords, SHAPE)
    built = _build(alloc, chords)

    assert built, alloc.summary()
    for a in alloc.accepted:
        pa, pb = built[(a.request.chain_a, a.request.chain_b)]
        assert round(far_closed_linking(pa, pb)) == a.windings
        assert min_separation(pa, pb) > 0.5


# ---------------------------------------------------------------------------
# Budgets
# ---------------------------------------------------------------------------

def test_short_chords_reject_rather_than_overlap():
    """When the room runs out, requests are refused, not stacked."""
    chords = _parallel_chords(4, length=9.0)
    alloc = allocate_contacts(
        [ContactRequest(1, 0, priority=3.0),
         ContactRequest(1, 2, priority=2.0),
         ContactRequest(1, 3, priority=1.0)],
        chords, SHAPE)
    assert alloc.rejected
    assert len(alloc.partners.get(1, [])) < 3


def test_max_partners_is_enforced():
    chords = _parallel_chords(4, length=120.0)
    alloc = allocate_contacts(
        [ContactRequest(1, 0, priority=3.0),
         ContactRequest(1, 2, priority=2.0),
         ContactRequest(1, 3, priority=1.0)],
        chords, SHAPE, max_partners=1)
    assert len(alloc.partners.get(1, [])) == 1
    assert any(r.reason == "partner budget full" for r in alloc.rejected)


def test_windings_are_reduced_not_refused_when_room_is_short():
    """A request too ambitious for the chord is granted at the count that
    fits, and the shortfall stays visible.

    Refusing outright would silently drop an entanglement the user asked
    for; serving the full count would silently compress the braid and cost
    clearance. Granting less and saying so keeps both facts.
    """
    chords = _parallel_chords(2, length=14.0)
    alloc = allocate_contacts([ContactRequest(0, 1, windings=6)],
                              chords, SHAPE)
    assert len(alloc.accepted) == 1
    granted = alloc.accepted[0]
    assert granted.windings < granted.request.windings
    assert granted.windings >= 1


def test_priority_decides_who_gets_the_room():
    chords = _parallel_chords(3, length=11.0)
    hi = allocate_contacts(
        [ContactRequest(0, 1, priority=10.0), ContactRequest(1, 2, priority=1.0)],
        chords, SHAPE)
    assert hi.accepted[0].request.chain_a == 0
    lo = allocate_contacts(
        [ContactRequest(0, 1, priority=1.0), ContactRequest(1, 2, priority=10.0)],
        chords, SHAPE)
    assert lo.accepted[0].request.chain_b == 2


# ---------------------------------------------------------------------------
# Obstruction: the hazard specific to long-range pairs
# ---------------------------------------------------------------------------

def test_second_neighbour_contact_blocked_by_the_chain_between():
    """A contact's axis sits midway between its partners, so on a lattice a
    second-neighbour pair puts its axis exactly where the first neighbour is.

    Left unchecked this is not a near miss. Measured on this arrangement: the
    1-3 braid drove chain 1 to within 0.90 of chain 2, produced an unwanted
    linking number of 1.00 with chain 2, and pushed the prescribed 1-2 braid
    from 2 windings to 3.
    """
    chords = _parallel_chords(4, length=90.0)
    alloc = allocate_contacts([ContactRequest(1, 3, windings=1)], chords, SHAPE)

    assert not alloc.accepted
    assert alloc.rejected[0].reason == "a third chain lies in the braid volume"


def test_the_same_pair_is_fine_once_the_obstruction_is_gone():
    """It is the intervening chain that blocks it, not the distance."""
    chords = _parallel_chords(4, length=90.0)
    del chords[2]                                   # remove the chain between
    alloc = allocate_contacts([ContactRequest(1, 3, windings=1)], chords, SHAPE)

    assert len(alloc.accepted) == 1
    paths = {c: compose_chain_path(c, alloc, chords, 800) for c in chords}
    a = alloc.accepted[0]
    assert round(far_closed_linking(paths[1], paths[3], a.contact)) == 1


def test_no_unwanted_links_across_a_whole_allocation():
    """The property that matters at system scale: every prescribed pair
    realises its count, and no unprescribed pair links at all."""
    chords = _parallel_chords(5, length=90.0)
    alloc = allocate_contacts(
        [ContactRequest(1, 0, windings=1, priority=4.0),
         ContactRequest(1, 2, windings=2, priority=3.0),
         ContactRequest(1, 3, windings=1, priority=2.0),
         ContactRequest(3, 4, windings=1, priority=1.0)],
        chords, SHAPE)
    paths = {c: compose_chain_path(c, alloc, chords, 900) for c in chords}

    prescribed = set()
    for a in alloc.accepted:
        r = a.request
        prescribed.add(frozenset((r.chain_a, r.chain_b)))
        assert round(far_closed_linking(
            paths[r.chain_a], paths[r.chain_b], a.contact)) == a.windings

    for i in chords:
        for j in chords:
            if i < j and frozenset((i, j)) not in prescribed:
                assert abs(far_closed_linking(paths[i], paths[j])) < 0.5, \
                    f"chains {i}-{j} linked without being asked to"


def test_obstruction_check_can_be_disabled():
    """Kept switchable: the check costs a pass over every chord, and a caller
    that has already filtered its candidates may not need it."""
    chords = _parallel_chords(4, length=90.0)
    alloc = allocate_contacts([ContactRequest(1, 3)], chords, SHAPE,
                              check_obstruction=False)
    assert len(alloc.accepted) == 1


# ---------------------------------------------------------------------------
# Bad input
# ---------------------------------------------------------------------------

def test_self_pair_is_rejected():
    chords = _parallel_chords(2)
    alloc = allocate_contacts([ContactRequest(0, 0)], chords, SHAPE)
    assert not alloc.accepted
    assert alloc.rejected[0].reason == "self-pair"


def test_unknown_chain_is_rejected():
    chords = _parallel_chords(2)
    alloc = allocate_contacts([ContactRequest(0, 99)], chords, SHAPE)
    assert not alloc.accepted
    assert alloc.rejected[0].reason == "chain has no chord"


def test_empty_request_set_is_fine():
    alloc = allocate_contacts([], _parallel_chords(2), SHAPE)
    assert not alloc.accepted and not alloc.rejected


def test_allocation_is_deterministic():
    chords = _parallel_chords(4, length=90.0)
    reqs = [ContactRequest(0, 1, priority=3.0),
            ContactRequest(1, 2, priority=2.0),
            ContactRequest(2, 3, priority=1.0)]
    first = allocate_contacts(reqs, chords, SHAPE)
    second = allocate_contacts(reqs, chords, SHAPE)
    assert ([a.request for a in first.accepted]
            == [a.request for a in second.accepted])
