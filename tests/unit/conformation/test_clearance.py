"""Paths drawn around the beads that are already there.

The point of `Clearance` is that a routed path must not start inside another
bead. That is not a tidiness concern: overlapping beads are the one condition
under which minimisation can push two chains through each other, and a chain
crossing is the only thing that can change a topology once it is built.
"""
import numpy as np
import pytest

from topon.conformation.paths import (
    Clearance,
    bridging_walk,
    loop_around,
    walk_through,
)

BOX = np.array([20.0, 20.0, 20.0])


@pytest.fixture
def melt():
    """A relaxed-melt stand-in: points on a lattice, so there are real gaps."""
    g = np.arange(0.5, 20.0, 1.4)
    return np.array([[x, y, z] for x in g for y in g for z in g])


def test_reports_distance_to_the_nearest_bead():
    c = Clearance(np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]]), radius=0.9)
    assert c.near([[2.0, 0.0, 0.0]])[0] == pytest.approx(2.0)
    assert c.near([[4.0, 0.0, 0.0]])[0] == pytest.approx(1.0)


def test_distance_is_minimum_image():
    c = Clearance(np.array([[0.5, 10.0, 10.0]]), BOX, radius=0.9)
    # Across the periodic face, not the long way round the box.
    assert c.near([[19.5, 10.0, 10.0]])[0] == pytest.approx(1.0)


def test_a_point_on_the_upper_face_is_accepted():
    # cKDTree's periodic box is half-open and refuses the whole tree over a
    # single point sitting exactly on the boundary.
    c = Clearance(np.array([[20.0, 20.0, 20.0], [1.0, 1.0, 1.0]]), BOX)
    assert np.isfinite(c.near([[10.0, 10.0, 10.0]])[0])


def test_nothing_in_the_box_means_nothing_to_avoid():
    c = Clearance(np.empty((0, 3)), BOX)
    assert c.near([[1.0, 2.0, 3.0]])[0] == np.inf
    assert c.ok([[1.0, 2.0, 3.0]])


def test_worst_is_the_tightest_contact_on_the_whole_path():
    c = Clearance(np.array([[0.0, 0.0, 0.0]]), radius=0.9)
    path = np.array([[3.0, 0.0, 0.0], [0.4, 0.0, 0.0], [5.0, 0.0, 0.0]])
    assert c.worst(path) == pytest.approx(0.4)
    assert not c.ok(path)


def test_avoiding_beats_not_avoiding(melt):
    """Aggregate, not a single draw.

    Choosing the first acceptable step is greedy, so on one seed it can walk
    into a pocket a blind draw happened to miss. The claim is about the
    typical path, and it is a large effect: measured drawing a chain through a
    real relaxed melt, the tightest contact went from 0.081 to 0.822 sigma.

    Interior beads only. The two junctions are fixed and shared, so whichever
    of them lies closest to a bead would otherwise set `worst` for both.
    """
    a, b = np.array([1.0, 1.0, 1.0]), np.array([15.0, 15.0, 15.0])
    c = Clearance(melt, BOX, radius=0.6)
    blind = [c.worst(bridging_walk(a, b, 40, 0.97,
                                   np.random.default_rng(s))[1:-1])
             for s in range(12)]
    aware = [c.worst(bridging_walk(a, b, 40, 0.97,
                                   np.random.default_rng(s), c)[1:-1])
             for s in range(12)]
    assert np.median(aware) > np.median(blind)


def test_avoiding_still_lands_on_the_far_junction(melt):
    a, b = np.array([1.0, 1.0, 1.0]), np.array([12.0, 9.0, 14.0])
    c = Clearance(melt, BOX, radius=0.6)
    p = bridging_walk(a, b, 40, 0.97, np.random.default_rng(5), c)
    assert p[-1] == pytest.approx(b)
    assert p[0] == pytest.approx(a)


def test_avoiding_does_not_rescale_any_bond(melt):
    c = Clearance(melt, BOX, radius=0.6)
    p = bridging_walk(np.zeros(3), np.array([9.0, 3.0, 2.0]), 30, 0.97,
                      np.random.default_rng(11), c)
    assert np.linalg.norm(np.diff(p, axis=0), axis=1) == pytest.approx(0.97)


def test_an_empty_clearance_changes_nothing_that_matters():
    """Not the same path, but the same kind of path.

    With obstacles the walk draws `tries` steps at once and keeps the first
    acceptable one, so it consumes the generator differently and the sequence
    is not the one a plain walk produces. What must not change is what the
    walk guarantees.
    """
    a, b = np.zeros(3), np.array([8.0, 2.0, 1.0])
    p = bridging_walk(a, b, 30, 0.97, np.random.default_rng(2),
                      Clearance(np.empty((0, 3))))
    assert p[0] == pytest.approx(a)
    assert p[-1] == pytest.approx(b)
    # The closing bond is whatever is left over when the walk arrives, so it
    # is the one that need not be exact.
    assert np.linalg.norm(np.diff(p, axis=0), axis=1) == pytest.approx(
        0.97, abs=1e-3)
    assert len(p) == 31


def test_walk_through_threads_clearance_into_every_leg(melt):
    """The legs improve; the waypoints are landed on either way.

    `walk_through` hits each waypoint exactly, by contract, so a waypoint
    dropped on an occupied site stays there however the legs are drawn. That
    is why `loop_around` has to do its own clearing rather than leave it to
    the walk.
    """
    c = Clearance(melt, BOX, radius=0.6)
    way = [np.array([5.8, 5.8, 5.8]), np.array([9.6, 9.6, 9.6])]
    a, b = np.array([1.0, 1.0, 1.0]), np.array([15.0, 15.0, 15.0])
    blind = [c.worst(walk_through(a, b, way, 60, 0.97,
                                  np.random.default_rng(s))[1:-1])
             for s in range(8)]
    aware = [c.worst(walk_through(a, b, way, 60, 0.97,
                                  np.random.default_rng(s), c)[1:-1])
             for s in range(8)]
    assert np.median(aware) > np.median(blind)


# ---------------------------------------------------------------------------
# The ring
# ---------------------------------------------------------------------------

def _target():
    z = np.linspace(0.0, 20.0, 81)
    return np.column_stack([np.full_like(z, 10.0), np.full_like(z, 10.0), z])


def test_ring_finds_more_room_than_a_rigid_rotation(melt):
    c = Clearance(melt, BOX, radius=0.9)
    t = _target()
    rigid = loop_around(t, 40, 2.0, 6, 0.3)
    per_point = loop_around(t, 40, 2.0, 6, 0.3, c)
    assert c.worst(per_point) > c.worst(rigid)


def test_ring_points_stay_in_their_own_sector(melt):
    """Each waypoint may move, but not past its neighbour.

    That is what keeps the winding: a path through the points in order goes
    once around only while their angular order is the one it was given.
    """
    c = Clearance(melt, BOX, radius=0.9)
    t = _target()
    n = 6
    ring = loop_around(t, 40, 2.0, n, 0.0, c)
    ang = np.unwrap(np.arctan2(ring[:, 1] - 10.0, ring[:, 0] - 10.0))
    assert np.all(np.diff(ang) > 0)
    assert ang[-1] - ang[0] < 2.0 * np.pi


def test_partial_span_makes_a_hook_not_a_loop():
    t = _target()
    half = loop_around(t, 40, 2.0, 6, 0.0, span=0.5)
    ang = np.unwrap(np.arctan2(half[:, 1] - 10.0, half[:, 0] - 10.0))
    assert ang[-1] - ang[0] == pytest.approx(np.pi, abs=1e-6)


def test_full_span_is_unchanged_by_the_span_knob():
    t = _target()
    assert loop_around(t, 40, 2.0, 6, 0.3) == pytest.approx(
        loop_around(t, 40, 2.0, 6, 0.3, span=1.0))
