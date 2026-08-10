"""Ranking entanglement candidates on a conformation rather than on crosslinks.

The property: a candidate whose two chains actually run alongside each other
outscores one whose crosslinks happen to be close, and a candidate whose
chains never meet scores zero and so leaves the draw.
"""

import numpy as np
import pytest

from topon.assignment.entanglements import compute_proximity_weights
from topon.conformation.paths import bridging_walk, straight


def _cand(a, b):
    """A candidate pair in the ((u, v, key), (u, v, key)) form."""
    return ((a[0], a[1], 0), (b[0], b[1], 0))


def test_a_pair_that_runs_alongside_outscores_one_that_crosses():
    """Two chains sharing a long stretch beat two that touch at a point."""
    n = 60
    together_a = straight((0, 0, 0), (30, 0, 0), n)
    together_b = straight((0, 1, 0), (30, 1, 0), n)          # parallel, 1 apart
    crossing_a = straight((0, 10, 0), (30, 10, 0), n)
    crossing_b = straight((15, 10, -15), (15, 10, 15), n)    # meets at a point

    paths = {frozenset((1, 2)): together_a, frozenset((3, 4)): together_b,
             frozenset((5, 6)): crossing_a, frozenset((7, 8)): crossing_b}
    w = compute_proximity_weights(
        [_cand((1, 2), (3, 4)), _cand((5, 6), (7, 8))], paths, cutoff=2.0)
    assert w[0] > w[1] > 0


def test_chains_that_never_meet_score_zero():
    """Zero removes the candidate from a weighted draw, which is the point.

    The pipeline ranks candidates on the distance between their crosslinks,
    and a pair can be close by that measure while its chains never come near
    each other anywhere. A kink placed there aims one chain at a partner that
    is not present.
    """
    a = straight((0, 0, 0), (30, 0, 0), 40)
    b = straight((0, 50, 0), (30, 50, 0), 40)
    paths = {frozenset((1, 2)): a, frozenset((3, 4)): b}
    w = compute_proximity_weights([_cand((1, 2), (3, 4))], paths, cutoff=2.0)
    assert w == [0.0]


def test_a_candidate_with_no_path_scores_zero_rather_than_raising():
    paths = {frozenset((1, 2)): straight((0, 0, 0), (10, 0, 0), 20)}
    w = compute_proximity_weights([_cand((1, 2), (3, 4))], paths, cutoff=2.0)
    assert w == [0.0]


def test_the_periodic_image_is_used_when_a_box_is_given():
    """A pair that is close across the boundary must count as close."""
    box = np.array([40.0, 40.0, 40.0])
    a = straight((0, 1.0, 0), (30, 1.0, 0), 40)
    b = straight((0, 39.0, 0), (30, 39.0, 0), 40)    # 2 apart through the face
    paths = {frozenset((1, 2)): a, frozenset((3, 4)): b}
    cand = [_cand((1, 2), (3, 4))]
    assert compute_proximity_weights(cand, paths, box=box, cutoff=3.0)[0] > 0
    assert compute_proximity_weights(cand, paths, box=None, cutoff=3.0) == [0.0]


def test_cutoff_selects_how_close_counts():
    a = straight((0, 0, 0), (30, 0, 0), 40)
    b = straight((0, 3.0, 0), (30, 3.0, 0), 40)
    paths = {frozenset((1, 2)): a, frozenset((3, 4)): b}
    cand = [_cand((1, 2), (3, 4))]
    assert compute_proximity_weights(cand, paths, cutoff=2.0) == [0.0]
    assert compute_proximity_weights(cand, paths, cutoff=4.0)[0] > 0


# ---------------------------------------------------------------------------
# The conformation the ranking needs
# ---------------------------------------------------------------------------

def test_bridging_walk_lands_on_its_far_junction():
    rng = np.random.default_rng(0)
    p = bridging_walk((0, 0, 0), (12, 0, 0), 40, bond=0.97, rng=rng)
    assert len(p) == 41
    assert np.allclose(p[0], [0, 0, 0])
    assert np.allclose(p[-1], [12, 0, 0])


def test_bridging_walk_keeps_its_bond_length():
    """Every bond, not just on average: the walk must not be rescaled."""
    rng = np.random.default_rng(0)
    for chord in (5.0, 20.0, 35.0):
        p = bridging_walk((0, 0, 0), (chord, 0, 0), 40, bond=0.97, rng=rng)
        d = np.linalg.norm(np.diff(p, axis=0), axis=1)
        assert d.max() <= 0.97 + 1e-9
        assert d[:-1].min() > 0.9        # the closing bond may be short


def test_bridging_walk_wanders_where_a_straight_line_does_not():
    """The reason this exists: a coiled chain leaves its chord.

    Whether two chains meet is decided by where the chains go, not by where
    their crosslinks sit, and a straight chain cannot express that.
    """
    rng = np.random.default_rng(0)
    a, b = np.array([0.0, 0, 0]), np.array([12.0, 0, 0])
    p = bridging_walk(a, b, 40, bond=0.97, rng=rng)
    chord = b - a
    t = ((p - a) @ chord) / (chord @ chord)
    off = np.linalg.norm((p - a) - t[:, None] * chord, axis=1)
    assert off.max() > 2.0


def test_a_fully_extended_walk_has_nowhere_to_wander():
    """Not a defect: there is one way to span a chord of n*bond."""
    rng = np.random.default_rng(0)
    p = bridging_walk((0, 0, 0), (40 * 0.97, 0, 0), 40, bond=0.97, rng=rng)
    assert np.abs(p[:, 1]).max() < 1e-6
    assert np.abs(p[:, 2]).max() < 1e-6
