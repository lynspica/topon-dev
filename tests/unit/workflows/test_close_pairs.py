"""Which chain pairs can carry an entanglement.

Restricting the measurement to chains that come into contact is exact, not a
sample: two chains that never come close cannot be entangled. It is also what
makes a system-wide density measurable at all, since the whole-system Z1+
export is refused once a chain is longer than the periodic cell and has a size
ceiling besides.
"""
import numpy as np
import pytest

from tests.workflows.entangle_density import close_pairs

BOX = np.array([20.0, 20.0, 20.0])


def line(start, direction, n=40, step=0.5):
    t = np.arange(n)[:, None] * step
    return np.asarray(start, float) + t * np.asarray(direction, float)


def test_touching_chains_are_found():
    paths = {0: line([5, 5, 2], [0, 0, 1]),
             1: line([5, 6, 2], [0, 0, 1])}
    assert close_pairs(paths, BOX, cutoff=3.0) == [(0, 1)]


def test_distant_chains_are_not():
    paths = {0: line([2, 2, 2], [0, 0, 1]),
             1: line([15, 15, 2], [0, 0, 1])}
    assert close_pairs(paths, BOX, cutoff=3.0) == []


def test_the_cutoff_is_what_decides():
    paths = {0: line([5, 5, 2], [0, 0, 1]),
             1: line([5, 9, 2], [0, 0, 1])}
    assert close_pairs(paths, BOX, cutoff=3.0) == []
    assert close_pairs(paths, BOX, cutoff=5.0) == [(0, 1)]


def test_contact_across_the_periodic_boundary_counts():
    """A pair touching only through the wall is still a pair.

    Missing these would under-report the density, and by construction they are
    the ones a naive distance would never see.
    """
    paths = {0: line([0.5, 5, 2], [0, 0, 1]),
             1: line([19.5, 5, 2], [0, 0, 1])}
    assert close_pairs(paths, BOX, cutoff=2.0) == [(0, 1)]


def test_a_chain_is_never_paired_with_itself():
    paths = {0: line([5, 5, 2], [0, 0, 1], n=60)}
    assert close_pairs(paths, BOX, cutoff=5.0) == []


def test_every_pair_appears_once_and_in_order():
    paths = {2: line([5, 5, 2], [0, 0, 1]),
             0: line([5, 6, 2], [0, 0, 1]),
             1: line([5, 7, 2], [0, 0, 1])}
    got = close_pairs(paths, BOX, cutoff=3.0)
    assert got == sorted(set(got))
    assert all(a < b for a, b in got)


def test_scales_to_a_realistic_chain_count():
    """The pairwise form is quadratic twice over and is minutes here."""
    rng = np.random.default_rng(0)
    paths = {k: line(rng.uniform(0, 20, 3), [0, 0, 1], n=20)
             for k in range(400)}
    got = close_pairs(paths, BOX, cutoff=1.5)
    assert isinstance(got, list)
    # Every reported pair really is in contact.
    for a, b in got[:20]:
        d = paths[a][:, None, :] - paths[b][None, :, :]
        d -= BOX * np.round(d / BOX)
        assert np.linalg.norm(d, axis=2).min() <= 1.5 + 1e-9


# ---------------------------------------------------------------------------
# The watch set
# ---------------------------------------------------------------------------

def test_every_designed_pair_is_watched_and_counts_for_itself():
    """The designed pairs are where the increment is meant to appear.

    Sampling them would put noise exactly where the signal is.
    """
    from tests.workflows.entangle_density import watch_set
    paths = {k: line([2 + k, 5, 2], [0, 0, 1]) for k in range(8)}
    designed = [(0, 1), (2, 3)]
    pairs, scale = watch_set(paths, BOX, designed, cutoff=3.0, sample=2,
                             rng=np.random.default_rng(0))
    for q in designed:
        assert q in pairs
        assert scale[q] == 1.0


def test_the_sample_stands_for_the_pairs_left_out():
    from tests.workflows.entangle_density import watch_set
    paths = {k: line([1 + 0.7 * k, 5, 2], [0, 0, 1]) for k in range(20)}
    pairs, scale = watch_set(paths, BOX, [], cutoff=3.0, sample=5,
                             rng=np.random.default_rng(1))
    assert len(pairs) == 5
    rest = len(close_pairs(paths, BOX, 3.0))
    assert sum(scale.values()) == pytest.approx(rest)


def test_a_small_system_is_measured_in_full():
    from tests.workflows.entangle_density import watch_set
    paths = {0: line([5, 5, 2], [0, 0, 1]), 1: line([5, 6, 2], [0, 0, 1])}
    pairs, scale = watch_set(paths, BOX, [], cutoff=3.0, sample=100,
                             rng=np.random.default_rng(2))
    assert pairs == [(0, 1)]
    assert scale[(0, 1)] == 1.0


def test_a_designed_pair_out_of_contact_is_still_watched():
    """It is out of contact now; routing is about to put it in contact."""
    from tests.workflows.entangle_density import watch_set
    paths = {0: line([2, 2, 2], [0, 0, 1]), 1: line([15, 15, 2], [0, 0, 1])}
    pairs, _s = watch_set(paths, BOX, [(0, 1)], cutoff=3.0, sample=10,
                          rng=np.random.default_rng(3))
    assert (0, 1) in pairs
