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
