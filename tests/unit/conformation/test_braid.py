"""Waypoint braid: the winding count must be exactly what was asked for.

The acceptance measure throughout is ``far_closed_linking`` -- close each open
path through a detour that runs far from the system, then take the Gauss
linking number of the two resulting loops. That is an integer topological
invariant meaning "does this path wind around its partner". A pair that pulls
apart freely measures 0.

The far detour matters. Closing with each chain's own straight chord is the
more obvious choice, and it is what ``chord_closed_linking`` does, but the
chords are then part of the loops: once a braid sweeps near the partner's
chord the integral diverges. Measured on a single-winding braid, chord
closure reads 1.38 at gap 3, 3.47 at gap 2 and 6.21 at gap 1, while the far
closure reads 1.002 throughout and the real chain-chain clearance never
changes.

These tests use bare chord pairs, no lattice and no periodic images, so a
failure here is a failure of the geometry and nothing else.
"""

import numpy as np
import pytest

from topon.conformation.entanglement import (
    BraidShape,
    braid_pair,
    braid_path,
    far_closed_linking,
    closest_approach,
    make_contact,
    min_separation,
    plan_braid,
)

# A pair long enough that the axial budget never binds, so these tests
# measure the braid rather than the room available to it.
LONG = dict(a0=[0, 0, 0], a1=[40, 0, 0], b0=[0, 3, 0], b1=[40, 3, 0])


def _pair(e, n_beads=400, shape=None, **chords):
    c = {**LONG, **chords}
    return braid_pair(c["a0"], c["a1"], c["b0"], c["b1"], e,
                      shape=shape, n_beads=n_beads)


# ---------------------------------------------------------------------------
# The headline property
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("e", [1, 2, 3, 4, 5])
def test_winding_count_is_exactly_as_requested(e):
    """e turns in, e turns out. This is the whole point of the construction."""
    pa, pb, e_max = _pair(e)
    assert e <= e_max, "test pair should have room; adjust LONG, not the code"
    assert round(far_closed_linking(pa, pb)) == e


def test_no_braid_gives_zero_linking():
    """Null control. Straight chains must measure exactly 0.

    Without this, a nonzero reading elsewhere proves nothing -- two nearby
    chains wind around each other by chance often enough that the raw linked
    fraction is meaningless on its own.
    """
    a = np.linspace([0, 0, 0], [40, 0, 0], 400)
    b = np.linspace([0, 3, 0], [40, 3, 0], 400)
    assert abs(far_closed_linking(a, b)) < 1e-6


def test_partners_stay_apart():
    """The braid must not press the chains together to achieve its winding."""
    for e in (1, 2, 3):
        pa, pb, _ = _pair(e)
        assert min_separation(pa, pb) > 0.8


def test_clearance_does_not_decay_with_winding_count():
    """More turns must cost axial room, not clearance.

    Clearance is set by the pitch; if it fell as e rose, the braid would be
    silently compressing instead of extending.
    """
    seps = [min_separation(*_pair(e)[:2]) for e in (1, 2, 3, 4)]
    assert max(seps) - min(seps) < 0.1


# ---------------------------------------------------------------------------
# The three cases the legacy mid-strand kink cannot handle
# ---------------------------------------------------------------------------

def test_off_midpoint_contact_still_hooks():
    """Legacy kink peaks at t=0.5 always; hook rate falls to 18% past offset
    0.25. The shared frame places the braid at the real contact instead."""
    pa, pb, e_max = _pair(2, b0=[10, 3, 0], b1=[50, 3, 0])
    assert e_max >= 2
    assert round(far_closed_linking(pa, pb)) == 2


@pytest.mark.parametrize("b0,b1", [
    ([20, 3, -20], [20, 3, 20]),        # 90 degrees
    ([20, 3, -18], [32, 3, 18]),        # oblique
])
def test_skew_chords_hook(b0, b1):
    """The legacy kink assumes the partner lies along a supplied orientation
    vector; a skew partner defeats it. The shared axis is the mean direction,
    so skew is ordinary."""
    pa, pb, e_max = _pair(2, b0=b0, b1=b1)
    assert e_max >= 2
    assert round(far_closed_linking(pa, pb)) == 2


@pytest.mark.parametrize("length", [20.0, 12.0, 8.0])
def test_unequal_partner_lengths_hook_and_stay_apart(length):
    """Each legacy chain scales its reach by its own length, so a long chain
    and a short one aim at different places and miss.

    This also guards the axial-keying bug specifically: when phase was keyed
    to each chain's own parameter, a 40-vs-8 pair still *linked* but closed
    to 0.05 separation, because the partners reached opposite phase at
    different axial positions. Assert both properties.
    """
    mid = 20.0
    pa, pb, e_max = _pair(
        2, b0=[mid - length / 2, 3, 0], b1=[mid + length / 2, 3, 0])
    if e_max >= 2:
        assert round(far_closed_linking(pa, pb)) == 2
    assert min_separation(pa, pb) > 0.5


# ---------------------------------------------------------------------------
# Budget reporting
# ---------------------------------------------------------------------------

def test_contact_at_a_junction_reports_no_room():
    """A pair meeting at a chord end has nowhere to braid, and says so.

    This is the simple-cubic nearest-neighbour case: every non-parallel pair
    meets at a junction, where a hook would be held by the crosslink rather
    than by topology. Reporting zero room is the honest answer, not a bug.
    """
    _, _, e_max = braid_pair([0, 0, 0], [40, 0, 0],
                             [38, 3, 0], [78, 3, 0], 2, n_beads=200)
    assert e_max == 0


def test_e_max_grows_with_chord_length():
    """Axial room is the binding budget, so a longer chord carries more."""
    rooms = []
    for L in (12.0, 24.0, 48.0):
        _, _, e_max = braid_pair([0, 0, 0], [L, 0, 0],
                                 [0, 3, 0], [L, 3, 0], 1, n_beads=200)
        rooms.append(e_max)
    assert rooms == sorted(rooms)
    assert rooms[-1] > rooms[0]


def test_over_budget_request_is_reported_not_hidden():
    """Asking beyond e_max compresses the braid; the caller must be able to
    see that rather than silently receive a squashed one."""
    _, _, e_max = braid_pair([0, 0, 0], [10, 0, 0],
                             [0, 3, 0], [10, 3, 0], 8, n_beads=200)
    assert e_max < 8


@pytest.mark.parametrize("gap", [0.8, 1.0, 1.5, 2.0, 3.0, 6.0])
def test_tight_gaps_still_give_the_right_winding(gap):
    """The braid must shrink to fit a narrow gap rather than over-reach.

    Left unclamped at gap 1.5 the braid swings past the partner's chord: the
    partners have effectively swapped sides, so nothing separates them during
    minimisation. Fitting the radius to the gap keeps each chain on its own
    side, and the winding is still exactly what was asked for.
    """
    pa, pb, _ = braid_pair([0, 0, 0], [90, 0, 0],
                           [0, gap, 0], [90, gap, 0], 1, n_beads=800)
    c = make_contact([0, 0, 0], [90, 0, 0], [0, gap, 0], [90, gap, 0])
    assert round(far_closed_linking(pa, pb, c)) == 1
    assert pa[:, 1].max() < gap, "chain A reached past B's chord"
    assert pb[:, 1].min() > 0.0, "chain B reached past A's chord"


def test_fit_to_gap_leaves_a_roomy_braid_alone():
    shape = BraidShape(n_radius=0.9)
    assert shape.fit_to_gap(10.0) is shape
    tight = shape.fit_to_gap(1.0)
    assert tight.n_radius < shape.n_radius
    assert tight.m_radius < shape.m_radius       # ellipse keeps its proportions
    assert tight.pitch == shape.pitch            # clearance comes from pitch


def test_chord_closure_diverges_where_far_closure_does_not():
    """Pin the limitation, so the weaker measure is not reached for by
    mistake. If chord closure ever becomes accurate at tight gaps this test
    fails and the docstring should be revisited."""
    from topon.conformation.entanglement import chord_closed_linking

    pa, pb, _ = braid_pair([0, 0, 0], [90, 0, 0],
                           [0, 1, 0], [90, 1, 0], 1, n_beads=800)
    c = make_contact([0, 0, 0], [90, 0, 0], [0, 1, 0], [90, 1, 0])
    assert round(far_closed_linking(pa, pb, c)) == 1
    assert abs(chord_closed_linking(pa, pb)) > 2.0


def test_span_matches_shape_arithmetic():
    shape = BraidShape(pitch=3.0, ramp=1.5)
    assert shape.span(2) == pytest.approx(2 * 3.0 + 2 * 1.5)
    assert shape.turns_within(shape.span(3) / 2.0) == 3


# ---------------------------------------------------------------------------
# Frame and contact
# ---------------------------------------------------------------------------

def test_contact_frame_is_orthonormal():
    c = make_contact([0, 0, 0], [40, 0, 0], [20, 3, -20], [20, 3, 20])
    for v in (c.axis, c.toward, c.across):
        assert np.linalg.norm(v) == pytest.approx(1.0)
    assert abs(c.axis @ c.toward) < 1e-9
    assert abs(c.axis @ c.across) < 1e-9
    assert abs(c.toward @ c.across) < 1e-9


def test_toward_points_from_a_to_b():
    """The frame's `toward` must aim at the partner, or the sides invert."""
    c = make_contact([0, 0, 0], [40, 0, 0], [0, 5, 0], [40, 5, 0])
    assert c.toward @ np.array([0.0, 1.0, 0.0]) > 0.9
    assert c.gap == pytest.approx(5.0)


def test_antiparallel_chords_give_a_real_axis():
    """When the chords run opposite ways the mean direction would cancel;
    the builder flips one first."""
    c = make_contact([0, 0, 0], [40, 0, 0], [40, 3, 0], [0, 3, 0])
    assert np.linalg.norm(c.axis) == pytest.approx(1.0)
    assert abs(c.axis @ np.array([1.0, 0.0, 0.0])) > 0.99


def test_closest_approach_is_clamped_to_the_segments():
    s, t = closest_approach([0, 0, 0], [10, 0, 0], [50, 3, 0], [60, 3, 0])
    assert 0.0 <= s <= 1.0 and 0.0 <= t <= 1.0
    assert s == pytest.approx(1.0)
    assert t == pytest.approx(0.0)


def test_explicit_contact_position_is_honoured():
    """Placing a contact deliberately is what allows several along one chain."""
    c = make_contact([0, 0, 0], [40, 0, 0], [0, 3, 0], [40, 3, 0],
                     s_a=0.25, s_b=0.25)
    assert c.origin[0] == pytest.approx(10.0)


# ---------------------------------------------------------------------------
# Path shape
# ---------------------------------------------------------------------------

def test_path_endpoints_are_the_junctions():
    """A braid must not move the crosslinks it hangs between."""
    pa, pb, _ = _pair(3)
    assert np.allclose(pa[0], [0, 0, 0], atol=1e-9)
    assert np.allclose(pa[-1], [40, 0, 0], atol=1e-9)
    assert np.allclose(pb[0], [0, 3, 0], atol=1e-9)
    assert np.allclose(pb[-1], [40, 3, 0], atol=1e-9)


def test_path_returns_to_the_chord_outside_the_braid():
    """The detour is local; the rest of the chain is untouched."""
    pa, _, _ = _pair(1)
    ts = np.linspace(0.0, 1.0, len(pa))
    chord = np.outer(ts, np.array([40.0, 0.0, 0.0]))
    off = np.linalg.norm(pa - chord, axis=1)
    assert off[0] < 1e-9 and off[-1] < 1e-9
    assert off.max() > 0.5                      # it did braid somewhere


def test_bead_count_is_respected():
    for n in (2, 17, 64):
        pa, pb, _ = _pair(1, n_beads=n)
        assert len(pa) == n and len(pb) == n


def test_zero_half_span_falls_back_to_the_chord():
    """Degenerate planning must give a straight chain, not a crash."""
    c = make_contact([0, 0, 0], [10, 0, 0], [0, 3, 0], [10, 3, 0])
    p = braid_path([0, 0, 0], [10, 0, 0], c, 2, -1, 50, half_span=0.0)
    assert np.allclose(p[0], [0, 0, 0])
    assert np.allclose(p[-1], [10, 0, 0])
    assert np.linalg.norm(p - np.linspace([0, 0, 0], [10, 0, 0], 50)) < 1e-9


def test_thin_across_axis_stays_thin():
    """The ellipse must not become circular.

    A circular braid swings a half-gap sideways, reaching into neighbouring
    strand corridors on a lattice and creating entanglements nobody asked
    for, at several times the contour cost and no clearance gain.
    """
    shape = BraidShape(n_radius=0.9, m_radius=1.5)
    pa, _, _ = _pair(2, shape=shape)
    assert np.abs(pa[:, 2]).max() <= shape.m_radius + 1e-6


def test_partners_sit_on_opposite_sides_of_the_axis():
    """Anti-phase is the mechanism, so check it directly.

    Measured geometrically rather than topologically: routing both partners
    on the same side does not give a clean "unlinked" reading to assert
    against, it gives two chains occupying the same space, and the linking
    number of a degenerate overlap is meaningless. The honest test is that
    the realised paths sit on opposite sides of the shared axis.
    """
    pa, pb, _ = _pair(2)
    c = make_contact(LONG["a0"], LONG["a1"], LONG["b0"], LONG["b1"])

    def side_of(path):
        rel = path - c.origin
        along = (rel @ c.axis)[:, None] * c.axis
        return (rel - along) @ c.toward          # signed, along `toward`

    # Take the extreme excursion of each: they must have opposite sign.
    ea, eb = side_of(pa), side_of(pb)
    assert ea.min() < -0.5, "A never swings away from the partner"
    assert eb.max() > 0.5, "B never swings away from the partner"
    # And at the axial midpoint they are on opposite sides of the axis.
    ka = int(np.argmin(np.abs((pa - c.origin) @ c.axis)))
    kb = int(np.argmin(np.abs((pb - c.origin) @ c.axis)))
    assert np.sign(ea[ka]) == -np.sign(eb[kb])
