"""Per-axis periodic boundaries in the Python generator.

The C searcher has always honoured `p_dims` per axis; the Python
builders wrapped unconditionally and ignored the `periodicity` config
value entirely, so a config asking for `"110"` silently got `"111"`.

An open axis omits its wrap-around bonds. The site set is unchanged --
only connectivity differs -- so the lattice grows a free surface there
and the sites on it lose coordination.

The default is fully periodic, which is exactly what the builders did
before, so nothing that omits `periodicity` moves.
"""

import random

import networkx as nx
import numpy as np
import pytest

from topon.config.schema import GeneratorConfig
from topon.topology.generator_python import PythonTopologyGenerator
from topon.topology.generator_python_diamond import create_diamond_lattice

DIMS = (4, 4, 4)


class _Cfg:
    def __init__(self, lattice_type, periodicity="111", mix=None, dims=DIMS):
        self.lattice_type = lattice_type
        self.lattice_size = dims
        self.max_functionality = 64
        self.degree_distribution = ""
        self.periodicity = periodicity
        self.mix_fractions = mix or {"SC": 1.0, "BCC": 0.0, "FCC": 0.0}
        self.mix_cutoff = 1.0


def build(lattice_type, periodicity="111", mix=None, dims=DIMS, seed=0):
    random.seed(seed)
    gen = PythonTopologyGenerator(_Cfg(lattice_type, periodicity, mix, dims))
    return gen._create_lattice(dims, lattice_type)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "value,expected",
    [
        ("111", (True, True, True)),
        ("110", (True, True, False)),
        ("011", (False, True, True)),
        ("000", (False, False, False)),
        (None, (True, True, True)),
        (True, (True, True, True)),
        (False, (False, False, False)),
        ((1, 0, 1), (True, False, True)),
        ([True, True, False], (True, True, False)),
    ],
)
def test_periodicity_parsing(value, expected):
    assert PythonTopologyGenerator._parse_periodicity(value) == expected


@pytest.mark.parametrize("bad", ["11", "1111", "abc", 5, (1, 0)])
def test_unparseable_periodicity_falls_back_to_periodic(bad):
    """Fall back loudly to fully periodic rather than to an open box.

    Guessing "open" on a malformed value would silently introduce free
    surfaces, which changes the physics; guessing "periodic" reproduces
    the behaviour every existing config already relies on.
    """
    assert PythonTopologyGenerator._parse_periodicity(bad) == (True, True, True)


# ---------------------------------------------------------------------------
# The default must not move anything
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("lattice_type", ["SC", "BCC", "FCC"])
def test_default_is_fully_periodic_and_unchanged(lattice_type):
    """A config with no `periodicity` must build exactly what it used to."""
    class _NoPeriodicity:
        lattice_size = DIMS
        max_functionality = 64
        degree_distribution = ""
        mix_fractions = {"SC": 1.0, "BCC": 0.0, "FCC": 0.0}
        mix_cutoff = 1.0

    cfg = _NoPeriodicity()
    cfg.lattice_type = lattice_type
    silent = PythonTopologyGenerator(cfg)._create_lattice(DIMS, lattice_type)
    explicit = build(lattice_type, "111")

    assert {frozenset(e) for e in silent.edges()} == {frozenset(e) for e in explicit.edges()}
    assert silent.number_of_nodes() == explicit.number_of_nodes()


# ---------------------------------------------------------------------------
# Effect on connectivity
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "lattice_type,expected",
    [
        # 4x4x4. Closing an axis removes one N^2 layer of wrap bonds for
        # SC; the centred lattices lose more because several offsets
        # cross the same face.
        ("SC", {"111": 192, "110": 176, "100": 160, "000": 144}),
        ("BCC", {"111": 512, "110": 448, "100": 392, "000": 343}),
        ("FCC", {"111": 1536, "110": 1408, "100": 1288, "000": 1176}),
    ],
)
def test_opening_axes_removes_wrap_bonds(lattice_type, expected):
    for periodicity, n_edges in expected.items():
        g = build(lattice_type, periodicity)
        assert g.number_of_edges() == n_edges, (
            f"{lattice_type} periodicity={periodicity}"
        )
        # Sites never change; only connectivity does.
        assert g.number_of_nodes() == build(lattice_type, "111").number_of_nodes()


def test_sc_open_edge_count_matches_the_closed_form():
    """3*N^2*(N-1) bonds with every axis open, against 3*N^3 closed."""
    for n in (3, 4, 5):
        dims = (n, n, n)
        assert build("SC", "000", dims=dims).number_of_edges() == 3 * n * n * (n - 1)
        assert build("SC", "111", dims=dims).number_of_edges() == 3 * n ** 3


def test_open_axis_creates_a_free_surface():
    """Sites on an opened face lose coordination; interior sites keep it."""
    g = build("SC", "110")           # z open
    pos = nx.get_node_attributes(g, "pos")
    surface = [n for n, p in pos.items() if p[2] in (0.0, 3.0)]
    interior = [n for n, p in pos.items() if p[2] in (1.0, 2.0)]

    assert {g.degree(n) for n in surface} == {5}
    assert {g.degree(n) for n in interior} == {6}


def test_no_edge_crosses_an_open_face():
    """Every bond must be short in the open direction, i.e. no wrap."""
    g = build("SC", "100")           # y and z open
    pos = {n: np.asarray(p, float) for n, p in nx.get_node_attributes(g, "pos").items()}
    for u, v in g.edges():
        d = np.abs(pos[u] - pos[v])
        assert d[1] <= 1.0 and d[2] <= 1.0, f"edge {u}-{v} wrapped an open axis"


def test_periodicity_is_recorded_on_the_graph():
    for periodicity, expected in (("111", (True, True, True)),
                                  ("101", (True, False, True))):
        assert build("SC", periodicity).graph["periodicity"] == expected


# ---------------------------------------------------------------------------
# Mixed lattices
# ---------------------------------------------------------------------------

def test_mixture_honours_periodicity():
    """The MIX neighbour search wraps only the periodic axes."""
    mix = {"SC": 0.2, "BCC": 0.4, "FCC": 0.4}
    closed = build("MIX", "111", mix, seed=3)
    open_z = build("MIX", "110", mix, seed=3)

    # Same draw, so the same sites; fewer bonds once z stops wrapping.
    assert closed.number_of_nodes() == open_z.number_of_nodes()
    assert open_z.number_of_edges() < closed.number_of_edges()

    box = np.asarray(open_z.graph["box"], float)
    pos = {n: np.asarray(p, float) for n, p in nx.get_node_attributes(open_z, "pos").items()}
    for u, v in open_z.edges():
        assert abs(pos[u][2] - pos[v][2]) <= open_z.graph["box"][2] / 2, (
            "a MIX bond wrapped the open z axis"
        )


# ---------------------------------------------------------------------------
# Diamond
# ---------------------------------------------------------------------------

def test_diamond_is_four_coordinated_when_periodic():
    g = create_diamond_lattice(3, 3, 3)
    assert g.number_of_nodes() == 8 * 27
    assert g.number_of_edges() == 2 * g.number_of_nodes()
    assert {d for _, d in g.degree()} == {4}


def test_diamond_open_axis_reduces_coordination():
    closed = create_diamond_lattice(4, 4, 4)
    open_z = create_diamond_lattice(4, 4, 4, (True, True, False))

    assert closed.number_of_nodes() == open_z.number_of_nodes()
    assert open_z.number_of_edges() < closed.number_of_edges()
    assert min(d for _, d in open_z.degree()) < 4


def test_diamond_reachable_through_the_main_generator():
    """`lattice_type: "Diamond"` must work on the config path too.

    The Diamond logic stays in its own module; `_create_lattice` only
    dispatches to it, so both generators accept the same config.
    """
    cfg = GeneratorConfig(lattice_type="Diamond", lattice_size="3x3x3",
                          max_functionality=4, degree_distribution="")
    g = PythonTopologyGenerator(cfg)._create_lattice((3, 3, 3), "Diamond")
    assert g.number_of_nodes() == 216
    assert {d for _, d in g.degree()} == {4}
    assert g.graph["box"] == (3.0, 3.0, 3.0)


def test_config_accepts_diamond():
    assert GeneratorConfig(lattice_type="Diamond").lattice_type == "Diamond"
