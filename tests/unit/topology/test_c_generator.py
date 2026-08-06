"""Shared-surface tests for the vendored C generator.

``topon/topology/csrc/generator.c`` is the standalone searcher for long
runs; ``topon/topology/generator_python.py`` is the quick in-process path
and the pipeline default. They are independent programs, not a library
and a wrapper. What they share, and what these tests police, is the
lattice construction and the ``.nodes`` / ``.edges`` format: a study that
switches between them must get the same kind of network out.

The C generator seeds with ``srand(time(NULL))`` and draws from
``rand()``, so it cannot reproduce a given Python draw. These tests
therefore check what must match regardless of the stream: lattice site
counts, edge counts, coordinates, the ``# BOX`` header, and that the C
can sculpt the configurations Python can.

That last one is load-bearing. The archive holds several near-identical
``generator_serial_debug11.c`` files, and the newest by timestamp is an
experimental variant that sculpts almost nothing. Picking a source by
date alone is not safe, so
``test_c_sculpts_the_configs_python_sculpts`` fails on it.

Every test compiles the C source on the fly and skips when no compiler is
available, so the suite still runs on machines without one.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from topon.config.schema import GeneratorConfig
from topon.topology.generator import format_lattice_arg
from topon.topology.loader import load_graph, read_box_header

ROOT = Path(__file__).resolve().parents[3]
CSRC = ROOT / "topon/topology/csrc/generator.c"

pytestmark = pytest.mark.skipif(
    shutil.which("gcc") is None and shutil.which("cc") is None,
    reason="no C compiler on PATH",
)


@pytest.fixture(scope="module")
def generator_exe(tmp_path_factory):
    """Compile the vendored generator once for the module."""
    if not CSRC.exists():
        pytest.skip(f"C source not found: {CSRC}")
    cc = shutil.which("gcc") or shutil.which("cc")
    out_dir = tmp_path_factory.mktemp("csrc")
    exe = out_dir / ("generator.exe" if sys.platform == "win32" else "generator")
    proc = subprocess.run(
        [cc, "-O2", "-o", str(exe), str(CSRC), "-lm"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, (
        f"C generator failed to compile:\n{proc.stdout}\n{proc.stderr}"
    )
    return exe


def run_c(exe, workdir, lattice_arg, dims="4x4x4", max_func=64,
          trials=20, degree_dist="0:0,1:0", periodicity="111"):
    """Run the C generator and return (nodes_path, edges_path) or None.

    ``max_func`` defaults above any degree these lattices reach (a mixture
    tops out near 20), so sculpting has nothing to prune and returns the
    base lattice directly. That is deliberate: these tests are about what
    the two generators *build*, and the C sculptor seeds from the clock,
    so leaving pruning in the loop would make them flaky for reasons that
    have nothing to do with parity.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [str(exe), dims, periodicity, str(max_func), str(trials), "1",
         degree_dist, "0", lattice_arg],
        cwd=workdir, capture_output=True, text=True, timeout=600,
    )
    nodes = sorted((workdir / "output").glob("*.nodes")) if (workdir / "output").exists() else []
    edges = sorted((workdir / "output").glob("*.edges")) if (workdir / "output").exists() else []
    if not nodes or not edges:
        return None, proc
    return (nodes[0], edges[0]), proc


# ---------------------------------------------------------------------------
# The .nodes format contract
# ---------------------------------------------------------------------------

def test_c_emits_a_box_header_the_python_loader_reads(generator_exe, tmp_path):
    """The C generator has to record the cell, same spelling as Python.

    Without it a BCC/FCC .nodes file reloads through the positional
    estimate, which overshoots the cell by half and misplaces a third of
    the edges.
    """
    paths, proc = run_c(generator_exe, tmp_path / "sc", "SC", max_func=6)
    assert paths is not None, f"C generator produced nothing:\n{proc.stdout[-1500:]}"
    nodes_p, edges_p = paths

    assert read_box_header(nodes_p) == (4.0, 4.0, 4.0)
    assert nodes_p.read_text().splitlines()[0] == "# BOX 4 4 4"

    G, dims = load_graph(nodes_path=nodes_p, edges_path=edges_p)
    assert np.allclose(dims, [4.0, 4.0, 4.0])
    assert G.graph["box"] == (4.0, 4.0, 4.0)


def test_c_sc_output_matches_the_python_sc_lattice(generator_exe, tmp_path):
    """Same sites and the same count of them, for an unsculpted lattice."""
    paths, proc = run_c(generator_exe, tmp_path / "sc", "SC", max_func=6)
    assert paths is not None, f"C generator produced nothing:\n{proc.stdout[-1500:]}"
    G, _ = load_graph(nodes_path=paths[0], edges_path=paths[1])

    assert G.number_of_nodes() == 64
    assert G.number_of_edges() == 192
    positions = {tuple(d["pos"]) for _, d in G.nodes(data=True)}
    expected = {(float(x), float(y), float(z))
                for x in range(4) for y in range(4) for z in range(4)}
    assert positions == expected


# ---------------------------------------------------------------------------
# Diamond and periodicity: exact parity, both sides deterministic
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "lattice_type,periodicity",
    [(lat, per)
     for lat in ("SC", "BCC", "FCC", "Diamond")
     for per in ("111", "110", "100", "000")],
)
def test_c_and_python_build_identical_lattices(generator_exe, tmp_path,
                                               lattice_type, periodicity):
    """Same sites AND the same edge set, for every lattice and boundary.

    Lattice construction is deterministic on both sides (only MIX draws
    randomly), so this is an exact comparison rather than a statistical
    one. It covers the two features that used to exist on one side only:
    Diamond was Python-only, per-axis periodicity C-only.
    """
    import random

    from topon.topology.generator_python import PythonTopologyGenerator
    from topon.topology.generator_python_diamond import create_diamond_lattice

    paths, proc = run_c(generator_exe, tmp_path / f"{lattice_type}_{periodicity}",
                        lattice_type, periodicity=periodicity, trials=5,
                        degree_dist="")
    assert paths is not None, (
        f"C produced nothing for {lattice_type} periodicity={periodicity}:\n"
        f"{proc.stdout[-800:]}{proc.stderr[-800:]}"
    )
    G_c, _ = load_graph(nodes_path=paths[0], edges_path=paths[1])

    axes = tuple(c == "1" for c in periodicity)
    if lattice_type == "Diamond":
        G_py = create_diamond_lattice(4, 4, 4, axes)
    else:
        class _Cfg:
            lattice_size = (4, 4, 4)
            max_functionality = 64
            degree_distribution = ""
            mix_fractions = {"SC": 1.0, "BCC": 0.0, "FCC": 0.0}
            mix_cutoff = 1.0
        cfg = _Cfg()
        cfg.lattice_type = lattice_type
        cfg.periodicity = periodicity
        random.seed(0)
        G_py = PythonTopologyGenerator(cfg)._create_lattice((4, 4, 4), lattice_type)

    assert G_c.number_of_nodes() == G_py.number_of_nodes()
    assert ({frozenset(e) for e in G_c.edges()}
            == {frozenset(e) for e in G_py.edges()}), (
        f"{lattice_type} periodicity={periodicity}: edge sets differ "
        f"(C {G_c.number_of_edges()}, Python {G_py.number_of_edges()})"
    )


def test_c_diamond_is_four_coordinated(generator_exe, tmp_path):
    """Diamond's whole point: 4-regular without any sculpting."""
    paths, proc = run_c(generator_exe, tmp_path / "dia", "Diamond",
                        dims="3x3x3", max_func=4, trials=5, degree_dist="")
    assert paths is not None, f"C produced nothing:\n{proc.stdout[-800:]}"
    G, dims = load_graph(nodes_path=paths[0], edges_path=paths[1])

    assert G.number_of_nodes() == 8 * 27
    assert {d for _, d in G.degree()} == {4}
    assert np.allclose(dims, [3.0, 3.0, 3.0])

    pos = {n: np.asarray(d["pos"], float) for n, d in G.nodes(data=True)}
    box = np.asarray(dims, float)
    lengths = np.asarray([
        np.linalg.norm((pos[u] - pos[v]) - box * np.round((pos[u] - pos[v]) / box))
        for u, v in G.edges()
    ])
    assert np.allclose(lengths, np.sqrt(3) / 4)


def test_c_accepts_both_diamond_spellings(generator_exe, tmp_path):
    for spelling in ("Diamond", "DIAMOND"):
        paths, proc = run_c(generator_exe, tmp_path / spelling, spelling,
                            dims="2x2x2", max_func=4, trials=3, degree_dist="")
        assert paths is not None, f"{spelling} rejected:\n{proc.stderr[-400:]}"


def test_c_open_axis_drops_only_wrap_bonds(generator_exe, tmp_path):
    """Opening an axis leaves the sites alone and removes wrap bonds."""
    closed, _ = run_c(generator_exe, tmp_path / "closed", "SC",
                      max_func=6, trials=5, degree_dist="")
    open_z, _ = run_c(generator_exe, tmp_path / "open", "SC",
                      max_func=6, trials=5, degree_dist="", periodicity="110")
    assert closed is not None and open_z is not None

    Gc, _ = load_graph(nodes_path=closed[0], edges_path=closed[1])
    Go, _ = load_graph(nodes_path=open_z[0], edges_path=open_z[1])

    assert Gc.number_of_nodes() == Go.number_of_nodes() == 64
    assert Gc.number_of_edges() == 192
    assert Go.number_of_edges() == 176        # one 4x4 layer of z-wraps gone


# ---------------------------------------------------------------------------
# MIX parity
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "lattice_arg,n_nodes,n_edges",
    [
        # The pure corners are deterministic, so exact counts apply and
        # must equal what the Python builder produces for the same mix.
        ("MIX:1,0,0", 64, 192),      # == SC
        ("MIX:0,1,0", 128, 896),     # BCC sites, degree 14 under the 1.0 cutoff
        ("MIX:0,0,1", 256, 2304),    # FCC sites, degree 18
    ],
)
def test_c_mix_pure_corners_match_python(generator_exe, tmp_path,
                                         lattice_arg, n_nodes, n_edges):
    """C and Python must agree exactly where the mixture is deterministic."""
    import random

    from topon.topology.generator_python import PythonTopologyGenerator

    paths, proc = run_c(generator_exe, tmp_path / lattice_arg.replace(":", "_"),
                        lattice_arg)
    assert paths is not None, f"C generator produced nothing:\n{proc.stdout[-1500:]}"
    G_c, _ = load_graph(nodes_path=paths[0], edges_path=paths[1])

    frac = dict(zip(("SC", "BCC", "FCC"),
                    (float(v) for v in lattice_arg.split(":")[1].split(","))))

    class _Cfg:
        lattice_type = "MIX"
        lattice_size = (4, 4, 4)
        max_functionality = 64
        degree_distribution = ""
        mix_fractions = frac
        mix_cutoff = 1.0

    random.seed(0)
    G_py = PythonTopologyGenerator(_Cfg())._create_lattice((4, 4, 4), "MIX")

    assert G_c.number_of_nodes() == n_nodes
    assert G_c.number_of_edges() == n_edges
    assert G_py.number_of_nodes() == n_nodes
    assert G_py.number_of_edges() == n_edges

    py_pos = {tuple(d["pos"]) for _, d in G_py.nodes(data=True)}
    c_pos = {tuple(d["pos"]) for _, d in G_c.nodes(data=True)}
    assert c_pos == py_pos


def test_c_mixture_site_count_follows_the_same_formula(generator_exe, tmp_path):
    """N*(1 + f_bcc + 3*f_fcc), the same expectation the Python builder has.

    Random draws differ between the two RNGs, so this is a tolerance
    check on the mean rather than an equality.
    """
    frac = (0.2, 0.4, 0.4)
    expected = 64 * (1 + frac[1] + 3 * frac[2])

    counts = []
    for i in range(5):
        paths, proc = run_c(generator_exe, tmp_path / f"mix{i}",
                            "MIX:0.2,0.4,0.4")
        assert paths is not None, f"C generator produced nothing:\n{proc.stdout[-1500:]}"
        G, _ = load_graph(nodes_path=paths[0], edges_path=paths[1])
        counts.append(G.number_of_nodes())

    assert abs(np.mean(counts) - expected) / expected < 0.10


@pytest.mark.parametrize(
    "dims,max_func,degree_dist",
    [
        ("4x4x4", 4, "0:0,1:0"),
        ("4x4x4", 6, "0:0,1:0"),
        ("5x5x5", 4, "0:0,1:0"),
        ("5x5x5", 3, "0:0,1:0"),
        ("4x4x4", 4, "0:2,1:0"),
        ("5x5x5", 4, "e:200"),
    ],
)
def test_c_sculpts_the_configs_python_sculpts(generator_exe, tmp_path,
                                              dims, max_func, degree_dist):
    """The C searcher must succeed wherever the Python one does.

    This is the check that catches a wrong source being vendored. The
    archive holds a later variant (md5 83d7f9d3, 2026-02-27) which swaps
    the per-degree count check in is_move_safe for a cumulative one; it
    fails five of these six, every case where max_func sits below the
    lattice coordination. It was vendored first on the strength of being
    newest and the mistake was invisible until someone ran a real
    sculpting job, because the other tests all use a max_func high enough
    that nothing needs pruning.

    Kept to small lattices and a low trial count so it stays fast; these
    all converge within a few trials on the correct source.
    """
    paths, proc = run_c(generator_exe, tmp_path / "sculpt", "SC",
                        dims=dims, max_func=max_func, trials=500,
                        degree_dist=degree_dist)
    assert paths is not None, (
        f"C failed to sculpt {dims} max_func={max_func} "
        f"degree_distribution={degree_dist!r}, which the Python generator "
        f"handles. Check which source is vendored.\n{proc.stdout[-800:]}"
    )
    G, _ = load_graph(nodes_path=paths[0], edges_path=paths[1])
    assert max(d for _, d in G.degree()) <= max_func


def test_c_survives_the_densest_mixture(generator_exe, tmp_path):
    """Guard the site-buffer allocation against a high-occupancy draw.

    ``create_mixed_lattice`` sizes its scratch array from the worst case,
    5 sites per cell (corner + body + 3 faces). Sizing it at 4, the FCC
    site count, looks safe on the mean but overruns: at f_bcc=0.1,
    f_fcc=0.9 roughly 0.2% of 4x4x4 draws exceed 4*N sites, which is a
    heap overflow rather than a wrong answer. Repeated because the failure
    is probabilistic and the C generator seeds from the clock.
    """
    for i in range(12):
        paths, proc = run_c(generator_exe, tmp_path / f"dense{i}",
                            "MIX:0.0,0.1,0.9", trials=5)
        assert paths is not None, (
            f"run {i} produced nothing (exit {proc.returncode}):\n"
            f"{proc.stdout[-1000:]}\n{proc.stderr[-1000:]}"
        )
        G, _ = load_graph(nodes_path=paths[0], edges_path=paths[1])
        # 1 corner + 0.1 body + 3*0.9 faces = 3.8 per cell on average.
        assert 64 <= G.number_of_nodes() <= 5 * 64


def test_c_mix_edges_respect_the_cutoff(generator_exe, tmp_path):
    """Every edge within the cutoff, and the same shells Python produces."""
    paths, proc = run_c(generator_exe, tmp_path / "mixcut", "MIX:0.2,0.4,0.4,1.0")
    assert paths is not None, f"C generator produced nothing:\n{proc.stdout[-1500:]}"
    G, dims = load_graph(nodes_path=paths[0], edges_path=paths[1])

    pos = {n: np.asarray(d["pos"], float) for n, d in G.nodes(data=True)}
    lengths = np.asarray([
        np.linalg.norm((pos[u] - pos[v]) - dims * np.round((pos[u] - pos[v]) / dims))
        for u, v in G.edges()
    ])
    assert lengths.max() <= 1.0 + 1e-9
    # The same four shells the Python builder produces: the 0.5 body-to-face
    # contact, FCC's 0.707, BCC's 0.866, and SC's 1.0. A subset is allowed
    # because a sparse draw may not realise every one of them.
    assert set(np.round(lengths, 4)) <= {0.5, 0.7071, 0.866, 1.0}
    assert 1.0 in set(np.round(lengths, 4)), "corner-corner shell missing"


# ---------------------------------------------------------------------------
# Argument formatting and rejection
# ---------------------------------------------------------------------------

def test_format_lattice_arg_passes_pure_types_through():
    for lt in ("SC", "BCC", "FCC"):
        assert format_lattice_arg(GeneratorConfig(lattice_type=lt)) == lt


def test_format_lattice_arg_encodes_mix():
    cfg = GeneratorConfig(lattice_type="MIX",
                          mix_fractions={"SC": 0.2, "BCC": 0.4, "FCC": 0.4})
    assert format_lattice_arg(cfg) == "MIX:0.2,0.4,0.4,1"


@pytest.mark.parametrize(
    "bad_arg,expect",
    [
        ("MIX:0.5,0.2,0.0", "sum to 1"),
        ("MIX:0.5,0.5", "three fractions"),
        ("HCP", "Invalid lattice type"),
        ("MIX:0.2,0.4,0.4,0", "cutoff must be positive"),
        # A prefix match on "MIX" would swallow these and quietly build a
        # pure-SC lattice, so a typo would cost a whole study.
        ("MIXED", "Invalid lattice type"),
        ("MIXTURE:0.2,0.4,0.4", "Invalid lattice type"),
    ],
)
def test_c_rejects_bad_lattice_args(generator_exe, tmp_path, bad_arg, expect):
    """Bad input must fail loudly at startup, not produce a silent lattice."""
    paths, proc = run_c(generator_exe, tmp_path / "bad", bad_arg, trials=1)
    assert paths is None
    assert proc.returncode != 0
    assert expect in proc.stderr


@pytest.mark.parametrize("degree_dist", ["0:0,1:0,7:5", "0:0,1:0,6:100"])
def test_c_rejects_degree_targets_above_max_func(generator_exe, tmp_path,
                                                 degree_dist):
    """A target no node can reach must fail, not silently pass.

    Sculpting enforces ``max_func``, so a node can never finish above it.
    The completion check only scans up to ``max_func``, so before this
    guard such a target was parsed, stored, then never looked at: the run
    printed "SUCCESS: Target distribution met!" over a network with none
    of the requested nodes. Both cases below asked for degrees above
    ``max_func=4`` and got a network back that had no such node.
    """
    paths, proc = run_c(generator_exe, tmp_path / "unreach", "SC",
                        dims="5x5x5", max_func=4, trials=20,
                        degree_dist=degree_dist)
    assert paths is None, "unreachable target produced a network"
    assert proc.returncode != 0
    assert "unreachable" in proc.stderr


@pytest.mark.parametrize("degree_dist", ["0:0,1:0", "0:0,1:0,4:20", "0:0,1:0,7:0"])
def test_c_still_accepts_reachable_targets(generator_exe, tmp_path, degree_dist):
    """The guard must not catch legitimate requests.

    ``7:0`` is the interesting one: it forbids degree-7 nodes rather than
    demanding them, so it is satisfiable even though 7 exceeds max_func.
    """
    paths, proc = run_c(generator_exe, tmp_path / f"ok{abs(hash(degree_dist))}",
                        "SC", dims="5x5x5", max_func=4, trials=500,
                        degree_dist=degree_dist)
    assert paths is not None, (
        f"reachable target {degree_dist!r} was rejected:\n{proc.stderr[-500:]}"
    )


def test_c_runs_started_together_are_independent(generator_exe, tmp_path):
    """Back-to-back runs must not produce the same network.

    `srand(time(NULL))` alone advances once a second, so a script looping
    this executable to collect N networks got byte-identical output from
    every run that started within the same second. Verified before the
    fix: three back-to-back MIX runs produced the same file. The pid is
    now mixed into the seed.

    MIX is the probe because its site draw makes any shared stream
    obvious in the node count; the pure lattices would look identical
    regardless.
    """
    counts, digests = [], []
    for i in range(4):
        paths, proc = run_c(generator_exe, tmp_path / f"indep{i}",
                            "MIX:0.2,0.4,0.4", trials=3, degree_dist="")
        assert paths is not None, f"run {i} produced nothing:\n{proc.stdout[-600:]}"
        counts.append(len(
            [ln for ln in paths[0].read_text().splitlines()
             if ln and not ln.startswith("#")]
        ))
        digests.append(paths[0].read_text())

    assert len(set(digests)) > 1, (
        f"all {len(digests)} runs produced identical output (node counts "
        f"{counts}); the seed is not varying between runs"
    )


def test_c_seed_env_makes_runs_reproducible(generator_exe, tmp_path):
    """TOPON_SEED pins the stream, the way seeding `random` does in Python."""
    outputs = []
    for i in range(2):
        workdir = tmp_path / f"seeded{i}"
        workdir.mkdir(parents=True, exist_ok=True)
        proc = subprocess.run(
            [str(generator_exe), "4x4x4", "111", "64", "3", "1", "", "0",
             "MIX:0.2,0.4,0.4"],
            cwd=workdir, capture_output=True, text=True, timeout=600,
            env={**os.environ, "TOPON_SEED": "42"},
        )
        nodes = sorted((workdir / "output").glob("*.nodes"))
        assert nodes, f"seeded run {i} produced nothing:\n{proc.stdout[-600:]}"
        outputs.append(nodes[0].read_text())

    assert outputs[0] == outputs[1], "TOPON_SEED did not make the run reproducible"


def test_c_mix_negative_fraction_rejected(generator_exe, tmp_path):
    paths, proc = run_c(generator_exe, tmp_path / "neg", "MIX:1.5,-0.5,0.0", trials=1)
    assert paths is None
    assert "non-negative" in proc.stderr
