"""Parity tests between the vendored C generator and the Python one.

``topon/topology/csrc/generator.c`` and
``topon/topology/generator_python.py`` implement the same algorithm. The
Python one is the pipeline default; the C one is opt-in via
``generator.exe_path`` and used for large runs on HPC. They have to agree
on what they build and on the file format they emit, or a study that
switches paths silently changes meaning.

The C generator seeds with ``srand(time(NULL))`` and draws from
``rand()``, so it cannot reproduce a given Python draw. These tests
therefore check the things that must match regardless of the stream:
lattice site counts, edge counts, coordinates, and the ``.nodes`` header.

Every test compiles the C source on the fly and skips when no compiler is
available, so the suite still runs on machines without one.
"""

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
          trials=20, degree_dist="0:0,1:0"):
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
        [str(exe), dims, "111", str(max_func), str(trials), "1",
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
    ],
)
def test_c_rejects_bad_lattice_args(generator_exe, tmp_path, bad_arg, expect):
    """Bad input must fail loudly at startup, not produce a silent lattice."""
    paths, proc = run_c(generator_exe, tmp_path / "bad", bad_arg, trials=1)
    assert paths is None
    assert proc.returncode != 0
    assert expect in proc.stderr


def test_c_mix_negative_fraction_rejected(generator_exe, tmp_path):
    paths, proc = run_c(generator_exe, tmp_path / "neg", "MIX:1.5,-0.5,0.0", trials=1)
    assert paths is None
    assert "non-negative" in proc.stderr
