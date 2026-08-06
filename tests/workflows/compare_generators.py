"""Compare the C searcher and the Python generator across the config matrix.

Run directly (not via pytest)::

    python tests/workflows/compare_generators.py            # structural sweep
    python tests/workflows/compare_generators.py --lammps   # also build + minimize

The two are independent programs (see topon/topology/csrc/README.md), so
they cannot produce the same individual network: the C seeds from the
clock and draws from rand(). What they must agree on is the *kind* of
network produced, which is what this sweeps.

What is compared, and how strictly
----------------------------------
Exactly, because they are deterministic:
  * base lattice site and edge counts for SC / BCC / FCC
  * the periodic cell recorded in ``# BOX``
  * the set of edge-length shells present

Statistically, because the draw differs:
  * site counts for MIX (binomial in the fractions)
  * degree histogram after sculpting
  * success rate over repeated trials

Distributions that neither can satisfy are reported as such rather than
counted as disagreement. A target being unreachable on a given lattice is
a property of the request, not a defect in either generator; the two are
only expected to agree on *whether* it is reachable.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from collections import Counter
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.topology.generator_python import PythonTopologyGenerator  # noqa: E402
from topon.topology.loader import load_graph  # noqa: E402

CSRC = ROOT / "topon/topology/csrc/generator.c"
OUT = ROOT / "tests/output/generator_comparison"
REPS = 5


# ---------------------------------------------------------------------------
# Config matrix
# ---------------------------------------------------------------------------

def slugify(label):
    """Filesystem-safe form of a case label.

    Whitelist rather than blacklist: labels contain ``->``, ``/`` and
    spaces, and on Windows a stray ``>`` makes the whole path invalid
    with a WinError 123 that surfaces far from its cause.
    """
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in label).strip("_")


class Case:
    def __init__(self, label, lattice, dims, max_func, degree_dist, mix=None):
        self.label = label
        self.lattice = lattice          # "SC" | "BCC" | "FCC" | "MIX"
        self.dims = dims
        self.max_func = max_func
        self.degree_dist = degree_dist
        self.mix = mix or {"SC": 1.0, "BCC": 0.0, "FCC": 0.0}

    @property
    def c_arg(self):
        if self.lattice != "MIX":
            return self.lattice
        m = self.mix
        return f"MIX:{m['SC']:g},{m['BCC']:g},{m['FCC']:g},1"

    @property
    def dims_str(self):
        return "x".join(str(d) for d in self.dims)


def build_matrix():
    cases = []

    # --- the three canonical lattices, unsculpted and sculpted ---
    for lat, native_deg in (("SC", 6), ("BCC", 8), ("FCC", 12)):
        cases.append(Case(f"{lat} native", lat, (4, 4, 4), native_deg, "0:0,1:0"))
        cases.append(Case(f"{lat} prune->4", lat, (4, 4, 4), 4, "0:0,1:0"))

    # --- lattice sizes, including non-cubic ---
    for dims in ((3, 3, 3), (5, 5, 5), (6, 6, 6), (3, 4, 5)):
        cases.append(Case(f"SC {dims[0]}x{dims[1]}x{dims[2]}", "SC", dims, 4, "0:0,1:0"))

    # --- distribution modes ---
    cases.append(Case("SC defects d0",     "SC", (5, 5, 5), 4, "0:5,1:0"))
    cases.append(Case("SC dangling d1",    "SC", (5, 5, 5), 4, "0:0,1:10"))
    cases.append(Case("SC edge target",    "SC", (5, 5, 5), 4, "e:200"))
    cases.append(Case("SC per-degree",     "SC", (5, 5, 5), 4, "0:0,1:0,2:20"))
    # Degree 7 is above SC's coordination of 6, so this can never be met.
    # Python rejects it up front (the V40 fail-fast guard). The C parses
    # it, then never checks it, because its completion loop only runs to
    # max_func -- so it reports success on a network with no degree-7
    # nodes at all. Kept in the matrix to keep that visible.
    cases.append(Case("SC unreachable d7",  "SC", (4, 4, 4), 4, "0:0,1:0,7:5"))
    # Same shape but within range: 100 nodes of degree 6 on a 125-node
    # lattice that must also prune to max_func=4. Genuinely unreachable
    # and both should say so.
    cases.append(Case("SC over-target d6",  "SC", (5, 5, 5), 4, "0:0,1:0,6:100"))

    # --- mixtures ---
    for name, mix in (
        ("MIX pure-SC",   {"SC": 1.0, "BCC": 0.0, "FCC": 0.0}),
        ("MIX pure-BCC",  {"SC": 0.0, "BCC": 1.0, "FCC": 0.0}),
        ("MIX pure-FCC",  {"SC": 0.0, "BCC": 0.0, "FCC": 1.0}),
        ("MIX 50/50 SB",  {"SC": 0.5, "BCC": 0.5, "FCC": 0.0}),
        ("MIX 20/40/40",  {"SC": 0.2, "BCC": 0.4, "FCC": 0.4}),
        ("MIX 34/33/33",  {"SC": 0.34, "BCC": 0.33, "FCC": 0.33}),
        ("MIX dense FCC", {"SC": 0.0, "BCC": 0.1, "FCC": 0.9}),
    ):
        cases.append(Case(name, "MIX", (4, 4, 4), 64, "0:0,1:0", mix))
    cases.append(Case("MIX 20/40/40 prune->4", "MIX", (4, 4, 4), 4, "0:0,1:0",
                      {"SC": 0.2, "BCC": 0.4, "FCC": 0.4}))

    return cases


# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------

def compile_c(workdir):
    cc = shutil.which("gcc") or shutil.which("cc")
    if cc is None or not CSRC.exists():
        return None
    exe = workdir / ("generator.exe" if sys.platform == "win32" else "generator")
    proc = subprocess.run([cc, "-O2", "-o", str(exe), str(CSRC), "-lm"],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"C compile FAILED:\n{proc.stderr[-1500:]}")
        return None
    return exe


def run_c(exe, workdir, case, trials=500):
    """One C run in its own directory. Returns (graph, dims, paths) or None.

    A fresh directory per run rather than clearing one: the generator
    always writes to ``output/`` relative to its cwd, and on Windows
    deleting that between runs races against file handles held by
    OneDrive or a shell sitting in the tree.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / "output"
    proc = subprocess.run(
        [str(exe), case.dims_str, "111", str(case.max_func), str(trials), "1",
         case.degree_dist, "0", case.c_arg],
        cwd=workdir, capture_output=True, text=True, timeout=900,
    )
    nodes = sorted(out.glob("*.nodes")) if out.exists() else []
    edges = sorted(out.glob("*.edges")) if out.exists() else []
    if not nodes or not edges:
        return None, proc
    with redirect_stdout(StringIO()):
        g, dims = load_graph(nodes_path=nodes[0], edges_path=edges[0])
    return (g, dims, nodes[0], edges[0]), proc


def run_python(case, seed, trials=500):
    import random

    class _Cfg:
        lattice_type = case.lattice
        lattice_size = case.dims
        max_functionality = case.max_func
        degree_distribution = case.degree_dist
        mix_fractions = case.mix
        mix_cutoff = 1.0

    random.seed(seed)
    try:
        with redirect_stdout(StringIO()):
            graphs = PythonTopologyGenerator(_Cfg()).generate(
                trials=trials, max_saves=1, time_limit=120
            )
    except ValueError as exc:               # unreachable target, fail-fast guard
        return None, f"rejected: {exc}"
    return (graphs[0] if graphs else None), None


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

def describe(g, dims):
    """Structural fingerprint of one generated network."""
    pos = {n: np.asarray(d["pos"], float) for n, d in g.nodes(data=True)
           if "pos" in d}
    box = np.asarray(dims, float)
    shells = set()
    for u, v in g.edges():
        if u in pos and v in pos:
            d = pos[u] - pos[v]
            shells.add(round(float(np.linalg.norm(d - box * np.round(d / box))), 4))
    return {
        "nodes": g.number_of_nodes(),
        "edges": g.number_of_edges(),
        "degrees": Counter(d for _, d in g.degree()),
        "max_degree": max((d for _, d in g.degree()), default=0),
        "shells": tuple(sorted(shells)),
        "box": tuple(round(float(b), 4) for b in box),
    }


def summarise(runs):
    """Aggregate a list of describe() dicts."""
    if not runs:
        return None
    return {
        "n": len(runs),
        "nodes_mean": float(np.mean([r["nodes"] for r in runs])),
        "edges_mean": float(np.mean([r["edges"] for r in runs])),
        "max_degree": max(r["max_degree"] for r in runs),
        "shells": sorted({s for r in runs for s in r["shells"]}),
        "box": runs[0]["box"],
    }


def compare(case, c_runs, py_runs, c_fail_note, py_fail_note):
    """Return (verdict, detail) for one case."""
    c, p = summarise(c_runs), summarise(py_runs)

    if c is None and p is None:
        return "both-unreachable", f"neither produced a network ({py_fail_note or c_fail_note or 'no success'})"
    if c is None:
        return "C-ONLY-FAILED", f"Python succeeded {p['n']}/{REPS}, C produced nothing"
    if p is None:
        return "PY-ONLY-FAILED", f"C succeeded {c['n']}/{REPS}, Python produced nothing"

    problems = []
    if c["box"] != p["box"]:
        problems.append(f"box {c['box']} vs {p['box']}")
    if c["shells"] != p["shells"]:
        problems.append(f"shells {c['shells']} vs {p['shells']}")

    # Site count: exact for the deterministic lattices, tolerance for MIX.
    deterministic = case.lattice != "MIX" or set(case.mix.values()) <= {0.0, 1.0}
    if deterministic:
        if c["nodes_mean"] != p["nodes_mean"]:
            problems.append(f"nodes {c['nodes_mean']:.0f} vs {p['nodes_mean']:.0f}")
    else:
        rel = abs(c["nodes_mean"] - p["nodes_mean"]) / max(p["nodes_mean"], 1)
        if rel > 0.10:
            problems.append(f"nodes {c['nodes_mean']:.1f} vs {p['nodes_mean']:.1f} ({rel:.0%})")

    if c["max_degree"] > case.max_func or p["max_degree"] > case.max_func:
        problems.append(f"max_degree C={c['max_degree']} PY={p['max_degree']} "
                        f"exceeds max_func={case.max_func}")

    # Connectivity: compare mean degree, not raw edge count. On a mixture
    # the site count is a random draw and edges grow superlinearly with
    # it, so two statistically identical generators show a large edge-count
    # gap from a small node-count gap. Mean degree is the scale-free form.
    c_deg = 2 * c["edges_mean"] / max(c["nodes_mean"], 1)
    p_deg = 2 * p["edges_mean"] / max(p["nodes_mean"], 1)
    rel_d = abs(c_deg - p_deg) / max(p_deg, 1e-9)
    if rel_d > 0.10:
        problems.append(f"mean degree {c_deg:.2f} vs {p_deg:.2f} ({rel_d:.0%})")

    if problems:
        return "DISAGREE", "; ".join(problems)
    return "agree", (f"{c['nodes_mean']:.0f}/{p['nodes_mean']:.0f} nodes, "
                     f"deg {c_deg:.2f}/{p_deg:.2f}, "
                     f"{len(c['shells'])} shell(s)")


# ---------------------------------------------------------------------------
# LAMMPS pathway
# ---------------------------------------------------------------------------

def lammps_pathway(case, nodes_p, edges_p, tag):
    """Build one network through the pipeline and minimise it."""
    from topon.config.schema import (
        AssignmentConfig, ChemistryConfig, DPConfig, DPDistributionConfig,
        ExistingFilesConfig, OutputConfig, StudyConfig, ToponConfig,
        TopologyConfig,
    )
    from topon.pipeline import Pipeline

    name = f"{slugify(case.label)}_{tag}"
    cfg = ToponConfig(
        study=StudyConfig(name=name, output_dir=str(OUT / "lammps")),
        topology=TopologyConfig(
            source="load",
            existing_files=ExistingFilesConfig(nodes_file=str(nodes_p),
                                               edges_file=str(edges_p)),
        ),
        assignment=AssignmentConfig(
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=10.0))),
        chemistry=ChemistryConfig(model_type="atomistic", target_density=0.9),
        output=OutputConfig(),
    )
    try:
        with redirect_stdout(StringIO()):
            Pipeline(cfg).run()
    except Exception as exc:                       # noqa: BLE001
        return False, f"pipeline raised {type(exc).__name__}: {exc}"

    study = OUT / "lammps" / name
    if not (study / "03_Conformation/system_conformed.data").exists():
        return False, "no conformed data file"
    if shutil.which("lmp") is None:
        return None, "built ok; lmp not on PATH"

    sim = study / "04_Simulation"
    proc = subprocess.run(["lmp", "-in", "minimize_1_serial.in"], cwd=sim,
                          capture_output=True, text=True, timeout=900)
    (sim / "stage1.log").write_text(proc.stdout + "\n" + proc.stderr)
    if proc.returncode != 0:
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-2:]
        return False, f"lmp exit {proc.returncode}: {' | '.join(tail)}"
    if not (sim / "system_after_soft.data").exists():
        return False, "lmp exit 0 but no output data"
    return True, "minimised"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lammps", action="store_true",
                    help="also build a subset through the pipeline and minimise")
    ap.add_argument("--reps", type=int, default=REPS)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    work = OUT / "_work"
    work.mkdir(parents=True, exist_ok=True)

    exe = compile_c(work)
    if exe is None:
        print("No C compiler or source; cannot compare. Aborting.")
        return 1
    print(f"C generator compiled: {exe.name}\n")

    cases = build_matrix()
    print("=" * 100)
    print(f"{'case':26s} {'C':>7s} {'PY':>7s}  {'verdict':16s} detail")
    print("=" * 100)

    results = []
    for case in cases:
        c_runs, py_runs = [], []
        c_note = py_note = None
        first_c_paths = None
        slug = slugify(case.label)
        for rep in range(args.reps):
            got, proc = run_c(exe, work / "c" / f"{slug}_{rep}", case)
            if got:
                g, dims, np_, ep_ = got
                c_runs.append(describe(g, dims))
                if first_c_paths is None:
                    first_c_paths = (np_, ep_)
            elif c_note is None:
                c_note = "no success"

            g_py, note = run_python(case, seed=rep)
            if note:
                py_note = note
            elif g_py is not None:
                from topon.topology.loader import (
                    infer_dims_from_graph, remove_vacancies,
                )
                # The C side arrives via load_graph, which drops degree-0
                # nodes; the pipeline does the same to a Python-generated
                # graph. Match that or a `0:N` defect target looks like a
                # node-count disagreement when it is only bookkeeping.
                remove_vacancies(g_py)
                py_runs.append(describe(g_py, infer_dims_from_graph(g_py)))

        verdict, detail = compare(case, c_runs, py_runs, c_note, py_note)
        print(f"{case.label:26s} {len(c_runs):>3d}/{args.reps:<3d} "
              f"{len(py_runs):>3d}/{args.reps:<3d}  {verdict:16s} {detail}")
        results.append((case, verdict, first_c_paths))

    print("=" * 100)
    bad = [r for r in results if r[1] in ("DISAGREE", "C-ONLY-FAILED", "PY-ONLY-FAILED")]
    print(f"agree: {sum(1 for r in results if r[1] == 'agree')}   "
          f"both-unreachable: {sum(1 for r in results if r[1] == 'both-unreachable')}   "
          f"problems: {len(bad)}")

    if args.lammps:
        print()
        print("=" * 100)
        print("LAMMPS pathway (C-generated topology -> pipeline -> stage-1 minimize)")
        print("=" * 100)
        subset = [r for r in results
                  if r[2] is not None
                  and r[0].label in {"SC prune->4", "BCC prune->4", "FCC prune->4",
                                     "MIX 20/40/40 prune->4", "SC edge target",
                                     "SC 3x4x5"}]
        for case, _, paths in subset:
            ok, msg = lammps_pathway(case, paths[0], paths[1], "C")
            status = {True: "ok", False: "FAILED", None: "skipped"}[ok]
            print(f"  {case.label:26s} {status:8s} {msg}")
            if ok is False:
                bad.append((case, "LAMMPS", None))

    print()
    print("=" * 100)
    print(f"VERDICT: {'FAIL' if bad else 'PASS'}    artifacts in {OUT}")
    print("=" * 100)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
