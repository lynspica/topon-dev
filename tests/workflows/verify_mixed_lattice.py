"""End-to-end check that a mixed lattice builds and runs.

Run directly (not via pytest)::

    python tests/workflows/verify_mixed_lattice.py

The unit tests in ``tests/unit/topology/test_mixed_lattice.py`` cover the
point set and the edge rule. This script answers the question they cannot:
does a mixed lattice survive the rest of the pipeline and minimise in
LAMMPS?

That matters because mixing widens the spread of junction separations. A
BCC body centre and an FCC face centre can sit 0.5 cells apart against the
1.0 of the simple-cubic shell, and DP is assigned independently of edge
length, so strands of the same DP get built at bond lengths differing by
up to 2x. Whether that is tolerable is an empirical question about FENE
strain, so this script measures it rather than asserting it.

Compares a pure SC run against a 0.2/0.4/0.4 mixture at matched cell count
and DP. Artifacts land in ``tests/output/mixed_lattice/``.
"""
from __future__ import annotations

import random
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from topon.config.schema import (  # noqa: E402
    AssignmentConfig,
    ChemistryConfig,
    DPConfig,
    DPDistributionConfig,
    ExistingFilesConfig,
    OutputConfig,
    StudyConfig,
    ToponConfig,
    TopologyConfig,
)
from topon.pipeline import Pipeline  # noqa: E402
from topon.topology.generator_python import PythonTopologyGenerator  # noqa: E402
from topon.topology.loader import save_nodes_edges  # noqa: E402

OUT = ROOT / "tests/output/mixed_lattice"
DIMS = (4, 4, 4)
SEED = 42

CASES = [
    ("sc_baseline", {"SC": 1.0, "BCC": 0.0, "FCC": 0.0}),
    ("mix_20_40_40", {"SC": 0.2, "BCC": 0.4, "FCC": 0.4}),
]


class _Cfg:
    def __init__(self, mix, degree_distribution="0:0,1:0", max_func=4):
        self.lattice_type = "MIX"
        self.lattice_size = DIMS
        self.max_functionality = max_func
        self.degree_distribution = degree_distribution
        self.mix_fractions = mix
        self.mix_cutoff = 1.0


def lattice_edge_lengths(G):
    box = np.asarray(G.graph["box"], dtype=float)
    pos = {n: np.asarray(d["pos"], dtype=float) for n, d in G.nodes(data=True)}
    return np.asarray([
        np.linalg.norm((pos[u] - pos[v]) - box * np.round((pos[u] - pos[v]) / box))
        for u, v in G.edges()
    ])


def parse_data_file(path):
    """Return (box_lengths, {atom_id: xyz}, [(a, b), ...]) from a data file."""
    box = np.zeros(3)
    atoms, bonds = {}, []
    section = None
    for raw in Path(path).read_text().splitlines():
        line = raw.split("#")[0].strip()
        if not line:
            continue
        for i, tag in enumerate(("xlo xhi", "ylo yhi", "zlo zhi")):
            if line.endswith(tag):
                lo, hi = (float(v) for v in line.split()[:2])
                box[i] = hi - lo
        w = line.split()
        if w[0] in ("Atoms", "Bonds", "Velocities", "Masses",
                    "Angles", "Dihedrals", "Impropers"):
            section = w[0]
            continue
        if section == "Atoms" and len(w) >= 7:
            atoms[int(w[0])] = np.array([float(v) for v in w[4:7]])
        elif section == "Bonds" and len(w) >= 4:
            bonds.append((int(w[2]), int(w[3])))
    return box, atoms, bonds


def bond_lengths(path):
    box, atoms, bonds = parse_data_file(path)
    out = []
    for a, b in bonds:
        if a in atoms and b in atoms:
            d = atoms[a] - atoms[b]
            out.append(np.linalg.norm(d - box * np.round(d / box)))
    return np.asarray(out)


def build_config(name, nodes_p, edges_p):
    return ToponConfig(
        study=StudyConfig(name=name, output_dir=str(OUT)),
        topology=TopologyConfig(
            source="load",
            existing_files=ExistingFilesConfig(
                nodes_file=str(nodes_p), edges_file=str(edges_p)
            ),
        ),
        assignment=AssignmentConfig(
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=10.0)),
        ),
        chemistry=ChemistryConfig(model_type="atomistic", target_density=0.9),
        output=OutputConfig(),
    )


def run_lammps(sim_dir):
    script = sim_dir / "minimize_1_serial.in"
    if not script.exists():
        return False, "stage-1 script was not written"
    if shutil.which("lmp") is None:
        return None, "lmp not on PATH, skipped"
    proc = subprocess.run(
        ["lmp", "-in", script.name], cwd=sim_dir,
        capture_output=True, text=True, timeout=900,
    )
    (sim_dir / "stage1.log").write_text(proc.stdout + "\n" + proc.stderr)
    if proc.returncode != 0:
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-3:]
        return False, f"exit {proc.returncode}: {' | '.join(tail)}"
    if not (sim_dir / "system_after_soft.data").exists():
        return False, "exit 0 but no system_after_soft.data"
    return True, "completed"


def run_case(name, mix):
    print(f"--- {name}  {mix} ---")
    random.seed(SEED)
    gen = PythonTopologyGenerator(_Cfg(mix))
    graphs = gen.generate(trials=8000, max_saves=1)
    if not graphs:
        print("    FAIL: sculpting produced no graph")
        return None
    G = graphs[0]

    lat = lattice_edge_lengths(G)
    shells = sorted(set(np.round(lat, 4)))
    print(f"    sculpted         : {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges, box {G.graph['box']}")
    print(f"    lattice shells   : {shells}")
    print(f"    junction spread  : {lat.max() / lat.min():.2f}x "
          f"({lat.min():.4f} to {lat.max():.4f} cells)")

    topo = OUT / "topology"
    nodes_p, edges_p = topo / f"{name}.nodes", topo / f"{name}.edges"
    save_nodes_edges(G, nodes_p, edges_p)

    Pipeline(build_config(name, nodes_p, edges_p)).run()

    conformed = OUT / name / "03_Conformation/system_conformed.data"
    if not conformed.exists():
        print(f"    FAIL: no conformed output at {conformed}")
        return None

    bl = bond_lengths(conformed)
    ok, msg = run_lammps(OUT / name / "04_Simulation")
    # Report the upper tail, not max/min. The data file mixes backbone,
    # pendant and hydrogen bonds, whose lengths differ for reasons that
    # have nothing to do with the lattice, so a max/min ratio is dominated
    # by the shortest hydrogen bond and says nothing useful.
    p99 = float(np.percentile(bl, 99))
    print(f"    bonds            : {len(bl)}, mean {bl.mean():.3f} A, "
          f"p99 {p99:.3f} A, max {bl.max():.3f} A")
    print(f"    LAMMPS stage 1   : {msg}")
    print()
    return {"name": name, "lattice": lat, "bonds": bl, "p99": p99,
            "shells": shells, "lammps_ok": ok}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("=" * 78)
    print("Mixed lattice end-to-end: pipeline + LAMMPS stage-1 minimize")
    print("=" * 78)

    results = [r for r in (run_case(n, m) for n, m in CASES) if r]
    if len(results) != len(CASES):
        print("VERDICT: FAIL (a case did not build)")
        return 1

    print("-" * 78)
    print(f"{'case':16s} {'shells':>7s} {'junction':>10s} "
          f"{'p99 bond':>10s} {'max bond':>10s} {'LAMMPS':>10s}")
    for r in results:
        status = {True: "ok", False: "FAILED", None: "skipped"}[r["lammps_ok"]]
        print(f"{r['name']:16s} {len(r['shells']):7d} "
              f"{r['lattice'].max() / r['lattice'].min():9.2f}x "
              f"{r['p99']:10.3f} {r['bonds'].max():10.3f} {status:>10s}")
    print("-" * 78)
    print("More shells is the point of mixing: it spreads strand end-to-end")
    print("distances over several values instead of one. The cost shows in")
    print("the bond tail, since DP does not track edge length, so a strand")
    print("spanning a 1.0-cell gap gets the same monomer count as one")
    print("spanning 0.5.")

    failed = [r["name"] for r in results if r["lammps_ok"] is False]
    print()
    print("=" * 78)
    verdict = "FAIL" if failed else "PASS"
    print(f"VERDICT: {verdict}"
          + (f"  (LAMMPS failed: {', '.join(failed)})" if failed else "")
          + f"    artifacts in {OUT}")
    print("=" * 78)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
