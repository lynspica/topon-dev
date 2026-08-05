"""End-to-end verification of the recorded-periodic-cell fix.

Run directly (not via pytest)::

    python tests/workflows/verify_lattice_box.py

What this checks
----------------
``infer_dims_from_graph`` used to derive the periodic cell from node
positions as ``max - min + 1``. That is exact for simple cubic, whose
sites sit on integer coordinates with unit spacing, but wrong for every
lattice with fractional basis sites. Generators now record the true cell
in ``G.graph["box"]`` and ``.nodes`` files carry it in a ``# BOX`` header.

The old behaviour is still reachable through a supported code path: a
``.nodes`` file written *without* the header falls back to the estimate.
So this script writes the same BCC network twice, once with the header
and once without, and runs both through the real pipeline. No
monkeypatching is involved; the contrast comes entirely from whether the
box was recorded.

Part 1  lattice-level audit of SC / BCC / FCC / Diamond
Part 2  full pipeline on both .nodes variants, then LAMMPS stage-1 minimize

Artifacts land in ``tests/output/lattice_box_fix/``.
"""
from __future__ import annotations

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
    GeneratorConfig,
    OutputConfig,
    StudyConfig,
    ToponConfig,
    TopologyConfig,
)
from topon.pipeline import Pipeline  # noqa: E402
from topon.topology.generator_python import PythonTopologyGenerator  # noqa: E402
from topon.topology.generator_python_diamond import (  # noqa: E402
    create_diamond_lattice,
)
from topon.topology.loader import (  # noqa: E402
    infer_dims_from_graph,
    save_nodes_edges,
)

OUT = ROOT / "tests/output/lattice_box_fix"

# Nearest-neighbour separation of each lattice, in units of the cubic
# cell. These are the values every edge of a full lattice must take.
EXPECTED_NN = {
    "SC": 1.0,
    "BCC": np.sqrt(3) / 2,
    "FCC": np.sqrt(2) / 2,
    "Diamond": np.sqrt(3) / 4,
}


class _Cfg:
    def __init__(self, lattice_type, lattice_size, degree_distribution="", max_func=6):
        self.lattice_type = lattice_type
        self.lattice_size = lattice_size
        self.max_functionality = max_func
        self.degree_distribution = degree_distribution


def edge_lengths(G, box):
    """Minimum-image edge lengths under ``box``."""
    box = np.asarray(box, dtype=float)
    pos = {n: np.asarray(d["pos"], dtype=float) for n, d in G.nodes(data=True)}
    return np.asarray([
        np.linalg.norm((pos[u] - pos[v]) - box * np.round((pos[u] - pos[v]) / box))
        for u, v in G.edges()
    ])


# ---------------------------------------------------------------------------
# Part 1: lattice-level audit
# ---------------------------------------------------------------------------

def audit_lattices():
    print("=" * 78)
    print("PART 1  Periodic cell and edge lengths, recorded box vs old estimate")
    print("=" * 78)
    print(f"{'lattice':9s} {'recorded':>10s} {'estimate':>10s} "
          f"{'true NN':>9s} {'bad edges (est.)':>18s} {'worst len':>10s}")
    print("-" * 78)

    lattices = {
        lt: PythonTopologyGenerator(_Cfg(lt, (4, 4, 4)))._create_lattice((4, 4, 4), lt)
        for lt in ("SC", "BCC", "FCC")
    }
    lattices["Diamond"] = create_diamond_lattice(4, 4, 4)

    ok = True
    for name, G in lattices.items():
        recorded = infer_dims_from_graph(G)
        stripped = G.copy()
        del stripped.graph["box"]
        estimate = infer_dims_from_graph(stripped)

        nn = EXPECTED_NN[name]
        len_rec = edge_lengths(G, recorded)
        len_est = edge_lengths(G, estimate)
        bad = int(np.sum(~np.isclose(len_est, nn)))

        print(f"{name:9s} {recorded[0]:10.2f} {estimate[0]:10.2f} {nn:9.4f} "
              f"{bad:6d} / {len(len_est):<9d} {len_est.max():10.4f}")

        if not np.allclose(len_rec, nn):
            print(f"    FAIL: recorded box does not give uniform {nn:.4f} bonds")
            ok = False

    print("-" * 78)
    print("All edges of a full lattice are nearest-neighbour bonds under the")
    print("recorded box. Under the old estimate the flagged edges resolve to")
    print("the wrong periodic replica and measure twice their true length.")
    print()
    return ok, lattices


# ---------------------------------------------------------------------------
# Part 2: full pipeline + LAMMPS
# ---------------------------------------------------------------------------

def parse_data_file(path):
    """Return (box_lengths, atom_xyz_by_id, bond_pairs) from a LAMMPS data file."""
    lines = Path(path).read_text().splitlines()
    box = np.zeros(3)
    atoms, bonds = {}, []
    section = None
    for raw in lines:
        line = raw.split("#")[0].strip()
        if not line:
            continue
        for i, tag in enumerate(("xlo xhi", "ylo yhi", "zlo zhi")):
            if line.endswith(tag):
                lo, hi = (float(v) for v in line.split()[:2])
                box[i] = hi - lo
        if line.split()[0] in ("Atoms", "Bonds", "Velocities", "Masses",
                               "Angles", "Dihedrals", "Impropers"):
            section = line.split()[0]
            continue
        parts = line.split()
        if section == "Atoms" and len(parts) >= 7:
            # id mol type charge x y z  (full style)
            atoms[int(parts[0])] = np.array([float(v) for v in parts[4:7]])
        elif section == "Bonds" and len(parts) >= 4:
            bonds.append((int(parts[2]), int(parts[3])))
    return box, atoms, bonds


def bond_length_stats(path):
    box, atoms, bonds = parse_data_file(path)
    lengths = []
    for a, b in bonds:
        if a not in atoms or b not in atoms:
            continue
        d = atoms[a] - atoms[b]
        lengths.append(np.linalg.norm(d - box * np.round(d / box)))
    return np.asarray(lengths)


def build_config(name, nodes_p, edges_p, out_dir):
    return ToponConfig(
        study=StudyConfig(name=name, output_dir=str(out_dir)),
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
    """Run the stage-1 soft minimize. Returns (ok, message)."""
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
    return True, "stage-1 minimize completed"


def pipeline_case(label, nodes_p, edges_p):
    print(f"--- {label} ---")
    study_dir = OUT / label
    Pipeline(build_config(label, nodes_p, edges_p, OUT)).run()

    # Stage 4 writes every atom at the origin; the real coordinates are
    # applied in stage 5, so the conformed file is the one to measure.
    data = study_dir / "03_Conformation/system_conformed.data"
    if not data.exists():
        print(f"    FAIL: no conformed output at {data}")
        return None

    box, _, _ = parse_data_file(data)
    lengths = bond_length_stats(data)
    stretched = int(np.sum(lengths > 3.0))
    print(f"    box              : {box[0]:.3f} A")
    print(f"    bonds            : {len(lengths)}")
    print(f"    bond length      : mean {lengths.mean():.3f}  "
          f"max {lengths.max():.3f} A")
    print(f"    bonds over 3 A   : {stretched}   "
          f"(chains left stretched across the box)")

    ok, msg = run_lammps(study_dir / "04_Simulation")
    print(f"    LAMMPS stage 1   : {msg}")
    print()
    return {"lengths": lengths, "lammps_ok": ok, "lammps_msg": msg}


def run_pipeline_comparison():
    print("=" * 78)
    print("PART 2  Same BCC network through the pipeline, with and without")
    print("        the recorded box, then LAMMPS stage-1 minimize")
    print("=" * 78)

    gen = PythonTopologyGenerator(
        _Cfg("BCC", (4, 4, 4), degree_distribution="0:0,1:0", max_func=4)
    )
    graphs = gen.generate(trials=4000, max_saves=1)
    if not graphs:
        print("FAIL: sculpting produced no BCC graph")
        return False
    G = graphs[0]
    print(f"BCC 4x4x4 sculpted: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges, box {G.graph['box']}\n")

    topo = OUT / "topology"
    topo.mkdir(parents=True, exist_ok=True)

    with_box = topo / "bcc_with_box.nodes"
    save_nodes_edges(G, with_box, topo / "bcc_with_box.edges")

    # The old behaviour: strip the box so no header is written and the
    # loader falls back to the max-min+1 estimate. Passing box=None would
    # not do it, since that means "take it from the graph".
    no_box = topo / "bcc_no_box.nodes"
    stripped = G.copy()
    stripped.graph.pop("box", None)
    save_nodes_edges(stripped, no_box, topo / "bcc_no_box.edges")

    old = pipeline_case("bcc_no_box_OLD", no_box, topo / "bcc_no_box.edges")
    new = pipeline_case("bcc_with_box_FIXED", with_box, topo / "bcc_with_box.edges")

    if old is None or new is None:
        return False

    print("-" * 78)
    print(f"{'':22s} {'max bond (A)':>14s} {'bonds > 3 A':>13s} {'LAMMPS':>12s}")
    for lbl, r in (("old (estimated box)", old), ("fixed (recorded box)", new)):
        status = {True: "ok", False: "FAILED", None: "skipped"}[r["lammps_ok"]]
        print(f"{lbl:22s} {r['lengths'].max():14.3f} "
              f"{int(np.sum(r['lengths'] > 3.0)):13d} {status:>12s}")
    print("-" * 78)

    # The fixed run must leave no chain stretched across the box.
    n_stretched = int(np.sum(new["lengths"] > 3.0))
    if n_stretched:
        print(f"FAIL: {n_stretched} bonds over 3 A in the fixed run")
        return False
    return new["lammps_ok"] is not False


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ok_lattice, _ = audit_lattices()
    ok_pipeline = run_pipeline_comparison()
    print()
    print("=" * 78)
    verdict = "PASS" if (ok_lattice and ok_pipeline) else "FAIL"
    print(f"VERDICT: {verdict}    artifacts in {OUT}")
    print("=" * 78)
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
