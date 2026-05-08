"""v41: sweep water bead type x water density x protein variant.

Generates one MARTINI 3 protein-network system per cell, runs LAMMPS stages
1-3 (hierarchical soft-push -> LJ ramp -> NVT @ 310 K), and records the
final equilibrated state. Outputs land in tests/output/v41/<label>/.

Matrix (7 cells):
    --- nat_pro varying water content ---
    cell 0: nat_pro,  dry  (water=0)            -- collapse probe
    cell 1: nat_pro,  W medium (4 beads/nm^3)
    cell 2: nat_pro,  W bulk   (8 beads/nm^3)

    --- nat_pro varying water bead size at same H2O number-density ---
    cell 3: nat_pro,  SW (3 H2O/bead) at 5.33/nm^3  ~= W=4 H2O density
    cell 4: nat_pro,  TW (2 H2O/bead) at 8.0/nm^3   ~= W=4 H2O density

    --- different protein resolution at fixed water condition ---
    cell 5: high_pro, W=4   (540 atoms/chain)
    cell 6: no_pro,   W=4   (468 atoms/chain)

After all runs, prints a table of per-cell metrics:
    n_atoms, gel_point_conv, stage3 final PE/T/ebond/evdwl, wall time.

Usage:
    python tests/workflows/run_v41_matrix.py
"""
from __future__ import annotations

import shutil
import subprocess
import time
from pathlib import Path

from topon.protein_network import workflow

# Site-specific LAMMPS executable.
LMP = r"C:/Users/ahmet/AppData/Local/LAMMPS 64-bit 2Apr2025/bin/lmp.exe"

# (label, block_seq, water_type, density, n_chains, n_repeats, n_na, n_cl)
MATRIX: list[tuple] = [
    # nat_pro x water content
    ("dry_nat",       "GGRPSDSYGAPGGGN", "W",  0.0,  8, 6, 0,  0),
    ("W4_nat",        "GGRPSDSYGAPGGGN", "W",  4.0,  8, 6, 16, 16),
    ("W8_nat",        "GGRPSDSYGAPGGGN", "W",  8.0,  8, 6, 16, 16),
    # nat_pro x water bead size (same H2O number density as W=4 -> 16 H2O/nm^3)
    ("SW5p33_nat",    "GGRPSDSYGAPGGGN", "SW", 5.33, 8, 6, 16, 16),
    ("TW8_nat",       "GGRPSDSYGAPGGGN", "TW", 8.0,  8, 6, 16, 16),
    # protein variant x default water
    ("W4_highpro",    "GPRPSDSYGAPGPGN", "W",  4.0,  8, 6, 16, 16),
    ("W4_nopro",      "GGRGSDSYGAGGGGN", "W",  4.0,  8, 6, 16, 16),
]

OUTPUT_ROOT = Path("tests/output/v41")


def parse_thermo_last(log_path: Path) -> dict[str, float]:
    """Extract the last thermo line from a LAMMPS log. Returns dict by column name."""
    cols: list[str] = []
    last_line: str | None = None
    in_block = False
    with log_path.open() as f:
        for line in f:
            s = line.rstrip()
            if "Step" in s and "PotEng" in s:
                cols = s.split()
                in_block = True
                continue
            if in_block and s.startswith("Loop time"):
                in_block = False
                continue
            if in_block and s.strip() and not s.lstrip().startswith(("WARNING", "Generated")):
                parts = s.split()
                # rows are numeric; ignore lines that aren't (e.g. unfix messages)
                try:
                    [float(x) for x in parts]
                except ValueError:
                    continue
                last_line = s
    if last_line is None or not cols:
        return {}
    parts = last_line.split()
    return {c: float(v) for c, v in zip(cols, parts) if len(parts) >= len(cols)}


def cell_volume_initial(data_path: Path) -> float:
    box = [0.0, 0.0, 0.0]
    with data_path.open() as f:
        for line in f:
            s = line.strip()
            if "xlo xhi" in s: box[0] = float(s.split()[1])
            elif "ylo yhi" in s: box[1] = float(s.split()[1])
            elif "zlo zhi" in s: box[2] = float(s.split()[1])
            elif s.startswith(("Atoms", "Masses")):
                break
    return box[0] * box[1] * box[2]


def run_cell(cell: tuple) -> dict:
    label, block_seq, water_type, density, n_chains, n_repeats, n_na, n_cl = cell
    out_dir = OUTPUT_ROOT / label
    # Skip already-completed cells (resumable matrix).
    eq_data = out_dir / "system_equilibrated.data"
    if eq_data.exists():
        relax = out_dir / "relaxation"
        final = parse_thermo_last(relax / "stage3.log") if (relax / "stage3.log").exists() else {}
        return {
            "label": label, "status": "OK_CACHED",
            "PE": final.get("PotEng"),
            "T": final.get("Temp"),
            "ebond": final.get("E_bond"),
            "eangle": final.get("E_angle"),
            "edihed": final.get("E_dihed"),
            "evdwl": final.get("E_vdwl"),
            "ecoul": final.get("E_coul"),
            "press": final.get("Press"),
            "vol_nm3": cell_volume_initial(out_dir / f"{label}.data") / 1000.0,
            "n_atoms": _read_n_atoms(out_dir / f"{label}.data"),
            "gel_conv": _read_gel_conv(out_dir / f"{label}_topology.json"),
            "wall_gen_s": 0.0, "wall_lmp_s": 0.0,
            "n_atoms_per_nm3": _read_n_atoms(out_dir / f"{label}.data") / max(1.0, cell_volume_initial(out_dir / f"{label}.data") / 1000.0),
        }
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== {label} === block={block_seq} water={water_type}@{density}/nm^3 "
          f"chains={n_chains} repeats={n_repeats} ions={n_na}NA+{n_cl}CL")
    t0 = time.time()
    try:
        artifacts = workflow.run_protein_network(
            block_seq=block_seq,
            n_repeats=n_repeats, n_chains=n_chains,
            output_dir=out_dir,
            snapshot_label="gel_point",
            segs_per_block=2, target_packing=0.45,
            equil_steps=20000,
            n_extra_snapshots=1, snapshot_delta_conv=0.05,
            seed=42,
            water_density_w_per_nm3=density,
            water_bead_type=water_type,
            water_exclusion_ang=4.0,
            n_na_ions=n_na, n_cl_ions=n_cl, ion_exclusion_ang=4.0,
            base_name=label,
            use_itp_template=True,
            verbose=False,
        )
    except Exception as exc:
        return {"label": label, "status": "GEN_FAILED", "error": str(exc)}
    t_gen = time.time() - t0

    # Count atoms / box from data file
    data_path = artifacts["data"]
    vol_ang3 = cell_volume_initial(data_path)
    n_atoms = 0
    with data_path.open() as f:
        for line in f:
            s = line.strip()
            if s.endswith(" atoms"):
                n_atoms = int(s.split()[0]); break

    # Run LAMMPS stages 1-3. Per-stage timeout is generous since the LJ ramp
    # in stage 2 scales with atom count (W=8 W/nm^3 has ~50k atoms in this box).
    relax = artifacts["stage1"].parent
    stage_times: dict[str, float] = {}
    for stage in ("stage1", "stage2", "stage3"):
        in_path = artifacts[stage]
        log_path = relax / f"{stage}.log"
        ts = time.time()
        try:
            r = subprocess.run(
                [LMP, "-in", in_path.name, "-log", log_path.name],
                cwd=str(relax),
                capture_output=True, text=True, timeout=1800,
            )
        except subprocess.TimeoutExpired:
            stage_times[stage] = time.time() - ts
            return {
                "label": label, "status": f"{stage}_TIMEOUT",
                "wall_lmp_s": round(sum(stage_times.values()), 1),
            }
        stage_times[stage] = time.time() - ts
        if r.returncode != 0:
            return {
                "label": label, "status": f"{stage}_FAILED",
                "rc": r.returncode,
                "stderr_tail": "\n".join(r.stderr.splitlines()[-12:]),
            }

    final_thermo = parse_thermo_last(relax / "stage3.log")
    return {
        "label": label,
        "status": "OK",
        "n_atoms": n_atoms,
        "vol_nm3": vol_ang3 / 1000.0,
        "n_atoms_per_nm3": n_atoms / (vol_ang3 / 1000.0),
        "gel_conv": _read_gel_conv(out_dir / f"{label}_topology.json"),
        "PE": final_thermo.get("PotEng"),
        "T": final_thermo.get("Temp"),
        "ebond": final_thermo.get("E_bond"),
        "eangle": final_thermo.get("E_angle"),
        "edihed": final_thermo.get("E_dihed"),
        "evdwl": final_thermo.get("E_vdwl"),
        "ecoul": final_thermo.get("E_coul"),
        "press": final_thermo.get("Press"),
        "wall_gen_s": round(t_gen, 1),
        "wall_lmp_s": round(sum(stage_times.values()), 1),
    }


def _read_gel_conv(topo_path: Path) -> float | None:
    import json
    if not topo_path.exists(): return None
    data = json.loads(topo_path.read_text())
    for snap in data["snapshots"]:
        if snap["label"] == "gel_point":
            return float(snap["conv"])
    return None


def _read_n_atoms(data_path: Path) -> int:
    if not data_path.exists(): return 0
    with data_path.open() as f:
        for line in f:
            s = line.strip()
            if s.endswith(" atoms"):
                return int(s.split()[0])
    return 0


def main() -> int:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []
    for cell in MATRIX:
        results.append(run_cell(cell))

    # Print summary
    print("\n" + "=" * 110)
    print(f"{'label':<14} {'status':<14} {'n_atoms':>8} {'vol_nm3':>9} "
          f"{'gel':>6} {'PE':>10} {'T':>7} {'ebond':>9} {'evdwl':>11} "
          f"{'press':>9} {'gen_s':>6} {'lmp_s':>6}")
    print("-" * 110)
    for r in results:
        if r["status"] not in ("OK", "OK_CACHED"):
            print(f"{r['label']:<14} {r['status']:<14}  ... see {OUTPUT_ROOT / r['label']}")
            continue
        print(f"{r['label']:<14} {r['status']:<14} {r['n_atoms']:>8d} "
              f"{r['vol_nm3']:>9.1f} {r.get('gel_conv', 0) or 0:>6.3f} "
              f"{r['PE'] or 0:>10.0f} {r['T'] or 0:>7.1f} {r['ebond'] or 0:>9.0f} "
              f"{r['evdwl'] or 0:>11.0f} {r['press'] or 0:>9.1f} "
              f"{r['wall_gen_s']:>6.1f} {r['wall_lmp_s']:>6.1f}")
    print("=" * 110)

    return 0 if all(r["status"] in ("OK", "OK_CACHED") for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
