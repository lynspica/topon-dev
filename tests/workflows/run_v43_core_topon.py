"""v43: core topon CG + atomistic network smoke test at DP=25.

Runs the existing CG and atomistic combined workflows (entanglements + grafts
on a 5x5x5 lattice) with DP overridden to 25, and confirms LAMMPS minimization
sequences run cleanly. Outputs land in:

    tests/output/v43/
    ├── config_cg.json                  # patched DP=25 copy
    ├── config_atomistic.json           # patched DP=25 copy
    ├── cg_combined/                    # CG workflow output
    │   └── 04_Simulation/              # LAMMPS data + scripts + logs
    └── atomistic_combined/             # atomistic workflow output
        └── 04_Simulation/

The two workflows live in `tests/workflows/generate_cg_combined.py` and
`generate_atomistic_combined.py` and read their config from `examples/`. We
patch the module-level `CONFIG_FILE` constant before calling each
`run_workflow()` so we don't disturb the canonical example configs.

Run:
    python tests/workflows/run_v43_core_topon.py
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PKG))

V43_DP = 25
V43_DIR = PKG / "tests" / "output" / "v43"


def patched_config(src_path: Path, dest_path: Path, dp: int) -> Path:
    """Copy a JSON config and override degree_of_polymerization."""
    cfg = json.loads(src_path.read_text())
    cfg["chemistry"]["degree_of_polymerization"] = dp
    # Disable auto LAMMPS run so we have full control over which scripts
    # are executed and can capture per-stage exit codes.
    cfg.setdefault("execution", {})["auto_run"] = False
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    dest_path.write_text(json.dumps(cfg, indent=4))
    return dest_path


def run_lammps_sequence(sim_dir: Path, scripts: list[str]) -> dict:
    """Run a sequence of LAMMPS .in scripts inside `sim_dir`. Returns a
    dict keyed by script name with rc + wall-time + last log lines."""
    out: dict[str, dict] = {}
    for script in scripts:
        in_path = sim_dir / script
        log_name = script.replace(".in", ".log")
        log_path = sim_dir / log_name
        if not in_path.exists():
            out[script] = {"status": "MISSING_INPUT", "rc": None}
            continue
        ts = time.time()
        try:
            r = subprocess.run(
                ["lmp", "-in", script, "-log", log_name],
                cwd=str(sim_dir),
                capture_output=True, text=True, timeout=900,
            )
        except subprocess.TimeoutExpired:
            out[script] = {"status": "TIMEOUT", "rc": None,
                           "wall_s": round(time.time() - ts, 1)}
            return out  # don't continue past a timeout
        wall = round(time.time() - ts, 1)
        if r.returncode != 0:
            tail = "\n".join(r.stderr.splitlines()[-10:])
            out[script] = {"status": "FAILED", "rc": r.returncode,
                           "wall_s": wall, "stderr_tail": tail}
            return out
        # Pull the last thermo line + Total wall time
        last_thermo = None
        wall_lammps = None
        if log_path.exists():
            with log_path.open() as f:
                for line in f:
                    s = line.rstrip()
                    if s.startswith("Total wall time"):
                        wall_lammps = s
                    elif s.strip() and s[0].isdigit() and len(s.split()) > 4:
                        try:
                            [float(x) for x in s.split()]
                            last_thermo = s
                        except ValueError:
                            pass
        out[script] = {
            "status": "OK", "rc": 0, "wall_s": wall,
            "wall_lammps": wall_lammps,
            "last_thermo": last_thermo,
        }
    return out


def run_cg(out_root: Path) -> dict:
    """Run the CG combined workflow at DP=25."""
    print("\n" + "=" * 70)
    print(f"=== CG (DP={V43_DP})")
    print("=" * 70)
    src = PKG / "examples" / "config_cg_combined.json"
    cfg = patched_config(src, V43_DIR / "config_cg.json", V43_DP)

    import tests.workflows.generate_cg_combined as cg_mod  # type: ignore
    cg_mod.CONFIG_FILE = cfg

    t0 = time.time()
    try:
        cg_mod.run_workflow(out_root)
    except Exception as exc:
        return {"stage": "GEN_FAILED", "error": str(exc),
                "trace": traceback.format_exc()}
    t_gen = round(time.time() - t0, 1)

    sim_dir = Path(out_root) / "cg_combined" / "04_Simulation"
    if not sim_dir.exists():
        return {"stage": "GEN_FAILED",
                "error": f"04_Simulation not found at {sim_dir}"}

    scripts = ["minimize_1_serial.in",
               "minimize_2_parallel.in",
               "minimize_3_parallel.in"]
    print(f"\n--- LAMMPS for CG: {sim_dir} ---")
    lmp = run_lammps_sequence(sim_dir, scripts)
    return {"stage": "OK", "wall_gen_s": t_gen,
            "sim_dir": str(sim_dir), "lammps": lmp}


def run_atomistic(out_root: Path) -> dict:
    """Run the atomistic combined workflow at DP=25."""
    print("\n" + "=" * 70)
    print(f"=== ATOMISTIC (DP={V43_DP})")
    print("=" * 70)
    src = PKG / "examples" / "config_atomistic_combined.json"
    cfg = patched_config(src, V43_DIR / "config_atomistic.json", V43_DP)

    import tests.workflows.generate_atomistic_combined as atom_mod  # type: ignore
    atom_mod.CONFIG_FILE = cfg

    t0 = time.time()
    try:
        atom_mod.run_workflow(out_root)
    except Exception as exc:
        return {"stage": "GEN_FAILED", "error": str(exc),
                "trace": traceback.format_exc()}
    t_gen = round(time.time() - t0, 1)

    sim_dir = Path(out_root) / "atomistic_combined" / "04_Simulation"
    if not sim_dir.exists():
        return {"stage": "GEN_FAILED",
                "error": f"04_Simulation not found at {sim_dir}"}

    scripts = ["minimize_1_serial.in",
               "minimize_2_parallel.in",
               "minimize_3_parallel.in"]
    print(f"\n--- LAMMPS for atomistic: {sim_dir} ---")
    lmp = run_lammps_sequence(sim_dir, scripts)
    return {"stage": "OK", "wall_gen_s": t_gen,
            "sim_dir": str(sim_dir), "lammps": lmp}


def print_lmp_summary(label: str, result: dict) -> None:
    print(f"\n  [{label}] generation: {result.get('stage')} "
          f"(gen wall {result.get('wall_gen_s', '-')} s)")
    lmp = result.get("lammps") or {}
    if not lmp:
        print(f"    no LAMMPS results: {result.get('error', '-')}")
        return
    for script, info in lmp.items():
        line = f"    {script:<30} {info['status']:<10}"
        if info.get("wall_s") is not None:
            line += f"  wall={info['wall_s']}s"
        if info.get("last_thermo"):
            line += f"  thermo: {info['last_thermo'].strip()[:80]}"
        if info.get("stderr_tail"):
            line += f"\n      stderr: {info['stderr_tail']}"
        print(line)


def main() -> int:
    if V43_DIR.exists():
        shutil.rmtree(V43_DIR)
    V43_DIR.mkdir(parents=True)

    print(f"v43 core topon test at DP={V43_DP}")
    print(f"output: {V43_DIR}")

    cg_result = run_cg(V43_DIR)
    atom_result = run_atomistic(V43_DIR)

    print("\n" + "=" * 70)
    print("=== v43 SUMMARY")
    print("=" * 70)
    print_lmp_summary("CG", cg_result)
    print_lmp_summary("atomistic", atom_result)

    cg_ok = cg_result.get("stage") == "OK" and all(
        info.get("status") == "OK" for info in (cg_result.get("lammps") or {}).values()
    )
    atom_ok = atom_result.get("stage") == "OK" and all(
        info.get("status") == "OK" for info in (atom_result.get("lammps") or {}).values()
    )
    return 0 if (cg_ok and atom_ok) else 1


if __name__ == "__main__":
    raise SystemExit(main())
