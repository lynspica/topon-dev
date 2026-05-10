"""Smoke test: JSON config -> load_config_full -> Pipeline + LAMMPS stage 1.

Loads a real JSON config that includes the historic extras
(`conformation`, `simulation`, `execution` — sections not covered by the
ToponConfig schema), exercises `topon.config.load_config_full` to split
them into (validated config, raw dict), then runs Pipeline with both.
Pins the P0-A fix that added the schema-vs-raw split at load time.

Workarounds for still-open issues: none — this is the first smoke test
that exercises the full CLI-equivalent path (`topon generate <config>`)
end-to-end.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from topon.config import load_config_full
from topon.pipeline import Pipeline


pytestmark = [pytest.mark.requires_lammps]


REPO_ROOT = Path(__file__).resolve().parents[2]
SAMPLE_NODES = REPO_ROOT / "tests" / "sample_graphs" / "network_N5x5x5_trial3.nodes"
SAMPLE_EDGES = REPO_ROOT / "tests" / "sample_graphs" / "network_N5x5x5_trial3.edges"
FIXTURE_JSON = REPO_ROOT / "tests" / "smoke" / "fixtures" / "json_load_smoke.json"


def test_json_load_pipeline_runs_and_lammps_minimizes(tmp_path: Path) -> None:
    """Read a JSON config with extras, run Pipeline+LAMMPS via load_config_full."""
    assert FIXTURE_JSON.exists(), f"fixture missing: {FIXTURE_JSON}"
    assert SAMPLE_NODES.exists(), f"sample graph missing: {SAMPLE_NODES}"

    # Read fixture, fill in absolute paths for this test's tmp_path.
    fixture = json.loads(FIXTURE_JSON.read_text(encoding="utf-8"))
    fixture["study"]["output_dir"] = str(tmp_path)
    fixture["topology"]["existing_files"]["nodes_file"] = str(SAMPLE_NODES)
    fixture["topology"]["existing_files"]["edges_file"] = str(SAMPLE_EDGES)
    test_config_path = tmp_path / "smoke.json"
    test_config_path.write_text(json.dumps(fixture, indent=2), encoding="utf-8")

    # Load via the new split-aware loader. Schema sections validate; the
    # conformation / simulation / execution sections come back in raw_cfg.
    config, raw_cfg = load_config_full(test_config_path)
    assert sorted(raw_cfg.keys()) == ["conformation", "execution", "simulation"], (
        f"expected raw extras to include conformation/simulation/execution, "
        f"got {sorted(raw_cfg.keys())}"
    )
    # Note: validate_config has its own defaults-related warnings that
    # aren't relevant to the load+run smoke check; the CLI's `topon
    # validate` runs that pass independently.

    # Run Pipeline with both config and raw extras.
    Pipeline(config, raw_config=raw_cfg).run()

    study_dir = tmp_path / "smoke_json"
    assert (study_dir / "02_Chemistry" / "system.data").exists()

    sim_dir = study_dir / "04_Simulation"
    minimize_in = sim_dir / "minimize_1_serial.in"
    assert minimize_in.exists()

    result = subprocess.run(
        ["lmp", "-in", "minimize_1_serial.in"],
        cwd=sim_dir,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert result.returncode == 0, (
        f"LAMMPS stage-1 minimize failed (exit {result.returncode}):\n"
        f"--- stdout ---\n{result.stdout[-2000:]}\n"
        f"--- stderr ---\n{result.stderr[-2000:]}"
    )
    assert (sim_dir / "system_after_soft.data").exists()
