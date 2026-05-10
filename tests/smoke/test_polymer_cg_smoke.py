"""Smoke test: small CG (Kremer-Grest) network through Pipeline + LAMMPS stage 1.

Mirror of test_polymer_atomistic_smoke.py for the CG path. Loads the
5x5x5 sample graph, runs the full 6-stage Pipeline at DP=10 with
`model_type="coarse_grained"`, then invokes LAMMPS to run the stage-1
minimize script. Asserts the pipeline emits the expected files and
LAMMPS exits 0.

Why CG smoke matters: P0-C used to silently mis-route any
`coarse_grained` system through the atomistic writer branch (the writer
only knew the legacy `"cg"` literal, not the schema's `"coarse_grained"`).
This test pins the fix.

Workarounds for still-open Pipeline issues (see
internal/DEVELOPMENT_INTERNAL.md sec.1):
- Uses `topology.source="load"` (P0-B blocks `"generate"`).
- Builds `ToponConfig` programmatically (P0-A blocks `load_config`).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from topon.config.schema import (
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
from topon.pipeline import Pipeline


pytestmark = [pytest.mark.requires_lammps]


REPO_ROOT = Path(__file__).resolve().parents[2]
SAMPLE_NODES = REPO_ROOT / "tests" / "sample_graphs" / "network_N5x5x5_trial3.nodes"
SAMPLE_EDGES = REPO_ROOT / "tests" / "sample_graphs" / "network_N5x5x5_trial3.edges"


@pytest.fixture
def smoke_cg_config(tmp_path: Path) -> ToponConfig:
    """Minimal-scale CG config: 5x5x5 SC sample graph, DP=10."""
    return ToponConfig(
        study=StudyConfig(name="smoke_cg", output_dir=str(tmp_path)),
        topology=TopologyConfig(
            source="load",
            existing_files=ExistingFilesConfig(
                nodes_file=str(SAMPLE_NODES),
                edges_file=str(SAMPLE_EDGES),
            ),
        ),
        assignment=AssignmentConfig(
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=10.0)),
        ),
        chemistry=ChemistryConfig(model_type="coarse_grained", target_density=0.85),
        output=OutputConfig(),
    )


def test_cg_pipeline_runs_and_lammps_minimizes(
    smoke_cg_config: ToponConfig, tmp_path: Path
) -> None:
    """Generate small CG network and run LAMMPS stage-1 minimize."""
    assert SAMPLE_NODES.exists(), f"sample graph missing: {SAMPLE_NODES}"
    assert SAMPLE_EDGES.exists(), f"sample graph missing: {SAMPLE_EDGES}"

    Pipeline(smoke_cg_config).run()

    study_dir = tmp_path / "smoke_cg"
    assert (study_dir / "02_Chemistry" / "system.data").exists(), (
        "stage 4 chemistry output missing"
    )

    sim_dir = study_dir / "04_Simulation"
    minimize_in = sim_dir / "minimize_1_serial.in"
    assert minimize_in.exists(), "stage 6 LAMMPS input script missing"

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
    assert (sim_dir / "system_after_soft.data").exists(), (
        "LAMMPS exited 0 but did not write stage-1 output data file"
    )
