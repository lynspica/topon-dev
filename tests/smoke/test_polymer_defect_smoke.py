"""Smoke test: atomistic network with primary-loop defects + LAMMPS stage 1.

5 primary loops (parallel edges) injected on a 5x5x5 atomistic network
at DP=5. Runs Pipeline through stage-1 LAMMPS minimize. Pins the defect
injection + valence-protection paths (V21.2).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from topon.config.schema import (
    AssignmentConfig,
    ChemistryConfig,
    DefectsConfig,
    DPConfig,
    DPDistributionConfig,
    ExistingFilesConfig,
    OutputConfig,
    StudyConfig,
    TargetConfig,
    ToponConfig,
    TopologyConfig,
)
from topon.pipeline import Pipeline


pytestmark = [pytest.mark.requires_lammps]


REPO_ROOT = Path(__file__).resolve().parents[2]
SAMPLE_NODES = REPO_ROOT / "tests" / "sample_graphs" / "network_N5x5x5_trial3.nodes"
SAMPLE_EDGES = REPO_ROOT / "tests" / "sample_graphs" / "network_N5x5x5_trial3.edges"


@pytest.fixture
def smoke_defect_config(tmp_path: Path) -> ToponConfig:
    """5x5x5 atomistic with 5 primary-loop defects at DP=5."""
    return ToponConfig(
        study=StudyConfig(name="smoke_defect", output_dir=str(tmp_path)),
        topology=TopologyConfig(
            source="load",
            existing_files=ExistingFilesConfig(
                nodes_file=str(SAMPLE_NODES),
                edges_file=str(SAMPLE_EDGES),
            ),
        ),
        assignment=AssignmentConfig(
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=5.0)),
            defects=DefectsConfig(
                primary_loops=TargetConfig(
                    enabled=True, target=5, target_type="count"
                ),
            ),
        ),
        chemistry=ChemistryConfig(model_type="atomistic", target_density=0.9),
        output=OutputConfig(),
    )


def test_defect_pipeline_runs_and_lammps_minimizes(
    smoke_defect_config: ToponConfig, tmp_path: Path
) -> None:
    """Build atomistic defected network and run LAMMPS stage-1 minimize."""
    Pipeline(smoke_defect_config).run()

    study_dir = tmp_path / "smoke_defect"
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
