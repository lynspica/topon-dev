"""Smoke test: CG network with entanglements + LAMMPS stage 1.

5 Gaussian-Kink entanglements on a 5x5x5 CG bead-spring network at DP=10.
Runs Pipeline through stage-1 LAMMPS minimize. Pins the entanglement
geometry placement (V15-V21).
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
    EntanglementsConfig,
    ExistingFilesConfig,
    KinkParams,
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
def smoke_entanglement_config(tmp_path: Path) -> ToponConfig:
    """5x5x5 CG with 5 Gaussian-Kink entanglements at DP=10."""
    return ToponConfig(
        study=StudyConfig(name="smoke_entanglement", output_dir=str(tmp_path)),
        topology=TopologyConfig(
            source="load",
            existing_files=ExistingFilesConfig(
                nodes_file=str(SAMPLE_NODES),
                edges_file=str(SAMPLE_EDGES),
            ),
        ),
        assignment=AssignmentConfig(
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=10.0)),
            entanglements=EntanglementsConfig(
                enabled=True,
                target=5,
                target_type="count",
                kink_params=KinkParams(overshoot=0.2, z_amp=0.5, sigma=0.15),
            ),
        ),
        chemistry=ChemistryConfig(model_type="coarse_grained", target_density=0.85),
        output=OutputConfig(),
    )


def test_entanglement_pipeline_runs_and_lammps_minimizes(
    smoke_entanglement_config: ToponConfig, tmp_path: Path
) -> None:
    """Build CG entangled network and run LAMMPS stage-1 minimize."""
    Pipeline(smoke_entanglement_config).run()

    study_dir = tmp_path / "smoke_entanglement"
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
