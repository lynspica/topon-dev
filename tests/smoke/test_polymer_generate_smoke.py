"""Smoke test: Pipeline with `topology.source="generate"` + LAMMPS stage 1.

Generates a tiny 4x4x4 SC topology in-process via the pure-Python generator
(no C binary required), then runs the full Pipeline at DP=10 atomistic and
invokes LAMMPS to run the stage-1 minimize. Pins the P0-B fix that added
the Python-vs-C dispatch to `Pipeline._generate_topology`.

Workarounds for still-open Pipeline issues:
- Builds `ToponConfig` programmatically (P0-A blocks `load_config`).

Note: the 4x4x4 lattice is small enough that the Python generator usually
finds a valid graph within ~5-50 trials. `max_trials=2000` is a generous
budget for the trivial constraint we use here.
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
    GeneratorConfig,
    OutputConfig,
    StudyConfig,
    ToponConfig,
    TopologyConfig,
)
from topon.pipeline import Pipeline


pytestmark = [pytest.mark.requires_lammps]


@pytest.fixture
def smoke_generate_config(tmp_path: Path) -> ToponConfig:
    """Tiny generated 4x4x4 SC topology, atomistic, DP=10."""
    return ToponConfig(
        study=StudyConfig(name="smoke_generate", output_dir=str(tmp_path)),
        topology=TopologyConfig(
            source="generate",
            generator=GeneratorConfig(
                exe_path=None,                    # forces the Python path
                lattice_size="4x4x4",
                lattice_type="SC",
                periodicity="111",
                max_functionality=4,
                degree_distribution="0:0,1:0",   # trivial: every node has degree >= 2
                max_trials=2000,
                max_saves=1,
            ),
        ),
        assignment=AssignmentConfig(
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=10.0)),
        ),
        chemistry=ChemistryConfig(model_type="atomistic", target_density=0.9),
        output=OutputConfig(),
    )


def test_generate_pipeline_runs_and_lammps_minimizes(
    smoke_generate_config: ToponConfig, tmp_path: Path
) -> None:
    """Generate a 4x4x4 topology in-process and run LAMMPS stage-1 minimize."""
    Pipeline(smoke_generate_config).run()

    study_dir = tmp_path / "smoke_generate"
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
