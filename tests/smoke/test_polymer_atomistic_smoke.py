"""Smoke test: small atomistic network through Pipeline + LAMMPS stage 1.

Loads a pre-generated 5x5x5 SC topology from tests/sample_graphs/, runs the
full 6-stage Pipeline at DP=5 atomistic (DREIDING + PDMS), then invokes
LAMMPS to run the stage-1 minimize script. Asserts the pipeline emits
the expected files and LAMMPS exits 0.

CURRENTLY XFAIL — the polymer Pipeline path is blocked by three
pre-existing bugs (P0-B / P0-C / P0-D in
internal/DEVELOPMENT_INTERNAL.md sec.1). This test is wired and ready;
remove the xfail marker once the upstream bugs are fixed.

P0-B blocks `topology.source="generate"` (run_generator signature mismatch
+ no Python-only path). Worked around here with `source="load"`.
P0-C blocks `chemistry.model_type="coarse_grained"` (literal mismatch in
the writer). Worked around here with `"atomistic"`.
P0-D is the active blocker: a `'>=' tuple vs int` comparison crashes the
chemistry/conformation handoff after Stage 4 finishes writing displacement
files. Root cause not yet traced; the smoke test reproduces it
deterministically with a 5x5x5 + DP=5 atomistic config.

Why a smoke test matters: pure-Python unit tests can't catch regressions
where the pipeline emits a syntactically valid LAMMPS file that LAMMPS
nevertheless rejects. This is the cheapest end-to-end check that the
Python -> LAMMPS handoff still works.
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


pytestmark = [
    pytest.mark.requires_lammps,
    pytest.mark.xfail(
        reason=(
            "Pipeline blocked by P0-D (TypeError: '>=' tuple vs int after Stage 4) "
            "+ P0-B + P0-C; see internal/DEVELOPMENT_INTERNAL.md sec.1"
        ),
        strict=False,
    ),
]


REPO_ROOT = Path(__file__).resolve().parents[2]
SAMPLE_NODES = REPO_ROOT / "tests" / "sample_graphs" / "network_N5x5x5_trial3.nodes"
SAMPLE_EDGES = REPO_ROOT / "tests" / "sample_graphs" / "network_N5x5x5_trial3.edges"


@pytest.fixture
def smoke_atomistic_config(tmp_path: Path) -> ToponConfig:
    """Minimal-scale atomistic config: 5x5x5 SC sample graph, DP=5."""
    return ToponConfig(
        study=StudyConfig(name="smoke_atomistic", output_dir=str(tmp_path)),
        topology=TopologyConfig(
            source="load",
            existing_files=ExistingFilesConfig(
                nodes_file=str(SAMPLE_NODES),
                edges_file=str(SAMPLE_EDGES),
            ),
        ),
        assignment=AssignmentConfig(
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=5.0)),
        ),
        chemistry=ChemistryConfig(model_type="atomistic", target_density=0.9),
        output=OutputConfig(),
    )


def test_atomistic_pipeline_runs_and_lammps_minimizes(
    smoke_atomistic_config: ToponConfig, tmp_path: Path
) -> None:
    """Generate small atomistic network and run LAMMPS stage-1 minimize."""
    assert SAMPLE_NODES.exists(), f"sample graph missing: {SAMPLE_NODES}"
    assert SAMPLE_EDGES.exists(), f"sample graph missing: {SAMPLE_EDGES}"

    Pipeline(smoke_atomistic_config).run()

    study_dir = tmp_path / "smoke_atomistic"
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
