"""Smoke test: POSS at chain end-caps (atomistic) + LAMMPS stage 1.

Loads the 5x5x5 sample graph, maps degree-1 (chain-cap) nodes to
POSS_AM0270, leaves degree-2/3 as bare Si, runs the full Pipeline at
DP=5 atomistic, and runs LAMMPS stage-1 minimize.

This mirrors the documented original workflow (legacy
`generate_atomistic_poss.py`): POSS at dangling chain ends, where each
POSS has exactly one chain attached. The internal-junction case (POSS
at degree-4 nodes) is broken on the periodic boundary — see P0-H in
internal/DEVELOPMENT_INTERNAL.md — and is intentionally NOT exercised
here.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from topon.config.schema import (
    AssignmentConfig,
    ChemistryConfig,
    DegreeNodeTypeConfig,
    DPConfig,
    DPDistributionConfig,
    EdgeChemistryConfig,
    EdgeTypesConfig,
    ExistingFilesConfig,
    MonomerConfig,
    NodeMoleculeConfig,
    NodeTypesConfig,
    OutputConfig,
    StudyConfig,
    ToponConfig,
    TopologyConfig,
    UniformEdgeConfig,
)
from topon.pipeline import Pipeline


pytestmark = [pytest.mark.requires_lammps]


REPO_ROOT = Path(__file__).resolve().parents[2]
SAMPLE_NODES = REPO_ROOT / "tests" / "sample_graphs" / "network_N5x5x5_trial3.nodes"
SAMPLE_EDGES = REPO_ROOT / "tests" / "sample_graphs" / "network_N5x5x5_trial3.edges"


@pytest.fixture
def smoke_poss_config(tmp_path: Path) -> ToponConfig:
    """5x5x5 atomistic with degree-1 (chain-cap) nodes -> POSS_AM0270."""
    return ToponConfig(
        study=StudyConfig(name="smoke_poss", output_dir=str(tmp_path)),
        topology=TopologyConfig(
            source="load",
            existing_files=ExistingFilesConfig(
                nodes_file=str(SAMPLE_NODES),
                edges_file=str(SAMPLE_EDGES),
            ),
        ),
        assignment=AssignmentConfig(
            node_types=NodeTypesConfig(
                method="degree",
                # POSS at degree-1 caps; bare Si elsewhere. Mirrors the
                # legacy generate_atomistic_poss.py pattern.
                degree=DegreeNodeTypeConfig(
                    mapping={"1": "POSS", "2": "A", "3": "A", "4": "A"}
                ),
            ),
            edge_types=EdgeTypesConfig(
                method="uniform",
                uniform=UniformEdgeConfig(type="A"),
            ),
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=5.0)),
        ),
        chemistry=ChemistryConfig(
            model_type="atomistic",
            target_density=1.1,
            node_type_map={
                "A": NodeMoleculeConfig(molecule="Si"),
                "POSS": NodeMoleculeConfig(
                    molecule="POSS_AM0270", is_end_cap=True
                ),
            },
            edge_type_map={"A": EdgeChemistryConfig(monomer="PDMS")},
            monomers={
                "PDMS": MonomerConfig(
                    smiles="[Si](C)(C)O", chain_head="Si", chain_tail="O"
                ),
            },
        ),
        output=OutputConfig(),
    )


def test_poss_pipeline_runs_and_lammps_minimizes(
    smoke_poss_config: ToponConfig, tmp_path: Path
) -> None:
    """Build POSS-junction atomistic system and run LAMMPS stage-1 minimize."""
    Pipeline(smoke_poss_config).run()

    study_dir = tmp_path / "smoke_poss"
    assert (study_dir / "02_Chemistry" / "system.data").exists()

    sim_dir = study_dir / "04_Simulation"
    minimize_in = sim_dir / "minimize_1_serial.in"
    assert minimize_in.exists()

    result = subprocess.run(
        ["lmp", "-in", "minimize_1_serial.in"],
        cwd=sim_dir,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"LAMMPS stage-1 minimize failed (exit {result.returncode}):\n"
        f"--- stdout ---\n{result.stdout[-2000:]}\n"
        f"--- stderr ---\n{result.stderr[-2000:]}"
    )
    assert (sim_dir / "system_after_soft.data").exists()
