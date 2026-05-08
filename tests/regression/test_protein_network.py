"""Regression: MARTINI protein-network generator (sequence -> LAMMPS).

Reference golden at: tests/output/v33_protein_network_resilin_dry/

Asserts:
  - Two consecutive runs with the same seed produce byte-identical files
    (data, settings, groups, in).
  - Regenerated output matches the frozen golden byte-for-byte (data, settings,
    groups, in). Topology JSON is compared structurally because it carries a
    `created` timestamp.

To regenerate the golden after an intentional change, run with the env var
`UPDATE_PN_GOLDEN=1`.

Workflow under test:
    topon.protein_network.workflow.run_protein_network(
        block_seq="GGRPSDSYGAPGGGN", n_repeats=2, n_chains=4,
        equil_steps=200, seed=42, water_density_w_per_nm3=0.0)
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent.parent
GOLDEN_DIR = ROOT / "tests" / "output" / "v33_protein_network_resilin_dry"
UPDATE = os.environ.get("UPDATE_PN_GOLDEN") == "1"

ARTIFACT_FILES = (
    "protein_network.data",
    "protein_network.in.settings",
    "protein_network.in.groups",
    "relaxation/protein_network_stage1.in",
    "relaxation/protein_network_stage2.in",
    "relaxation/protein_network_stage3.in",
)


def _sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _run_workflow(out_dir: Path) -> dict:
    from topon.protein_network.workflow import run_protein_network
    return run_protein_network(
        block_seq="GGRPSDSYGAPGGGN",
        n_repeats=2,
        n_chains=4,
        output_dir=out_dir,
        equil_steps=200,
        seed=42,
        water_density_w_per_nm3=0.0,
        verbose=False,
    )


def test_workflow_is_deterministic_for_same_seed(tmp_path: Path):
    a = tmp_path / "run_a"
    b = tmp_path / "run_b"
    _run_workflow(a)
    _run_workflow(b)
    for name in ARTIFACT_FILES:
        assert _sha256(a / name) == _sha256(b / name), (
            f"{name} differs between two runs of the same seed"
        )


def test_workflow_matches_frozen_golden(tmp_path: Path):
    out = tmp_path / "run"
    _run_workflow(out)

    if UPDATE or not (GOLDEN_DIR / ARTIFACT_FILES[0]).exists():
        GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
        for name in ARTIFACT_FILES:
            dst = GOLDEN_DIR / name
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(out / name, dst)
        topo_src = out / "protein_network_topology.json"
        if topo_src.exists():
            shutil.copy2(topo_src, GOLDEN_DIR / topo_src.name)
        if UPDATE:
            pytest.skip(f"Updated golden at {GOLDEN_DIR}; re-run to verify")

    for name in ARTIFACT_FILES:
        new_hash = _sha256(out / name)
        ref_hash = _sha256(GOLDEN_DIR / name)
        assert new_hash == ref_hash, (
            f"{name} drifted from frozen golden.\n"
            f"  new:   {out / name}\n"
            f"  ref:   {GOLDEN_DIR / name}\n"
            f"To accept this change deliberately, set UPDATE_PN_GOLDEN=1 and re-run."
        )


def test_topology_json_is_structurally_consistent(tmp_path: Path):
    """Topology JSON has a `created` timestamp -> compare keys/values minus that."""
    out = tmp_path / "run"
    _run_workflow(out)
    new = json.loads((out / "protein_network_topology.json").read_text())
    if not (GOLDEN_DIR / "protein_network_topology.json").exists():
        pytest.skip("topology JSON golden not present; run test_workflow_matches_frozen_golden first")
    ref = json.loads((GOLDEN_DIR / "protein_network_topology.json").read_text())
    new.pop("created", None)
    ref.pop("created", None)
    assert new == ref, "topology JSON drifted (config or snapshots changed)"
