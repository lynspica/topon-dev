"""Smoke test: CHARMM atomistic protein-network builder + LAMMPS stage 1.

Generates a small BFM topology (8 chains x 8 repeats), builds the dry
atomistic system through `topon.protein_network.charmm.build_systems`,
and runs LAMMPS stage 1 (soft overlap removal). Stages 2 and 3 are not
exercised here — they take O(minutes) per stage even at this size and
are out of scope for a smoke test.

The CHARMM36m PRM/RTF/CMAP files bundled with the package are used; if
that data dir gets relocated, this test will fail at FF-load time.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from topon.protein_network.bfm import generate_topology
from topon.protein_network.charmm.topology_io import save_topology


pytestmark = [pytest.mark.requires_lammps]


def test_charmm_build_and_stage1_minimize(tmp_path: Path) -> None:
    """Build a small CHARMM36m system and run LAMMPS stage 1."""
    # 1. Topology
    topo = generate_topology(
        n_chains=8, n_repeats=8, segs_per_block=2,
        equil_steps=10000, n_extra_snapshots=0, snapshot_delta_conv=0.05,
        seed=42, verbose=False,
    )
    topo_path = tmp_path / "topo.json"
    save_topology(topo, str(topo_path))

    # 2. Build LAMMPS files via the CLI entry point (covers the real path
    # users will run, including PRM/RTF/CMAP data resolution).
    sys_root = tmp_path / "sys"
    proc = subprocess.run(
        [
            "python", "-m", "topon.protein_network.charmm.build_systems",
            "--topology", str(topo_path),
            "--snapshot", "0",
            "--n_repeats", "8",
            "--water_contents", "0",
            "--output", str(sys_root),
        ],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, (
        f"build_systems failed:\n--- stdout ---\n{proc.stdout[-2000:]}\n"
        f"--- stderr ---\n{proc.stderr[-1000:]}"
    )

    # 3. Sanity-check outputs
    w0 = sys_root / "w0"
    assert (w0 / "protein_network.data").exists()
    assert (w0 / "protein_network.in.settings").exists()
    assert (w0 / "protein_network.in.groups").exists()
    assert (w0 / "charmm36m.cmap").exists()
    relax = w0 / "relaxation"
    assert (relax / "protein_network_stage1.in").exists()
    assert (relax / "protein_network_stage2.in").exists()
    assert (relax / "protein_network_stage3.in").exists()

    # 4. Run LAMMPS stage 1
    result = subprocess.run(
        ["lmp", "-in", "protein_network_stage1.in"],
        cwd=relax, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, (
        f"LAMMPS stage-1 failed (exit {result.returncode}):\n"
        f"--- stdout ---\n{result.stdout[-2000:]}\n"
        f"--- stderr ---\n{result.stderr[-1000:]}"
    )
    # stage1 writes system_after_soft.data into the relaxation/ dir
    assert (relax / "system_after_soft.data").exists()
