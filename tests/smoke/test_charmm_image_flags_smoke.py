"""Smoke test: CHARMM atomistic writer emits MPI-safe image flags.

Builds a small GGRPSDSYGAPGGGN (resilin consensus) system through the real
`build_systems` CLI and verifies the 2026-05-19 image-flag fix:

  * Atoms rows are 10-column ``id mol type q x y z ix iy iz``.
  * Some atoms carry non-zero image flags (i.e. chains really do wrap the
    box, so the test actually exercises the fix).
  * Every emitted bond is minimum-image once unwrapped via its image flags
    — no phantom ~box-length bonds that break parallel-MPI ghost shells
    (the `bond atoms missing` failure mode).
  * Net charge is integer and ~zero.

The build + invariant checks run without LAMMPS. A final gated check
(`requires_lammps`) confirms LAMMPS actually reads the 10-column file and
runs stage 1.

Why this matters: the CHARMM writer used to emit 7-column wrapped coords
with no image flags. For a percolated/crosslinked network run under MPI,
a bond whose atoms wrap to opposite box faces looks ~box-long and the
ghost-shell can't reconstruct it. This test guards the port of the
priority-MST image-flag pass (already proven on the MARTINI writer) into
the CHARMM path.
"""
from __future__ import annotations

import math
import re
import subprocess
import sys
from pathlib import Path

import pytest

from topon.protein_network.bfm import generate_topology
from topon.protein_network.charmm.topology_io import save_topology

BLOCK_SEQ = "GGRPSDSYGAPGGGN"
N_REPEATS = 8


@pytest.fixture(scope="module")
def built_data(tmp_path_factory) -> Path:
    """Build a small GGRPSDSYGAPGGGN CHARMM system once (dry, w0)."""
    tmp = tmp_path_factory.mktemp("charmm_imgflags")
    topo = generate_topology(
        n_chains=8, n_repeats=N_REPEATS, segs_per_block=2,
        equil_steps=10_000, n_extra_snapshots=0, seed=42, verbose=False,
    )
    topo_path = tmp / "topo.json"
    save_topology(topo, str(topo_path))
    out = tmp / "sys"
    proc = subprocess.run(
        [sys.executable, "-m", "topon.protein_network.charmm.build_systems",
         "--topology", str(topo_path), "--snapshot", "0",
         "--block_seq", BLOCK_SEQ, "--n_repeats", str(N_REPEATS),
         "--water_contents", "0", "--output", str(out)],
        capture_output=True, text=True, timeout=300,
    )
    assert proc.returncode == 0, (
        f"build_systems failed:\n{proc.stdout[-1500:]}\n{proc.stderr[-1500:]}"
    )
    data = out / "w0" / "protein_network.data"
    assert data.exists()
    return data


def _parse_atoms_bonds(text: str):
    L = [float(re.search(rf"(\S+) (\S+) {ax}lo {ax}hi", text).group(2))
         for ax in "xyz"]
    am = re.search(r"Atoms.*?\n\n(.*?)\n\n", text, flags=re.S)
    rows = [l for l in am.group(1).splitlines()
            if l.split() and l.split()[0].isdigit()]
    pos, img, q = {}, {}, 0.0
    for l in rows:
        t = l.split("#")[0].split()
        aid = int(t[0])
        q += float(t[3])
        pos[aid] = (float(t[4]), float(t[5]), float(t[6]))
        if len(t) >= 10:
            img[aid] = (int(t[7]), int(t[8]), int(t[9]))
    bm = re.search(r"Bonds.*?\n\n(.*?)\n\n", text, flags=re.S)
    bonds = [(int(t[2]), int(t[3]))
             for l in bm.group(1).splitlines()
             if len(t := l.split()) >= 4 and t[0].isdigit()]
    return L, pos, img, q, rows, bonds


def test_atoms_section_is_10_column_with_image_flags(built_data):
    text = built_data.read_text(encoding="utf-8")
    _, _, img, _, rows, _ = _parse_atoms_bonds(text)
    ncol = len(rows[0].split("#")[0].split())
    assert ncol == 10, f"expected 10-column Atoms (with ix iy iz), got {ncol}"
    # The fix is only meaningful if chains actually wrap the box.
    assert any(any(c != 0 for c in flag) for flag in img.values()), (
        "no non-zero image flags — test not exercising the wrap case; "
        "shrink the box / raise n_repeats"
    )


def test_every_bond_is_minimum_image_under_image_flags(built_data):
    text = built_data.read_text(encoding="utf-8")
    L, pos, img, _, _, bonds = _parse_atoms_bonds(text)
    max_r = 0.0
    for a, b in bonds:
        du = [(pos[a][i] + img[a][i] * L[i]) - (pos[b][i] + img[b][i] * L[i])
              for i in range(3)]
        max_r = max(max_r, math.sqrt(sum(x * x for x in du)))
    # Every real CHARMM bond is < ~1.6 A; allow generous 6 A for the crude
    # initial placement. A phantom box-length bond (the bug) would be ~80 A.
    assert max_r < 6.0, (
        f"max unwrapped bond = {max_r:.1f} A — image flags are NOT making "
        f"bonds minimum-image (phantom box-length bond present)"
    )


def test_build_is_charge_neutral(built_data):
    text = built_data.read_text(encoding="utf-8")
    _, _, _, q, _, _ = _parse_atoms_bonds(text)
    assert abs(q - round(q)) < 1e-3, f"net charge {q:+.4f} not integer"
    assert abs(q) < 1e-3, f"net charge {q:+.4f} not neutral"


@pytest.mark.requires_lammps
def test_lammps_reads_10col_and_runs_stage1(built_data):
    relax = built_data.parent / "relaxation"
    r = subprocess.run(
        ["lmp", "-in", "protein_network_stage1.in"],
        cwd=relax, capture_output=True, text=True, timeout=300,
    )
    assert r.returncode == 0, (
        f"LAMMPS stage 1 failed on 10-column data:\n"
        f"{r.stdout[-1500:]}\n{r.stderr[-1000:]}"
    )
    assert (relax / "system_after_soft.data").exists()
