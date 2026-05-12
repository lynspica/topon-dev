"""Build a CHARMM36m atomistic LAMMPS system from a BFM topology JSON.

Same shape as `bfm_to_martini/run.py` but using the atomistic CHARMM
chemistry stage (`topon.protein_network.charmm.build_systems`) instead
of the MARTINI builder.

A typical workflow:
  1. Generate topology JSONs with `bfm_gel_point_sweep/run.py`.
  2. Pick one that hit gel.
  3. Run this script with that JSON path + your sequence + a list of
     water contents (0/35/55/65/75 wt% is the common dietary swing).

Bundled CHARMM36m PRM/RTF/CMAP files in `topon/protein_network/charmm/data/`
are used by default. Override via the CLI flags inside `build_systems`
if you need a different force field.

Per water content, the output is a self-contained directory ready for
the 3-stage relaxation (~5 s + ~3 min + ~30 s per water content for
the small example here):

    <output>/sys/wXX/
        protein_network.data          # ~11 k atoms dry, ~19 k at 35 wt%
        protein_network.in.settings
        protein_network.in.groups
        charmm36m.cmap                # backbone CMAP correction grid
        relaxation/
            protein_network_stage1/2/3.in

Usage:
    python examples/workflows/bfm_to_charmm/run.py
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from topon.protein_network.bfm import generate_topology
from topon.protein_network.charmm.topology_io import save_topology


OUTPUT_DIR = Path(__file__).parent / "output"

# ---- knobs --------------------------------------------------------------

BLOCK_SEQ = "GGRPSDSYGAPGGGN"   # resilin repeat unit
N_CHAINS = 8
N_REPEATS = 8
SEGS_PER_BLOCK = 2
EQUIL_STEPS = 10_000
SEED = 42

WATER_CONTENTS = [0, 35]        # wt%; add 55, 65, 75 for a fuller swing
SALT_CONC = 0.15                # mol/L NaCl background
TARGET_DENSITY = 0.85           # g/cm^3 for auto lattice scaling
SNAPSHOT = 0                    # 0 = first snapshot (gel_point if reached)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"--- BFM -> CHARMM36m: {BLOCK_SEQ} x {N_REPEATS} x {N_CHAINS} chains ---")
    print(f"  Water contents: {WATER_CONTENTS} wt%")
    print(f"  Output: {OUTPUT_DIR.resolve()}")

    print("\n[1] Generating BFM topology...")
    topology = generate_topology(
        n_chains=N_CHAINS,
        n_repeats=N_REPEATS,
        segs_per_block=SEGS_PER_BLOCK,
        equil_steps=EQUIL_STEPS,
        n_extra_snapshots=2,
        seed=SEED,
        verbose=False,
    )
    print(f"  -> {len(topology['snapshots'])} snapshots: "
          f"{[s['label'] for s in topology['snapshots']]}")

    topo_path = OUTPUT_DIR / "topo.json"
    save_topology(topology, str(topo_path))
    print(f"  topology saved to {topo_path.relative_to(Path.cwd())}")

    print("\n[2] Building CHARMM36m atomistic LAMMPS systems...")
    sys_dir = OUTPUT_DIR / "sys"
    sys_dir.mkdir(exist_ok=True)
    cmd = [
        sys.executable, "-m", "topon.protein_network.charmm.build_systems",
        "--topology", str(topo_path),
        "--snapshot", str(SNAPSHOT),
        "--block_seq", BLOCK_SEQ,
        "--n_repeats", str(N_REPEATS),
        "--water_contents", ",".join(str(w) for w in WATER_CONTENTS),
        "--salt_conc", str(SALT_CONC),
        "--target_density", str(TARGET_DENSITY),
        "--output", str(sys_dir),
    ]
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        sys.exit(f"build_systems failed (rc={rc})")

    print(f"\nDone. Per-water-content output under {sys_dir.relative_to(Path.cwd())}/")
    for wc in WATER_CONTENTS:
        relax = sys_dir / f"w{int(wc)}" / "relaxation"
        if relax.exists():
            print(f"\n  cd {relax.relative_to(Path.cwd())}")
            print(f"  lmp -in protein_network_stage1.in   # ~5 s soft min")
            print(f"  lmp -in protein_network_stage2.in   # ~3 min LJ ramp")
            print(f"  lmp -in protein_network_stage3.in   # ~30 s tight min + NVT/NPT")


if __name__ == "__main__":
    main()
