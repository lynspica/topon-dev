"""End-to-end runner for the CHARMM atomistic demo.

Reads `config.json` next to this file, runs the BFM topology generator, then
calls `topon.protein_network.charmm.build_systems` to produce LAMMPS files.
LAMMPS itself is intentionally not invoked — see README.md for the three
stage commands. This script just gets you to the point where `lmp -in
protein_network_stage1.in` works.

Usage::

    python examples/demos/protein/charmm/run.py [--output runs/charmm_demo]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from topon.protein_network.bfm import generate_topology
from topon.protein_network.charmm.topology_io import save_topology


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", default="runs/charmm_demo",
        help="Output root directory",
    )
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    cfg = json.loads((here / "config.json").read_text())
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"  CHARMM atomistic demo  ->  {out.resolve()}")
    print("=" * 60)

    print("\n[1/2] Generating BFM topology ...")
    topo = generate_topology(verbose=False, **cfg["topology"])
    topo_path = out / "topo.json"
    save_topology(topo, str(topo_path))
    print(f"  {len(topo['snapshots'])} snapshot(s): "
          f"{', '.join(s['label'] for s in topo['snapshots'])}")

    print("\n[2/2] Building atomistic LAMMPS systems ...")
    chem = cfg["chemistry"]
    cmd = [
        sys.executable, "-m", "topon.protein_network.charmm.build_systems",
        "--topology", str(topo_path),
        "--snapshot", str(chem["snapshot"]),
        "--block_seq", chem["block_seq"],
        "--n_repeats", str(cfg["topology"]["n_repeats"]),
        "--water_contents", ",".join(str(w) for w in chem["water_contents"]),
        "--salt_conc", str(chem["salt_conc_M"]),
        "--target_density", str(chem["target_density_g_cm3"]),
        "--output", str(out / "sys"),
    ]
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        sys.exit(f"build_systems failed with exit code {rc}")

    print(f"\nDone. To run LAMMPS:")
    for wc in chem["water_contents"]:
        relax = out / "sys" / f"w{int(wc)}" / "relaxation"
        print(f"  cd {relax}")
        print(f"  lmp -in protein_network_stage1.in")
        print(f"  lmp -in protein_network_stage2.in")
        print(f"  lmp -in protein_network_stage3.in")


if __name__ == "__main__":
    main()
