"""MARTINI 3 protein-network demo — small dry resilin reference.

Builds a 4 chains x 18 repeats (= 270 residues, vendored nat_pro ITP)
MARTINI 3 protein network through the full topon pipeline:
    BFM topology -> MARTINI chemistry -> LAMMPS data + 3 stage scripts.

This is the smallest functioning resilin cell: 4 chains is enough for
the BFM lattice to find Y-Y adjacencies at the v41/v42-default packing
(0.45), and gel happens at low conversion. ~3 seconds wall + write.

Outputs (under expected_output/):
    protein_network_topology.json   BFM snapshots (gel_point, post_gel_*)
    protein_network.data            LAMMPS atom_style full data file
    protein_network.in.settings     pair / bond / angle / dihedral coeffs
    protein_network.in.groups       protein / water / ions groups
    relaxation/
      protein_network_stage1.in     soft-push overlap removal
      protein_network_stage2.in     LJ epsilon ramp
      protein_network_stage3.in     tight CG min + brief NVT/NPT @ 310 K

After stage 3, system_equilibrated.data lands in the parent folder.

Usage:
    python examples/demos/protein/martini/run.py
"""
from __future__ import annotations

import time
from pathlib import Path

from topon.protein_network.workflow import run_protein_network


OUTPUT_DIR = Path(__file__).parent / "expected_output"


def main() -> None:
    print("--- MARTINI 3 dry resilin (8 chains, 18 repeats = 270 residues) ---")
    print("    (v41/v42 reference; gels at conv=0.125, fully percolated)")
    t0 = time.perf_counter()
    paths = run_protein_network(
        block_seq="GGRPSDSYGAPGGGN",
        n_chains=8,
        n_repeats=18,  # vendored nat_pro ITP locks this to 18
        output_dir=str(OUTPUT_DIR),
        snapshot_label="gel_point",
        segs_per_block=2,
        target_packing=0.45,
        equil_steps=20_000,
        n_extra_snapshots=2,
        snapshot_delta_conv=0.05,
        water_density_w_per_nm3=0.0,   # dry
        seed=42,
        hierarchical_stage1=True,
        verbose=False,
    )
    dt = time.perf_counter() - t0
    print(f"  built in {dt:.1f} s")

    print("\nArtifacts:")
    for k, p in paths.items():
        print(f"  {k:25s} {Path(p).resolve()}")

    relax = Path(paths["stage1"]).parent
    print(f"\nTo run LAMMPS:")
    print(f"  cd {relax}")
    print(f"  lmp -in protein_network_stage1.in   # ~3 s soft overlap removal")
    print(f"  lmp -in protein_network_stage2.in   # ~30 s LJ epsilon ramp")
    print(f"  lmp -in protein_network_stage3.in   # ~5 s CG min + brief NVT/NPT")


if __name__ == "__main__":
    main()
