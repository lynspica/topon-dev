"""Build a MARTINI 3 coarse-grained LAMMPS system from a BFM topology JSON.

A typical workflow:
  1. Generate (or sweep) topology JSONs with `bfm_gel_point_sweep/run.py`.
  2. Pick a topology JSON that hit gel at a reasonable conversion.
  3. Run this script with that JSON path + your sequence + water content.

The MARTINI builder lives in `topon.protein_network.workflow.run_protein_network`
(also reachable as the `python -m topon.protein_network generate` CLI).
This wrapper just shows the parameter surface and tees the call from
Python so you can drop it into a larger pipeline (e.g. sweep water
content for the same topology, or build several sequences against one
gel snapshot).

Output directory layout:
    <output>/
        protein_network_topology.json   # original BFM snapshots (re-saved)
        protein_network.data            # LAMMPS data file (atom_style full)
        protein_network.in.settings     # pair / bond / angle / dihedral coeffs
        protein_network.in.groups       # protein / water / ions groups
        relaxation/
            protein_network_stage1.in   # soft overlap removal
            protein_network_stage2.in   # LJ epsilon ramp
            protein_network_stage3.in   # tight CG min + brief NVT/NPT @ 310 K

After stage 3, `system_equilibrated.data` is the production-ready file.

Usage:
    python examples/workflows/bfm_to_martini/run.py
"""
from __future__ import annotations

from pathlib import Path

from topon.protein_network.bfm import generate_topology
from topon.protein_network.workflow import run_protein_network


# ---- knobs --------------------------------------------------------------

OUTPUT_DIR = Path(__file__).parent / "output"
BLOCK_SEQ = "GGRPSDSYGAPGGGN"   # resilin repeat unit (8 of these = 120 residues)
N_REPEATS = 6
N_CHAINS = 4
SEGS_PER_BLOCK = 2
EQUIL_STEPS = 5_000
TARGET_PACKING = 0.45
WATER_CONTENT = 4               # water density bucket (0-7); 0=dry
SEED = 42
SNAPSHOT_LABEL = "gel_point"    # or "post_gel_1", "post_gel_2", etc.


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"--- BFM -> MARTINI: {BLOCK_SEQ} x {N_REPEATS} x {N_CHAINS} chains ---")
    print(f"  Output: {OUTPUT_DIR.resolve()}")

    # 1. Generate the BFM topology (or load an existing one — see comment
    #    block below for the "load from disk" variant).
    print("\n[1] Generating BFM topology...")
    topology = generate_topology(
        n_chains=N_CHAINS,
        n_repeats=N_REPEATS,
        segs_per_block=SEGS_PER_BLOCK,
        target_packing=TARGET_PACKING,
        equil_steps=EQUIL_STEPS,
        seed=SEED,
        verbose=False,
    )
    print(f"  -> {len(topology['snapshots'])} snapshots: "
          f"{[s['label'] for s in topology['snapshots']]}")

    # 2. Run the MARTINI builder (chemistry + LAMMPS file emission).
    print("\n[2] Building MARTINI system + LAMMPS files...")
    paths = run_protein_network(
        block_seq=BLOCK_SEQ,
        n_repeats=N_REPEATS,
        n_chains=N_CHAINS,
        output_dir=str(OUTPUT_DIR),
        snapshot_label=SNAPSHOT_LABEL,
        segs_per_block=SEGS_PER_BLOCK,
        target_packing=TARGET_PACKING,
        equil_steps=EQUIL_STEPS,
        water_density=WATER_CONTENT,
        seed=SEED,
    )

    print(f"\nDone. Files:")
    for key, path in paths.items():
        print(f"  {key:25s} {Path(path).relative_to(Path.cwd())}")

    print(f"\nNext steps:")
    print(f"  cd {Path(paths['stage1']).parent}")
    print(f"  lmp -in protein_network_stage1.in   # soft overlap removal")
    print(f"  lmp -in protein_network_stage2.in   # LJ epsilon ramp")
    print(f"  lmp -in protein_network_stage3.in   # tight CG min + NVT/NPT")


# ---- variant: load an existing topology JSON --------------------------------
# If you already have a topology from bfm_gel_point_sweep/run.py, replace
# the generate_topology() call above with a no-op and pass the path:
#
#     from topon.protein_network.charmm.topology_io import load_topology
#     topology = load_topology("../bfm_gel_point_sweep/output/baseline.json")
#     # then pass --snapshot-from-disk to run_protein_network (or fork
#     # run_protein_network's internals — see topon/protein_network/workflow.py).
# -----------------------------------------------------------------------------


if __name__ == "__main__":
    main()
