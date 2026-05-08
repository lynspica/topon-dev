"""argparse CLI for the MARTINI protein-network generator.

Run via ``python -m topon.protein_network`` (which dispatches to this module).

Examples
--------
* Single run, dry::

    python -m topon.protein_network generate --block-seq GGRPSDSYGAPGGGN \\
        --n-repeats 6 --n-chains 4 --output runs/resilin_dry --seed 42

* Sweep multiple water contents (mirroring topro's wXX/ layout)::

    python -m topon.protein_network sweep --block-seq GGRPSDSYGAPGGGN \\
        --n-repeats 6 --n-chains 4 --water-densities 0,4,8,10 --output runs/resilin

* Two-stage: just generate the BFM topology JSON for inspection::

    python -m topon.protein_network topology --n-chains 16 --n-repeats 12 \\
        --output topo.json
"""
from __future__ import annotations

import argparse
from pathlib import Path

from . import bfm, topology_io, workflow


def _add_common_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--block-seq", default="GGRPSDSYGAPGGGN",
                   help="One-letter repeat block (default: resilin GGRPSDSYGAPGGGN).")
    p.add_argument("--n-repeats", type=int, default=6, help="Number of repeats per chain.")
    p.add_argument("--n-chains", type=int, default=4, help="Number of chains.")
    p.add_argument("--segs-per-block", type=int, default=2, help="BFM segments per repeat.")
    p.add_argument("--equil-steps", type=int, default=5_000,
                   help="Monte Carlo equilibration steps (0 = skip).")
    p.add_argument("--target-packing", type=float, default=0.45)
    p.add_argument("--min-intrachain-sep", type=int, default=2)
    p.add_argument("--lattice-scale-ang", type=float, default=None,
                   help="Angstroms per BFM lattice unit (default: auto-scaled "
                        "so BB-BB equilibrium length = MARTINI 3.6 A).")
    p.add_argument("--sc-jitter-ang", type=float, default=1.5,
                   help="Sidechain bead random offset magnitude (A).")
    p.add_argument("--snapshot-label", default="gel_point",
                   help="Which BFM snapshot to build chemistry from.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--quiet", action="store_true")
    p.add_argument("--hierarchical-stage1", action="store_true",
                   help="Use core-topon-style progressive freeze/unfreeze in stage 1 "
                        "(safer for the BFM-derived topology; mirrors "
                        "topon/writers/lammps_inputs.py:write_serial_soft_minimization).")


def _cmd_generate(args: argparse.Namespace) -> int:
    paths = workflow.run_protein_network(
        block_seq=args.block_seq,
        n_repeats=args.n_repeats,
        n_chains=args.n_chains,
        output_dir=args.output,
        snapshot_label=args.snapshot_label,
        segs_per_block=args.segs_per_block,
        equil_steps=args.equil_steps,
        target_packing=args.target_packing,
        min_intrachain_sep=args.min_intrachain_sep,
        lattice_scale_ang=args.lattice_scale_ang,
        sc_jitter_ang=args.sc_jitter_ang,
        water_density_w_per_nm3=args.water_density,
        water_exclusion_ang=args.water_exclusion,
        water_bead_type=getattr(args, "water_bead", "W"),
        n_na_ions=getattr(args, "n_na_ions", 0),
        n_cl_ions=getattr(args, "n_cl_ions", 0),
        seed=args.seed,
        hierarchical_stage1=getattr(args, "hierarchical_stage1", False),
        verbose=not args.quiet,
    )
    print("Files written:")
    for kind, p in paths.items():
        print(f"  [{kind}] {p}")
    return 0


def _cmd_sweep(args: argparse.Namespace) -> int:
    base_out = Path(args.output)
    base_out.mkdir(parents=True, exist_ok=True)
    densities = [float(x) for x in args.water_densities.split(",")]
    for d in densities:
        # Subfolder name encodes density in W/nm^3 (integer if whole, else decimal).
        label = f"w{int(d)}" if d == int(d) else f"w{d:g}".replace(".", "p")
        sub = base_out / label
        print(f"[sweep] water density {d:.2f} W/nm^3 -> {sub}")
        workflow.run_protein_network(
            block_seq=args.block_seq,
            n_repeats=args.n_repeats,
            n_chains=args.n_chains,
            output_dir=sub,
            snapshot_label=args.snapshot_label,
            segs_per_block=args.segs_per_block,
            equil_steps=args.equil_steps,
            target_packing=args.target_packing,
            min_intrachain_sep=args.min_intrachain_sep,
            lattice_scale_ang=args.lattice_scale_ang,
            sc_jitter_ang=args.sc_jitter_ang,
            water_density_w_per_nm3=d,
            water_exclusion_ang=args.water_exclusion,
            water_bead_type=getattr(args, "water_bead", "W"),
            n_na_ions=getattr(args, "n_na_ions", 0),
            n_cl_ions=getattr(args, "n_cl_ions", 0),
            seed=args.seed,
            hierarchical_stage1=getattr(args, "hierarchical_stage1", False),
            verbose=not args.quiet,
        )
    return 0


def _cmd_topology(args: argparse.Namespace) -> int:
    topo = bfm.generate_topology(
        n_chains=args.n_chains,
        n_repeats=args.n_repeats,
        segs_per_block=args.segs_per_block,
        equil_steps=args.equil_steps,
        target_packing=args.target_packing,
        min_intrachain_sep=args.min_intrachain_sep,
        seed=args.seed,
        hierarchical_stage1=getattr(args, "hierarchical_stage1", False),
        verbose=not args.quiet,
    )
    topology_io.save_topology(topo, str(args.output))
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="topon.protein_network",
        description="MARTINI 3 protein-network generator (sequence -> LAMMPS).",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("generate", help="Run topology + chemistry + LAMMPS write.")
    _add_common_args(g)
    g.add_argument("--output", required=True, help="Output directory.")
    g.add_argument("--water-density", type=float, default=0.0,
                   help="Water beads per nm^3 (0 = dry; ~10 for bulk MARTINI water = W bead default).")
    g.add_argument("--water-exclusion", type=float, default=4.0,
                   help="Min protein-water distance (A).")
    g.add_argument("--water-bead", default="W", choices=["W", "SW", "TW"],
                   help="MARTINI water bead type: W = 4 H2O/bead (default, bulk water), "
                        "SW = 3 H2O/bead (small, for confined water), "
                        "TW = 2 H2O/bead (tiny, for very tight pockets).")
    g.add_argument("--n-na-ions", type=int, default=0, help="NA+ ions to pack.")
    g.add_argument("--n-cl-ions", type=int, default=0, help="CL- ions to pack.")
    g.set_defaults(func=_cmd_generate)

    s = sub.add_parser("sweep", help="Run a water-content sweep into wXX/ subdirs.")
    _add_common_args(s)
    s.add_argument("--output", required=True, help="Base output directory.")
    s.add_argument("--water-densities", default="0,4,8",
                   help="Comma-separated water-bead densities to sweep (default: 0,4,8).")
    s.add_argument("--water-exclusion", type=float, default=4.0)
    s.add_argument("--water-bead", default="W", choices=["W", "SW", "TW"],
                   help="MARTINI water bead type per sweep point (W=4 H2O, SW=3, TW=2).")
    s.add_argument("--n-na-ions", type=int, default=0, help="NA+ ions to pack per density.")
    s.add_argument("--n-cl-ions", type=int, default=0, help="CL- ions to pack per density.")
    s.set_defaults(func=_cmd_sweep)

    t = sub.add_parser("topology", help="Generate just the BFM topology JSON.")
    _add_common_args(t)
    t.add_argument("--output", required=True, help="Output JSON path.")
    t.set_defaults(func=_cmd_topology)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
