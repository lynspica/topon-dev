#!/usr/bin/env python
"""Build CHARMM36m atomistic LAMMPS systems from a BFM protein-network topology.

Forked from `legacy/subprojects/protein_network/topro/scripts/build_systems.py`.
Two integration changes from the legacy CLI:

1. Default PRM/RTF/CMAP files now resolve to the bundled
   `topon/protein_network/charmm/data/` directory rather than an
   absolute Windows path.
2. Topology JSON is interchangeable with `topon.protein_network.bfm`
   output — `python -m topon.protein_network.bfm` (legacy: gen_topology.py)
   produces the same schema, so a generate step is no longer required if
   you already have a topology file.

Output layout per water content::

    <output_dir>/
        w0/   protein_network.data
              protein_network.in.settings
              protein_network_stage1/2/3.in
        w35/  ...
        w55/  ...

CLI mirrors topro's `scripts/build_systems.py`. See the legacy README at
`subprojects/protein_network/topro/README.md` (gitignored) for the full
parameter reference.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

from .charmm_ff import CHARMMForceField
from .builder import (
    build_protein_system,
    add_water_and_ions,
    compute_lattice_scale,
    find_cmap_crossterms,
)
from .lammps_writer import (
    find_angles,
    find_dihedrals,
    build_type_maps,
    write_lammps_data,
    write_lammps_settings,
    write_lammps_groups,
    write_lammps_input,
)
from .topology_io import load_topology, get_snapshot, list_snapshots
from ..sequence import build_full_sequence, get_node_residue_mapping


_DATA = Path(__file__).resolve().parent / "data"
_DEFAULT_PRM = _DATA / "par_all36m_prot_C2L.prm"
_DEFAULT_RTF = _DATA / "top_all36m_prot_C2L.rtf"
_DEFAULT_CMAP = _DATA / "charmm36m.cmap"


def parse_args():
    p = argparse.ArgumentParser(
        description="Build CHARMM36m atomistic LAMMPS systems from a BFM topology.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--topology", "-t", required=True,
                   help="Path to topology JSON file (from topon.protein_network.bfm)")
    p.add_argument("--snapshot", default="gel_point",
                   help="Snapshot label or integer index to use")
    p.add_argument("--block_seq", default="GGRPSDSYGAPGGGN",
                   help="One-letter repeat-block sequence")
    p.add_argument("--n_repeats", type=int, default=None,
                   help="Override number of repeats (default: from topology)")
    p.add_argument("--charmm_prm", default=str(_DEFAULT_PRM),
                   help="CHARMM PRM file path")
    p.add_argument("--charmm_rtf", default=str(_DEFAULT_RTF),
                   help="CHARMM RTF file path")
    p.add_argument("--charmm_cmap", default=str(_DEFAULT_CMAP),
                   help="CHARMM CMAP file path (copied next to the .data files)")
    p.add_argument("--water_contents", default="0,35,55,65,75",
                   help="Comma-separated weight-percent water values")
    p.add_argument("--salt_conc", type=float, default=0.15,
                   help="NaCl background concentration in mol/L")
    p.add_argument("--lattice_scale", type=float, default=None,
                   help="Override lattice scale in A/BFM unit (default: auto)")
    p.add_argument("--target_density", type=float, default=0.85,
                   help="Target initial density in g/cm^3 for box sizing")
    p.add_argument("--output", "-o", default="output",
                   help="Root output directory")
    p.add_argument("--prefix", default="protein_network",
                   help="File name prefix for LAMMPS output files")
    return p.parse_args()


def main():
    args = parse_args()

    print("=" * 60)
    print("  topon.protein_network.charmm — Atomistic system builder")
    print("=" * 60)

    print(f"\n[1] Loading topology: {args.topology}")
    topo = load_topology(args.topology)
    list_snapshots(topo)

    snap_label = args.snapshot
    try:
        snap_label = int(snap_label)
    except ValueError:
        pass
    snapshot = get_snapshot(topo, snap_label)
    print(f"\n    Using snapshot: '{snapshot['label']}'  conv={snapshot['conv']:.4f}")

    cfg = topo["config"]
    n_repeats = args.n_repeats or cfg["n_repeats"]
    segs_per_block = cfg["segs_per_block"]
    y_offset = cfg["y_offset_in_block"]

    print(f"\n[2] Parsing CHARMM36m force field ...")
    ff = CHARMMForceField(args.charmm_prm, args.charmm_rtf)
    print(f"    {len(ff.masses)} atom types | {len(ff.bonds_prm)} bonds | "
          f"{len(ff.angles_prm)} angles | {len(ff.dihedrals_prm)} dihedrals")
    for patch in ["NTER", "GLYP", "CTER", "DITY"]:
        tag = "[OK]" if patch in ff.patches else "[MISSING]"
        print(f"    {tag} patch {patch}")

    print(f"\n[3] Building sequence ({n_repeats}x '{args.block_seq}') ...")
    full_seq = build_full_sequence(args.block_seq, n_repeats)
    node_to_res = get_node_residue_mapping(
        n_repeats, segs_per_block,
        y_offset_in_block=y_offset,
        block_seq=args.block_seq,
    )
    print(f"    {len(full_seq)} residues | {len(node_to_res)} node mappings")

    print(f"\n[4] Building dry protein system ...")
    Nx = snapshot["Nx"]
    dry_scale = args.lattice_scale or 15.0
    dry_atoms, dry_bonds, dry_impropers, dry_box, xlinks = build_protein_system(
        ff, snapshot, full_seq, node_to_res, lattice_scale=dry_scale
    )
    print(f"    Atoms: {len(dry_atoms)} | Bonds: {len(dry_bonds)} | "
          f"Crosslinks: {len(xlinks)}")

    from .builder import _estimate_protein_mass
    protein_mass = _estimate_protein_mass(dry_atoms)
    total_charge = round(sum(a.charge for a in dry_atoms))
    print(f"    Protein mass: {protein_mass:.0f} Da  Net charge: {total_charge:+d} e")

    water_contents = [float(w) for w in args.water_contents.split(",")]

    for wc in water_contents:
        label = f"w{int(wc)}"
        out_dir = Path(args.output) / label
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = str(out_dir / args.prefix)

        print(f"\n{'-' * 55}")
        print(f"  Water content: {wc:.0f} wt%  ->  {out_dir}")
        print(f"{'-' * 55}")

        if args.lattice_scale is not None:
            scale = args.lattice_scale
        else:
            scale = compute_lattice_scale(
                Nx, protein_mass, wc, target_density=args.target_density
            )
        print(f"  Lattice scale: {scale:.3f} A/unit | box ~= {Nx*scale:.1f}^3 A")

        atoms, bonds, impropers, box, _ = build_protein_system(
            ff, snapshot, full_seq, node_to_res, lattice_scale=scale
        )
        add_water_and_ions(
            atoms, bonds, box,
            water_content_pct=wc,
            salt_conc_M=args.salt_conc,
        )

        print("  Computing angles and dihedrals ...")
        atom_idx_set = set(a.idx for a in atoms)
        angles = find_angles(bonds, atom_idx_set)
        dihedrals = find_dihedrals(bonds, atom_idx_set)
        print(f"  Angles: {len(angles)} | Dihedrals: {len(dihedrals)}")

        crossterms = find_cmap_crossterms(atoms)
        print(f"  CMAP crossterms: {len(crossterms)}")

        (atom_type_map, bond_type_map, angle_type_map,
         dihedral_type_map, improper_type_map) = build_type_maps(
            atoms, bonds, angles, dihedrals, impropers
        )

        cmap_filename = None
        if os.path.exists(args.charmm_cmap):
            cmap_dst = out_dir / os.path.basename(args.charmm_cmap)
            shutil.copy2(args.charmm_cmap, cmap_dst)
            cmap_filename = cmap_dst.name
        else:
            print(f"  [WARN] CMAP file not found: {args.charmm_cmap}")

        data_file = f"{prefix}.data"
        settings_file = f"{prefix}.in.settings"
        groups_file = f"{prefix}.in.groups"

        write_lammps_data(
            data_file, atoms, bonds, angles, dihedrals, impropers,
            atom_type_map, bond_type_map, angle_type_map,
            dihedral_type_map, improper_type_map, box, ff,
            crossterms=crossterms,
        )
        write_lammps_settings(
            settings_file, ff,
            atom_type_map, bond_type_map, angle_type_map,
            dihedral_type_map, improper_type_map,
        )
        write_lammps_groups(
            groups_file, atoms,
            atom_type_map, bond_type_map, angle_type_map,
        )
        write_lammps_input(
            prefix,
            os.path.basename(data_file),
            os.path.basename(settings_file),
            box,
            groups_file=os.path.basename(groups_file),
            cmap_file=cmap_filename,
        )

        final_charge = sum(a.charge for a in atoms)
        print(f"  [OK] {len(atoms)} atoms | final charge = {final_charge:.4f} e")
        print(f"  [OK] {data_file}")

    print(f"\n{'=' * 60}")
    print("  All systems built successfully.")
    print(f"  Output root: {Path(args.output).resolve()}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    sys.exit(main())
