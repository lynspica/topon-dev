"""End-to-end MARTINI protein-network workflow.

Pulls together topology, chemistry, water packing, and LAMMPS writing in one
function. Mirrors the shape of topro's `build_systems.py` but for a single
output directory per call (no built-in water-content sweep -- the caller can
loop). Two-stage (separate topology JSON) is available by calling the lower-
level modules directly: `bfm.generate_topology` -> `topology_io.save_topology`
-> `build_protein_system` -> `write_lammps`.
"""
from __future__ import annotations

from pathlib import Path

from . import bfm, builder, ions, lammps_writer, residues, sequence, topology_io, water
from . import itp_template, template_builder
from .martini_ff import MartiniLibrary


def run_protein_network(
    block_seq: str = "GGRPSDSYGAPGGGN",
    n_repeats: int = 6,
    n_chains: int = 4,
    *,
    output_dir: str | Path = "protein_network_run",
    snapshot_label: str = "gel_point",
    snapshot_fallback_index: int = -1,
    segs_per_block: int = 2,
    target_packing: float = 0.45,
    equil_steps: int = 5_000,
    n_extra_snapshots: int = 4,
    snapshot_delta_conv: float = 0.05,
    min_intrachain_sep: int = 2,
    lattice_scale_ang: float | None = None,
    crosslink_method: str = "adjacent",
    max_crosslink_distance_ang: float | None = None,
    pre_gel_conversions: list[float] | None = None,
    sc_jitter_ang: float = 1.5,
    water_density_w_per_nm3: float = 0.0,
    water_content_pct: float | None = None,
    target_density_g_cm3: float = 0.85,
    water_exclusion_ang: float = 4.0,
    water_bead_type: str = "W",
    n_na_ions: int = 0,
    n_cl_ions: int = 0,
    ion_exclusion_ang: float = 4.0,
    seed: int = 42,
    base_name: str = "protein_network",
    save_topology_json: bool = True,
    hierarchical_stage1: bool = True,
    use_itp_template: bool = True,
    chain_itp_path: str | Path | None = None,
    library: MartiniLibrary | None = None,
    verbose: bool = True,
) -> dict[str, Path]:
    """Generate a MARTINI protein network and write LAMMPS input files.

    Returns a dict mapping artifact name -> path:
      ``topology_json``, ``data``, ``settings``, ``groups``,
      ``stage1``, ``stage2``, ``stage3`` (the three relaxation input scripts
      under ``relaxation/``).

    The workflow is deterministic given ``seed``: the same arguments produce
    bitwise-identical files across runs.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if library is None:
        library = MartiniLibrary.from_package_data()

    # 1a. Try to load the ITP template first -- it locks the chain length, so
    # if available we must use its n_repeats for the BFM topology.
    template = None
    if use_itp_template:
        try:
            if chain_itp_path:
                template = itp_template.load_chain_template(chain_itp_path)
            else:
                template = itp_template.load_vendored_template(block_seq)
            if verbose:
                print(f"[workflow] using ITP template: {template.name}, "
                      f"{template.n_atoms} atoms, {len(template.angles)} angles, "
                      f"{len(template.dihedrals_proper)} propers, "
                      f"{len(template.dihedrals_improper)} impropers")
        except KeyError as e:
            if verbose:
                print(f"[workflow] no vendored ITP for {block_seq!r}; "
                      f"falling back to canonical extraction. ({e})")

    bfm_n_repeats = n_repeats
    if template is not None:
        actual_n_repeats = template.n_residues // len(block_seq)
        if actual_n_repeats != n_repeats:
            if verbose:
                print(f"[workflow] overriding n_repeats={n_repeats} -> "
                      f"{actual_n_repeats} to match ITP template "
                      f"({template.n_residues} residues / {len(block_seq)} per block)")
            bfm_n_repeats = actual_n_repeats

    # 1b. BFM topology (using template-corrected n_repeats if applicable).
    topology = bfm.generate_topology(
        n_chains=n_chains,
        n_repeats=bfm_n_repeats,
        segs_per_block=segs_per_block,
        target_packing=target_packing,
        equil_steps=equil_steps,
        n_extra_snapshots=n_extra_snapshots,
        snapshot_delta_conv=snapshot_delta_conv,
        min_intrachain_sep=min_intrachain_sep,
        seed=seed,
        verbose=verbose,
        lattice_scale_ang=lattice_scale_ang,
        max_crosslink_distance_ang=max_crosslink_distance_ang,
        crosslink_method=crosslink_method,
        pre_gel_conversions=pre_gel_conversions,
    )

    artifacts: dict[str, Path] = {}
    if save_topology_json:
        topo_path = out / f"{base_name}_topology.json"
        topology_io.save_topology(topology, str(topo_path))
        artifacts["topology_json"] = topo_path

    # 2. Pick the snapshot we'll build chemistry on
    try:
        snapshot = topology_io.get_snapshot(topology, snapshot_label)
    except KeyError:
        snapshot = topology["snapshots"][snapshot_fallback_index]
        if verbose:
            print(
                f"[workflow] requested snapshot {snapshot_label!r} not found; "
                f"falling back to index {snapshot_fallback_index} "
                f"(label={snapshot['label']!r})"
            )

    # Optional: size the box for a target WEIGHT-PERCENT water content,
    # mirroring the CHARMM builder's `water_content_pct` semantics (water mass
    # as % of total mass: water_mass = protein_mass * pct/(100-pct)). When set,
    # this OVERRIDES lattice_scale_ang so the box holds protein + water at
    # `target_density_g_cm3`, and fixes the number of W beads packed so the
    # system actually has that weight fraction (the MARTINI density knob
    # `water_density_w_per_nm3` is ignored on this path). NPT then fine-tunes
    # the density; we don't force it.
    n_water_target: int | None = None
    if water_content_pct is not None:
        if not (0.0 <= water_content_pct < 100.0):
            raise ValueError(
                f"water_content_pct must be in [0, 100); got {water_content_pct}"
            )
        Nx = snapshot["Nx"]
        # Protein mass is geometry-independent: sum of per-chain bead masses.
        if template is not None:
            per_chain_mass = sum(library.get_mass(a.bead_type) for a in template.atoms)
        else:
            _seq3 = sequence.build_full_sequence(block_seq, n_repeats)
            per_chain_mass = sum(
                library.get_mass(bt)
                for rn in _seq3
                for (_an, bt, _q) in residues.RESIDUES[rn]["beads"]
            )
        protein_mass = per_chain_mass * len(snapshot["chains"])
        w_bead_mass = library.get_mass(water_bead_type)
        if water_content_pct > 0.0:
            water_mass = protein_mass * water_content_pct / (100.0 - water_content_pct)
            n_water_target = int(round(water_mass / w_bead_mass))
        else:
            water_mass = 0.0
            n_water_target = 0
        total_mass = protein_mass + water_mass
        # V[A^3] = M[Da] * 1.66054 (A^3/Da per g/cm^3) / target_density
        V_ang3 = total_mass * 1.66054 / target_density_g_cm3
        lattice_scale_ang = (V_ang3 ** (1.0 / 3.0)) / Nx
        if verbose:
            print(
                f"[workflow] water_content_pct={water_content_pct:.1f} wt%: "
                f"protein_mass={protein_mass:.3e} amu -> "
                f"n_water_target={n_water_target} {water_bead_type} beads; "
                f"target_density={target_density_g_cm3} g/cm^3 -> "
                f"lattice_scale={lattice_scale_ang:.3f} A/unit "
                f"(box {lattice_scale_ang * Nx:.1f} A)"
            )

    if template is not None:
        # ITP-template: replicate the polyply chain topology per BFM chain.
        sys_ = template_builder.build_protein_system_from_template(
            snapshot, template, library,
            block_seq=block_seq, lattice_scale_ang=lattice_scale_ang,
            sc_jitter_ang=sc_jitter_ang, seed=seed,
        )
    else:
        seq3 = sequence.build_full_sequence(block_seq, n_repeats)
        sys_ = builder.build_protein_system(
            snapshot, seq3, library,
            block_seq=block_seq, lattice_scale_ang=lattice_scale_ang,
            sc_jitter_ang=sc_jitter_ang, seed=seed,
        )
    if verbose:
        print(
            f"[workflow] built {sys_.n_atoms()} beads, {len(sys_.bonds)} bonds, "
            f"{len(sys_.constraints)} constraints, {len(sys_.angles)} angles, "
            f"{len(sys_.dihedrals)} dihedrals, total charge "
            f"{sys_.total_charge():+.4f}"
        )

    # 4a. Optional ion packing FIRST (so the rare ions get prime real-estate
    # before water fills the box). Matches typical GROMACS workflow where
    # `genion` adds ions before water-equilibration.
    if n_na_ions > 0 or n_cl_ions > 0:
        nna, ncl = ions.pack_ions(
            sys_, library,
            n_na=n_na_ions, n_cl=n_cl_ions,
            exclusion_radius_ang=ion_exclusion_ang,
            seed=seed + 2,
        )
        if verbose:
            print(f"[workflow] packed {nna} NA + {ncl} CL ions")

    # 4b. Water packing.
    if n_water_target is not None:
        # WEIGHT-PERCENT path: pack a fixed COUNT of W beads, uniformly spread
        # (shuffle) across the density-sized box, so the system's water mass is
        # `water_content_pct` % of the total. Candidate grid density is set
        # generously above the per-box target so exclusions near protein don't
        # starve the count.
        n_w = 0
        if n_water_target > 0:
            box_vol_nm3 = (lattice_scale_ang * snapshot["Nx"] / 10.0) ** 3
            cand_density = max(water.DEFAULT_W_DENSITY_NM3,
                               2.0 * n_water_target / max(box_vol_nm3, 1e-9))
            n_w = water.pack_water(
                sys_, library,
                density_w_per_nm3=cand_density,
                exclusion_radius_ang=water_exclusion_ang,
                seed=seed + 1, bead_type=water_bead_type,
                max_beads=n_water_target, shuffle=True,
            )
        if verbose:
            w_mass = library.get_mass(water_bead_type) * n_w
            p_mass = sum(library.get_mass(b.bead_type)
                         for b in sys_.beads if b.bead_type != water_bead_type)
            achieved = 100.0 * w_mass / max(w_mass + p_mass, 1e-9)
            print(f"[workflow] packed {n_w}/{n_water_target} {water_bead_type} "
                  f"beads -> {achieved:.1f} wt% water "
                  f"(target {water_content_pct:.1f} wt%)")
    elif water_density_w_per_nm3 > 0:
        # DENSITY path (legacy): fill the BFM box at a fixed number density.
        n_w = water.pack_water(
            sys_, library,
            density_w_per_nm3=water_density_w_per_nm3,
            exclusion_radius_ang=water_exclusion_ang,
            seed=seed + 1,
            bead_type=water_bead_type,
        )
        if verbose:
            ratio = water.WATER_BEAD_TYPES.get(water_bead_type, 4)
            print(f"[workflow] packed {n_w} {water_bead_type} water beads "
                  f"({ratio} H2O/bead, ~{n_w * ratio} H2O equivalent)")

    # 5. LAMMPS files
    paths = lammps_writer.write_lammps(
        sys_, library, out, base_name=base_name,
        hierarchical_stage1=hierarchical_stage1,
    )
    artifacts.update(paths)
    if verbose:
        print(f"[workflow] wrote {len(artifacts)} files to {out}")
    return artifacts
