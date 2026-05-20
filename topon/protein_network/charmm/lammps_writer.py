"""
topro.lammps.writer — Write LAMMPS data, settings, and 3-stage input files.

Public API
----------
find_angles(bonds, atom_idx_set)       → list of (i,j,k)
find_dihedrals(bonds, atom_idx_set)    → list of (i,j,k,l)
write_lammps_data(...)                 → old_to_new dict
write_lammps_settings(...)
write_lammps_groups(...)               → writes .in.groups include file
write_lammps_input(base_name, ...)     → writes stage1/2/3 .in files
build_type_maps(atoms, bonds, angles, dihedrals, impropers)
"""

import os
from collections import defaultdict


# ── Topology traversal ────────────────────────────────────────────────────────

def find_angles(bonds, atom_idx_set=None):
    """Find all 1-2-3 angle triplets from the bond list."""
    adj = defaultdict(set)
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)

    angles = []
    for center in adj:
        nbrs = sorted(adj[center])
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                angles.append((nbrs[i], center, nbrs[j]))
    return angles


def find_dihedrals(bonds, atom_idx_set=None):
    """Find all 1-2-3-4 dihedral quadruplets from the bond list."""
    adj = defaultdict(set)
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)

    dihedrals = []
    for a, b in bonds:
        for i in adj[a]:
            if i == b:
                continue
            for j in adj[b]:
                if j == a or j == i:
                    continue
                dihedrals.append((i, a, b, j))
    return dihedrals


# ── Type map builder ──────────────────────────────────────────────────────────

def build_type_maps(atoms, bonds, angles, dihedrals, impropers, ff=None):
    """
    Build integer type maps for all interaction kinds.

    ``dihedral_type_map`` maps a canonical atom-type quad (and its reverse)
    to a LIST of LAMMPS dihedral type IDs — one per CHARMM Fourier term.
    CHARMM proper dihedrals are sums of cosine terms (up to ~4 per quad,
    each with its own multiplicity); `dihedral_style charmm` represents
    each term as a separate dihedral type with the atom quad listed once
    per term in the Dihedrals section. ``ff`` is needed to count the terms;
    if omitted, every quad collapses to a single type (legacy behaviour,
    which silently dropped all but the first term).

    Returns
    -------
    atom_type_map, bond_type_map, angle_type_map,
    dihedral_type_map (quad -> [type_id, ...]), improper_type_map
    """
    idx_to_atom = {a.idx: a for a in atoms}

    # Atom types
    unique_atypes = sorted(set(a.atype for a in atoms))
    atom_type_map = {t: i + 1 for i, t in enumerate(unique_atypes)}

    # Bond types
    btype_set = set()
    for a1, a2 in bonds:
        at1 = idx_to_atom.get(a1)
        at2 = idx_to_atom.get(a2)
        if at1 and at2:
            btype_set.add(tuple(sorted([at1.atype, at2.atype])))
    bond_type_map = {k: i + 1 for i, k in enumerate(sorted(btype_set))}

    # Angle types (canonical = lexicographically smaller of (t1,t2,t3) vs (t3,t2,t1))
    atype_set = set()
    for a1, a2, a3 in angles:
        at1 = idx_to_atom.get(a1)
        at2 = idx_to_atom.get(a2)
        at3 = idx_to_atom.get(a3)
        if at1 and at2 and at3:
            k = (at1.atype, at2.atype, at3.atype)
            atype_set.add(min(k, (k[2], k[1], k[0])))
    angle_type_map = {k: i + 1 for i, k in enumerate(sorted(atype_set))}
    # Add reverse entries pointing to same type ID
    for k, v in list(angle_type_map.items()):
        angle_type_map[(k[2], k[1], k[0])] = v

    # Dihedral types — one LAMMPS type per CHARMM Fourier term per quad.
    dtype_set = set()
    for a1, a2, a3, a4 in dihedrals:
        at1 = idx_to_atom.get(a1)
        at2 = idx_to_atom.get(a2)
        at3 = idx_to_atom.get(a3)
        at4 = idx_to_atom.get(a4)
        if all([at1, at2, at3, at4]):
            k = (at1.atype, at2.atype, at3.atype, at4.atype)
            dtype_set.add(min(k, (k[3], k[2], k[1], k[0])))
    dihedral_type_map = {}
    next_tid = 1
    for k in sorted(dtype_set):
        n_terms = 1
        if ff is not None:
            prm = ff.lookup_dihedral(*k)
            if prm:
                n_terms = len(prm)
        tids = list(range(next_tid, next_tid + n_terms))
        next_tid += n_terms
        dihedral_type_map[k] = tids
        dihedral_type_map[(k[3], k[2], k[1], k[0])] = tids

    # Improper types
    itype_set = set()
    for a1, a2, a3, a4 in impropers:
        at1 = idx_to_atom.get(a1)
        at2 = idx_to_atom.get(a2)
        at3 = idx_to_atom.get(a3)
        at4 = idx_to_atom.get(a4)
        if all([at1, at2, at3, at4]):
            itype_set.add((at1.atype, at2.atype, at3.atype, at4.atype))
    improper_type_map = {k: i + 1 for i, k in enumerate(sorted(itype_set))}

    return atom_type_map, bond_type_map, angle_type_map, dihedral_type_map, improper_type_map


# ── LAMMPS data file ──────────────────────────────────────────────────────────

def write_lammps_data(filename, atoms, bonds, angles, dihedrals, impropers,
                      atom_type_map, bond_type_map, angle_type_map,
                      dihedral_type_map, improper_type_map, box, ff,
                      crossterms=None, image_flags=None):
    """
    Write a LAMMPS data file in 'full' atom style.

    All atom IDs are renumbered contiguously.  Bond/angle/dihedral/improper
    entries are pre-collected into lists so header counts are exact and IDs
    are always sequential even when some interactions are skipped.

    ``image_flags`` (optional): mapping ``atom.idx -> (ix, iy, iz)``. When
    given, each Atoms row is written with 10 columns
    ``(id mol type q x y z ix iy iz)`` instead of the legacy 7. This is
    required for percolated/crosslinked networks run under MPI: without
    image flags, a bond whose two atoms wrap to opposite box faces appears
    ~box-long in wrapped coords, which breaks parallel ghost-shell
    construction (``bond atoms missing``). Image flags are assigned by a
    priority-MST over the bond graph (see ``build_systems``); with them,
    every emitted bond is minimum-image. Omit (None) for the legacy
    7-column behaviour.

    Returns
    -------
    old_to_new : dict  (original atom.idx → new sequential 1-based ID)
    """
    old_to_new = {a.idx: new_id for new_id, a in enumerate(atoms, start=1)}

    def remap(idx):
        return old_to_new.get(idx)

    def atom_at(idx):
        nid = old_to_new.get(idx)
        return atoms[nid - 1] if nid else None

    # ── Pre-collect all interactions ──────────────────────────────────────────
    bond_rows = []
    for a1, a2 in bonds:
        na1, na2 = remap(a1), remap(a2)
        if na1 is None or na2 is None:
            continue
        at1, at2 = atoms[na1 - 1].atype, atoms[na2 - 1].atype
        bt = bond_type_map.get(tuple(sorted([at1, at2])), 1)
        bond_rows.append((bt, na1, na2))

    angle_rows = []
    for a1, a2, a3 in angles:
        na1, na2, na3 = remap(a1), remap(a2), remap(a3)
        if None in (na1, na2, na3):
            continue
        at1, at2, at3 = (atoms[na1-1].atype, atoms[na2-1].atype, atoms[na3-1].atype)
        k = (at1, at2, at3)
        at = angle_type_map.get(k) or angle_type_map.get((at3, at2, at1), 1)
        angle_rows.append((at, na1, na2, na3))

    dih_rows = []
    for a1, a2, a3, a4 in dihedrals:
        na1, na2, na3, na4 = remap(a1), remap(a2), remap(a3), remap(a4)
        if None in (na1, na2, na3, na4):
            continue
        at1, at2 = atoms[na1-1].atype, atoms[na2-1].atype
        at3, at4 = atoms[na3-1].atype, atoms[na4-1].atype
        k = (at1, at2, at3, at4)
        dts = (dihedral_type_map.get(k)
               or dihedral_type_map.get((at4, at3, at2, at1))
               or [1])
        # One Dihedrals row per CHARMM Fourier term of this quad.
        for dt in dts:
            dih_rows.append((dt, na1, na2, na3, na4))

    impr_rows = []
    for a1, a2, a3, a4 in impropers:
        na1, na2, na3, na4 = remap(a1), remap(a2), remap(a3), remap(a4)
        if None in (na1, na2, na3, na4):
            continue
        at1, at2 = atoms[na1-1].atype, atoms[na2-1].atype
        at3, at4 = atoms[na3-1].atype, atoms[na4-1].atype
        it = improper_type_map.get((at1, at2, at3, at4), 1)
        impr_rows.append((it, na1, na2, na3, na4))

    # ── Write file ────────────────────────────────────────────────────────────
    with open(filename, "w", encoding='utf-8') as f:
        f.write("LAMMPS Protein Network — CHARMM36/36m All-Atom (topro)\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(bond_rows)} bonds\n")
        f.write(f"{len(angle_rows)} angles\n")
        f.write(f"{len(dih_rows)} dihedrals\n")
        f.write(f"{len(impr_rows)} impropers\n")
        if crossterms:
            f.write(f"{len(crossterms)} crossterms\n")
        f.write("\n")
        f.write(f"{len(set(atom_type_map.values()))} atom types\n")
        f.write(f"{len(set(bond_type_map.values()))} bond types\n")
        f.write(f"{len(set(angle_type_map.values()))} angle types\n")
        n_dih_types = len({tid for tids in dihedral_type_map.values() for tid in tids})
        f.write(f"{n_dih_types} dihedral types\n")
        f.write(f"{len(set(improper_type_map.values()))} improper types\n\n")
        f.write(f"0.0 {box[0]:.4f} xlo xhi\n")
        f.write(f"0.0 {box[1]:.4f} ylo yhi\n")
        f.write(f"0.0 {box[2]:.4f} zlo zhi\n\n")

        f.write("Masses\n\n")
        for atype, tid in sorted(atom_type_map.items(), key=lambda x: x[1]):
            mass = ff.masses.get(atype, 12.011)
            f.write(f"{tid} {mass:.5f} # {atype}\n")

        f.write("\nAtoms # full\n\n")
        for a in atoms:
            nid = old_to_new[a.idx]
            tid = atom_type_map[a.atype]
            pos = a.pos % box
            if image_flags is not None:
                ix, iy, iz = image_flags.get(a.idx, (0, 0, 0))
                f.write(f"{nid} {a.mol_id} {tid} {a.charge:.6f} "
                        f"{pos[0]:.4f} {pos[1]:.4f} {pos[2]:.4f} "
                        f"{ix} {iy} {iz} "
                        f"# {a.res_name} {a.name}\n")
            else:
                f.write(f"{nid} {a.mol_id} {tid} {a.charge:.6f} "
                        f"{pos[0]:.4f} {pos[1]:.4f} {pos[2]:.4f} "
                        f"# {a.res_name} {a.name}\n")

        f.write("\nBonds\n\n")
        for i, (bt, na1, na2) in enumerate(bond_rows, 1):
            f.write(f"{i} {bt} {na1} {na2}\n")

        f.write("\nAngles\n\n")
        for i, (at, na1, na2, na3) in enumerate(angle_rows, 1):
            f.write(f"{i} {at} {na1} {na2} {na3}\n")

        f.write("\nDihedrals\n\n")
        for i, (dt, na1, na2, na3, na4) in enumerate(dih_rows, 1):
            f.write(f"{i} {dt} {na1} {na2} {na3} {na4}\n")

        f.write("\nImpropers\n\n")
        for i, (it, na1, na2, na3, na4) in enumerate(impr_rows, 1):
            f.write(f"{i} {it} {na1} {na2} {na3} {na4}\n")

        if crossterms:
            f.write("\nCMAP\n\n")
            ct_id = 0
            for (cmap_type, c_prev, n_i, ca_i, c_i, n_next) in crossterms:
                a1 = old_to_new.get(c_prev)
                a2 = old_to_new.get(n_i)
                a3 = old_to_new.get(ca_i)
                a4 = old_to_new.get(c_i)
                a5 = old_to_new.get(n_next)
                if None in (a1, a2, a3, a4, a5):
                    continue
                ct_id += 1
                # 5 unique backbone atoms: C(i-1) N(i) CA(i) C(i) N(i+1)
                f.write(f"{ct_id} {cmap_type} {a1} {a2} {a3} {a4} {a5}\n")

    return old_to_new


# ── LAMMPS settings file ──────────────────────────────────────────────────────

def write_lammps_settings(filename, ff, atom_type_map, bond_type_map,
                          angle_type_map, dihedral_type_map, improper_type_map):
    """Write LAMMPS force-field coefficients file."""
    with open(filename, "w", encoding='utf-8') as f:
        f.write("# CHARMM36/36m Force Field Coefficients (topro)\n\n")

        f.write("# Pair Coeffs  epsilon [kcal/mol]  sigma [Å]\n")
        for atype, tid in sorted(atom_type_map.items(), key=lambda x: x[1]):
            eps, rmin2 = ff.vdw_prm.get(atype, (-0.01, 1.0))
            eps_pos = abs(eps)
            rmin = 2.0 * abs(rmin2)
            sigma = rmin / (2.0 ** (1.0 / 6.0))
            f.write(f"pair_coeff {tid} {tid} {eps_pos:.6f} {sigma:.6f} # {atype}\n")

        f.write("\n# Bond Coeffs  K [kcal/mol/Å²]  r0 [Å]\n")
        for bkey, bt in sorted(bond_type_map.items(), key=lambda x: x[1]):
            prm = ff.lookup_bond(bkey[0], bkey[1])
            if prm:
                kb, b0 = prm
                f.write(f"bond_coeff {bt} {kb:.4f} {b0:.4f} # {bkey[0]}-{bkey[1]}\n")
            else:
                f.write(f"bond_coeff {bt} 300.0 1.5 # {bkey[0]}-{bkey[1]} (DEFAULT)\n")

        f.write("\n# Angle Coeffs  Ktheta  theta0  Kub  S0  (charmm style)\n")
        for akey, at in sorted(angle_type_map.items(), key=lambda x: x[1]):
            prm = ff.lookup_angle(akey[0], akey[1], akey[2])
            if prm:
                kth, th0, kub, s0 = prm
                f.write(f"angle_coeff {at} {kth:.4f} {th0:.4f} {kub:.4f} {s0:.4f} "
                        f"# {akey[0]}-{akey[1]}-{akey[2]}\n")
            else:
                f.write(f"angle_coeff {at} 50.0 109.5 0.0 0.0 "
                        f"# {akey[0]}-{akey[1]}-{akey[2]} (DEFAULT)\n")

        f.write("\n# Dihedral Coeffs  K  n  d  weight  (charmm style)\n")
        # dihedral_type_map maps each quad (and its reverse) to a LIST of
        # type IDs, one per CHARMM Fourier term. Emit one coeff line per
        # term; dedup the reverse-key alias via the first term's id.
        emitted: set[int] = set()
        for dkey, dts in sorted(dihedral_type_map.items(), key=lambda x: x[1][0]):
            if dts[0] in emitted:
                continue
            emitted.add(dts[0])
            label = f"{dkey[0]}-{dkey[1]}-{dkey[2]}-{dkey[3]}"
            prm = ff.lookup_dihedral(dkey[0], dkey[1], dkey[2], dkey[3])
            if prm:
                for (k, n, d), dt in zip(prm, dts):
                    f.write(f"dihedral_coeff {dt} {k:.4f} {n} {int(d)} 0.0 "
                            f"# {label}\n")
            else:
                for dt in dts:
                    f.write(f"dihedral_coeff {dt} 0.0 1 0 0.0 "
                            f"# {label} (DEFAULT)\n")

        f.write("\n# Improper Coeffs  K  chi0  (harmonic)\n")
        for ikey, it in sorted(improper_type_map.items(), key=lambda x: x[1]):
            prm = ff.lookup_improper(ikey[0], ikey[1], ikey[2], ikey[3])
            if prm:
                kpsi, psi0 = prm
                f.write(f"improper_coeff {it} {kpsi:.4f} {int(psi0)} "
                        f"# {ikey[0]}-{ikey[1]}-{ikey[2]}-{ikey[3]}\n")
            else:
                f.write(f"improper_coeff {it} 20.0 0 "
                        f"# {ikey[0]}-{ikey[1]}-{ikey[2]}-{ikey[3]} (DEFAULT)\n")


# ── LAMMPS groups include file ────────────────────────────────────────────────

def write_lammps_groups(filename, atoms, atom_type_map, bond_type_map, angle_type_map):
    """
    Write a LAMMPS include file that defines atom groups and SHAKE type variables.

    Groups defined
    --------------
    water    — OT + HT atom types
    ions     — SOD, CLA, CAL, ZN, FE3P (any present)
    protein  — everything that is not water or ions
    chain1 … chainN — one group per protein molecule ID

    Variables defined
    -----------------
    water_bond_type   — integer type ID for the OT-HT bond (for SHAKE)
    water_angle_type  — integer type ID for the HT-OT-HT angle (for SHAKE)
    """
    _WATER_TYPES = ("OT", "HT")
    _ION_TYPES   = ("SOD", "CLA", "CAL", "ZN", "FE3P")

    water_tids = sorted(atom_type_map[t] for t in _WATER_TYPES if t in atom_type_map)
    ion_tids   = sorted(atom_type_map[t] for t in _ION_TYPES   if t in atom_type_map)
    wi_tids    = set(water_tids) | set(ion_tids)

    # Protein molecule IDs: atoms whose type is NOT water or ion
    protein_mol_ids = sorted(set(
        a.mol_id for a in atoms if atom_type_map.get(a.atype, 0) not in wi_tids
    ))

    # Bond type ID for OT-HT (used for SHAKE)
    water_bond_tid  = bond_type_map.get(("HT", "OT")) or bond_type_map.get(("OT", "HT"))
    # Angle type ID for HT-OT-HT (used for SHAKE)
    water_angle_tid = angle_type_map.get(("HT", "OT", "HT"))

    with open(filename, "w", encoding="ascii") as f:
        f.write("# Group and SHAKE-type definitions -- auto-generated by topro\n")
        f.write("# Include this file after read_data in every LAMMPS input script.\n\n")

        f.write("# --- Water and ion groups ---\n")
        if water_tids:
            f.write(f"group           water   type {' '.join(str(t) for t in water_tids)}"
                    f"  # OT HT (TIP3P)\n")
        else:
            f.write("# (no water molecules in this system)\n")
        if ion_tids:
            ion_names = [t for t in _ION_TYPES if t in atom_type_map]
            f.write(f"group           ions    type {' '.join(str(t) for t in ion_tids)}"
                    f"  # {' '.join(ion_names)}\n")
        else:
            f.write("# (no ions in this system)\n")

        # Protein group: subtract known non-protein groups (or alias `all` when dry)
        subtract_list = (["water"] if water_tids else []) + (["ions"] if ion_tids else [])
        if subtract_list:
            f.write(f"group           protein subtract all {' '.join(subtract_list)}\n\n")
        else:
            f.write("group           protein union all\n\n")

        f.write("# --- Per-chain groups ---\n")
        for mid in protein_mol_ids:
            f.write(f"group           chain{mid:02d}  molecule {mid}\n")

        f.write("\n# --- SHAKE type IDs for water (fix shake water ...) ---\n")
        if water_bond_tid is not None:
            f.write(f"variable        water_bond_type  equal {water_bond_tid}"
                    f"  # OT-HT bond\n")
        if water_angle_tid is not None:
            f.write(f"variable        water_angle_type equal {water_angle_tid}"
                    f"  # HT-OT-HT angle\n")
        if water_bond_tid is not None:
            f.write(f"# SHAKE usage: fix s water shake 1e-4 100 0"
                    f" b ${{water_bond_type}} a ${{water_angle_type}}\n")

    return os.path.basename(filename)


# ── LAMMPS 3-stage input scripts ──────────────────────────────────────────────

def write_lammps_input(base_name, data_file, settings_file, box,
                       groups_file=None, cmap_file=None):
    """
    Write three LAMMPS input scripts for staged minimisation + equilibration.

    Files created:  {base_name}_stage1.in  (soft overlap removal)
                    {base_name}_stage2.in  (LJ epsilon ramp)
                    {base_name}_stage3.in  (tight minimise + NVT + NPT)

    Parameters
    ----------
    groups_file : str or None
        Basename of the .in.groups include file (from write_lammps_groups).
        When provided, each stage script includes it after read_data.
    cmap_file : str or None
        Basename of the CMAP correction file (e.g. charmm36m.cmap).
        When provided, fix cmap is activated before read_data in each stage.
    """
    # Stage scripts go in a relaxation/ subfolder; parent files referenced via ../
    out_dir      = os.path.dirname(os.path.abspath(base_name))
    relax_dir    = os.path.join(out_dir, "relaxation")
    os.makedirs(relax_dir, exist_ok=True)
    stage_prefix = os.path.join(relax_dir, os.path.basename(base_name))

    d = "../" + os.path.basename(data_file)
    s = "../" + os.path.basename(settings_file)
    g = ("../" + os.path.basename(groups_file)) if groups_file else None
    grp_line = f"include         {g}\n" if g else ""

    # CMAP lines: fix must be declared before read_data; read_data needs suffix
    if cmap_file:
        cm = "../" + os.path.basename(cmap_file)
        cmap_pre  = (f"fix             cmap all cmap {cm}\n"
                     f"fix_modify      cmap energy yes\n")
        cmap_rdat = " fix cmap crossterm CMAP"
    else:
        cmap_pre  = ""
        cmap_rdat = ""

    # ── Stage 1: Soft minimisation ────────────────────────────────────────────
    with open(f"{stage_prefix}_stage1.in", "w", encoding='utf-8') as f:
        f.write(f"""\
# Stage 1 — Serial Soft Minimisation (topro / CHARMM36/36m)
# Resolves all atomic overlaps using pair_style soft.
# Run from this relaxation/ directory.

units           real
atom_style      full
boundary        p p p
bond_style      harmonic
angle_style     charmm
dihedral_style  charmm
improper_style  harmonic
pair_style      lj/charmm/coul/long 10.0 12.0
kspace_style    pppm 1.0e-4

{cmap_pre}read_data       {d}{cmap_rdat}
include         {s}
{grp_line}
special_bonds   charmm
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

# ── Switch to soft potential (disable long-range electrostatics) ──────────────
kspace_style    none
pair_style      soft 1.0
pair_coeff      * * 0.0
variable        prefactor equal ramp(0,60)

thermo          100
thermo_style    custom step pe ke etotal evdwl ecoul epair ebond eangle edihed eimp press vol

# Stage A: ramped soft push
fix             soft_push all adapt 1 pair soft a * * v_prefactor
min_style       cg
minimize        1.0e-4 1.0e-6 1000 10000
unfix           soft_push
write_data      min_stage_A.data

# Stage B: nve/limit dynamics under soft potential
reset_timestep  0
timestep        1.0
fix             soft_push all adapt 1 pair soft a * * v_prefactor
fix             nve_limit all nve/limit 0.1
run             1000
unfix           nve_limit
unfix           soft_push

# Stage C: final soft minimisation
fix             soft_push all adapt 1 pair soft a * * v_prefactor
minimize        1.0e-4 1.0e-6 1000 10000
unfix           soft_push

write_data      system_after_soft.data
write_restart   1.restart
""")

    # ── Stage 2: LJ epsilon ramp ──────────────────────────────────────────────
    with open(f"{stage_prefix}_stage2.in", "w", encoding='utf-8') as f:
        f.write(f"""\
# Stage 2 — Parallel LJ Ramp (topro / CHARMM36/36m)
# Gradually scales LJ epsilon from 0.001 → 1.0 under nve/limit.
# Run from this relaxation/ directory.

units           real
atom_style      full
boundary        p p p
bond_style      harmonic
angle_style     charmm
dihedral_style  charmm
improper_style  harmonic
pair_style      soft 1.0

{cmap_pre}read_data       system_after_soft.data{cmap_rdat}
{grp_line}
neigh_modify    one 10000
comm_modify     mode single cutoff 14.0

# Pre-minimise under soft potential
pair_style      soft 1.0
pair_coeff      * * 1.0
min_style       cg
minimize        1.0e-4 1.0e-6 1000 10000

# Switch to real LJ/Coul (lj/cut/coul/long supports fix adapt epsilon)
pair_style      lj/cut/coul/long 12.0
kspace_style    pppm 1.0e-4
include         {s}
special_bonds   charmm

variable        scale equal ramp(0.001,1.0)
timestep        1.0

thermo          1000
thermo_style    custom step pe ke etotal evdwl ecoul epair ebond eangle edihed eimp press vol temp

fix             1 all adapt 1 pair lj/cut/coul/long epsilon * * v_scale scale yes
fix             fxnve all nve/limit 0.1
run             20000
unfix           fxnve
unfix           1

write_data      system_ramped.data
write_restart   2.restart
""")

    # ── Stage 3: Tight minimise + NVT + NPT ───────────────────────────────────
    with open(f"{stage_prefix}_stage3.in", "w", encoding='utf-8') as f:
        f.write(f"""\
# Stage 3 — Tight Minimisation + Short Equilibration (topro / CHARMM36/36m)
# Run from this relaxation/ directory.
# Output: ../system_equilibrated.data (in the parent wXX/ directory)

units           real
atom_style      full
boundary        p p p
bond_style      harmonic
angle_style     charmm
dihedral_style  charmm
improper_style  harmonic
pair_style      lj/charmm/coul/long 10.0 12.0

{cmap_pre}read_data       system_ramped.data{cmap_rdat}
{grp_line}
neigh_modify    one 10000
pair_style      lj/charmm/coul/long 10.0 12.0
kspace_style    pppm 1.0e-4
include         {s}
special_bonds   charmm

thermo          100
thermo_style    custom step pe ke etotal evdwl ecoul epair ebond eangle edihed eimp press vol temp

# Tight minimisation
min_style       cg
minimize        1.0e-6 1.0e-8 100000 1000000
write_data      system_minimized_final.data

# Short NVT (10 000 steps x 1 fs)
reset_timestep  0
variable        T equal 300
velocity        all create ${{T}} 12345
timestep        1.0
fix             1 all nvt temp ${{T}} ${{T}} 100.0
run             10000
unfix           1
write_data      after_nvt.data

# Short NPT (10 000 steps x 1 fs)
fix             1 all npt temp ${{T}} ${{T}} 100.0 iso 1.0 1.0 1000.0
run             10000
unfix           1

write_data      ../system_equilibrated.data
print "All stages complete."
""")

    print(f"  Stage scripts: relaxation/{os.path.basename(base_name)}_stage1/2/3.in")
