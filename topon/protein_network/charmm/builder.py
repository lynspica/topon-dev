"""
topro.protein.builder — Build an atomistic protein system from a BFM snapshot.

Key functions
-------------
build_protein_system(ff, snapshot, full_sequence, node_to_res, lattice_scale)
    Returns dry atoms, bonds, improper_defs, box, crosslink_bonds.

add_water_and_ions(atoms, bonds, box, water_content_pct, salt_conc_M)
    Solvates the system in-place.  Handles charge neutralisation automatically.

compute_lattice_scale(Nx, protein_mass_Da, water_content_pct, target_density)
    Returns the Å/BFM-unit scale for a given water content and target density.
"""

import random
from collections import defaultdict

import numpy as np

from .charmm_ff import CHARMMForceField


# ── Atom data class ───────────────────────────────────────────────────────────

class Atom:
    __slots__ = ["idx", "name", "atype", "charge", "res_name",
                 "res_id", "chain_id", "pos", "mol_id"]

    def __init__(self, idx, name, atype, charge, res_name,
                 res_id, chain_id, pos, mol_id):
        self.idx = idx
        self.name = name
        self.atype = atype
        self.charge = charge
        self.res_name = res_name
        self.res_id = res_id
        self.chain_id = chain_id
        self.pos = np.asarray(pos, dtype=float)
        self.mol_id = mol_id


# ── Box / density helpers ─────────────────────────────────────────────────────

def compute_lattice_scale(Nx, protein_mass_Da, water_content_pct=0.0,
                          target_density=0.85):
    """
    Compute lattice_scale (Å per BFM unit) such that the simulation box
    accommodates protein + water at the specified target density.

    Parameters
    ----------
    Nx : int           Lattice edge length (cubic).
    protein_mass_Da : float  Total dry protein mass in Daltons.
    water_content_pct : float  Weight percent water (0–100).
    target_density : float  g/cm³.  Default 0.85 (slightly loose initial box).

    Returns
    -------
    float : lattice_scale in Å/BFM unit
    """
    total_mass = protein_mass_Da
    if water_content_pct > 0.0:
        water_mass = protein_mass_Da * water_content_pct / (100.0 - water_content_pct)
        total_mass += water_mass

    # V [Å³] = M [Da] × (1.66054 Å³/Da / target_density)
    V = total_mass * 1.66054 / target_density
    L_box = V ** (1.0 / 3.0)
    return L_box / Nx


def _estimate_protein_mass(atoms):
    """Rough protein mass estimation from atom types (Da)."""
    mass = 0.0
    for a in atoms:
        t = a.atype
        if t.startswith("C"):
            mass += 12.0
        elif t.startswith("O"):
            mass += 16.0
        elif t.startswith("N"):
            mass += 14.0
        elif t.startswith("H"):
            mass += 1.0
        elif t.startswith("S"):
            mass += 32.0
        else:
            mass += 10.0
    return mass


# ── Main system builder ────────────────────────────────────────────────────────

def build_protein_system(ff, snapshot, full_sequence, node_to_res,
                         lattice_scale=15.0):
    """
    Build the full atomistic (dry) protein system from a BFM topology snapshot.

    Parameters
    ----------
    ff : CHARMMForceField
    snapshot : dict
        Single snapshot from a topology file.
        Must have: chains, Nx, Ny, Nz, crosslinker_positions, reactions.
    full_sequence : list of str
        3-letter residue names for the whole chain
        (length = n_residues = len(node_to_res possible residues)).
    node_to_res : dict
        chain_node_index → residue_index mapping (from sequence.get_node_residue_mapping).
    lattice_scale : float
        Å per BFM lattice unit.  Default 15.0.

    Returns
    -------
    atoms, bonds, improper_defs, box, crosslink_bonds
    """
    chains = snapshot["chains"]
    Nx, Ny, Nz = snapshot["Nx"], snapshot["Ny"], snapshot["Nz"]
    xlink_positions = set(snapshot["crosslinker_positions"])
    n_residues = len(full_sequence)

    box = np.array([Nx, Ny, Nz], dtype=float) * lattice_scale

    def lattice_pos(flat_idx):
        z = flat_idx // (Nx * Ny)
        rem = flat_idx % (Nx * Ny)
        y = rem // Nx
        x = rem % Nx
        return np.array([x, y, z], dtype=float) * lattice_scale

    atoms = []
    bonds = []
    improper_defs = []
    crosslink_bonds = []

    atom_counter = 0
    global_res_counter = 0

    # (chain_id, node_idx_in_chain, atom_name) → global atom idx
    node_atom_lookup = {}

    for chain_id, chain_nodes in enumerate(chains):
        n_chain_nodes = len(chain_nodes)

        # Build lattice positions for each node
        node_positions = {ni: lattice_pos(fi) for ni, fi in enumerate(chain_nodes)}

        # Interpolate residue coordinates along the chain
        residue_positions = {}
        sorted_nodes = sorted(node_to_res.keys())
        for seg_i in range(len(sorted_nodes) - 1):
            ni_start = sorted_nodes[seg_i]
            ni_end = sorted_nodes[seg_i + 1]

            if ni_start >= n_chain_nodes or ni_end >= n_chain_nodes:
                continue

            r_start = node_to_res[ni_start]
            r_end = node_to_res[ni_end]
            p_start = node_positions[ni_start]
            p_end = node_positions[ni_end]

            # Minimum image displacement
            diff = p_end - p_start
            diff -= box * np.round(diff / box)

            n_seg_res = r_end - r_start + 1
            for k, r in enumerate(range(r_start, r_end + 1)):
                if r not in residue_positions:
                    frac = k / max(n_seg_res - 1, 1)
                    residue_positions[r] = p_start + frac * diff

        # Instantiate atoms residue by residue
        chain_atom_ids = {}   # (res_idx, atom_name) → global_atom_idx
        prev_c_idx = None

        for res_idx in range(n_residues):
            res_name = full_sequence[res_idx]
            res_tmpl = ff.residues.get(res_name)
            if not res_tmpl:
                continue

            global_res_counter += 1
            center = residue_positions.get(res_idx, np.zeros(3))

            atom_list = dict(res_tmpl["atoms"])
            bond_list = list(res_tmpl["bonds"])
            impr_list = list(res_tmpl.get("impropers", []))
            deletes = []

            # Terminal patches
            if res_idx == 0:
                patch_name = "GLYP" if res_name == "GLY" else "NTER"
                _apply_patch(ff, patch_name, atom_list, bond_list, impr_list, deletes)

            if res_idx == n_residues - 1:
                _apply_patch(ff, "CTER", atom_list, bond_list, impr_list, deletes)

            for d in deletes:
                atom_list.pop(d, None)

            for atom_name, (atype, charge) in atom_list.items():
                if atom_name.startswith("+") or atom_name.startswith("-"):
                    continue
                atom_counter += 1
                offset = np.zeros(3)
                if atom_name != "CA":
                    rng = np.random.default_rng(seed=atom_counter)
                    offset = rng.uniform(-0.3, 0.3, 3)

                a = Atom(
                    idx=atom_counter,
                    name=atom_name,
                    atype=atype,
                    charge=charge,
                    res_name=res_name,
                    res_id=global_res_counter,
                    chain_id=chain_id,
                    pos=center + offset,
                    mol_id=chain_id + 1,
                )
                atoms.append(a)
                chain_atom_ids[(res_idx, atom_name)] = atom_counter

            # Intra-residue bonds
            for a1, a2 in bond_list:
                if a1.startswith("+") or a1.startswith("-"):
                    continue
                if a2.startswith("+") or a2.startswith("-"):
                    continue
                id1 = chain_atom_ids.get((res_idx, a1))
                id2 = chain_atom_ids.get((res_idx, a2))
                if id1 and id2:
                    bonds.append((id1, id2))

            # Intra-residue impropers
            for quad in impr_list:
                names = list(quad)
                if any(n.startswith("+") or n.startswith("-") for n in names):
                    continue
                ids = [chain_atom_ids.get((res_idx, n)) for n in names]
                if all(ids):
                    improper_defs.append(tuple(ids))

            # Peptide bond to previous residue
            if res_idx > 0 and prev_c_idx is not None:
                n_idx = chain_atom_ids.get((res_idx, "N"))
                if n_idx:
                    bonds.append((prev_c_idx, n_idx))

            prev_c_idx = chain_atom_ids.get((res_idx, "C"))

            # Track Y atoms for crosslink stitching
            for ni, ri in node_to_res.items():
                if ri == res_idx and ni in xlink_positions:
                    for aname in atom_list:
                        aid = chain_atom_ids.get((res_idx, aname))
                        if aid:
                            node_atom_lookup[(chain_id, ni, aname)] = aid

    # ── Apply crosslinks ──────────────────────────────────────────────────────
    node_usage = defaultdict(list)
    for ci, chain in enumerate(chains):
        for ni, fi in enumerate(chain):
            if ni in xlink_positions:
                node_usage[fi].append((ci, ni))

    dity_patch = ff.patches.get("DITY", {})
    dity_atoms = dity_patch.get("atoms", {})

    dity_charge_map = {}
    for aname, (atype, charge) in dity_atoms.items():
        clean = aname[1:] if (len(aname) > 1 and aname[0] in ("1", "2")) else aname
        dity_charge_map[clean] = (atype, charge)

    crosslinked_residues = set()
    for fi, usages in node_usage.items():
        if len(usages) <= 1:
            continue
        for i in range(len(usages)):
            for j in range(i + 1, len(usages)):
                ci1, ni1 = usages[i]
                ci2, ni2 = usages[j]
                ce2_1 = node_atom_lookup.get((ci1, ni1, "CE2"))
                ce2_2 = node_atom_lookup.get((ci2, ni2, "CE2"))
                if ce2_1 and ce2_2:
                    bonds.append((ce2_1, ce2_2))
                    crosslink_bonds.append((ce2_1, ce2_2))
                    crosslinked_residues.add((ci1, ni1))
                    crosslinked_residues.add((ci2, ni2))

    # Apply DITY patch charges/types and delete HE2
    atoms_to_remove = set()
    for ci, ni in crosslinked_residues:
        for atom in atoms:
            if atom.chain_id == ci:
                aid_check = node_atom_lookup.get((ci, ni, atom.name))
                if aid_check == atom.idx:
                    if atom.name in dity_charge_map:
                        atom.atype, atom.charge = dity_charge_map[atom.name]
                    if atom.name == "HE2":
                        atoms_to_remove.add(atom.idx)

    if atoms_to_remove:
        atoms = [a for a in atoms if a.idx not in atoms_to_remove]
        bonds = [(a, b) for a, b in bonds
                 if a not in atoms_to_remove and b not in atoms_to_remove]

    return atoms, bonds, improper_defs, box, crosslink_bonds


def _apply_patch(ff, patch_name, atom_list, bond_list, impr_list, deletes):
    patch = ff.patches.get(patch_name, {})
    if not patch:
        return
    for d in patch.get("deletes", []):
        deletes.append(d)
    for aname, (atype, charge) in patch.get("atoms", {}).items():
        atom_list[aname] = (atype, charge)
    bond_list.extend(patch.get("bonds", []))
    impr_list.extend(patch.get("impropers", []))


# ── Solvation ─────────────────────────────────────────────────────────────────

def add_water_and_ions(atoms, bonds, box, water_content_pct=35.0,
                       salt_conc_M=0.15, cation_type="SOD", anion_type="CLA"):
    """
    Add TIP3P water and NaCl ions to the system (in-place).

    Water amount   → from weight-percent relative to protein mass.
    NaCl (background) → 0.15 M based on water volume.
    Charge neutralisation → extra cations or anions to zero the net charge.

    Parameters
    ----------
    atoms : list of Atom   (modified in-place)
    bonds : list of (int,int)  (modified in-place)
    box   : np.array([Lx,Ly,Lz])
    water_content_pct : float  weight percent water.  0 = no water.
    salt_conc_M : float  background NaCl concentration in mol/L.  Default 0.15.
    cation_type : str  CHARMM atom type for cation.  Default 'SOD'.
    anion_type  : str  CHARMM atom type for anion.   Default 'CLA'.
    """
    # Protein mass & charge
    protein_mass = _estimate_protein_mass(atoms)
    net_charge = round(sum(a.charge for a in atoms))

    # Water count from weight fraction
    num_waters = 0
    if water_content_pct > 0.0:
        water_mass = protein_mass * water_content_pct / (100.0 - water_content_pct)
        num_waters = int(water_mass / 18.015)

    # Background NaCl: n_NaCl = C [mol/L] × N_A × V_water [L]
    # V_water [L] = n_waters × 18.015 [g/mol] / (6.022e23 [/mol] × 1.0 [g/cm³]) × 1e-3
    n_nacl = 0
    if num_waters > 0 and salt_conc_M > 0.0:
        V_water_L = num_waters * 18.015 / (6.022e23 * 1.0) * 1e-3 * 1e24
        # V_water_L: num_waters * 18.015 / (6.022e23) mL → × 1e-3 for L
        V_water_L = num_waters * 18.015 / (6.022e23 * 1000.0)
        n_nacl = int(round(salt_conc_M * 6.022e23 * V_water_L))

    # Charge neutralisation
    # net_charge > 0: protein positive → extra Cl-
    # net_charge < 0: protein negative → extra Na+
    n_cations = n_nacl + max(0, -net_charge)
    n_anions  = n_nacl + max(0,  net_charge)

    print(f"\n  [Solvation]")
    print(f"    Protein mass : {protein_mass:.0f} Da | net charge = {net_charge:+d} e")
    print(f"    Water        : {num_waters} molecules ({water_content_pct:.0f} wt%)")
    print(f"    NaCl (bg)    : {n_nacl} pairs at {salt_conc_M:.3f} M")
    print(f"    Cations ({cation_type}) : {n_cations}")
    print(f"    Anions  ({anion_type})  : {n_anions}")

    start_idx = max(a.idx for a in atoms) + 1 if atoms else 1
    max_mol = max(a.mol_id for a in atoms) if atoms else 0

    def rand_pos(margin=1.5):
        return (
            random.uniform(margin, box[0] - margin),
            random.uniform(margin, box[1] - margin),
            random.uniform(margin, box[2] - margin),
        )

    # Cations
    for _ in range(n_cations):
        max_mol += 1
        atoms.append(Atom(start_idx, cation_type, cation_type, 1.0,
                          cation_type, max_mol, max_mol, rand_pos(), max_mol))
        start_idx += 1

    # Anions
    for _ in range(n_anions):
        max_mol += 1
        atoms.append(Atom(start_idx, anion_type, anion_type, -1.0,
                          anion_type, max_mol, max_mol, rand_pos(), max_mol))
        start_idx += 1

    # TIP3P water  (O-H bond = 0.9572 Å, H-O-H = 104.52°)
    H2_dx = 0.9572 * np.cos(np.radians(104.52))   # ≈ -0.2399
    H2_dy = 0.9572 * np.sin(np.radians(104.52))   # ≈  0.9266
    for _ in range(num_waters):
        max_mol += 1
        ox, oy, oz = rand_pos()
        iO  = start_idx
        iH1 = start_idx + 1
        iH2 = start_idx + 2
        atoms.append(Atom(iO,  "OH2", "OT", -0.834, "TIP3", max_mol, max_mol,
                          [ox,           oy,           oz], max_mol))
        atoms.append(Atom(iH1, "H1",  "HT",  0.417, "TIP3", max_mol, max_mol,
                          [ox + 0.9572,  oy,           oz], max_mol))
        atoms.append(Atom(iH2, "H2",  "HT",  0.417, "TIP3", max_mol, max_mol,
                          [ox + H2_dx,   oy + H2_dy,   oz], max_mol))
        bonds.append((iO, iH1))
        bonds.append((iO, iH2))
        start_idx += 3


# ── CMAP backbone crossterms ──────────────────────────────────────────────────

def find_cmap_crossterms(atoms):
    """
    Identify backbone CMAP crossterms for all non-terminal protein residues.

    For each internal residue i in a chain, the crossterm spans five backbone
    atoms: C(i-1) - N(i) - CA(i) - C(i) - N(i+1), defining the phi/psi pair.

    CMAP type (matching charmm36m.cmap / charmm36.cmap):
        1 = regular residue
        2 = regular residue before PRO
        3 = PRO
        4 = PRO before PRO
        5 = GLY
        6 = GLY before PRO

    Water, ion, and other non-protein atoms (identified by atype) are skipped.

    Returns
    -------
    list of (cmap_type, c_prev_idx, n_i_idx, ca_i_idx, c_i_idx, n_next_idx)
        All indices are original atom .idx values (pre-renumbering).
    """
    _SKIP_ATYPES = {"OT", "HT", "SOD", "CLA", "CAL", "ZN", "FE3P"}

    from collections import defaultdict

    # Group atoms by (chain_id, res_id) — protein only
    # chain_res[chain_id][res_id] = {atom_name: atom_idx}
    # chain_resname[chain_id][res_id] = res_name
    chain_res     = defaultdict(lambda: defaultdict(dict))
    chain_resname = defaultdict(dict)

    for a in atoms:
        if a.atype in _SKIP_ATYPES:
            continue
        chain_res[a.chain_id][a.res_id][a.name] = a.idx
        chain_resname[a.chain_id][a.res_id] = a.res_name

    crossterms = []

    for chain_id in sorted(chain_res.keys()):
        res_ids = sorted(chain_res[chain_id].keys())
        n_res   = len(res_ids)

        # Skip terminal residues (need i-1 and i+1)
        for pos in range(1, n_res - 1):
            prev_rid = res_ids[pos - 1]
            curr_rid = res_ids[pos]
            next_rid = res_ids[pos + 1]

            c_prev = chain_res[chain_id][prev_rid].get("C")
            n_i    = chain_res[chain_id][curr_rid].get("N")
            ca_i   = chain_res[chain_id][curr_rid].get("CA")
            c_i    = chain_res[chain_id][curr_rid].get("C")
            n_next = chain_res[chain_id][next_rid].get("N")

            if None in (c_prev, n_i, ca_i, c_i, n_next):
                continue

            curr_name = chain_resname[chain_id][curr_rid]
            next_name = chain_resname[chain_id][next_rid]

            before_pro = (next_name == "PRO")
            if curr_name == "GLY":
                cmap_type = 6 if before_pro else 5
            elif curr_name == "PRO":
                cmap_type = 4 if before_pro else 3
            else:
                cmap_type = 2 if before_pro else 1

            crossterms.append((cmap_type, c_prev, n_i, ca_i, c_i, n_next))

    return crossterms
