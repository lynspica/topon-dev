"""
topon.workflows.cg_network
==========================
Canonical workflow: Coarse-Grained (Kremer-Grest) polymer network.

Four-stage pipeline
-------------------
1. Topology   — load network graph → NetworkGraph
2. Chemistry  — assign bead types, build RDKit mol, write LAMMPS data + displacements
3. Conformation — apply displacements, resolve overlaps
4. Simulation — write LAMMPS input scripts (minimize / equilibrate)

Supports optional assignments:
  - Entanglements (kinked backbones, v15.2 midpoint logic)
  - Grafts (side chains attached to backbone beads, v20 dynamic scaling)

Usage (Python)
--------------
    from topon.workflows.cg_network import run

    run(
        nodes_path="tests/sample_graphs/network_N5x5x5_trial3.nodes",
        edges_path="tests/sample_graphs/network_N5x5x5_trial3.edges",
        config_path="examples/config_cg_combined.json",
        experimental_path="examples/experimental_test.json",
        output_dir="output/cg_run",
    )

Usage (CLI)
-----------
    python -m topon.workflows.cg_network \\
        --nodes tests/sample_graphs/network_N5x5x5_trial3.nodes \\
        --edges tests/sample_graphs/network_N5x5x5_trial3.edges \\
        --config examples/config_cg_combined.json \\
        --output output/cg_run
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem

# Stage 1 — topology interface
from topon.topology.network import load as load_network

# Stage 2 — chemistry (writers + displacement utilities)
from topon.writers import CGWriter, LammpsInputGenerator
from topon.utils import write_lammps_displacement_file
from topon.assignment.attributor import EntanglementsConfig
from topon.assignment.entanglements import select_entanglements
from topon.utils.network_helpers import calculate_entangled_kink

# Stage 3 — conformation
from topon.conformation import ConformationManager

# Stage 4 — simulation runner (optional, only if auto_run=True)
from topon.simulation import SimulationRunner


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run(
    nodes_path: str | Path,
    edges_path: str | Path,
    config_path: str | Path,
    experimental_path: str | Path,
    output_dir: str | Path,
    gpickle_path: str | Path | None = None,
    seed: int | None = None,
) -> Path:
    """
    Run the full CG network workflow.

    Parameters
    ----------
    nodes_path, edges_path : .nodes / .edges topology files
    gpickle_path : alternative to nodes/edges (takes precedence)
    config_path : JSON config (chemistry, assignment, conformation, simulation)
    experimental_path : JSON experimental config (KG dynamics parameters)
    output_dir : root output directory (study subfolder created inside)

    Returns
    -------
    Path to the study root directory (contains 02_Chemistry/, 03_Conformation/,
    04_Simulation/).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    config = json.loads(Path(config_path).read_text())
    experimental = json.loads(Path(experimental_path).read_text())

    DP = config["chemistry"]["degree_of_polymerization"]
    DENSITY = config["chemistry"]["bead_density"]
    ASSIGNMENT = config.get("assignment", {})

    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    # =========================================================================
    # Stage 1: Topology
    # =========================================================================
    print("=" * 60)
    print("[Stage 1] Loading topology...")

    network = load_network(
        nodes_path=gpickle_path or nodes_path,
        edges_path=None if gpickle_path else edges_path,
        gpickle_path=gpickle_path,
    )
    # Unpack to (G, dims) for the legacy code used in stages 2-3
    G, dims = network.to_legacy()

    # --- Optional: Entanglements ---
    ent_conf_dict = ASSIGNMENT.get("entanglements", {})
    if ent_conf_dict.get("enabled"):
        print(f"  Applying entanglements (target={ent_conf_dict.get('target')})...")
        ent_config = EntanglementsConfig(**ent_conf_dict)
        select_entanglements(G, ent_config, dims)

    # --- Optional: Graft config ---
    grafts_cfg = ASSIGNMENT.get("grafts", {})
    grafts_enabled = grafts_cfg.get("enabled", False)
    per_edge_grafts = grafts_cfg.get("per_edge_type", {})
    extension_factor = experimental.get("cg", {}).get("graft_extension_factor", 0.5)

    print(f"  {network}")
    print(f"  Entanglements: {ent_conf_dict.get('enabled', False)}, "
          f"Grafts: {grafts_enabled}")

    # =========================================================================
    # Stage 2: Chemistry — bead typing, mol building, writers
    # =========================================================================
    print("[Stage 2] Chemistry embedding...")

    study_name = config.get("study", {}).get("name", "cg_network")
    root = output_dir / study_name
    chem_dir = root / "02_Chemistry"
    chem_dir.mkdir(parents=True, exist_ok=True)

    mol = Chem.RWMol()
    node_map = {}       # graph node id → mol atom idx
    edge_atom_map = {}  # edge index → backbone atom idx list
    graft_atom_map = {} # edge index → [(backbone_pos_k, graft_atom_idx_list)]

    # Nodes (junction beads)
    for node in sorted(G.nodes()):
        idx = mol.AddAtom(Chem.Atom("Si"))
        mol.GetAtomWithIdx(idx).SetProp("bead_type", "J")
        node_map[node] = idx

    edges = list(G.edges(data=True))
    total_grafts = 0

    for i, (u, v, data) in enumerate(edges):
        dp_val = data.get("dp", DP)
        edge_type = data.get("type", "A")

        graft_conf = per_edge_grafts.get(edge_type)
        should_graft = grafts_enabled and graft_conf is not None
        if should_graft:
            graft_density = graft_conf.get("graft_density", 0.0)
            side_dp = graft_conf.get("side_chain_dp", 5)

        chain_atoms = []
        edge_grafts = []
        prev = node_map[u]

        for k in range(dp_val):
            # Backbone bead
            idx = mol.AddAtom(Chem.Atom("Si"))
            mol.GetAtomWithIdx(idx).SetProp("bead_type", "A")
            chain_atoms.append(idx)
            mol.AddBond(prev, idx, Chem.BondType.SINGLE)

            # Graft side chain
            if should_graft and random.random() < graft_density:
                g_prev = idx
                graft_chain = []
                for _ in range(side_dp):
                    g_idx = mol.AddAtom(Chem.Atom("Si"))
                    mol.GetAtomWithIdx(g_idx).SetProp("bead_type", "G")
                    graft_chain.append(g_idx)
                    mol.AddBond(g_prev, g_idx, Chem.BondType.SINGLE)
                    g_prev = g_idx
                edge_grafts.append((k, graft_chain))
                total_grafts += 1

            prev = idx

        mol.AddBond(prev, node_map[v], Chem.BondType.SINGLE)
        edge_atom_map[i] = chain_atoms
        graft_atom_map[i] = edge_grafts

    print(f"  Grafts added: {total_grafts}")
    entangled_count = sum(1 for _, _, d in edges if d.get("entangled_with") is not None)
    print(f"  Entangled edges: {entangled_count}")

    # Box scale from density
    n_beads = mol.GetNumAtoms()
    scale = ((n_beads / DENSITY) / np.prod(dims)) ** (1 / 3.0)
    sx = sy = sz = scale

    # Write LAMMPS data file
    include_angles = config.get("simulation", {}).get("include_angles", True)
    CGWriter(mol, str(chem_dir / "system.data"), include_angles=include_angles).write()

    # Settings file (CG uses inline coeffs in input scripts, not a separate settings file)
    (chem_dir / "system.in.settings").write_text("# Dummy settings\n")

    # Displacement files: nodes
    node_coords = {
        idx: G.nodes[node].get("pos", (0, 0, 0))
        for node, idx in node_map.items()
    }
    write_lammps_displacement_file(
        node_coords, sx, sy, sz,
        str(chem_dir / "system_nodes.displace"), "nodes"
    )

    # Displacement files: backbone beads + grafts
    chain_coords = {}
    graft_coords = {}

    for i, atoms in edge_atom_map.items():
        u, v, data = edges[i]
        pos_u = np.array(G.nodes[u].get("pos", (0, 0, 0)))
        pos_v = np.array(G.nodes[v].get("pos", (0, 0, 0)))
        vec = pos_v - pos_u
        mic = vec - dims * np.round(vec / dims)
        edge_len = np.linalg.norm(mic)

        unit_vec = mic / (edge_len + 1e-9)
        rand_vec = np.random.randn(3)
        perp = np.cross(unit_vec, rand_vec)
        if np.linalg.norm(perp) < 1e-6:
            perp = np.cross(unit_vec, np.array([1.0, 0.0, 0.0]))
        perp_unit = perp / np.linalg.norm(perp)

        entangled_partner_key = data.get("entangled_with")
        backbone_xyz = []

        if entangled_partner_key is not None:
            p_u, p_v = entangled_partner_key[0], entangled_partner_key[1]
            p_pos_u = np.array(G.nodes[p_u]["pos"])
            p_pos_v = np.array(G.nodes[p_v]["pos"])
            p_vec = p_pos_v - p_pos_u
            p_mic = p_vec - dims * np.round(p_vec / dims)

            my_mid = pos_u + 0.5 * mic
            p_mid = p_pos_u + 0.5 * p_mic
            delta = p_mid - my_mid
            delta -= dims * np.round(delta / dims)
            p_mid_wrapped = my_mid + delta

            orient_vec = p_mid_wrapped - my_mid
            if np.linalg.norm(orient_vec) < 0.01:
                orient_vec = perp_unit

            kink_dict = calculate_entangled_kink(
                start_pos=np.zeros(3),
                end_pos=mic,
                num_atoms=len(atoms) + 2,
                orientation_vec=orient_vec,
                z_phase=1.0,
            )
            full_path = [kink_dict[k] for k in sorted(kink_dict.keys())]
            backbone_xyz = [pos_u + np.array(pt) for pt in full_path[1:-1]]
        else:
            for j in range(len(atoms)):
                frac = (j + 1) / (len(atoms) + 1)
                backbone_xyz.append(pos_u + frac * mic)

        for j, atom_idx in enumerate(atoms):
            chain_coords[atom_idx] = tuple(backbone_xyz[j])

        for k, g_atoms in graft_atom_map.get(i, []):
            anchor_pos = backbone_xyz[k]
            graft_dp = len(g_atoms)
            eff_factor = min(extension_factor, graft_dp / len(atoms))
            graft_vec = perp_unit * (edge_len * eff_factor)
            for m, g_idx in enumerate(g_atoms):
                g_frac = (m + 1) / len(g_atoms)
                graft_coords[g_idx] = tuple(anchor_pos + g_frac * graft_vec)

    write_lammps_displacement_file(
        chain_coords, sx, sy, sz,
        str(chem_dir / "system_beads.displace"), "beads"
    )
    write_lammps_displacement_file(
        graft_coords, sx, sy, sz,
        str(chem_dir / "system_grafts.displace"), "grafts"
    )

    (chem_dir / "system.groups").write_text(
        "# Groups\n"
        f"group nodes id {' '.join(str(idx+1) for idx in sorted(node_map.values()))}\n"
        "group beads subtract all nodes\n"
    )

    # =========================================================================
    # Stage 3: Conformation — apply displacements, resolve overlaps
    # =========================================================================
    print("[Stage 3] Conformation embedding...")

    cm = ConformationManager(str(output_dir), study_name)
    conformed, roles = cm.apply_displacements("system.data")
    noisy = cm.apply_noise(conformed, magnitude=1e-4)
    cm.resolve_overlaps(
        noisy, roles,
        cutoff=config["conformation"]["overlap_cutoff"],
        max_iters=config["conformation"].get("overlap_max_iters", 20),
    )

    # =========================================================================
    # Stage 4: Simulation — write LAMMPS input scripts
    # =========================================================================
    print("[Stage 4] Writing LAMMPS input scripts...")

    gen = LammpsInputGenerator(
        str(output_dir), study_name,
        config=config.get("simulation", {}),
        experimental=experimental,
    )
    gen.write_serial_soft_minimization(
        settings_file="system.in.settings", model_type="cg"
    )
    gen.write_parallel_production(
        settings_file="system.in.settings", model_type="cg"
    )

    print(f"CG network workflow complete -> {root}")

    # Optional: run LAMMPS
    if config.get("execution", {}).get("auto_run", False):
        exec_cfg = config.get("execution", {})
        runner = SimulationRunner(
            sim_dir=root / "04_Simulation",
            executable=exec_cfg.get("executable", "lmp"),
            n_procs=exec_cfg.get("n_procs", 1),
            use_mpi=False,
        )
        runner.run_sequence([
            "minimize_1_serial.in",
            "minimize_2_parallel.in",
            "minimize_3_parallel.in",
        ])

    return root


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="CG polymer network workflow (Kremer-Grest)"
    )
    p.add_argument("--nodes", required=True, help=".nodes topology file")
    p.add_argument("--edges", required=True, help=".edges topology file")
    p.add_argument("--config", required=True, help="JSON config file")
    p.add_argument("--experimental", required=True, help="JSON experimental config")
    p.add_argument("--output", required=True, help="Output directory")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducible graft/entanglement placement")
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    run(
        nodes_path=args.nodes,
        edges_path=args.edges,
        config_path=args.config,
        experimental_path=args.experimental,
        output_dir=args.output,
        seed=args.seed,
    )
