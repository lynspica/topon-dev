"""Post-process v41 matrix outputs to check for protein-network collapse.

For each cell, reads system_equilibrated.data and computes:
  * radius of gyration (Rg) of all protein backbone beads (combined)
  * mean nearest-neighbor distance between BB beads on different chains
  * fraction of protein atoms within 6 A of another protein atom on a
    different chain (a high fraction indicates inter-chain contact /
    aggregation, expected in dry case if it collapses)

Run after `run_v41_matrix.py`. Reads from tests/output/v41/<label>/.
"""
from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path

import numpy as np

V41_ROOT = Path("tests/output/v41")


def load_atoms(data_path: Path) -> tuple[list[dict], list[float]]:
    """Read positions/mol/type from a LAMMPS data file. Works with both the
    initial file (7-col + name comment) and LAMMPS-written equilibrated file
    (10-col with image flags, no comment)."""
    box = [0.0, 0.0, 0.0]
    atoms: list[dict] = []
    section = None
    with data_path.open() as f:
        for line in f:
            s = line.strip()
            if "xlo xhi" in s: box[0] = float(s.split()[1])
            elif "ylo yhi" in s: box[1] = float(s.split()[1])
            elif "zlo zhi" in s: box[2] = float(s.split()[1])
            elif s.startswith("Atoms"): section = "A"; continue
            elif s.startswith(("Velocities", "Bonds", "Masses", "Angles",
                               "Dihedrals", "Impropers", "Pair Coeffs",
                               "Bond Coeffs", "Improper Coeffs", "Angle Coeffs",
                               "Dihedral Coeffs")):
                section = None; continue
            if section == "A" and s and not s.startswith("#"):
                d, _, c = s.partition("#"); p = d.split()
                if p:
                    name = c.strip().split()[0] if c.strip() else ""
                    atoms.append({
                        "id": int(p[0]),
                        "mol": int(p[1]),
                        "type": int(p[2]),
                        "q": float(p[3]),
                        "pos": np.array([float(p[4]), float(p[5]), float(p[6])]),
                        "name": name,    # e.g. "TYR/BB"; empty for LAMMPS-written files
                    })
    return atoms, box


def load_atom_names(initial_path: Path) -> dict[int, str]:
    """Read atom_id -> name mapping from the initial data file (has comments)."""
    out: dict[int, str] = {}
    section = None
    with initial_path.open() as f:
        for line in f:
            s = line.strip()
            if s.startswith("Atoms"): section = "A"; continue
            elif s.startswith(("Velocities", "Bonds", "Masses", "Angles",
                               "Dihedrals", "Impropers", "Pair Coeffs",
                               "Bond Coeffs")):
                section = None; continue
            if section == "A" and s and not s.startswith("#"):
                d, _, c = s.partition("#"); p = d.split()
                if p and c.strip():
                    out[int(p[0])] = c.strip().split()[0]
    return out


def is_protein_bb(name: str) -> bool:
    return "/BB" in name


def is_protein(name: str) -> bool:
    # not water, not ion
    return name and not name.startswith(("W/", "SW/", "TW/", "ION/"))


def analyze_cell(label: str) -> dict:
    data = V41_ROOT / label / "system_equilibrated.data"
    initial = V41_ROOT / label / f"{label}.data"
    if not data.exists():
        return {"label": label, "status": "missing"}
    atoms, box = load_atoms(data)
    box_arr = np.array(box)

    # Names live in the initial file's comments; cross-reference by atom_id.
    name_by_id = load_atom_names(initial) if initial.exists() else {}
    for a in atoms:
        if not a["name"] and a["id"] in name_by_id:
            a["name"] = name_by_id[a["id"]]

    bb = [a for a in atoms if is_protein_bb(a["name"])]
    protein = [a for a in atoms if is_protein(a["name"])]
    if not bb or not protein:
        return {"label": label, "status": "no_protein"}

    # Rg of all BB (over the whole protein "cloud", min-image-folded around CoM)
    bb_pos = np.array([b["pos"] for b in bb])
    com = bb_pos.mean(axis=0)
    diff = bb_pos - com
    diff -= box_arr * np.round(diff / box_arr)
    rg = float(np.sqrt((diff ** 2).sum(axis=1).mean()))

    # Inter-chain BB-BB nearest-neighbor distances (one per BB, to the closest
    # BB on a DIFFERENT chain). Smaller = more collapsed.
    bb_by_mol: dict[int, list[int]] = defaultdict(list)
    for i, b in enumerate(bb):
        bb_by_mol[b["mol"]].append(i)
    if len(bb_by_mol) < 2:
        return {"label": label, "status": "single_chain"}

    near_dists: list[float] = []
    for i, b in enumerate(bb):
        own = b["mol"]
        others = np.array([bb_pos[j] for mid, idxs in bb_by_mol.items()
                           if mid != own for j in idxs])
        if len(others) == 0: continue
        d = others - bb_pos[i]
        d -= box_arr * np.round(d / box_arr)
        near_dists.append(float(np.sqrt((d ** 2).sum(axis=1).min())))
    mean_inter = float(np.mean(near_dists))
    median_inter = float(np.median(near_dists))

    # Fraction of protein atoms within 6 A of another chain
    # (cell-list for speed)
    cell_size = 6.0
    nx = max(1, int(math.ceil(box_arr[0] / cell_size)))
    ny = max(1, int(math.ceil(box_arr[1] / cell_size)))
    nz = max(1, int(math.ceil(box_arr[2] / cell_size)))
    grid: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    prot_pos = np.array([a["pos"] for a in protein])
    prot_mol = np.array([a["mol"] for a in protein])
    for i, p in enumerate(prot_pos):
        cx = int(p[0] // cell_size) % nx
        cy = int(p[1] // cell_size) % ny
        cz = int(p[2] // cell_size) % nz
        grid[(cx, cy, cz)].append(i)

    near_count = 0
    excl_sq = 6.0 ** 2
    for i, p in enumerate(prot_pos):
        own = prot_mol[i]
        cx = int(p[0] // cell_size) % nx
        cy = int(p[1] // cell_size) % ny
        cz = int(p[2] // cell_size) % nz
        found = False
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    if found: break
                    key = ((cx + dx) % nx, (cy + dy) % ny, (cz + dz) % nz)
                    for j in grid.get(key, ()):
                        if prot_mol[j] == own: continue
                        d = prot_pos[j] - p
                        d -= box_arr * np.round(d / box_arr)
                        if (d * d).sum() < excl_sq:
                            found = True; break
        if found: near_count += 1
    near_frac = near_count / len(protein)

    return {
        "label": label,
        "status": "ok",
        "n_chains": len(bb_by_mol),
        "n_bb": len(bb),
        "Rg": rg,
        "inter_chain_BB_NN_mean": mean_inter,
        "inter_chain_BB_NN_med": median_inter,
        "frac_protein_within_6A_other_chain": near_frac,
        "box_nm": [float(b) / 10.0 for b in box],
    }


def main() -> None:
    labels = sorted([d.name for d in V41_ROOT.iterdir() if d.is_dir()])
    rows = [analyze_cell(lab) for lab in labels]
    print(f"{'label':<14} {'Rg(A)':>8} {'BB-NN mean':>11} {'BB-NN med':>10} "
          f"{'frac<6A':>9} {'n_chains':>8}")
    print("-" * 70)
    for r in rows:
        if r["status"] != "ok":
            print(f"{r['label']:<14} [{r['status']}]")
            continue
        print(f"{r['label']:<14} {r['Rg']:>8.2f} "
              f"{r['inter_chain_BB_NN_mean']:>11.2f} "
              f"{r['inter_chain_BB_NN_med']:>10.2f} "
              f"{r['frac_protein_within_6A_other_chain']:>9.3f} "
              f"{r['n_chains']:>8d}")


if __name__ == "__main__":
    main()
