"""Rescale the logo network into Kremer-Grest units and write a runnable
LAMMPS input for minimisation + short equilibration.

The logo .data is a geometric artifact: bead spacing 0.2 lattice units, one
bond type, placeholder masses, no pair/bond coeffs. KG (the model topon's
README advertises: FENE + WCA) wants bond r0 ~0.97 sigma, so scale coordinates
by ~4.85 and attach real coefficients.

Bead types are preserved (1/3 quiet, 2/4 accent) so the render can still colour
the wordmark while the physics treats every bead identically -- the letters have
NO special forces holding them; they relax like anything else. That is the point
of the trajectory.

Usage:  python make_kg.py topon_entangled.data
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCALE = 4.85          # 0.2 lattice units -> ~0.97 sigma bond length
SRC = Path(sys.argv[1] if len(sys.argv) > 1 else "topon_entangled.data")
OUT = Path("kg_logo.data")


def main() -> None:
    L = SRC.read_text().splitlines()

    # header
    n_atoms = int(next(s for s in L if s.strip().endswith(" atoms")).split()[0])
    n_bonds = int(next(s for s in L if s.strip().endswith(" bonds")).split()[0])
    n_at = int(next(s for s in L if s.strip().endswith(" atom types")).split()[0])
    box = []
    for s in L[:30]:
        for t in ("xlo", "ylo", "zlo"):
            if t in s:
                p = s.split(); box.append((float(p[0]) * SCALE, float(p[1]) * SCALE))

    i = [k for k, s in enumerate(L) if s.startswith("Atoms")][0] + 2
    atoms = []
    while i < len(L) and L[i].strip():
        p = L[i].split()
        atoms.append((int(p[0]), int(p[1]), int(p[2]),
                      float(p[4]) * SCALE, float(p[5]) * SCALE, float(p[6]) * SCALE))
        i += 1
    j = [k for k, s in enumerate(L) if s.strip() == "Bonds"][0] + 2
    bonds = []
    while j < len(L) and L[j].strip():
        p = L[j].split(); bonds.append((int(p[0]), int(p[2]), int(p[3]))); j += 1

    # z is only 2 layers -> pad generously so the slab can breathe out-of-plane
    zlo, zhi = box[2]
    zpad = 6.0
    box[2] = (zlo - zpad, zhi + zpad)

    O = [f"topon logo, Kremer-Grest units (scaled x{SCALE})", "",
         f"{len(atoms)} atoms", f"{len(bonds)} bonds",
         f"{n_at} atom types", "1 bond types", "",
         f"{box[0][0]:.4f} {box[0][1]:.4f} xlo xhi",
         f"{box[1][0]:.4f} {box[1][1]:.4f} ylo yhi",
         f"{box[2][0]:.4f} {box[2][1]:.4f} zlo zhi", "",
         "Masses", ""] + [f"{t} 1.0" for t in range(1, n_at + 1)] + ["", "Atoms # molecular", ""]
    for aid, mol, t, x, y, z in atoms:
        O.append(f"{aid} {mol} {t} {x:.4f} {y:.4f} {z:.4f}")
    O += ["", "Bonds", ""]
    for bid, a, b in bonds:
        O.append(f"{bid} 1 {a} {b}")
    OUT.write_text("\n".join(O) + "\n", encoding="utf-8")
    print(f"[OK] {len(atoms)} beads, {len(bonds)} bonds -> {OUT}")

    # --- LAMMPS input: soft push-off, then FENE/WCA minimise + NVT ----------
    types = " ".join(str(t) for t in range(1, n_at + 1))
    inp = f"""# topon logo: KG minimisation + short equilibration
# Every bead is identical physics -- the letters are only a COLOURING, nothing
# holds them in shape. Frames are dumped so the trajectory can be animated.
units           lj
atom_style      molecular
boundary        p p p
read_data       {OUT.name}

# --- soft push-off: the built lattice has overlaps at 0.97 sigma spacing ---
pair_style      soft 1.12246
pair_coeff      * * 0.0
bond_style      harmonic
bond_coeff      1 100.0 0.97
neighbor        0.4 bin
neigh_modify    every 1 delay 0 check yes

variable        pre equal ramp(0,30)
fix             push all adapt 1 pair soft a * * v_pre
dump            d1 all custom 40 traj.lammpstrj id type x y z
dump_modify     d1 sort id
thermo          50
thermo_style    custom step temp pe ebond epair press

min_style       cg
minimize        1e-4 1e-6 400 4000
unfix           push

# --- real KG: FENE + WCA ------------------------------------------------
pair_style      lj/cut 1.12246
pair_coeff      * * 1.0 1.0 1.12246
pair_modify     shift yes
bond_style      fene
bond_coeff      1 30.0 1.5 1.0 1.0
special_bonds   fene

minimize        1e-4 1e-6 600 6000

# --- short equilibration (Langevin NVT) ---------------------------------
velocity        all create 0.4 4928459 dist gaussian
fix             1 all nve/limit 0.05
fix             2 all langevin 0.4 0.4 1.0 90210
timestep        0.004
run             4000
unfix           1
unfix           2
write_data      kg_relaxed.data nocoeff
print           "trajectory complete"
"""
    Path("in.kg_logo").write_text(inp, encoding="utf-8")
    print("[OK] in.kg_logo")


if __name__ == "__main__":
    main()
