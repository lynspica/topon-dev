# Topro issues to fix later

Discovered while building the MARTINI port (`topon/protein_network/`). They don't break topro's typical use cases (CHARMM atomistic with `lattice_scale` ~ 3 A, where wrap effects are small), but they would fail at any larger scale.

## 1. `pos % box` wrap breaks bonds across the periodic boundary

**Where:** `legacy/subprojects/protein_network/topro/topro/lammps/writer.py:220`
```python
pos = a.pos % box
```

**Why it's wrong:** LAMMPS does NOT use minimum-image distance for bonded
interactions (`bond_style harmonic` etc). It uses literal `r_a - r_b` from the
unwrapped positions. When `pos % box` puts two bonded atoms on opposite sides
of the periodic box, LAMMPS sees a `box - epsilon` Angstrom bond and either
errors with "bond atoms missing" or applies an enormous restoring force that
launches atoms.

**Why topro hides it:** with `lattice_scale ~ 3 A` and CHARMM bond lengths
~ 1.5 A, BFM segments rarely reach the box edge, so wraparound is rare. Once
in a while a chain might cross the boundary, but the relaxation protocol's
soft-push absorbs the resulting stretched bond before damage.

**Fix:** also write LAMMPS image flags (`ix iy iz` columns 8-10 in `atom_style
full`). Compute as `floor(pos / box)` per atom and write `wrapped_pos + image`.
LAMMPS reconstructs the unwrapped position correctly and bonds work. See
`topon/protein_network/lammps_writer.py:174-188` for the MARTINI port.

## 2. Interpolation overwrites anchor positions with chain-walk values

**Where:** `legacy/subprojects/protein_network/topro/topro/protein/builder.py:149-172`
```python
for k, r in enumerate(range(r_start, r_end + 1)):
    if r not in residue_positions:
        frac = k / max(n_seg_res - 1, 1)
        residue_positions[r] = p_start + frac * diff
```

**Why it's wrong:** when the loop hits the segment's END residue (`r = r_end,
k = n_seg_res - 1`), `frac = 1.0` and the position becomes `p_start + diff`
where `diff` is the **minimum-image displacement**. So the anchor residue is
NOT placed at its literal `lattice_pos(end_node)` -- it's placed at
`p_start + min_image_diff(p_start, p_end)`, which can differ by integer
multiples of the box.

For two chains crosslinked at the SAME merged BFM lattice node, each chain's
walk reaches that anchor from its own direction, so the two atoms land at
DIFFERENT unwrapped positions even though they should share the literal
`lattice_pos`. With `lattice_scale ~ 3 A` (topro), the offset is invisible
(< few Angstroms). With `lattice_scale ~ 27 A` (MARTINI), the offset is the
full box edge.

**Fix:** in two passes. First pass pins every anchor residue to its literal
`lattice_pos(chain[ni])`. Second pass interpolates ONLY strictly between
adjacent anchors. See the `_interpolate_residue_positions` function in
`topon/protein_network/builder.py:64-114`.

## 3. (minor) `prefactor = ramp(0, 30)` during minimize doesn't actually ramp

**Where:** `legacy/.../topro/topro/lammps/writer.py` and many user-facing scripts.

**Why:** LAMMPS `ramp(lo, hi)` is defined over a **run** of N timesteps, not
a minimize. During CG minimize, the variable just returns `lo` (= 0), so the
soft-push prefactor stays 0 and the minimize converges immediately doing
nothing (especially noticeable when atoms are at saddle points like r=0 of
the soft potential). Topro hides this by following the minimize with a brief
`nve/limit` dynamics block (where `ramp` works correctly), which kicks atoms
off saddles.

**Fix:** use `variable prefactor equal step*30/run_length` or just hardcode
the ramp endpoints in separate minimize calls.
