# The gallery

Two arc **animations** (CG and atomistic, `anim/`) plus eleven **stills** for the
repo README, in three sections: the lattice → minimised → equilibrated arc, the
copolymer sequences, and entanglements + side chains. **Every panel sits on a
strict-sculpted, heterogeneous network** (junction functionality 2–6, mean 4.0),
not a perfect lattice — that is the whole point. Nothing was built for the
picture, and every number in the captions was measured back out of the file it
shows.

## Rebuilding

```bash
# 1. generate the systems (written outside the repo -- big and transient)
export TOPON_GALLERY_SYSTEMS=/some/scratch/systems
python gen_systems.py       # sculpt_250, copoly_*, entangled_grafted, atom_sculpt

# 2. render the stills
"C:/v/ovito/Scripts/python.exe" render_gallery.py     # or: <panel> [<panel> ...]

# 3. the three arc animations (lattice -> non-ideal melt). Each needs a dedicated
#    movie run, then one builder assembles the GIF + MP4:
#      cg / copoly  -> movie_cg.in   (copoly first needs `lmp minimize_1_serial.in`
#                                      for its 1.restart)
#      atom         -> movie_atom.in
cp movie_cg.in   $TOPON_GALLERY_SYSTEMS/sculpt_250/sculpt_250/04_Simulation/
cp movie_cg.in   $TOPON_GALLERY_SYSTEMS/copoly_block/copoly_block/04_Simulation/
cp movie_atom.in $TOPON_GALLERY_SYSTEMS/atom_sculpt/atom_sculpt/04_Simulation/
# ... run each deck with lmp in its 04_Simulation dir (stage 1 SERIAL where needed) ...
"C:/v/ovito/Scripts/python.exe" make_arc_movie.py     # cg copoly atom -> anim/*.{gif,mp4}
```

`gen_systems.py` sits beside this file — a thin driver over
`topon.config.load_config` + `topon.pipeline.Pipeline`. The entanglement panel
shows the lattice state and never touches LAMMPS.

**The animations are real MD trajectories, not morphs.** A LAMMPS `dump` writes
positions but not bonds, so the movie loads the `03_Conformation` data file for
bond topology and updates coordinates per frame (`LoadTrajectoryModifier`), fixed
camera, then builds a seamless boomerang (forward + reversed).

**Why a dedicated movie run.** Under the production decks the lattice→melt
transition completes inside a single `minimize`, so a fixed-interval dump lands
the whole jump between two frames (median bond 0.17 → 1.00 in one step, then 65
static frames). The movie decks inflate the lattice under a *ramped soft push*
with bounded dynamics (`nve/limit`) instead — the transition is driven by excluded
volume, not by the bonds (a symmetric compressed lattice is force-balanced) —
dumping every step, so the strands open smoothly.

**Even pacing needs the right progress axis.** `make_arc_movie.py` selects ~35
frames evenly in a monotone progress metric, not in time. For CG/copoly that is
**bond-orientation disorder** (`1 − ⟨max axis component / |bond|⟩`), which keeps
rising through the inflation *and* the subsequent coiling — bond length saturates
at ~0.92 while the strands are still coiling, so pacing by it still jump-cuts. For
atomistic (where methyl-H bonds give a nonzero disorder baseline) it uses **mean
minimum-image displacement from the lattice**. Result: largest per-frame visual
gap 0.011 (CG) / 0.018 (copoly) / 0.14 (atom), versus a single 0.83 cliff for a
naive fixed-interval dump.

**Per-resolution notes.** For CG/copoly the soft push inflates, then an LJ
epsilon-ramp + NVT tail coils the melt (harmonic bonds kept — FENE would snap the
freshly-inflated bonds). The atomistic soft push (real units, `movie_atom.in`)
already carries the whole arc on its own — bonds 0.44 → 1.12 Å *and* continued
coiling (mean displacement rising to 3.5 Å) — so it needs no separate tail;
coulomb/PPPM are dropped for the movie (geometry is set by bonds, angles and
excluded volume).

## Provenance

Everything is **generated fresh** — no `tests/output/` golden is used. The three
arc animations are run through LAMMPS locally (via the movie decks); the still
panels are `03_Conformation` (lattice) state only.

| section | system | LAMMPS |
|---|---|---|
| CG arc (GIF) | `sculpt_250` — 5×5×5 SC sculpted 375 → 250 edges | yes (`movie_cg.in`) |
| atomistic arc (GIF) | `atom_sculpt` — 3×3×3 PDMS sculpted 81 → 54 edges | yes (`movie_atom.in`) |
| copolymer arc (GIF) + stills | `copoly_*` — 4×4×4 sculpted 192 → 128 edges | GIF yes (`movie_cg.in`); stills no |
| entanglements (stills) | `entangled_grafted` — 5×5×5 sculpted to 250 edges, 12 entanglements + grafts | no |

The earlier `v21_*` goldens are **not** used. The atomistic one is a perfect grid;
the CG one is stale — its `03_Conformation` carries 200 bonds far past the median
(190 at 7.45 σ = √3 × its spacing, 10 at 10.53 σ), bead-ends pinned to lattice
nodes, which fresh output never shows. (That golden *loads* a 210-edge sculpted
graph and its config no longer exists, so it can't be re-run for a controlled
comparison; the solid evidence it predates a change is that the bead-spacing
formula moved `a/dp` → `a/(dp+1)`.)

`tests/output/` is gitignored (`.gitignore:48`) and read-only (CLAUDE.md), so the
gallery deliberately depends on nothing in it.

## What the panels show

| panel | |
|---|---|
| `cg_lattice` / `cg_minimised` / `cg_equilibrated` | 5 125 beads, 125 junctions (functionality 2–6, mean **4.0**: 15/24/41/36/9), 250 strands of DP 20, 18.2 σ box. Median bond **0.17 → 0.96 → 0.97 σ** |
| `atom_lattice` / `atom_minimised` / `atom_equilibrated` | 10 896 atoms, DREIDING PDMS, 3×3×3 sculpted to mean degree 4.0, 53.2 Å box. Median bond **0.44 → 1.09 → 1.12 Å**; NPT settles the box to 52.5 Å |
| `copoly_{block,random,alternating}` | 4×4×4 sculpted (mean f 4.0), 128 strands of DP 20, A/B 50:50 |
| `ent_full` / `ent_zoom` | 250 strands, 12 entanglements → **24** kinked, 259 grafts; closest entangled pair **0.39 σ** |

## Facts worth keeping

Each cost a wrong picture first.

**The as-built lattice is 5× compressed, and that is the point, not a bug.**
topon strings DP beads along a lattice edge, so their spacing is
`edge_length / (dp+1)` ≈ 0.198 σ — against a harmonic equilibrium of **0.97**
(`bond_coeff 1 466.1 0.97`). The chains are taut wires; MD lets them coil out to
their real length while the junctions hold the graph.

**The CG "minimisation" does nothing at all — the panel is why we know.**
`system_after_soft` renders identical to the as-built lattice, and it is not
because stage 1 pins the nodes (it does `unfix freeze_nodes` before the final
free minimise). Every `minimize` in the CG deck stops after **1 iteration** on
"energy tolerance": `lammps_inputs.py` hardcodes `etol=1.0e-4`, LAMMPS tests
|ΔE|/(|E₁|+|E₂|), and at E ≈ 3.9e6 the first line search already scores ~4e-7.
Measured: max displacement from as-built **2e-4 σ**, zero atoms past 1e-3. Even
the soft push-off is a no-op. **All** of the 0.198 → 0.963 relaxation is the
`run 20000` NVE ramp in stage 2 (which ends at T = 18.4, so `system_ramped` is
not a minimised state either — the middle panel is captioned "pushed off &
ramped", not "minimised"). The atomistic deck uses `etol=1e-8` and genuinely
minimises.

**A melt cannot be seen into; pick the state, then the bead size.** At physical
bead size (r ≈ 0.5 σ) a Kremer–Grest melt is an opaque solid and no slab or
camera angle helps. All three CG panels use a *fixed* r = 0.11 σ, so what the
viewer sees change is the chains, not the radius. Bond width (0.22) is set ≥ the
0.198 bead spacing so a strand renders as a solid rod, not a dotted line.

**Chain-walking must strip the grafts first.** A backbone bead carrying a side
chain has degree 3. A walk that refuses degree-3 beads stops dead there and
splits one chain into two fragments — which then look like two different chains
"0.40 σ apart" when they are *literally the same beads*. Peeling degree-1
non-junction beads removes the dead-end side chains and leaves a clean degree-2
backbone. The check on the ideal (unsculpted) 5×5×5: `walk_chains` returns
**375 chains of exactly 20 beads, no bead used twice** — exactly its 375 edges.
On the sculpted `sculpt_250` it returns **250**, matching the sculpt target.

**Entanglement geometry is real.** On the sculpted entanglement network, ask for
12 and exactly **24** strands bow > 1 σ off their axis while the other 226 measure
**0.00**; the no-entanglement control gives 0. Two per entanglement is guaranteed
by `assignment/entanglements.py:229-231`, which skips any candidate whose edges
are already used — *not* by the count, and not by geometry (the closest crossing
varies run to run; this realisation was chosen from six for its tight 0.39 σ
pair). Cite the code, not the arithmetic.

**`gradient` ignores the composition.** `w = max(0, 1 - |t - pivot| * n)` scales
the window by `n`, so overlap needs `2/n > 1/(n-1)` ⟺ **n > 2**: with two
monomers the ranges [0, 0.5] and [0.5, 1] merely touch and nothing ever blends.
`weights.append(w * f)` cannot rescue it — multiplying by the fraction cannot
change *which* weight is zero, so **every** requested composition comes out 50:50
(asked A=0.1 → delivered A=0.50 over 200 seeds). It coincides with `block`
exactly when fractions are equal *and* DP is even; at odd DP the midpoint bead
hits the `sum(weights)==0` fallback and is drawn at random, so it matches only
~half the time. No gradient panel is shipped.

**Type IDs do not mean the same thing across systems.** `atomistic_combined` is
`1=Si 2=O 3=C 4=H`; the fPDMS box is `1=C 2=C 3=H 4=Si 5=F 6=O`. A fixed palette
painted every *hydrogen* blue-as-nitrogen and *carbon* gold-as-silicon.
`element_palette()` keys off each file's own mass table. For CG it is worse: in
the copolymer systems the junctions and the A monomers are **both type 1**, so
types alone cannot separate them — `node_ids()` reads topon's own
`system.groups` instead.

**"The beads look too small" was 20% of the bonds being missing.** The as-built
lattice puts nodes at x=0 — exactly *on* the periodic face — and the conformation
stage jitters every bead by ±1e-4, so wrapping throws the ones that land at -1e-4
across to x≈L while their neighbours stay at x≈0. Measured on a full 5×5×5 as-built
lattice: **2 761 of 7 625 beads** within 0.01 of a face and **1 596 of 7 875 bonds
(20%)** spanning the box, each drawn as two stubs at opposite faces. Whole strands
render as dashes, which reads as a bead-size problem — and no radius or bond width
fixes it, because the bonds genuinely are not there. `half_cell_offset()` shifts
the lattice planes into the interior; on the full lattice what remains is **75**
bonds, exactly the genuinely periodic ones (25 crossing edges per axis × 3).
*Lesson: when a render looks like a style problem, measure before restyling.*

**Wrapping breaks bonds unless the shifts are rebuilt.** OVITO fixes each bond's
PBC shift when the file is read; any later position change leaves them describing
the old geometry, and every boundary-crossing bond is drawn as a ray across the
box. `rebuild_bond_pbc` must follow every wrap — and every recentre.

**To zoom on a kink, recentre — don't crop on minimum-image distance.** Cropping
by MIC distance drags in beads from the far side of the box, which render as a
detached island floating beside the kink. Sliding the contact point to the box
centre *before* wrapping keeps the neighbourhood contiguous.
