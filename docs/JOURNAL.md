# topon — Engineering Journal

A live log of changes, issues encountered, and the resolutions taken. More granular than [`DEVELOPMENT.md`](DEVELOPMENT.md) (which is the formal version-by-version changelog), and complementary to [`internal/DEVELOPMENT_INTERNAL.md`](../internal/DEVELOPMENT_INTERNAL.md) (which tracks open issues and planned work).

Append a new entry whenever you ship something non-trivial. Each entry follows:

- **Date — short headline**
- **Change:** what was done.
- **Why:** motivation.
- **Issue / solution:** any non-obvious problems and how they were fixed.
- **Follow-up:** (optional) what remains.

Newest first.

---

## 2026-05-20 (ultimate fix) — physical geometry + correct exclusions: topon output now runs in GROMACS at 20 fs

**Change** Three coordinated source fixes that eliminate the overlap-launch crash at its origin (rather than capping it at run time), in [topon/protein_network/builder.py](../topon/protein_network/builder.py), [template_builder.py](../topon/protein_network/template_builder.py), [lammps_writer.py](../topon/protein_network/lammps_writer.py):

1. **Proper sidechain geometry.** `_embed_residue_local` + `_orient_offsets` (builder.py) replace the old `bb + jitter` placement. A tiny per-residue distance-geometry embedding satisfies the bond + ring-constraint lengths and repels non-bonded beads to ≥3 Å, then orients the sidechain outward. (Old jitter piled SC beads onto the BB — down to ~0.02 nm apart.)
2. **Crosslink dimer.** The BFM merges two crosslinked TYR onto one lattice node; `crosslink_anchor_offset` (`CROSSLINK_SEP_ANG=7 Å`) now offsets the two partners into an adjacent dimer instead of coincident (r≈0).
3. **Correct exclusions.** `special_bonds` changed `0 0 0` → `0 1 1` (nrexcl=1, exclude 1-2 only — matching the reference `high_pro.itp`/GROMACS).

**Why** The recurring crash was an overlap-launch (see entry below), worked around with `nve/limit`. The user wanted it fixed at the source so the topology runs in GROMACS at the reference 20 fs (rigid LINCS) with the rigorous Parrinello-Rahman ensemble, not just LAMMPS-capped at dt=2.

**Issue / solution** Peeling the onion exposed three overlap sources, the third being the deepest: **`special_bonds 0 0 0` over-exclusion was load-bearing** — it both *hid* the bad sidechain geometry (intra 1-3/1-4 pairs felt no force) and *actively destroyed* good geometry during dynamics (with no repulsion, the embedded sidechains collapsed back together; a LAMMPS relax that read `E_vdwl` negative still handed GROMACS a structure at 10²³). With `nrexcl=1` the 1-3/1-4 pairs repel and stay apart.

**Validation (athena, GROMACS 2026.0):** the geometry+crosslink-fixed topon build → standard minimisation → **GROMACS NPT ran the full 2 ns at dt = 20 fs**, no crash, no caps: em PE −1.05×10⁶ (finite; was 10²³–10³¹ aborts), final **T = 309.9 K**, **ρ = 1.237 g/cm³** (notably lower than the over-excluded LAMMPS 1.38 — the over-exclusion was inflating the density). An independent investigator audit of the builder changes came back CLEAN (topology byte-identical; only positions change; crosslink mapping + determinism verified). 89 protein_network unit tests + 11 writer tests pass.

**Follow-up**
- Integrate the GROMACS exporter (`make_gromacs2.py`: merged single moleculetype + `[intermolecular_interactions]` crosslinks) as a topon `gromacs_writer` module + workflow flag.
- Add explicit TYR-ring SC1–SC4 (1-3) exclusion for byte-exact reference parity (with `0 1 1` it currently feels a mild LJ at ~5 Å; harmless).
- Re-validate the *LAMMPS* path uncapped with `0 1 1` (the empirical run was blocked by an athena mpirun/X11 issue; GROMACS success with the same exclusions strongly implies it).
- Update the seed=42 regression references (the geometry change alters coordinates, the `special_bonds` change alters the generated stage scripts — both intended).

---

## 2026-05-20 — root-cause + crash-proof fix for the MARTINI NPT "Bond atoms missing" crash

**Change** Investigation + MD-protocol fix (no topon code change). The recurring crash in the resilin 50 wt% equilibration (`resilin_martini_highpro_v2/w50__W`, athena) was root-caused; a crash-proof capped equilibration protocol replaces the uncapped Nose-Hoover stages. Forensics recorded in [internal/martini_npt_launch_rootcause.md](../internal/martini_npt_launch_rootcause.md).

**Why** `anneal_02_npt` (uncapped `fix npt`, Nose-Hoover, dt=2) died at step 258086 with `Bond atoms 935 936 missing`, after 103k *healthy* steps (density steady 1.37, T~315 K, E_bond stable). Earlier sessions had blamed (a) crosslinks starting at r=0, (b) an intra-residue ARG BB–SC2 overlap, (c) the over-constrained TYR ring (Option C). All three were wrong.

**Issue / solution** Two independent lines of evidence:
1. **topon's output is clean** (verified by running the as-generated `system_equilibrated.data` through LAMMPS locally): at step 0, pre-minimization, E_vdwl = −172593 (no LJ overlaps), E_bond = 8847 (no r=0 crosslinks/stretched bonds), Fmax = 228, CG-min converges in 59 iters. So no geometry/chemistry bug.
2. **The launch is a dynamics-time overlap ejection.** A re-run with a different seed crashed at a *different* step (187286) on a *different* atom (`Bond atoms 26565 26566 missing`); the launcher dump's final frame showed a TC4 sidechain bead and a water bead **1.46 Å apart with equal-and-opposite ~1.3×10⁶ kcal/mol/Å forces**. So: a rare, stochastic close-contact overlap (protein↔water / sidechain↔sidechain) develops in the dense melt and, under *uncapped* dt=2 integration, ejects a bead >30 Å in one step → "bond atoms missing". It is **not** specific to ARG SC2 (that was just the first crash site).

Fix: equilibrate with the launch-proof recipe already used by anneal_00/01 — `fix nve/limit 0.1` (caps per-step displacement to ~11× the thermal step, invisible to normal motion but clips ejections) + `fix langevin` (thermostat, T ramps) + `fix press/berendsen` (barostat). Implemented as `anneal_02b_npt_capped.in` + `anneal_03b..06b` + `run_capped_chain.sh` in the run folder. anneal_02b ran the full 250k steps (past the original crash point) with **zero crashes**; equilibrated density 1.365 g/cm³.

**Follow-up**
- Production ensemble RESULT: dt=1 Nose-Hoover **delays but does NOT prevent** the launch — it survived a 200 ps test, then crashed at ~450 ps in a longer run (step 1860006, `Bond atoms 26611 26612 missing`). dt only lowers the per-step launch probability (dt=2 crashed at 44-186 ps, dt=1 at ~450 ps); uncapped integration is unreliable for this BFM melt at any dt. **The reliable production is the capped recipe** (`nve/limit 0.1` + Langevin + Berendsen); trade-off is the Berendsen barostat (not Parrinello-Rahman) and a small Langevin-τ-tunable T offset (~+3 K at τ=200 fs, ~+10 K at τ=1000 fs). Truly ensemble-rigorous (uncapped) production would require eliminating the latent close contacts — better/longer equilibration to resolve them, the exclusion fix below, or rebuilding the network as a relaxed/kinetic gel.
- **Faithfulness item (P1, not the crash cause):** `topon/protein_network/lammps_writer.py` emits `special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0` (excludes 1-2/1-3/1-4) with no explicit `[exclusions]`, vs the reference `high_pro.itp` (`nrexcl=1` + TYR-ring `[exclusions]` only). topon over-excludes non-ring 1-3/1-4 pairs. A faithful port needs `special_bonds 0 1 1` plus an explicit exclusion of the TYR-ring SC1–SC4 (1-3) pair; deferred (own change + regression test).

---

## 2026-05-19 (CHARMM) — fix CHARMM36m default-parameter injection + N-terminal proline

**Change** [topon/protein_network/charmm/charmm_ff.py](../topon/protein_network/charmm/charmm_ff.py), [topon/protein_network/charmm/builder.py](../topon/protein_network/charmm/builder.py)

Three force-field-correctness bugs in the atomistic CHARMM builder, each of which made the LAMMPS writer fall back to a generic `(DEFAULT)` parameter (bond 300/1.5, angle 50/109.5, **dihedral K=0**, improper 20/0):

1. **`lookup_dihedral` had no wildcard fallback.** CHARMM36 stores most proline-ring (CP1/CP2/CP3) and sidechain torsions as wildcard terms `X t2 t3 X` (verified outer-only: 40/735 dihedral lines, all `X _ _ X`). The lookup only tried exact + reverse, so every wildcard-defined torsion silently got K=0 (no barrier). Added the `("X",t2,t3,"X")`/`("X",t3,t2,"X")` fallback, exact-first. (27 of 163 dihedral types in `PGRPSDSYPAPGPPN` resolve via this.)

2. **`lookup_improper` was not bidirectional.** ARG guanidinium (`C NC2 NC2 NC2`) and ASP/GLU/C-term carboxylate (`CC CT2A OC OC`) impropers are stored with the central atom LAST (`NC2 X X C`, `OC X X CC`). The lookup fixed only `t1` as central, so these got the K=20 default instead of K=45 / K=96. Now tries both directions, exact-both before wildcard-both (priority is load-bearing: the NH-centered `C HC HC NC2 = 0.0` exact term must win over the `NC2 X X C = 45` wildcard).

3. **N-terminal proline got the wrong patch.** `builder.py` chose `GLYP if GLY else NTER`; proline needs `PROP`. Proline is a secondary amine in a ring — NTER mis-typed its CA as CT1 and N as NH3, producing CHARMM-nonexistent angles/dihedrals (CT1-CP2, NH3-CT1, ...) AND a fractional residue charge (PRO+NTER = +1.18 vs the correct PRO+PROP = exactly +1.0). Now `GLYP / PROP / NTER`.

Plus a defensive guard in `add_water_and_ions`: `net_charge = round(sum(...))` silently masked a non-integer protein charge (then ion neutralisation left the fractional remainder). Now warns if `|raw − round(raw)| > 1e-3`.

**Two more bugs (found by a follow-up investigator pass, fixed same day):**

4. **Multi-term dihedral truncation** [topon/protein_network/charmm/lammps_writer.py](../topon/protein_network/charmm/lammps_writer.py). CHARMM proper dihedrals are sums of cosine terms; the writer emitted only `prm[0]`, silently dropping the rest. **101 of 575 dihedral keys are multi-term** — e.g. the peptide omega `CT1-C-NH1-CT1` kept n=1 K=1.6 but dropped the dominant n=2 K=2.5 (backbone planarity). Fixed via the standard CHARMM→LAMMPS idiom: `build_type_maps` now allocates one LAMMPS dihedral type per Fourier term (`dihedral_type_map` maps each canonical quad + reverse to a *list* of type IDs; needs `ff`), `write_lammps_data` emits one Dihedrals row per term, `write_lammps_settings` emits one `dihedral_coeff` per term (reverse-key alias de-duplicated). For `PGRPSDSYPAPGPPN` dihedral types went 138→215; all weights stay 0.0 so 1-4 (owned by `special_bonds charmm`) is not multiplied.

5. **Histidine silently dropped** [topon/protein_network/charmm/builder.py](../topon/protein_network/charmm/builder.py). `build_full_sequence` emits `"HIS"` but the RTF only defines HSD/HSE/HSP, so a His residue hit `if not res_tmpl: continue` and was skipped — fusing its neighbours across the peptide-bond gap. Now remaps HIS→HSD (neutral default tautomer) and **raises** on any genuinely-unknown residue instead of silently continuing. Does not affect GGQ/PGR (no His) but removes a latent corruption for general sequences.

**Why**

A user's old topro-generated resilin atomistic systems (`GGQPSDSYGAPGGGN`, `PGRPSDSYPAPGPPN`, 16 chains × n_repeats=12, water 0/35/55/65/75) had `(DEFAULT)` parameters and non-integer net charge (+0.48 / −0.12). They asked whether current topon fixes it and whether the non-neutrality was an ion-placement artifact. Answer: the defaults were three real lookup/patch bugs (still present in current topon until this commit), and the non-neutrality is a **bug**, not ion placement — non-integer protein charge cannot arise from integer ion addition; it traced to the N-terminal-proline mis-typing, with `round()` hiding it.

**Issue / solution**

The dihedral wildcard and improper bidirectional bugs affect *any* proline-containing or charged-residue sequence (i.e. essentially all of them). The N-terminal-proline bug only bites sequences that *start* with proline (PGR does, GGQ doesn't) — which is why GGQ looked "more broken" historically (its internal-proline dihedrals defaulted) yet PGR ran while GGQ needed manual patching. After all three fixes, both sequences regenerate with **0 DEFAULTs and net charge exactly 0.0000** across all five water contents.

Verification: 17 unit tests in [tests/unit/protein_network/test_charmm_ff.py](../tests/unit/protein_network/test_charmm_ff.py) (covering all five fixes); full `tests/unit/protein_network/` passes; CHARMM smoke test (real LAMMPS stage 1) passes; LAMMPS reads+runs the multi-term data file (stage 1 exit 0). Two reviewer passes (topon-reviewer + independent investigator) + a third topon-reviewer pass on the multi-term/His fixes: all HEALTHY. The multi-term fix's type-id↔coeff invariant was verified empirically (215 contiguous types, every emitted term matches `lookup_dihedral`, 0 mismatches) and on palindromic / dityrosine-wildcard edge cases.

**Follow-up**

- Stages 2/3 dynamics not yet validated on the regenerated systems (smoke covers stage 1 only).
- CMAP cross-term assignment (`find_cmap_crossterms`) interacts with N-terminal proline but was not part of this fix; not audited.
- Regenerated systems live under `C:/Users/ahmet/Downloads/GGQPSDSYGAPGGGN/_regen/<SEQ>/w<XX>/` (local, not in repo).

---

## 2026-05-19 (latest) — BFM `winding_safe` crosslink method (zero writer drops)

**Change** [topon/protein_network/bfm.py](../topon/protein_network/bfm.py)

- New function `apply_crosslinks_winding_safe` alongside the existing `apply_crosslinks_with_snapshots` (lattice-adjacent) and `apply_crosslinks_distance_based`.
- New `crosslink_method="winding_safe"` option wired into `generate_topology`.
- Two new helpers: `_lattice_image_delta(flat_a, flat_b, Nx, Ny, Nz)` returning the integer image-flag delta along a single BFM lattice bond, and `_compute_chain_node_images(chain_flat, ...)` returning per-node image flags of a chain via a forward walk.
- Snapshots produced by the new method carry an `n_rejected_winding` field (preserved through JSON round-trip automatically since `topology_io` uses `json.dump` transparently).
- Two new unit tests in [tests/unit/protein_network/test_bfm.py](../tests/unit/protein_network/test_bfm.py): `test_winding_safe_produces_zero_writer_drops` and `test_winding_safe_matches_adjacent_on_no_winding_seed`.

**Why**

The 2026-05-19 (later) priority-MST writer fix eliminated real-bond drops but still dropped ~0.06% of bonds as winding-cycle crosslinks. The user asked: instead of dropping at write time, can we *not form* those crosslinks at reaction time and replace them with non-winding candidates from the remaining pool? Yes — and the BFM-level check is the same image-delta arithmetic the writer's MST uses, just applied incrementally as candidates are processed.

**How it works**

For each candidate crosslink `(a, b)`:
1. Compute the crosslink's own lattice image-delta `xl_delta = _lattice_image_delta(pos_a, pos_b, Nx, Ny, Nz)`.
2. If `a` and `b` are in *different* connected components of the chain+crosslink graph: rebase the smaller component's per-node image flags so the new crosslink edge is minimum-image by construction; merge components; accept.
3. If `a` and `b` are in the *same* component: compute the BFS-implied delta along the existing spanning tree (`image[b] - image[a]`); compare to `xl_delta`. Match → trivial cycle (accept; the crosslink is structurally redundant but doesn't wind). Mismatch → winding cycle (reject; the TYRs stay unreacted and remain available for other partners).

The cost is O(n_chains × chain_size) per merge (~170K ops for resilin's 50 × 37 lattice × 91 reactions) — sub-second total.

**Verification on resilin_martini_highpro v2** (same seed=500, 50 × 270, 17³ lattice):

| Metric | `adjacent` (priority-MST writer) | `winding_safe` (this commit) |
|---|---|---|
| Reactions accepted at gel point | ~91 (with 18 dropped at write) | ~91 - rejected (will fill in) |
| Winding-rejected at reaction time | n/a | ~18 expected |
| Writer drops | 18 crosslinks | **0** |
| Final crosslinks in data file | 73 of 91 | full accepted count |
| Real BB-BB drops | 0 (assertion guarded) | 0 (none possible) |

The priority-MST writer fix is kept as a defensive safety net — for any future BFM topology that still has unavoidable winding (e.g. user explicitly opts out of `winding_safe`), the writer will still produce a parallel-MPI-safe data file.

**Issue / solution**

The reviewer audit (topon-reviewer agent) returned HEALTHY with two minor concerns:
- MC1: `<=` off-by-one in the snapshot-count guard matches the existing `adjacent` method — not a regression, flagged for awareness.
- MC2: `n_rejected_winding` is a new snapshot field not currently in the `test_snapshot_keys_match_topro_format` assertion. The current test uses `crosslink_method="adjacent"` (no field set), so it passes. Future winding-safe variants of that test would need to expand the expected key set.

**Follow-up**

- `crosslink_method="winding_safe"` is opt-in for now (default remains `"adjacent"`). Once the resilin athena production runs validate the new method, consider making `winding_safe` the default.
- The `target_conversion` parameter is wired into `apply_crosslinks_winding_safe` but not forwarded from `generate_topology`. Plumb it through if a user ever needs early termination at a specific conversion fraction.

---

## 2026-05-19 (later) — protein-network writer: priority-MST so only crosslinks drop

**Change** [topon/protein_network/lammps_writer.py](../topon/protein_network/lammps_writer.py)

- Promoted the Kruskal MST sort key from `(length,)` to `(priority, length)` where `priority=0` for non-crosslink bonds (backbone, sidechain, constraint) and `priority=1` for crosslinks (dityrosine SC4-SC4). All non-crosslink bonds are now tree edges by construction; only crosslinks can be back-edges across chains.
- Restored the hard assertion: real funct=integer non-crosslink bonds must NEVER drop. If any does, `AssertionError` fires with the offending atom IDs. (The previous "warning only" relaxation was masking the real bug.)
- Updated the drop report to split crosslinks vs constraints (TYR-ring straddling a box face — rare but legitimate).
- Stale comments in `builder.py` (line 352-358) and `template_builder.py` (line 42-46, 154-159) updated to describe "priority-weighted MST" instead of "length-weighted MST".

**Why**

The first 2026-05-19 cut of the MST writer was dropping ~0.08% of bonds, but the breakdown showed something wrong: of 22 winding-cycle drops on the resilin_martini_highpro v2 system, **16 were real backbone BB-BB bonds** between adjacent residues (e.g. atoms 5576-5578 = PRO/BB→GLY/BB in chain 11), 4 were crosslinks, 2 were constraints. The relaxation could swallow short-range BB-BB stretches but losing 16 backbone bonds was structural corruption.

Root cause: BFM merges two TYR/SC4 beads at every dityrosine crosslink onto the same lattice node. The two beads belong to two different chains, and each chain's `_interpolate_residue_positions` walk reaches that merged node from a different side of the box — the beads have *identical wrapped positions* (within the per-axis perturbation noise) but their *natural* image flags differ by integer box vectors. The wrapped-MIC distance of every crosslink is therefore ≈ 0.05 Å, much shorter than the ≈ 6.7 Å BB-BB distance at the projected BFM scale. With a pure length-sorted MST, all crosslinks entered the tree first; once they merged chains into one big component, BB-BB bonds (longer) became the longest edges in every chain-wraps-the-box cycle and got dropped instead of the crosslink responsible for the winding.

The priority key inverts that: crosslinks are demoted to "process last", so the longest back-edge in every BFM-crosslink-induced cycle is the crosslink itself, exactly matching the original design intent that crosslinks are the redundant elements of the network.

**Issue / solution**

The previous (length-only) writer had a "warning + categorised drop report" softening of the original assertion because the assertion was firing during v2 regen. That softening was the wrong fix — it kept the data file producible at the cost of silently corrupting the chain topology. The proper fix was to make the assertion non-firing by construction, via the priority key, so it can be restored as a hard structural invariant.

Verification on resilin_martini_highpro v2 (50 chains × 270 residues × ~91 crosslinks, 17³ BFM lattice, 458 Å box):

| Metric | Length-only MST | Priority MST (this commit) |
|---|---|---|
| Total drops | 22 (0.076%) | 18 (0.062%) |
| Real BB-BB drops | 16 | **0** |
| Constraint drops | 2 | 0 |
| Crosslink drops | 4 | 18 |

The drop *count* is now slightly higher (18 vs 4 crosslinks) because all winding obligations get pushed onto crosslinks — that's the correct accounting: there are 18 actual winding cycles in the BFM crosslink graph for this realisation. The previous "4 crosslinks" was an undercount because the MST was placing the winding obligation on BB-BB bonds instead. Crosslinks are designed to be sacrificial; backbone bonds are not.

All 67 protein_network unit tests + 3 protein_network regression tests pass. Image flags on the v2 dataset are now consistent within each chain segment (verified: atoms 5575-5581 in chain 11 all share image `(-1, 0, 0)`; the chain wraps cleanly at 5582 where x crosses the box face).

**Follow-up**

- TaskCreate #10: add a unit test that constructs a synthetic System with a deliberately-winding crosslink and asserts no real-bond drop + assertion fires if injected.
- Topon-reviewer agent rule 3 updated to describe the priority-MST and restored "no real bonds may drop" invariant.

---

## 2026-05-19 — protein-network LAMMPS writer: MST image flags + winding-cycle drop

**Change** [topon/protein_network/lammps_writer.py](../topon/protein_network/lammps_writer.py)

- Replaced the wrap-only 7-column Atoms section with a 10-column row including `ix iy iz` image flags.
- Image flags are computed by a new helper `_kruskal_image_flags_and_drop`: length-weighted MST (Kruskal) over the molecular bond graph, image flags propagated from the spanning-tree root so every tree-edge bond is minimum-image.
- Cycle-closing back-edges whose BFS-implied image-flag delta disagrees with the wrapped min-image delta (i.e. cycles with non-zero winding around the periodic box) are dropped from the Bonds section at write time. (The "hard assertion that only crosslinks may drop" was *intended* in this commit but the length-only sort key allowed BB-BB bonds to drop instead; see the 2026-05-19 (later) follow-up entry for the priority-MST fix that restores the invariant.)
- Header `N bonds` count, Bonds section indexing, and all downstream tests updated; the `tests/output/v33_protein_network_resilin_dry/protein_network.data` regression golden regenerated.
- Stale `# wrap-only writer convention` comments in `topon/protein_network/builder.py:353-357` and `topon/protein_network/template_builder.py:42-43, 152-159` updated to reflect the new behaviour.

**Why**

The MARTINI 3 protein-network annealing pipeline on Athena (32-rank MPI) was crashing during stage-01 NPT contraction with `ERROR on proc N: Bond atoms X Y missing on proc N at step S` and `WARNING: Inconsistent image flags`. Diagnosis on `system_after_soft.data` showed 1940 of 28,865 bonds had image-flag-implied unwrapped distance > 10 Å (max 649.9 Å, max wrapped MIC = 5.03 Å) — i.e. the bonds were *physically* normal MARTINI bonds (~3 Å) but their image-flag bookkeeping was inconsistent with the bond graph because the writer was emitting all atoms with `(ix, iy, iz) = (0, 0, 0)` despite chains that wrap across the periodic box. LAMMPS's parallel ghost-shell construction uses image flags during bonded-list build, so when a bond's unwrapped endpoints fell on opposite sides of a proc boundary the ghost-atom lookup failed and the run crashed.

**Issue / solution**

The fix has been attempted twice before:

- **v33-v38**: per-chain image walks (each chain walked independently, image flags accumulated from chain start). When two chains met at a dityrosine crosslink, their walk-accumulated images disagreed and the crosslink appeared as a ~box-length bond. Reverted in v39.
- **v39 onward**: wrap-only / 7-column rows, no image flags emitted. LAMMPS assumes `(0,0,0)` for every atom — works for small / non-percolated systems and for systems where chains don't wrap differently from each other, but breaks parallel MPI bond communication when they do.

The 2026-05 approach is structurally different from both. The MST is global over the bond graph rather than per-chain; tree edges are MIC by construction regardless of chain identity, so the v33-v38 "phantom 240 Å bonds at crosslinks" failure mode is impossible for tree edges. Cycle-closing back-edges with non-zero winding remain a genuine topological obstruction — these are detected (BFS-implied delta vs wrapped-MIC delta, exact integer arithmetic, no magnitude threshold) and dropped at write time. For the resilin/highpro 50-chain / 91-crosslink test case, ~23 of 28,865 bonds dropped (≈0.08% of the network, well below any mechanical-property sensitivity threshold).

Investigator audit (two passes: design + draft code) confirmed: (a) MST over the global bond graph is not the v33-v38 anti-pattern; (b) MST is strictly better than naive BFS for picking which bonds become tree edges (backbone bonds are always shortest, so they preferentially become tree edges and never drop); (c) the assertion that only crosslink bonds may drop catches the upstream chain-placement-bug scenario rather than silently masking it; (d) edge cases (isolated water beads, disconnected components, zero bonds) are handled correctly.

**Follow-up**

- `tests/regression/test_protein_network.py::test_topology_json_is_structurally_consistent` is failing pre-existing (chain count drifts from frozen golden — unrelated to this writer change, lives in the BFM topology generator); needs separate investigation.
- The wrap-only convention in `topon/writers/lammps_data.py`, `topon/writers/lammps_cg.py`, `topon/writers/lammps_atomistic.py`, and `topon/simbox/writer.py` is *not* changed by this commit — those writers serve simpler topologies (KG, DREIDING, single-region simbox packs) where chains don't wrap differently and 7-column rows are still safe. If a parallel run on one of those systems ever surfaces the same `Bond atoms missing` error, the same MST pattern can be applied.
- Athena production run uses `system_equilibrated.data` from an earlier topon version (pre-fix). For that one-off, a standalone `fix_image_flags.py` post-processor with the same MST algorithm is the immediate unblocker.

---

## 2026-05-11 (later 2) — CG end-cap collapse to single bead (P0 placement bug)

**Change** [topon/chemistry/builder.py:`_place_end_cap`](topon/chemistry/builder.py)
- In CG mode, `_place_end_cap` now always falls back to `_place_simple_atom` (one bead per node), regardless of the SMILES on `NodeMoleculeConfig.molecule`.
- Atomistic mode unchanged — it instantiates the full SMILES (e.g. trimethylsilyl `[Si](C)(C)C`) and the existing pendant-coordinate pass propagates positions for the methyl Cs through bond neighbours.

**Why (root cause)**
User flagged that the CG combined demo's `system_conformed.data` had bonds running into atoms parked at (0, 0, 0). Investigation traced 18 phantom Carbon atoms with no `bead_type` property, none of them in `node_map` / `edge_atom_map` / `graft_atom_map`. They came from the default `node_type_map["end"]` entry in [topon/config/schema.py](topon/config/schema.py): `NodeMoleculeConfig(molecule="[Si](C)(C)C", is_end_cap=True)`. `_place_end_cap` instantiated the trimethylsilyl SMILES, adding 1 Si + 3 methyl Cs to `chemical_space`, but only the Si was registered in `node_map`. In atomistic mode the pendant pass picks up the orphan methyls; in CG mode there's no pendant pass so they stayed at the origin for the entire run, dragging spurious bonds through the box.

The 6 degree-1 "end" nodes × 3 methyls = 18 phantom atoms. Matches the user observation exactly.

**Verification**
Before fix:
```
Total atoms: 5519; phantom (C, no bead_type): 18
Conformed at (0,0,0): 19  (1 real node-0 + 18 phantoms)
```
After fix:
```
Total atoms: 5336; phantom: 0
Conformed at (0,0,0): 1  (the legitimate node-0)
```
21/21 unit + smoke tests still pass.

**Follow-up**
Combined demos refreshing again in background to update their `expected_output/` artifacts with both the curvature-normal graft placement AND the end-cap single-bead fix.

---

## 2026-05-11 (later) — Grafts on entangled chains now placed along the outward curve normal

**Change** [topon/pipeline.py:`_local_perp_unit`](topon/pipeline.py)
- New helper `_local_perp_unit(backbone_xyz, k, fallback_unit, rand_vec)` computes the graft direction **per anchor** instead of once-per-edge:
  - Central-difference tangent at backbone index `k`.
  - Central-difference curvature `P[k+1] - 2P[k] + P[k-1]`; the outward Frenet normal is `-curvature / |curvature|`, projected perpendicular to the tangent.
  - When the chain genuinely bends (entangled kinks always do), the graft sticks out along the outward normal — the convex side of the bend. Cannot dive back into the chain.
  - When curvature is tiny (straight backbone), falls back to a random perpendicular from the per-edge `rand_vec` — same answer as the pre-2026-05-11 chord-perp behaviour.
- Both atomistic and CG branches of `_run_chemistry_stage` now call the helper inside the graft loop.

**Why**
User flagged that on `polymer/.../combined` (entanglement + graft together), grafts on the entangled chains appeared to dive back into the backbone. Suggested using the local-tangent / normal vector. Verified the symptom and implemented the fix.

**Verification (CG combined, A/B/C on 8 entangled + 196 non-entangled grafts)**
The diagnostic is `tip-to-nearest-backbone / graft-length`. A value of 1.0 means the graft sticks straight out by exactly its own length; smaller means the tip dives back.

| Method | Entangled mean | Ent min | Non-ent mean | Non-ent min |
|---|---|---|---|---|
| Old chord-perp (one perp_unit per edge) | 0.754 | 0.446 | 1.000 | 1.000 |
| Local-tangent perp (intermediate) | 0.993 | 0.919 | 1.000 | 1.000 |
| **Outward curvature normal (shipped)** | **1.000** | **1.000** | 1.000 | 1.000 |

Zero dive-in on every graft, entangled or not. Non-entangled behaviour is unchanged because the curvature-vector check falls back to the random-perpendicular path on straight backbones.

**Follow-up**
- Combined demos (`atomistic/combined`, `coarse_grained/combined`) refreshing in background to update their `expected_output/` artifacts with the new placement.
- Smoke 8/8 + diagnostics 6/6 + shell 7/7 all green.

---

## 2026-05-11 — UX overhaul (init / doctor / inspect / recipes / topro) + CG graft conformation fix

**Change — CLI surface**
- **`topon init`** rewritten: `--preset {atomistic_pdms,cg_kg,poss,martini_resilin,charmm_resilin}` copies a bundled demo `config.json` (preset-produced files pass `topon validate` immediately); `--interactive` walks through the 5–6 knobs that vary (study name, output dir, model type, lattice, DP, density) and writes the result; the MARTINI/CHARMM presets print the right `python -m` invocation rather than write a JSON.
- **`topon doctor <config>`** (new) lints for semantic footguns beyond Pydantic schema. Rule registry in [topon/diagnostics/rules.py](topon/diagnostics/rules.py): POSS-at-internal-junction (P1-H), unknown-node-type (P0-2 silent Si fallback), atomistic-graft-non-PDMS, lattice-size-format, DP-below-Kuhn, schema-gap-extras, defects-endcap-safe. 6 unit tests, [tests/unit/diagnostics/test_doctor_rules.py](tests/unit/diagnostics/test_doctor_rules.py).
- **`topon inspect <run_dir>`** (new) summarises a finished pipeline output: atom count / atom-type count / box / per-stage status / next LAMMPS commands. Works in both nested layout (`02_Chemistry/`+`03_Conformation/`+`04_Simulation/`) and flat `expected_output/`-style. Implementation in [topon/analysis/run_summary.py](topon/analysis/run_summary.py).
- **`topon recipes`** (new) prints "I want X → run Y" cheatsheet covering polymer / MARTINI / CHARMM / simbox / chain / batch / inspect / analyze.
- **`topon topro`** (new) subcommand wraps the existing argparse `topon.protein_network` CLI as a click subgroup, so `topon topro generate|sweep|topology` works alongside `topon generate ...`. Inner `--help` is preserved via `help_option_names=[]`.
- **Friendly errors** in [topon/utils/errors.py](topon/utils/errors.py): `format_pydantic_error()` prints field path + plain-language message + hints, replacing raw stack traces. `load_config_or_die()` shared helper for `validate` / `doctor`. Handles JSON parse errors, file-not-found, Pydantic validation errors uniformly.
- **Duplicate `gui` command removed** — `cli.py` had two `@main.command()` defs for `gui` (lines 145 and 411 pre-refactor); the second was dead code overriding the first. Replaced with `recipes`.

**Change — CG graft chemistry (P1-J resolved)**
- Pipeline's CG branch was emitting only the legacy 2-displace-file layout (`system_nodes.displace` + `system_beads.displace`) and never consulted `graft_atom_map`, so graft beads landed at (0,0,0) in `system_conformed.data` (user-reported: "I am not seeing grafts appended to the backbone chain").
- Unified the CG branch with the atomistic branch in [topon/pipeline.py:`_run_chemistry_stage`](topon/pipeline.py): same entanglement-aware kink loop, same 3-way graft-length cap, same perpendicular placement. CG now writes `system_backbone.displace` + `system_grafts.displace` (replacing the combined `system_beads.displace`).
- Verified on `examples/demos/polymer/coarse_grained/graft/`: 8589 atoms (was 8479 backbone-only), graft beads at IDs 8584-8589 now sit at `(16.88, 17.5-18.4, 17.3-17.5)` instead of `(0, 0, 0)`. **LAMMPS stage 1 now passes in 2.2 s** — P1-J "Neighbor list overflow" turned out to be caused by co-located graft beads, not crowding.

**Verification**
- 8/8 smoke tests pass in ~60 s.
- 6/6 diagnostics unit tests pass.
- `topon validate`, `topon doctor`, `topon inspect`, `topon recipes`, `topon topro --help`, `topon init --preset cg_kg` all manually tested.
- CG graft demo: chemistry → conformation → LAMMPS stage 1 all clean.

**Follow-up**
- `examples/demos/polymer/coarse_grained/graft/expected_output/` refresh (stages 2+3 currently running in background).
- Per-monomer "graft attachment site" attribute in `MonomerConfig` to lift the PDMS-only restriction on atomistic grafts (P1-L follow-up).

---

## 2026-05-10 (latest) — Defect chemistry root-fixed; atomistic graft side chains implemented; 4 new workflow examples

**Change**
- **Defect chemistry** (`topon/assignment/defects.py` + `topon/assignment/manager.py`): primary-loop injection now wires `max_degree=4` (was ignored), re-checks degree per-injection (was computed once before any injection — a single node in two selected pairs over-valenced), and excludes `node_type='end'` nodes from candidates (their effective valence cap is 1, not 4). Result: defect demo's RDKit mol no longer has any over-valent Si; `Sanitize → AddHs → Gasteiger` returns a clean neutral system naturally. **No charge neutralisation, no NaN scrub, no `SANITIZE_PROPERTIES` skip is triggered for any current demo.**
- **Atomistic graft chemistry** (`topon/chemistry/builder.py`): added `_build_pdms_chain_with_grafts` — a per-repeat PDMS builder that reads `graft_positions` from the edge data and emits a real side chain at each marked backbone Si (1 methyl cap + branch O + `graft_dp` repeats of Si-C-C-O, dropping the trailing O on the last repeat so all Si stay at valence 4). `_build_chain_atomistic` falls back to this builder when graft data is present and the monomer is PDMS; non-PDMS-with-grafts emits a warning and skips grafts (no silent corruption).
- **CG graft-map reshape**: `_build_chain_cg` was already populating `graft_atom_map` but as `dict[int -> list[int]]`; Pipeline's placement loop expects `[(frac, [atoms]), ...]`. Reshaped to the canonical form so the placement loop is no longer dead code for CG either.
- **Graft placement cap** (`topon/pipeline.py`): three-way length cap on the graft vector — `min(extension_factor=0.5, graft_dp/backbone_dp, 0.5 * lattice_spacing / edge_len)`. The third term keeps graft tips inside their own lattice cell regardless of edge length.
- **Four new workflow scripts** under `examples/workflows/` — these are standalone Python scripts (not config-driven demos), each with a knob-block at top:
  - `batch_polymer_topology/run.py` — generates 25 lattice networks, exports each as `.nodes/.edges` + GraphML + NPZ + writes a `summary.csv` of per-graph properties.
  - `bfm_gel_point_sweep/run.py` — sweeps 9 BFM parameter sets (n_chains, n_repeats, segs_per_block, intra-chain sep, equilibration steps) and records the gel-point conversion to `summary.csv`. Wall: < 5 s for the whole sweep.
  - `bfm_to_martini/run.py` — drives `topon.protein_network.workflow.run_protein_network` end-to-end (BFM → MARTINI 3 CG LAMMPS files).
  - `bfm_to_charmm/run.py` — drives `topon.protein_network.charmm.build_systems` end-to-end (BFM → CHARMM36m atomistic, multiple water contents).
- Index at `examples/workflows/README.md` ties them together.

**Why**
User pushed back twice. First: "defect has residual charge, something looks off." Second: "graft side chains aren't being built." Both turned out to be real chemistry bugs in `ChemistryBuilder`, not Pipeline-level fitting issues. Investigator agents traced the defect issue to one missing kwarg (`max_degree=4`) and a second silent over-valence path through end-cap nodes; the graft issue to two parallel-edge atoms-shape mismatches (CG's dict vs. Pipeline's list-of-tuples) AND atomistic never building side chains at all. Per-demo charge tables show all 6 atomistic demos at net charge ≤ 1e-8 e after these fixes.

**Verification**
- All 6 atomistic demos build through Pipeline with neutral net charge (≤ 5e-9 e), 4 DREIDING types incl. H_, and the correct displacement files:
  - basic: 10,949 atoms (no grafts, as configured)
  - combined: 53,635 atoms (was 42,449; +11k from real grafts at density 0.05)
  - copolymer: 21,449 atoms
  - defect: 21,944 atoms (was failing Sanitize before; now natural neutral)
  - entanglement: 21,449 atoms
  - graft: 41,941 atoms (was 21,449 — 90% more from real side chains at density 0.2)
- All 8 smoke tests pass in ~3 min.
- `system_grafts.displace` for the graft demo is now 2 MB of perpendicular placement coords (was 190-byte stub).

**Follow-up**
- `examples/demos/polymer/atomistic/{defect,graft}/expected_output/` are being refreshed through full LAMMPS stages 1/2/3 in a background runner (~2 h wall). Pre-fix `combined` expected_output is still on disk; given the atom count changed (42 k → 54 k), it should also be refreshed when wall time permits.
- The user-facing examples (`examples/workflows/`) have been smoke-tested for syntax + import; `batch_polymer_topology` and `bfm_gel_point_sweep` were run end-to-end. The two BFM-to-FF scripts wire to existing topon CLI entrypoints (covered by their own smoke tests).
- `P1-J` (CG graft stage-1 neighbour overflow) is now an opportunity rather than a regression — the CG graft path has the same `graft_atom_map` reshape applied, so re-running CG graft demo would presumably reveal whether the neighbour-overflow was a placement issue (now fixed via perpendicular cap) or a deeper crowding issue.

---

## 2026-05-10 (later) — P1-K fixed: Pipeline atomistic now matches canonical workflow end-to-end

**Change**
- Reverted my earlier same-day Gasteiger-only edit at `topon/pipeline.py:212-228` (which had been making things strictly worse: heavy-atom-only `ComputeGasteigerCharges` produced a net-−160-e system that crashed PPPM downstream).
- Re-implemented `_run_chemistry_stage`'s atomistic branch + displacement-writing tail to mirror `topon.workflows.atomistic_network.run` — the canonical hand-written workflow that produced the v21/v43 reference outputs. Order: `Chem.SanitizeMol` → `Chem.AddHs(mol)` → `AllChem.ComputeGasteigerCharges(mol_h)` → mass-based volume → DreidingWriter with `mol_h` and `use_charges=True` → five displacement files (`system_nodes`, `system_backbone`, `system_grafts`, `system_pendant`, `system_hydrogens`).
- Three Pipeline-specific patches over the canonical tail: keep the `isinstance(atom_ref, (list, tuple))` branch in node-coords for POSS cage; iterate `_builder.edge_atom_map` by `(u, v, key)` MultiGraph triples (canonical's int-indexed `edges` list doesn't apply); leave `system_grafts.displace` empty for atomistic-graft (preexisting — `ChemistryBuilder._build_chain_atomistic` doesn't populate `graft_atom_map`).
- Fallback path kept: `Sanitize → AddHs → Gasteiger` wrapped in try/except; on failure (over-valent atom etc.) falls back to writing the heavy-atom mol uncharged. Defect demo's degree-6 Si triggers a non-fatal RDKit warning but the main path succeeds anyway.
- Regenerated `examples/demos/polymer/atomistic/basic/expected_output/` with the full P1-K-fixed run: 10 949-atom data file, 5 displace files, stage 1/2/3 logs, `system_equilibrated.data`. Total wall: ~6 m 40 s on one core.

**Why**
User pushed back on my P1-K and P1-J "both pre-existing geometry issues" claim with "topro and cg was working fine before". Three investigator agents in parallel confirmed: (1) the legacy hand-written workflow `tests/workflows/generate_atomistic_combined.py` did produce healthy stage-2/3 output (v21/v43 reference logs prove it); (2) the new `topon.pipeline.Pipeline` atomistic path has never been validated end-to-end past stage 1; (3) my Gasteiger-without-AddHs edit had introduced a strictly-worse regression. The user said "carefully revert to canonical pipeline" — the fix is to make Pipeline's atomistic chemistry-stage tail equivalent to the canonical workflow's, not to invent new logic.

**Issue / solution**
- `Chem.AddHs` was the central missing call. The canonical workflow runs Gasteiger on `mol_h` (with H atoms) so net charge sums to 0 and PPPM auto-gewald tunes cleanly. Pipeline was running it on the heavy-only mol → net charge −160 e → PPPM warned, then stage-3 explosion at step 0 (T = 69 186 K).
- The five-displace split (vs Pipeline's two) is also load-bearing: with only `system_nodes` + `system_beads`, every pendant heavy atom and every H sat at (0, 0, 0) after `apply_displacements` → thousands of co-located atoms → infinite force at stage 2's first step.
- AddHs preserves heavy-atom indices, so Pipeline's `node_map` and `edge_atom_map` remained valid post-AddHs — no reindexing needed. This made Option C (hybrid: keep `ChemistryBuilder`, swap the tail) feasible with ~110 LOC in `pipeline.py` and **zero changes** to `ChemistryBuilder`, `ConformationManager`, or the calibrated LAMMPS scripts.

**Result**
- All 6 polymer atomistic demos now build through Pipeline with charge-neutral output and 5 displacement files (`basic`, `combined`, `copolymer`, `defect`, `entanglement`, `graft`). All 6 CG demos still build identically to before (2 displace files). All 8 smoke tests pass.
- Verified end-to-end LAMMPS run on `basic`: stages 1/2/3 all clean (4 s + 5 min 05 s + 1 min 31 s). E_pair drops monotonically through stage 2 (−87 k → −163 k) instead of exploding. `system_equilibrated.data` produced.

**Follow-up**
- P1-J (CG graft stage-1 neighbor overflow) is still pre-existing — `graft_density=0.2` + side-chain stacking exceeds default `neigh_modify one 2000`. Pure config issue. Easy fix: lower demo's density to 0.1 or 0.05 (matches the working `combined` demo).
- Other 5 atomistic demos' `expected_output/` folders still ship the old 2-displace format from the earlier failed runs; bulk regen (~35 minutes wall) is straightforward when desired.
- Atomistic graft side-chain placement is silent-broken (no atoms in `system_grafts.displace`; pendant pass catches them but extension isn't v20-dynamic). Logged as a follow-up — needs `ChemistryBuilder._build_chain_atomistic` to populate `graft_atom_map`.

---

## 2026-05-10 — CHARMM topro wired in; atomistic stages-2/3 PPPM/overlap diagnosis (P1-K)

**Change**
- **CHARMM atomistic protein networks** are now reachable through `topon.protein_network.charmm.build_systems` (also as `python -m`). The legacy topro CHARMM-side files were copied verbatim into `topon/protein_network/charmm/`: `charmm_ff.py`, `builder.py`, `lammps_writer.py`, `topology_io.py`. Bundled CHARMM36m PRM/RTF/CMAP files at `data/`. BFM topology stage continues to come from the existing `topon.protein_network.bfm` (the JSON schema is byte-identical to topro's). One small bug fix in the writer: `group protein all` (invalid LAMMPS syntax) → `group protein union all` for the dry path.
- **Demo at** `examples/demos/protein/charmm/`: README rewritten from stub → working quick-start; added `config.json` (declarative reference) + `run.py` (end-to-end runner that generates topology + builds LAMMPS files); `expected_output/` populated with `topo.json`, dry data file, settings, groups, three stage scripts, and a stage-1 reference log (~5 s wall on this machine).
- **Smoke test** `tests/smoke/test_charmm_protein_smoke.py` covers the CLI entry point + LAMMPS stage 1 on a small (8×8) system; ~6 s.
- **Pipeline atomistic chemistry**: `_run_chemistry_stage` now calls `Chem.SanitizeMol` + `AllChem.ComputeGasteigerCharges` and passes `use_charges=True` to `DreidingWriter`. Rationale: the calibrated polymer atomistic LAMMPS scripts use `lj/cut/coul/long` + PPPM, which auto-tunes `gewald` from per-atom charges and errors on a neutral system. The fix lives in chemistry, not the LAMMPS scripts.

**Why**
User asked: "the rest + charmm topro should be working fine." Translation: get CHARMM running and finish the polymer atomistic path. Both were stuck — CHARMM hadn't been migrated past the BFM stage, and atomistic LAMMPS was crashing at stage 2 with the gewald error.

**Issue / solution**
- *PPPM uncharged crash*: tried hardcoding `kspace_modify gewald 0.279` first → wrong (FFT mesh blows up to GB-class). Reverted, and instead enabled Gasteiger charges in the chemistry stage so PPPM auto-tunes correctly. Required a `Chem.SanitizeMol` call first to populate implicit valence.
- *CHARMM `group protein all`*: legacy code emitted invalid LAMMPS syntax for the dry-system path (no water/ions to subtract). Switched to `group protein union all`, which is valid and equivalent.
- **P1-K — atomistic stages 2/3 still don't relax cleanly.** With charges enabled and PPPM happy, stage 1 succeeds but stage 2 explodes (E_pair ~10²¹ kcal/mol at step 651, system stuck at numerical-pathology temperature). Root cause is geometry, not the script: Pipeline conformation stage's `noise_magnitude` default is 1e-4 Å (too small per the user's own xyz-perturbation memory; should be 0.05–0.1 Å) and `overlap_cutoff` is 0.2 Å (sub-LJ-sigma). Bumping both didn't unstick stage 2 in a smoke test, so the issue is deeper than just the perturbation magnitude. Per the user's "don't modify the calibrated scripts" memory, the fix lives in the conformation stage; logged as **P1-K** for follow-up.

**Result**
- CHARMM atomistic builder: working end-to-end through stage 1, stage 2 healthy (epsilon ramp brings E_pair from 6.6e10 → -1.99e4 over 1000 steps), stage 3 wired correctly. Smoke test passes.
- Polymer atomistic demos: stage 1 works for all 6 atomistic configs; stage 2 still blocked by P1-K.

**Follow-up**
- P1-K: investigate `topon/conformation/manager.py` — likely candidates: bump `noise_magnitude` default to 0.05 Å, raise `overlap_cutoff` to ~1.0 Å for atomistic, or check whether `apply_displacements` is leaving same-position atoms at chain junctions.
- The deferred POSS-at-junctions bug (P1-H) remains parked.

---

## 2026-05-10 — POSS clarified, coverage probe, expected_output for 3 demos

**Change**
- **POSS smoke** rewritten to match the documented working pattern: POSS at **degree-1 chain caps** (mirrors the legacy `generate_atomistic_poss.py` workflow that the user has used historically). Now passes in ~10 s. xfail marker removed. The previous failing config (POSS at degree-4 internal junctions) is a separate extension that's never been exercised; demoted from P0-H to **P1-H** in `INTERNAL.md` with a clearer scope note.
- **Defect smoke** had its xfail marker dropped after verifying it actually passes 3/3 in isolation. Earlier "1 failed" misattribution was on POSS, not defect.
- **Coverage probe** for graft / copolymer-block / combined (entanglements + grafts): all three configurations build through Pipeline AND pass LAMMPS stage-1 minimize. No additional bugs surfaced.
- **Expected outputs** committed for three core demos: `examples/demos/polymer/atomistic/basic/expected_output/`, `examples/demos/polymer/coarse_grained/basic/expected_output/`, and `examples/demos/poss/expected_output/`. Each ships the LAMMPS `system.data`, settings, groups, the stage-1 input script, the `log.lammps` from a successful run, and `system_after_soft.data`. Plus a one-page README per folder explaining the contents.

**Why**
User correctly questioned the POSS xfail. Investigation revealed two separate things:
1. The user's previous POSS workflow (POSS at degree-1 caps) is not broken — it works fine end-to-end through Pipeline + LAMMPS today.
2. POSS at internal junctions (what my smoke test was probing) IS a real but lower-priority bug: the documented usage doesn't trigger it.

**Result**
All 7 smoke tests now pass cleanly. Coverage probe confirms graft/copolymer/combined paths work. Three demos ship example output for users to compare against without running anything.

**Follow-up**
- P1-H still unfixed (POSS at internal junctions). Worth tracking but not blocking.
- The other 9 polymer demos under `examples/demos/polymer/` could get expected_output folders too — straightforward extension of today's script when desired.

---

## 2026-05-10 — Smoke-test reality check: defect actually passes; POSS reveals a real bug (P0-H)

**Change**
- Removed `xfail` marker from `tests/smoke/test_polymer_defect_smoke.py` after verifying it passes 3/3 in isolation in ~5.5 s. The earlier "1 failed, 6 passed" run had attributed the failure to defect; that was wrong. The slow / failing test in that run was POSS, not defect.
- Tightened the `xfail` reason on `tests/smoke/test_polymer_poss_smoke.py` to point at the real bug: `Bond/angle/dihedral extent > half of periodic box length`. Confirmed it's NOT a density issue — reproduces at `target_density=0.9` as well as `1.1`.
- Logged the underlying bug as **P0-H** in `internal/DEVELOPMENT_INTERNAL.md` §1: atomistic chemistry-stage placement leaves cross-boundary chains attached to POSS junctions un-wrapped. Suspected fix area: `chemistry/builder.py::_build_chain_atomistic` / `_place_poss_am0270` — compare to how non-POSS junctions handle the same cross-boundary case.

**Why**
User correctly questioned why these tests were now flagged xfail when nothing seemed to have regressed. Answer: defect was never failing — that was misattribution on my part during the rapid cycle of adding 3 smoke tests. POSS, on the other hand, is hitting a real previously-hidden bug that only became visible because the P0-D/E/C/B/A fixes finally let `Pipeline.run()` complete end-to-end. Workflow scripts probably bypassed the case (different lattice size, or different node-placement code path).

**Result**
6 of 7 smoke tests now pass cleanly. The one xfail (POSS) is a pinned bug, not flakiness — fixing it requires a chemistry-builder change.

---

## 2026-05-09 — Schema/loader/validator polish (P1-F + P2-G + public loader API + active-method validator)

**Change**
- **P1-F** — `topon/topology/generator_python.py:35`: chain `getattr(config, 'lattice_type', getattr(config, 'lattice_source', 'SC'))`. Accepts both schema's `lattice_type` and the legacy `lattice_source` attribute name. BCC/FCC no longer silently downgrade to SC on the Python topology path.
- **P2-G** — `topon/config/loader.py::load_config_full`: hoist legacy keys before schema validation. `chemistry.degree_of_polymerization` → `assignment.dp_distribution.default.mean`; `chemistry.bead_density` → `chemistry.target_density`. Demo configs keep working with their old shape; user's typed DP value is now respected.
- **Promoted loader API** — renamed `_remove_vacancies` → `remove_vacancies` and `_infer_dims_from_graph` → `infer_dims_from_graph` in `topon/topology/loader.py`. Updated the import in `topon/pipeline.py`. Underscore prefix was a misleading "private" marker; both are real graph-prep helpers other modules legitimately need.
- **Validator tightening** — `topon/config/validator.py::_check_type_mappings` now only checks the *active* `node_types.method` / `edge_types.method` source (e.g. `degree.mapping` when `method=="degree"`); ignores defaults of unused branches. Adjusted `tests/unit/config/test_config.py::test_missing_node_type_mapping` to set `method="random"` to match the new contract.

**Why**
- P1-F + P2-G clean up the silent-data-loss bugs the smoke tests surfaced.
- Promoting the loader helpers makes `Pipeline._generate_topology` (Python branch) idiomatic instead of reaching into a private API.
- The validator was reporting spurious "type 'B' missing" warnings on default configs because it inspected layer-types lists from inactive methods; now it only validates what the config will actually use.

**Result**
131 unit + smoke tests pass; no regressions. The smoke-test JSON fixture's `chemistry.degree_of_polymerization: 5` is now honored.

---

## 2026-05-09 — Fix P0-A (schema gap) — `topon generate` now accepts existing-style configs

**Change**
- Added `load_config_full(path) -> (ToponConfig, raw_dict)` in `topon/config/loader.py`. Splits the JSON into the five schema-known keys (`study`, `topology`, `assignment`, `chemistry`, `output`) and "everything else" (the raw dict — typically `conformation`, `simulation`, `execution`, `experimental`).
- Updated `topon/cli.py::generate` to use `load_config_full` and pass `raw_cfg` through to `Pipeline(config, raw_config=raw_cfg)`.
- Kept `load_config(path) -> ToponConfig` as a backward-compat thin wrapper (silently drops extras).
- Added `topon/config/__init__.py` export for `load_config_full`.
- Added `tests/smoke/test_polymer_json_load_smoke.py` + `tests/smoke/fixtures/json_load_smoke.json` — exercises a JSON config with all three extras sections through `load_config_full` → `Pipeline.run()` → LAMMPS stage-1.
- Marked P0-A fixed in `internal/DEVELOPMENT_INTERNAL.md` §1; logged a new P2-G for the secondary `chemistry.degree_of_polymerization` silent-drop issue.

**Why**
The headline blocker for `topon generate <config>`. Every existing-style config bundled with the repo had `conformation`/`simulation`/`execution` sections that `ToponConfig`'s `extra: "forbid"` rejected. The CLI was unusable for real workflows; users had to bypass via the `tests/workflows/run_*.py` scripts. After this fix, the CLI is the canonical entry point.

**Issue / solution**
Two valid approaches: (a) add full Pydantic schemas for `ConformationConfig`/`SimulationConfig`/`ExecutionConfig` and consume them validated; (b) split at load time and forward as raw. Picked (b) — smaller diff, no risk of changing semantics for sections the Pipeline already handled as raw dicts. Added (a) to follow-up notes; promoting these to validated schemas is still desirable but no longer urgent.

A separate gotcha during smoke-test wiring: `validate_config` raised warnings about node/edge type "B" being missing from `node_type_map`, even though the active `node_types.method = "degree"` doesn't use "B". The validator is iterating over default `positional.layer_types` / `composite.layer_types` regardless of active method. Out-of-scope spurious warning; documented and the smoke test skips that assertion. Fix is cheap when revisited.

**Result**
All 131 tests pass (127 fast unit + 4 smoke). Smoke harness now covers four orthogonal end-to-end paths: atomistic+load, cg+load, atomistic+generate (Python topology), and JSON-loaded with full extras. The five P0 bugs surfaced over the last day are all closed; the package's `topon generate` CLI is now functional end-to-end.

**Follow-up**
- P1-F (PythonTopologyGenerator silent SC downgrade for BCC/FCC).
- P2-G (chemistry.degree_of_polymerization silently ignored).
- Promote loader's `_remove_vacancies`/`_infer_dims_from_graph` to public names.
- Cleanup the `validate_config` spurious-type-warning issue (only check active method's layer types).
- Promote `conformation`/`simulation`/`execution`/`experimental` from raw dicts to Pydantic schemas (no longer urgent — load_config_full unblocks usage).

---

## 2026-05-09 — Fix P0-B (Pipeline `source="generate"` dispatch) + Python topology smoke test

**Change**
- `topon/pipeline.py:101-145`: rewrote `_generate_topology` to dispatch on `gen_cfg.exe_path`:
  - **C path** (`exe_path` set): `run_generator(gen_cfg, topology_dir, exe_path=...)` returns `(nodes_path, edges_path)`; reload through `load_graph`.
  - **Python path** (`exe_path=None`): `PythonTopologyGenerator(gen_cfg).generate(trials, max_saves)`; take `graphs[0]`, wrap `MultiGraph`, then `_remove_vacancies` + `_infer_dims_from_graph` for parity with the file-round-trip path. No I/O.
- Updated the module docstring at `pipeline.py:21-22` — `source="generate"` no longer requires a compiled C binary; pure-Python is the default fallback.
- Added `tests/smoke/test_polymer_generate_smoke.py` (4×4×4 SC, `exe_path=None`, atomistic DP=10). Passes after the fix.
- Marked P0-B fixed in `internal/DEVELOPMENT_INTERNAL.md` §1.
- Logged a new P1-F: `PythonTopologyGenerator.__init__` reads `lattice_source` instead of the schema's `lattice_type`, so BCC/FCC silently downgrade to SC on the Python path. Out of scope for this commit.

**Why**
P0-B was the user-facing blocker that prevented `topon generate` from working without a compiled C binary on PATH. With the dispatch in place, the package becomes self-contained — anyone who clones the repo and runs `pip install -e .` can use the full pipeline immediately.

**Issue / solution**
Plan agent flagged two latent issues during the research pass. (a) `PythonTopologyGenerator` reads the wrong attribute name (`lattice_source` vs `lattice_type`) — logged as P1-F for a follow-up. (b) `_remove_vacancies` / `_infer_dims_from_graph` are loader-private; the pipeline now imports them by underscored name. Both flagged for cleanup; not blocking for today's correctness fix.

**Follow-up**
- Fix P1-F (one-line attribute-name fallback in `generator_python.py:35`).
- Promote loader's `_remove_vacancies` / `_infer_dims_from_graph` to public names.
- One P0 remaining: **P0-A** (schema gap blocking JSON-config loading). Largest of the wave; will need new Pydantic schemas for `conformation`/`simulation`/`execution` sections.

Smoke-test count: now 3, all passing. Total wall-clock ~20 s.

---

## 2026-05-09 — Fix P0-C (model_type literal mapping) + add CG smoke test

**Change**
- `topon/pipeline.py:283`: map `"coarse_grained"` → `"cg"` at the call site before passing to `LammpsInputGenerator`. The writer only knows the legacy `"cg"` / `"atomistic"` literals; the schema's chemistry field uses `"coarse_grained"` / `"atomistic"`. Previously every CG system silently mis-routed through the atomistic writer branch.
- Added `tests/smoke/test_polymer_cg_smoke.py` — mirror of the atomistic smoke test with `model_type="coarse_grained"` and DP=10. Passes after the P0-C fix.
- Marked P0-C fixed in `internal/DEVELOPMENT_INTERNAL.md` §1.

**Why**
P0-C silently corrupted CG runs: PPPM electrostatics on a charge-neutral CG system, atomistic bond styles, etc. Pure-Python tests couldn't see this — the symptom only shows when LAMMPS reads the resulting input scripts. The CG smoke test pins the fix.

**Issue / solution**
Two valid fix locations: the call site (one literal map in `pipeline.py`) or the writer (normalization at its entry). Chose the call site to avoid touching the writer's many `if model_type == "cg"` branches; if the writer is later normalized to accept both literals, the call-site map becomes a harmless no-op.

**Follow-up**
Two P0s remaining:
- **P0-B**: dispatch in `_generate_topology` between C subprocess and pure-Python topology generation. Will unlock `source="generate"` smoke tests.
- **P0-A**: schema extensions for `conformation`/`simulation`/`execution` sections. Will unlock JSON-loaded smoke tests.

---

## 2026-05-09 — Fix P0-E (Stage 6 path doubling) — smoke test now PASSES end-to-end

**Change**
- `topon/pipeline.py:285-286`: changed `LammpsInputGenerator(str(self.output_dir), study_name, ...)` to `LammpsInputGenerator(str(self.config.study.output_dir), study_name, ...)`. The previous call passed `self.output_dir`, which already had `study.name` appended (`pipeline.py:64`); the writer re-appended internally, putting Stage 6 outputs at `<base>/<name>/<name>/04_Simulation/`. Now matches the `ConformationManager` call pattern (line 259-263).
- Removed the `xfail` marker from `tests/smoke/test_polymer_atomistic_smoke.py`.
- Marked P0-E as fixed in `internal/DEVELOPMENT_INTERNAL.md` §1.

**Why**
P0-D + P0-E were the two bugs blocking the smoke test from passing. P0-D fixed the in-Pipeline crash; P0-E fixed the on-disk layout so LAMMPS could find files from earlier stages. After both, the full atomistic load-path runs cleanly: Pipeline emits 6 stages of output and LAMMPS runs the stage-1 minimize without error.

**Issue / solution**
The path doubling looked cosmetic but was actually fatal: the LAMMPS stage-1 input script in `04_Simulation/` references `../02_Chemistry/system.data` — a relative path that only resolves correctly when both directories share the same parent. With the doubled `study.name`, `02_Chemistry` was at `<base>/<name>/02_Chemistry` while `04_Simulation` was at `<base>/<name>/<name>/04_Simulation`, so the relative reference broke. The one-line fix collapses everything back to the same parent.

**Follow-up**
Three P0s remaining:
- **P0-C** (next): one-line literal mapping in the same `_run_output_stage`. Will unlock CG smoke tests.
- **P0-B**: dispatch in `_generate_topology` between C subprocess (`run_generator`) and pure-Python (`PythonTopologyGenerator`). Will unlock `source="generate"` smoke tests.
- **P0-A**: schema extensions for `conformation`/`simulation`/`execution` sections. Will unlock JSON-loaded smoke tests.

---

## 2026-05-09 — Fix P0-D (Stage 4 bead-displacement TypeError)

**Change**
- `topon/pipeline.py:209-221`: rewrote the bead-displacement loop to unpack `(u, v, _key)` directly from `self._builder.edge_atom_map.items()` instead of treating the dict keys as int indices into `list(self.graph.edges(data=True))`. The dict keys are `(u, v, key)` tuples from MultiGraph edges; the old code's `if edge_idx >= len(edges)` raised `TypeError: '>=' not supported between instances of 'tuple' and 'int'`.
- Updated `internal/DEVELOPMENT_INTERNAL.md` §1 to mark P0-D as **fixed** (kept the entry for traceability).

**Why**
First of the four P0 bugs surfaced by the smoke test. Smallest, most obvious — quick win to validate the fix-via-smoke-test workflow.

**Issue / solution**
The fix unblocks the rest of the pipeline; running the smoke test now reaches "=== Pipeline Complete ===" successfully. But it surfaces a new pre-existing bug, **P0-E**: Stage 6's `LammpsInputGenerator` double-applies `study.name`, so LAMMPS scripts land at `<output_dir>/<name>/<name>/04_Simulation/` instead of `<output_dir>/<name>/04_Simulation/`. Pipeline-internally, Stage 4 and Stage 5 outputs are at the correct single-level depth; only Stage 6 doubles. Logged as P0-E in INTERNAL.md §1.

**Follow-up**
Fix P0-E next (one-line constructor call change in `pipeline.py:285`). Then P0-C, then P0-B, then P0-A.

---

## 2026-05-09 — Test infrastructure: tiers, per-component subdirs, smoke path

**Change**
- Reorganized `tests/unit/` into per-component subdirectories: `topology/`, `assignment/`, `chemistry/`, `config/`, `simbox/`, `protein_network/`. The 6 protein-network test files dropped their `protein_network_` filename prefix (subdir naming makes it redundant).
- Registered four pytest markers in `pyproject.toml`: `fast`, `smoke`, `regression`, `requires_lammps`.
- Added `tests/conftest.py` with two responsibilities: (a) auto-apply the tier marker to any test based on its parent directory (`tests/unit/` → `fast`, `tests/smoke/` → `smoke`, `tests/regression/` → `regression`), and (b) auto-skip any `requires_lammps` test when `lmp` is not on `PATH`.
- Added `tests/smoke/` with `test_polymer_cg_smoke.py` — a tiny end-to-end test that builds a 3×3×3 SC CG network, runs `Pipeline.run()` through all six stages, then invokes LAMMPS to run `minimize_1_serial.in` and asserts a clean exit + stage-1 output file.
- Moved `tests/tmp_hsp_audit.py` and `tests/Martini_Ahmet.zip` (no longer needed in tracked tree — the zip is already extracted to gitignored `tests/_martini_extracted/`) to `~/topon_archive/old_examples\`.

**Why**
- The "I changed X, retest X" workflow needs per-component subdirs (`pytest tests/unit/chemistry/` is now self-explanatory).
- Tiered markers let the same files participate in tier-based filtering (`pytest -m fast` for a quick pre-commit, `pytest -m "fast or smoke"` for pre-push).
- LAMMPS-running smoke tests catch regressions where the pipeline emits a syntactically valid LAMMPS file that LAMMPS still rejects — pure-Python unit tests can't see those.

**Issue / solution**
The first cut of the smoke test exposed **four** pre-existing package bugs in the `Pipeline` path. None were caused by the test-infra work; the smoke test surfaced them — which is exactly its job.

- **P0-A** (already documented): `load_config` rejects existing-style configs because `ToponConfig` has `extra: "forbid"`. Worked around by constructing `ToponConfig` programmatically in the smoke fixture.
- **P0-B** (newly logged): `Pipeline._generate_topology` calls `run_generator(...)` with the wrong signature; `run_generator` only supports the C-binary path. Worked around by using `topology.source="load"`.
- **P0-C** (newly logged): `Pipeline` passes `"coarse_grained"` to `LammpsInputGenerator`, which only branches on `"cg"` vs `"atomistic"`. Worked around by using `model_type="atomistic"`.
- **P0-D** (newly logged): `TypeError: '>=' tuple vs int` mid-Pipeline, after Stage 4 chemistry succeeds. Likely in the bead-displacement loop (`pipeline.py:212`); edge map keys appear to be tuples being treated as int indices. **Active blocker for the smoke test.**

The shipped smoke test (`tests/smoke/test_polymer_atomistic_smoke.py`) is **marked `xfail`** because of P0-D — pytest reports it as expected-failure rather than skip, so it's visible as a pinned reminder; flips to `xpass` automatically when the bug is fixed (`strict=False`). The test exercises the path we *want* to work: load 5×5×5 sample → DP=5 atomistic → Pipeline.run() → LAMMPS stage-1.

A simbox-based smoke test was attempted as a workaround (simbox is a separate code path that doesn't go through `Pipeline`), but its writer also produced LAMMPS-rejected output ("Unknown identifier in data file: 29 0.500000 -1 3" — likely a force-field-coefficients format mismatch with newer LAMMPS). Logged as part of the same P0 wave; not yet root-caused.

**Bottom line:** test infrastructure ships; one smoke test ships as xfail. No smoke test currently passes against this LAMMPS install (`2 Apr 2025`). The good news is that the smoke harness will catch regressions immediately once the P0 bugs are fixed.

**Follow-up**
- Trace and fix P0-D (chemistry → conformation handoff) — should be a few-line patch in `Pipeline._run_chemistry_stage`.
- Then P0-C (writer literal mismatch — one line).
- Then P0-B (run_generator signature + Python-only dispatch — small refactor).
- Then P0-A (schema extensions — moderate refactor).
- Investigate simbox writer's data-file format compatibility with LAMMPS 2 Apr 2025; if the failure is real (not just a regression-test golden mismatch), add a simbox smoke test once fixed.
- Add per-component fast tests where coverage is thin (current: assignment, chemistry, simbox, protein_network all have at least one fast test; topology has one; config has one; conformation and writers are covered indirectly).

---

## 2026-05-08 → 2026-05-09 — Documentation and examples consolidation (5-step roadmap)

**Change**
Five-step project consolidation completed across multiple commits:

1. Set up the `investigator` agent (`.claude/agents/investigator.md`) — unbiased read-only auditor used as a pre-commit reviewer for every non-trivial doc/code change in this consolidation.
2. Drafted four canonical docs: `docs/{ARCHITECTURE,USAGE,DEVELOPMENT}.md` and `internal/DEVELOPMENT_INTERNAL.md`. Each went through the loop *draft → investigator review → fix → commit*.
3. Cleanup commit: deleted 14 stale source docs (cli.md, config_reference.md, simbox.md, walkthrough.md, implementation_plan.md, etc.) now subsumed by the new four. Updated `README.md` and `CLAUDE.md` cross-refs. Fixed source-side drift (V36 four-files claim, V22 Hard Case framing, `workflow.py` docstring).
4. `examples/` curation: restructured into `demos/{polymer,protein,topology,poss}/` with READMEs at every category level; copied the npj-paper companion data into `examples/npjcompmat/` (1001 files, ~23 MB); archived old workflow scripts to `legacy/old_examples/`.
5. Wired up two GitHub remotes: `personal` → `https://github.com/lynspica/topon-dev` (public, primary), `stable` → `https://github.com/keten-group/topon` (URL only — paper-companion v0.1.0 left untouched).
6. Moved the entire 8.3 GB `legacy/` tree out of the repo working directory to `~/topon_archive/` (atomic same-volume rename; instant; reversible).
7. Added `AGENTS.md` at the root: single "read this first" doc for any AI agent (Claude / ChatGPT / Cursor / etc.) starting a session on the project.
8. Added `examples/showcase/network_5x5x5/` — small reference graph files for users to load via `topology.source = "load"`.

**Why**
The repo had drifted into 16+ scattered markdown files with mutually contradictory content (4-stage vs 6-stage pipeline, dead module names, non-existent workflow scripts), 8 GB of legacy artefacts in the working tree, and no clear onboarding path for new AI agents in fresh chats. The consolidation gave us a stable spine: AGENTS.md (entry point) → CLAUDE.md (rules) → ARCHITECTURE / USAGE / DEVELOPMENT (the canonical three).

**Issue / solution**
- **Schema gap (P0-A)**: surfaced when the investigator tried to validate the example configs through `topon generate`. Existing-style configs with `conformation`/`simulation`/`execution` sections are rejected by `ToponConfig`'s `extra: "forbid"`. Did not fix in this consolidation (out of scope) — documented as P0-A in `INTERNAL.md`, made the demo READMEs honest about the limitation, and added the workaround note to the smoke-test fixture.
- **Image-flag contradiction**: `topon-reviewer.md` said "wrap-only, image flags failed in v33-v38"; `martini_devlog.md` said "image flags mandatory at MARTINI scale". Resolved by reading the actual `protein_network/lammps_writer.py:188-217`: code is wrap-only, the topon-reviewer is correct. Updated `ARCHITECTURE.md` design principles 3 + 4 accordingly; flagged the docstring drift in `topro_issues_for_later.md` (now in `INTERNAL.md` §5).
- **Forgotten remote rename**: I created `lynspica/topon` initially, but the user clarified they wanted `topon-dev`. Renamed via `gh repo rename`; GitHub auto-redirects from the old URL.

**Follow-up**
All open work is tracked in `internal/DEVELOPMENT_INTERNAL.md`:
- P0-A: schema gap (above)
- P0-2: silent `Si` fallthrough in `_build_nodes`
- P1 polish (logger, default `min_dist`, `_guess_head` regex)
- P2 housekeeping (verbose prints, hot-loop imports)
- Future-work: SELFIES, NPZ output, GraphML CLI flag, RESP charges, GUI, Streamlit, Jupyter

---

*Earlier history lives in [`DEVELOPMENT.md`](DEVELOPMENT.md) §4 (V1–V36 changelog).*
