# MARTINI 3 protein-network — devlog & comparison

## What we built

`topon/protein_network/` — a peer of `topon/simbox/` and `topon/singlechain/`. Takes a residue sequence (e.g., the resilin consensus block `GGRPSDSYGAPGGGN`) plus chain count and water density, generates a coarse-grained MARTINI 3 polymer-network LAMMPS input set complete with relaxation and (optional) annealing protocols. Mirrors the **topro** CHARMM atomistic shape (BFM lattice topology → JSON snapshot → chemistry build → 3-stage relaxation) but emits MARTINI 3 beads.

## Topon family comparison

| | core topon (DREIDING/KG) | topro (CHARMM atomistic) | this (MARTINI 3) |
|---|---|---|---|
| Resolution | atomistic / CG | atomistic | coarse-grained |
| Input shape | lattice config (Pydantic) | block_seq × n_repeats CLI | block_seq × n_repeats CLI (same as topro) |
| Force field | DREIDING / Kremer-Grest | CHARMM36 / 36m + CMAP | MARTINI 3 (vendored, auto-extracted from `Martini_Ahmet/nat_pro.itp`) |
| Topology | direct lattice graph | BFM SAW lattice + Union-Find gel-point | BFM SAW lattice (forked verbatim from topro) |
| Crosslinks | none / Y-merge in cg | dityrosine via `CE2-CE2` + DITY patch | dityrosine via `SC4-SC4` (TN6 bead) |
| Water | n/a | TIP3P + ions, random + SHAKE | MARTINI W (single bead, no SHAKE) |
| Relaxation | hierarchical freeze/unfreeze (3-stage CG) | flat soft-push + nve/limit + tight min | both modes available (`--hierarchical-stage1`) |
| Annealing | none built-in | 14-stage T/P schedule | 14-stage MARTINI port available as templates |
| Output `.data` | `atom_style full` (no image flags) | `atom_style full` (no image flags) | `atom_style full` **with** image flags |

## Key technical decisions

* **Library lookup, not auto_martini.** The 8 amino acids in the resilin reference (GLY, ALA, ARG, PRO, SER, ASP, TYR, ASN) are auto-extracted from `tests/_martini_extracted/Martini_Ahmet/itp_files/nat_pro.itp` by `tools/extract_residues_from_itp.py`. No SMILES → bead conversion, no external polyply call. Re-run the extractor against a different polyply ITP to add residues.
* **Image flags are mandatory at MARTINI scale.** Topro's `pos % box` wrap (`writer.py:220`) works for them only because `lattice_scale ~ 3 Å`; at MARTINI's ~27 Å, segments cross the periodic boundary often enough that bonds explode without per-atom image flags. We compute `floor(pos / box)` and write `ix iy iz` in the Atoms section.
* **Drop oversized crosslinks.** The BFM merge-site semantics let two chains' SC4 atoms land at the same lattice node from different walks (different image cells). Crosslinks where the two atoms are >box/4 apart unwrapped are dropped with a warning (typically 7 of 13 in our test setups). Intra-chain bonds always survive.
* **Reaction-field electrostatics approximated.** GROMACS RF (`coulombtype = reaction-field, ε_r=15`) has no exact LAMMPS equivalent. We approximate with `pair_style lj/cut/coul/cut 12.0` + `dielectric 15.0`, which loses the RF correction term but uses stock LAMMPS only.
* **angle_style cosine/squared as surrogate for funct=10.** MARTINI 3 IDP backbone uses GROMACS restricted-bending angle `½ K (cos θ − cos θ₀)² / sin²θ`. Stock LAMMPS lacks the `1/sin²θ` factor; we use `cosine/squared` which matches the leading term and is good enough for IDP relaxation.

## Tested

* 130 unit + regression tests (pin DREIDING / KG / simbox golden byte-equivalence, plus new MARTINI tests).
* 6 small demo systems (nat_pro / high_pro × w0 / w1 / w4) end-to-end through stages 1+2+3 in both `flat` and `hierarchical` modes — topology preserved (atom count, bond count, crosslinks, percolation all identical) across initial → after_soft → ramped → equilibrated.
* Full 14-stage annealing port runs end-to-end on the dry small case (~1 min wall time, all restarts written, system_annealed.data emitted).
* Full resilin reference (50 chains × 18 repeats, matching `Martini_Ahmet.zip`) — generated, ready for production runs.

## Known limitations / deferred

* **NPT collapses dry boxes** (no solvent pressure). Stage 3's NPT block is commented out by default. For dry runs, NVT only.
* **Reaction-field correction term not applied.** Acceptable for screened hydrogel regimes; for quantitative electrostatics, USER-MISC `coul/diel` would be needed.
* **Restricted-bending angle approximation.** Energy difference is small for IDP backbone but worth flagging if the user runs folded proteins.
* **Residue coverage limited to 8 AA.** Adding more requires re-running the extractor against a polyply ITP that contains them.
* **Ions (NaCl) not packed.** The NaCl 0.15 M packing in topro/atomistic isn't replicated. The vendored ion ITP is in place; just needs a packer.
* **No virtual sites, elastic network, or Go contacts.** First-cut targets IDP proteins.

## Pre-generated outputs

* `tests/output/v33_protein_network_demo_{nat_pro,high_pro}_16chain/{w0,w1,w4}/` — flat-stage-1 demos
* `tests/output/v34_protein_network_demo_hier_{nat_pro,high_pro}_16chain/{w0,w1,w4}/` — hierarchical-stage-1 demos
* `tests/output/v35_resilin_full_zip_reference/{w0,w4}/` — full resilin reference (50 × 18, matching the user-supplied zip)
