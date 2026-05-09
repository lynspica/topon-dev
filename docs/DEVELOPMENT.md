# topon — Development

The history, methodology, and conventions for changes to the topon package.

For *what's planned next* (open issues, in-flight work, pitfalls collected during development), see [`internal/DEVELOPMENT_INTERNAL.md`](../internal/DEVELOPMENT_INTERNAL.md) — gitignored, local only.

For the package layout, see [ARCHITECTURE.md](ARCHITECTURE.md). For running it, [USAGE.md](USAGE.md).

---

## 1. Methodology — Questions → Research → Requirements → Roadmap

Every non-trivial change walks through these four steps. The point is to *make the change concrete before writing code*, so the diff is small and the regression boundary is obvious.

1. **Questions** — what is unclear? What assumption could be wrong? What would make this concrete? Write them down. Most "obvious" changes don't survive ten honest questions.
2. **Research** — read the relevant code and docs; cite `file:line`; check the regression suite to know what's currently locked. If the assumption from step 1 turns out wrong, restart.
3. **Requirements** — what the change *must* do, what it *must not* do, what the regression boundary is. Write the test before the code where possible.
4. **Roadmap** — sequence the change into reviewable steps. Each step is a separate commit. If you cannot describe a step in one sentence, split it.

Two project-specific anchors:

- **`CLAUDE.md`** at the repo root encodes the rules every change must follow (module boundaries, regression test requirements, output conventions, no-globals rule).
- **The investigator agent** (`.claude/agents/investigator.md`) does an unbiased review pass on any non-trivial change. The pattern is *draft → investigator review → fix → commit*. The `topon-reviewer` agent is the complement for code-level fixes.

The methodology has worked well in recent versions — V35 (peroxide fix) and V36 (MARTINI 3 protein network) were both shipped with this loop. The investigator pass catches stale module names, contradictory claims between code and docs, and silently dropped facts during consolidations.

---

## 2. Recent changes — current state

The latest five versions, in reverse chronological order, with full detail in §4 below.

| Version | Date | Summary |
|---|---|---|
| **V36** | 2026-05-07 | `topon/protein_network/` — MARTINI 3 protein-network generator (sequence → LAMMPS), the implementation behind **topro**. Peer of `simbox`/`singlechain`; does not plug into `Pipeline`. 65 unit-test assertions + 1 regression test added. |
| **V35** | 2026-05-01 | Critical peroxide-bond fix in `chemistry/builder.py`. Resolved §P0-1 of the solubility-downstream issue list. Reference outputs under `tests/output/solvent_effects/{PDMS,PTFPMS}/` will diff after this change. |
| **V34** | 2026-03 | Multi-solvent mixtures in `singlechain/`. HSP solubility module (`singlechain/solubility.py`). Better head/tail detection for non-Si polymers (EPDM, PIB, FKM, NBR). |
| **V33** | 2026-03 | Package hardening + single-chain. BCC/FCC lattices wired in Python generator. Grafts and copolymers in `AssignmentManager`. `topon simbox` / `topon analyze` / `topon chain` CLIs. `CLAUDE.md` written. |
| **V32** | 2026-03 | simbox v3 — full LAMMPS data + input script generation for crosslinking studies; POSS sweeps with 7 compositions. |

---

## 3. Closed milestones

Phases 1 through 13 from the original task tracker are largely complete. Items still open: a polyply integration placeholder (Phase 4), externalising LAMMPS parameters to `experimental.json` (Phase 6), the GUI / Streamlit / notebook stack (Phases 7 and 11), and broader documentation work — usage examples and API docs (Phase 10). The closed work covers:

- **Core package structure** (Phases 1–2) — module layout, `pyproject.toml`, Pydantic schema, config validation.
- **Topology and assignment** (Phase 3) — C generator subprocess wrapper, graph loaders, node/edge type assignment, DP distribution (Schulz-Zimm), entanglement selection, defect injection.
- **Chemistry** (Phase 4) — `ChemistryBuilder`, SMILES library, atomistic (DREIDING) and CG (Kremer-Grest) modes, copolymers, grafts, edge-type mapping. *Polyply integration placeholder remains open.*
- **Conformation** (Phase 5) — displacement-file generation, entanglement kink geometry, coordinate wrapping.
- **Writers** (Phase 6) — `DreidingWriter`, `CGWriter`, `ConformationManager`, `LammpsInputGenerator`. *Externalising LAMMPS parameters to `experimental.json` remains open.*
- **CLI** (Phase 7) — `topon generate`, `validate`, `init`, `simbox`, `analyze`, `chain`. *GUI/notebook items remain open.*
- **Python topology port** (Phase 8) — pure-Python port of the C generator (`generator_python.py`). Validated against C output. Solved the "Hard Case" benchmark in 0.23 s; 60 % success rate on a 20-case batch with ~0.3 s median success time.
- **Martini 3 experiment** (Phase 9, **archived**) — `chemistry.martini` and Polyply integration shelved at V24.2 due to parameter-mixing complexity. Material moved to `legacy/older_versions/martini_experiment_2026/`. The current MARTINI work in `protein_network/` (V36) takes a different path: library lookup from a polyply-generated ITP, no automated mixing.
- **Advanced features** (Phase 12) — GraphML export with entanglement edges, multi-entanglement, custom topology input, POSS junctions (full chemistry + sweeps).
- **simbox sub-system** (Phase 13) — `Molecule`, `BoxPacker`, `AssembledSystem`, `MoleculeLibrary`, writers, CLI, regression tests.

The 2026-03-23 code audit (eight issues, severity CRITICAL through MEDIUM) is fully resolved — see V33 / V35 changelog entries.

Open phases and planned next steps are tracked in [`internal/DEVELOPMENT_INTERNAL.md`](../internal/DEVELOPMENT_INTERNAL.md).

---

## 4. Changelog (V1 – V36)

Notable changes are documented in reverse chronological order.

### [V36] — 2026-05-07 — MARTINI 3 protein-network generator (topro)

#### Added
- **`topon/protein_network/`** — new peer package alongside `simbox/` and `singlechain/`. Sequence-string-to-LAMMPS MARTINI 3 generator for protein networks with stochastic dityrosine crosslinks. Mirrors the legacy CHARMM-atomistic shape (BFM cubic-lattice topology → JSON snapshot → chemistry build → LAMMPS write).
  - `bfm.py` — 6-neighbour cubic-lattice topology generator forked from `legacy/subprojects/protein_network/topro/topro/bfm/`. Self-avoiding-walk chain placement, MC equilibration (end / kink-crankshaft / reptation moves), Y–Y stochastic crosslinking with Union-Find gel-point detection.
  - `topology_io.py` — JSON I/O round-trip for topology snapshots; format byte-compatible with the legacy `topo_*.json`.
  - `sequence.py` — one-letter block + n_repeats expansion to 3-letter codes; BFM-node-to-residue mapping (anchors at chain ends and Y/TYR positions, NC nodes interpolated).
  - `martini_ff.py` — `MartiniLibrary` parser for `[ atomtypes ]` / `[ nonbond_params ]` / per-molecule blocks; geometric-mixing fallback for missing pairs; explicit unit-conversion helpers (GROMACS kJ/mol, nm → LAMMPS `units real` kcal/mol, Å).
  - `residues.py` — auto-extracted per-residue MARTINI 3 bead pattern + intra-residue bonded terms for the eight amino acids in the resilin reference (GLY, ALA, ARG, PRO, SER, ASP, TYR, ASN). Plus N/C terminal patches and inter-residue backbone parameters (BB-BB bonds with PRO variants, BBB restricted-bending angle, BBBB multi-term proper dihedrals).
  - `data/martini_v3_protein.itp` — master MARTINI 3 ITP pruned to the 14 bead types referenced by the resilin reference (16 MB → 8 KB, 0.05 % retained).
  - `data/martini_v3_water.itp` — copied verbatim from the MARTINI 3 solvents reference.
  - `system.py` — `Bead`, `Bond`, `Constraint`, `Angle`, `Dihedral`, `Exclusion`, `System` dataclasses shared between builder, water packer, and writer.
  - `builder.py` — `build_protein_system(snapshot, sequence_3letter, library)` walks each chain residue-by-residue, places BB beads at BFM-node-interpolated positions, jitters SC beads, applies terminal patches, emits intra- and inter-residue bonded terms, translates BFM `reactions` into SC4–SC4 dityrosine bonds. Uses minimum-image vectors across periodic boundaries.
  - `water.py` — voxel-grid W bead packer with cell-list-accelerated exclusion check around protein beads.
  - `lammps_writer.py` — emits the LAMMPS data + settings + groups files (`<base>.data` with `atom_style full` Bonds/Angles/Dihedrals/Impropers + Masses; `<base>.in.settings` with explicit pair_coeff for every bead-type pair plus bond/angle/dihedral/improper coeffs; `<base>.in.groups` for protein vs water) plus the three relaxation input scripts under `relaxation/` (`stage1.in` soft-push overlap removal; `stage2.in` LJ-epsilon ramp 0.001 → 1.0 via `nve/limit`; `stage3.in` tight CG min + brief NVT/NPT at 310 K). Stage scripts use `units real`, `lj/cut/coul/cut`, `dielectric=15`, no PPPM. Constraints emitted as ULTRA-stiff harmonic bonds matching the reference's `#ifdef FLEXIBLE` block.
  - `workflow.py` — `run_protein_network(...)` end-to-end orchestrator.
  - `cli.py` + `__main__.py` — argparse CLI (`python -m topon.protein_network {generate, sweep, topology}`).
- **`tools/extract_residues_from_itp.py`** — one-shot data-file generator. Reads a polyply-generated protein ITP and emits both `topon/protein_network/residues.py` and the pruned protein-FF ITP. Honours `#ifdef FLEXIBLE` / `#ifndef FLEXIBLE`. Re-running it after dropping a new polyply ITP regenerates the residue table without code edits.
- Topro user-facing reference doc — since consolidated into [`USAGE.md`](USAGE.md) §4.1 (2026-05-08).
- **`tests/output/v33_protein_network_resilin_dry/`** — frozen LAMMPS golden for the regression test.
- **`tests/regression/test_protein_network.py`** — byte-equivalence regression on a 4-chain × 2-repeat resilin system, `seed=42`.
- **`tests/unit/test_protein_network_*.py`** — six unit-test files (residues, FF, BFM/sequence, builder, water, writer) totalling 65 assertions.

#### Approximations (port of a GROMACS force field to LAMMPS)
- Reaction-field electrostatics approximated as `dielectric 15.0` + `pair_style lj/cut/coul/cut` (RF correction term omitted).
- GROMACS `funct=10` (restricted bending) emitted as `angle_style cosine/squared` (drops the `1/sin² θ` factor).
- `[ exclusions ]` collapsed into `special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0` (1-4 inclusive — more aggressive than the reference's `nrexcl=1`).

#### Notes
- Existing DREIDING + KG regression goldens (`tests/output/v21_*`, `tests/output/simbox_crosslink/v4/*`) are unchanged; the protein-network module is purely additive and shares no writer code with the existing pipeline.
- `topon/cli.py`, `topon/pipeline.py`, `topon/forcefield/{dreiding.py,kremer_grest.py}`, and all of `topon/writers/` are byte-identical (verified by Phase-0 SHA-256 baseline at `tests/.baseline_hashes.txt`).
- `topon.singlechain` and `topon.simbox` are not affected; the MARTINI work does **not** add a `martini` option to either.

### [V35] — 2026-05-01 — Peroxide fix

#### Fixed
- **CRITICAL** (`chemistry/builder.py::_create_chain_from_smiles`) — for O-terminal monomer SMILES (PDMS, PTFPMS), the chain SMILES had a stray trailing `[O]` placeholder appended after `(unit + linker) * dp`. Combined with the auto-bridge in `_build_chain_atomistic`, this produced a spurious peroxide `Si–O–O–Si` bond at every chain tail. Symptom in regression outputs: exactly one bond of type `O_3-O_3` per chain in `tests/output/solvent_effects/{PDMS,PTFPMS}/*/dp_*/system.data`, regardless of DP. Resolves §P0-1 of the solubility-downstream issue list.

#### Added
- **`chemistry/builder.py::_PEROXIDE_FIX_APPLIED`** module-level marker, `True`. Downstream projects (e.g. the `solubility` package's `_topon_local/` fork) can feature-detect upstream availability and retire their local patch.
- **`tests/verify_peroxide_fix.py`** — standalone check that builds a PDMS DP=10 chain via `ChemistryBuilder` and asserts no O–O bond in the resulting Mol.

#### Notes
- Reference LAMMPS data files under `tests/output/solvent_effects/PDMS/*` and `tests/output/solvent_effects/PTFPMS/*` will diff after this change: −1 O atom and −1 bond per chain. Regenerate when convenient.

### [V34] — 2026-03 — Solvent-effects study infrastructure

#### Added
- **Multi-solvent mixtures** (`singlechain/workflow.py`, `cli.py`) — `solvent_mixture` parameter accepts a list of `{"smiles": "...", "weight_fraction": float}` dicts. Molecule counts auto-calculated from target box size and weight fractions.
- **Chain centering** (`singlechain/workflow.py`) — polymer chain is centred at the box midpoint after packing for cleaner analysis.
- **Solubility-parameter module** (`singlechain/solubility.py`) — Hoftyzer–Van Krevelen group-contribution HSP estimation using RDKit SMARTS. Supports homopolymers, copolymers (weighted average), end-group correction, and Hansen distance R_a calculation.
- **Improved chain end-group detection** (`singlechain/workflow.py`) — `_guess_head` and `_guess_tail` now handle hydrocarbon monomers (EPDM, PIB), fluoropolymers (FKM), and nitrile-containing monomers (NBR).
- **`--solvent-mixture` CLI** (`cli.py`) — JSON-based multi-solvent specification for `topon chain`.

#### Changed
- `--solvent-smiles` defaults to `None` (internally falls back to toluene). `--n-solvent` auto-calculated when omitted.
- Chain workflow returns `solvent_species` metadata in the result dict.

### [V33] — 2026-03 — Package hardening + single-chain

#### Added
- **BCC and FCC lattices** (`topology/generator_python.py`) — full Python port of C `create_bcc_lattice` (8-neighbour) and `create_fcc_lattice` (12-neighbour), with high-resolution coordinate map and periodic BC. Both match C output exactly.
- **Graft assignment** (`assignment/manager.py`) — `_assign_grafts()` randomly selects backbone bead positions per edge using `graft_density`, writes `graft_positions`, `graft_dp`, `graft_monomer` edge attributes.
- **Copolymer assignment** (`assignment/manager.py`) — `_assign_copolymers()` generates per-bead monomer sequences via `generate_monomer_sequence()`, writes `monomer_sequence` edge attribute.
- **Gradient copolymer arrangement** (`chemistry/sequences.py`) — linear weight interpolation from head to tail; unknown arrangements emit `RuntimeWarning`.
- **CG graft / copolymer consumption** (`chemistry/builder.py`) — `_build_chain_cg()` reads `monomer_sequence` for per-bead `bead_type`; reads graft attributes to attach linear side-chain bead lists tracked in `builder.graft_atom_map`.
- **`topon simbox` CLI** (`cli.py`, `simbox/workflow.py`) — `run_workflow` moved from a test script into the package; CLI exposes `--n-epoxy`, `--n-amino`, `--n-poss`, `--density`, `--seed`.
- **`topon analyze` CLI** (`cli.py`) — wired to `analysis/report.py`; accepts `.gpickle`, `.nodes`, `.edges`; supports `--format text|json`.
- **Single-chain in solvent** (`singlechain/workflow.py`, `cli.py`) — new `topon chain` command; builds a single atomistic DREIDING polymer chain with optional grafts/copolymers and packs it with arbitrary solvent molecules.
- **`CLAUDE.md`** — root-level architecture rules.
- Complete CLI and config-reference documentation — since consolidated into [`USAGE.md`](USAGE.md) (2026-05-08).

#### Fixed
- **CRITICAL** (`chemistry/builder.py`) — SMILES chain fallback silently returned a 1-unit chain on failure; now raises `ValueError` with caller `RuntimeWarning`.
- **HIGH** (`forcefield/dreiding.py`) — removed debug `print`; parameter file now resolved from package `__file__` (no CWD dependency).
- **HIGH** (`forcefield/dreiding.py`) — wildcard entries were tuples, not dicts — caused `KeyError` at write time.
- **MEDIUM** (`assignment/manager.py`) — `max_entanglements` was `0` (falsy / misleading); changed to `None`.

#### Changed
- **Pipeline** (`pipeline.py`) — all six stages wired; `topon generate` CLI functional.
- **Docs** — retired five stale docs; `docs/planning/` removed.
- **Tests** — added POSS 50 % / 100 % regression tests. Full suite: 116 / 116 pass.

### [V32] — 2026-03 — simbox crosslink (v3)

- **simbox v3** — full LAMMPS data + input script generation for crosslinking studies.
  - `generate_all.py` and `setup_crosslink.py` workflow scripts.
  - POSS fraction sweeps with seven compositions (poss_0 through poss_6).

### [V31] — 2026-03 — POSS hydrogens

- **POSS hydrogens** — correct H atom placement and valence in AM0270 cage (`v31_fixed`); debug iterations (`v31_debug`, `v31_debug2`).

### [V30] — 2026-03 — POSS with explicit H

- AM0270 POSS cage conformer now includes explicit hydrogen atoms.

### [V29] — 2026-03 — POSS conformation

- Improved overlap resolution for POSS-containing systems.

### [V28] — 2026-03 — POSS atomistic + simbox v2

#### Added
- **POSS atomistic chemistry** — AM0270-POSS junction support in `ChemistryBuilder`.
  - `generate_atomistic_poss.py` workflow.
  - `run_poss_sweep.py` for parametric POSS-fraction sweeps.
  - `generate_poss_dataset.py` for dataset-scale POSS generation.
  - `_place_poss_am0270()` builder method (Si₈O₁₂ cage + aminopropyl + 7 isooctyl arms).
- **simbox v2** — general-purpose molecule packing sub-system.
  - `Molecule`, `BoxPacker`, `AssembledSystem`, `MoleculeLibrary` classes.
  - `generate_simbox_crosslink.py` workflow for Epoxy-PDMS / Amino-PDMS / POSS systems.
  - 15 iterative sub-versions (v2 through v2.15) refining crosslink setup — see §5 for the granular timeline.

### [V27] — 2026-02 — Custom topology input

- Load arbitrary user-defined graphs directly. Supports any `.nodes`/`.edges` or `.graphml` graph as topology source.

### [V26] — 2026-02 — Multi-entanglement + ABC block copolymer

- **Multi-entanglement** — multiple simultaneous entanglement pairs per system (3× density). Fast-mode variant with repulsive pair style.
- **ABC block copolymer entanglements** — three-component copolymer with entanglements.

### [V25] — 2026-02 — CG refinements + GraphML + ensemble

- **CG baseline refinements** — improved CG generation with repulsive pair style option.
- **GraphML workflow** — full round-trip topology storage. `GraphMLWriter.write_graphml()` exports a dual-graph representation with entanglement edges. `generate_cg_from_graphml.py` and `generate_atomistic_from_graphml.py` load topology from GraphML.
- **Ensemble & study workflows** — `generate_cg_ensemble.py` (1000 systems, 8-worker multiprocessing); `run_ensemble_study.py`, `run_study_v2.py`, `verify_ensemble_metrics.py`. *(The ensemble runner has since been removed; references in older docs do not resolve to current code.)*

### [V24.2] — 2026-01-23 — MARTINI rolled back

- Experimental MARTINI 3 code moved to `older_versions/martini_experiment_2026/`.
- Cleaned up partial implementation of `topon.chemistry.martini`.
- Removed `Auto_Martini` and `Polyply` integration hooks.
- Verified that core atomistic and Kremer-Grest workflows are fully functional.

### [V24.1] — 2026-01-22 — MARTINI experiment (archived)

- Attempted integration of MARTINI 3 CG model. `polyply_lib.py`, `generator.py`, `manager.py` in `chemistry/martini`.
- Verified heteropolymer parameter generation (Benzene / PDMS / FPDMS).
- Shelved due to parameter-mixing complexity (Polyply limitations) and prioritization of existing stable models.

### [V22.0] — 2026-01-20 — Python topology generator port

- Replicated C generator logic in pure Python (`topon.topology.generator_python`).
- Implemented "Strict Sculpting" and "Systematic Search" algorithms.
- Solved the "Hard Case" benchmark in 0.23 s. 60 % success rate on a 20-case batch (`network_candidates_SC_6x6x6_v2.txt`) with a 15 s timeout; ~0.3 s median success time. Validated against original C output.

### [V21.2] — 2026-01-21 — Defect injection

#### Added
- **Defect injection** — primary loops (parallel edges between same node pair). `inject_primary_loops()` in `topon.assignment.defects`. Integrated with `AssignmentManager`. Verification workflows for CG and atomistic.
- **Valence protection** — `max_degree` filter prevents hypervalence (d > 4) in atomistic Si networks.

#### Verified
- Baseline CG: 125 nodes, 210 edges (energy ~214).
- Defected CG: 215 edges, 5 loops (energy ~219).
- Defected atomistic: 11 k atoms, 5 loops; resolved RDKit valence errors by skipping fully-coordinated nodes; stable energy minimization (−38 k kcal/mol); exact atom conservation (+215 net atoms for 5 loops).

### [V21.1] — 2026-01-20 — Entanglement geometry fix

- Fixed zero-length bond issue by ensuring bead spacing N + 2.
- Corrected visualization of kinked chains.

### [V21.0] — 2026-01-20 — Combined features

- Simultaneous entanglements + grafts.
- Validated geometry transformation (linear → kinked → grafted).
- Tested with 5 entanglements + 0.2-density grafts.

### [V20.0] — 2026-01-20 — Dynamic graft scaling

- Graft geometry now scales with the ratio of graft DP to backbone DP: `Effective_Factor = min(Extension_Factor, Graft_DP / Backbone_DP)`. Prevents "spindly" grafts on short backbones.

### [V19.0] — 2026-01-20 — Scaled grafts

- Introduced `graft_extension_factor` (default 0.5).
- Graft length scales with edge length instead of absolute bond units.

### [V18.0] — 2026-01-20 — Grafts (side chains)

- Support for grafted side chains (e.g. siloxane on PDMS).
- CG: probabilistic attachment (type 'G').
- Atomistic: methyl substitution.

### [V17.0] — 2026-01-20 — Copolymers

- Support for block, random, alternating, and gradient sequences.
- `topon.chemistry.sequences` module.
- Verified A–B block copolymers in both CG and atomistic.

### [V16.0] — 2026-01-15 — Atomistic entanglements

- Ported Gaussian Kink geometry to DREIDING model.
- Verified stable overlap resolution for entangled atomistic chains.

### [V15.0] — 2026-01-15 — Entanglements (CG)

- Gaussian Kink implementation for topological knots.
- Selects spatially adjacent but topologically distant edges (~10 candidates).
- Verified kink geometry generation.

### [V14.0] — 2026-01-15 — Externalised parameters

- `run_steps` and dynamics settings moved to config / JSON.
- Created `experimental.json` for production settings.

### [V13.0] — 2026-01-15 — Pair-style configuration

- Added `pair_style` option: `repulsive` (WCA, 1.12 σ cutoff) vs `attractive` (LJ, 2.5 σ cutoff).
- Repulsive mode ~4× faster for equilibration.

### [V12] — 2026-01-15 — Automated execution

- `SimulationRunner` module for automated LAMMPS execution.
- `--run` flag for `generate_cg.py`.
- `execution` config section with `auto_run`, `executable`, `n_procs`.
- `GraphAttributor` module for graph re-attribution workflow.
- Default assignment files: `node_degree.json`, `edge_uniform.json`.
- `assign_graph.py` example workflow.
- `verify_workflows.py` comprehensive verification script.
- Separated `config_cg.json` / `config_atomistic.json`.
- Overlap resolution now configurable: CG (0.01, 10 iters), atomistic (0.1 Å, 20 iters).
- Updated node-degree mapping: `1 = end`, `2–8 = A`.

### [V11] — 2026-01-15 — Global angle toggle

- `include_angles` config flag — if `false`, no angles defined anywhere (data file, input scripts).
- `delete_bonds` only generated if `include_angles=true` AND `remove_cg_angles=true`.

### [V10] — 2026-01-15 — Configurable angle removal

- `remove_cg_angles` option in config.

### [V9] — 2026-01-15 — `delete_bonds` syntax

- Fixed `delete_bonds` command syntax: `delete_bonds all angle 1-1 remove`.

### [V8] — 2026-01-15

- Attempted `delete_bonds all angle * remove` (failed; superseded by V9).

### [V7] — 2026-01-15

- Added `remove` keyword to `delete_bonds` command.

### [V6] — 2026-01-15

- Explicit `bond_coeff` and `angle_coeff` in `minimize_3_parallel.in` Phase A.

### [V5] — 2026-01-15 — Gradual FENE switch

- Phase A: harmonic pre-minimization.
- Phase B: `delete_bonds` to remove angles, switch to FENE.
- Phase C: dynamics.

### [V4] — 2026-01-15

- Removed annealing/equilibration script generation for CG (min-only workflow).

### [V3] — 2026-01-15

- CG minimization scripts renamed to `_parallel.in` convention.
- Added explicit `theta0=180` to angle coefficients.
- Assigned 'J' bead type to junctions, 'A' to chain atoms.

### [V2] — 2026-01-14

- CG angle coefficients now include `theta0`.
- Atomistic DP reduced to 5 (~10 k atoms).

### [V1] — 2026-01-14 — Initial topon package

- Initial topon package structure.
- Core modules: topology, assignment, chemistry, conformation, writers.
- CLI with `generate`, `validate`, `init` commands.
- Pydantic config schema.
- CG and atomistic workflow examples.

---

## 5. simbox sub-version history

The simbox sub-system passed through 15 numbered iterations between V28 (introduction) and V32 (full LAMMPS generation), then was stabilised at v4 in V33. Most sub-versions are small refinements to packing, atom typing, or crosslink-template generation; this is the granular log for archaeological purposes.

| Version | Key change |
|---|---|
| v2.0 | Initial simbox packing with Epoxy-PDMS + Amino-PDMS |
| v2.1 – v2.9 | Refinements to packing algorithm, overlap detection, atom typing |
| v2.10 | `setup_crosslink.py` — manual crosslink-setup helper |
| v2.11 | `gen_pre_template.py`, `check_template_match.py` — bond/react templates |
| v2.11 (poss_only) | `patch_data.py` — POSS-only system test |
| v2.12 – v2.14 | Improved settings and atom-type consistency |
| v2.15 | `verify_consistency.py` — data-file validation |
| v3 | `generate_all.py` — full POSS-fraction sweep (poss_0 – poss_6) |
| v4 | `generate_simbox_crosslink.py` (canonical reference); fixed writer header bug (`max` not `len` for type counts); added `ff_coeffs.in` generation; `UniversalTypeMapper` enforces stable type IDs (N_3=4, H_=5) so `fix bond/react` templates remain compatible across compositions |

The current canonical entry point is `topon.simbox.workflow.run_workflow` (V33), exposed as the `topon simbox` CLI.

---

## 6. Pointers

- *What's planned next; pitfalls collected during development; private notes:* [`internal/DEVELOPMENT_INTERNAL.md`](../internal/DEVELOPMENT_INTERNAL.md) (gitignored, local only).
- *Package layout, module responsibilities, design principles:* [ARCHITECTURE.md](ARCHITECTURE.md).
- *How to run topon, CLI flags, config schema:* [USAGE.md](USAGE.md).
- *Rules every change must follow:* [`CLAUDE.md`](../CLAUDE.md) at the repo root.
