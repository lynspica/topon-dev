# Topon Changelog

All notable changes to this project are documented here.

---

## [V36] - 2026-05-07 (MARTINI 3 protein-network generator)

### Added
- **`topon/protein_network/`** - new peer package alongside `simbox/` and `singlechain/`. Sequence-string-to-LAMMPS MARTINI 3 generator for protein networks with stochastic dityrosine crosslinks, mirroring the topro CHARMM-atomistic shape (BFM cubic-lattice topology -> JSON snapshot -> chemistry build -> LAMMPS write).
  - `bfm.py` - 6-neighbour cubic-lattice topology generator forked from `legacy/subprojects/protein_network/topro/topro/bfm/`. Self-avoiding-walk chain placement, MC equilibration (end / kink-crankshaft / reptation moves), Y-Y stochastic crosslinking with Union-Find gel-point detection. Snapshot JSON format byte-compatible with topro's `topo_*.json`.
  - `topology_io.py` - JSON I/O round-trip for topology snapshots.
  - `sequence.py` - one-letter block + n_repeats expansion to 3-letter codes; BFM-node-to-residue mapping (anchors at chain ends and Y/TYR positions, NC nodes interpolated).
  - `martini_ff.py` - `MartiniLibrary` parser for `[ atomtypes ]` / `[ nonbond_params ]` / per-molecule blocks; geometric-mixing fallback for missing pairs; explicit unit-conversion helpers from GROMACS (kJ/mol, nm) to LAMMPS `units real` (kcal/mol, A).
  - `residues.py` - **auto-extracted** per-residue MARTINI 3 bead pattern + intra-residue bonded terms (bonds, constraints, angles, dihedrals, impropers, exclusions) for the 8 amino acids in the resilin reference: GLY, ALA, ARG, PRO, SER, ASP, TYR, ASN. Plus N/C terminal patches (Q5+1, Q5-1) and inter-residue backbone parameters (BB-BB bonds with PRO variants, BBB restricted-bending angle, BBBB multi-term proper dihedrals tagged BBBB / GGGX / GGXG / GXGG / XGGG).
  - `data/martini_v3_protein.itp` - master MARTINI 3 ITP pruned to the 14 bead types referenced by the resilin reference (16 MB -> 8 KB; 0.05% retained).
  - `data/martini_v3_water.itp` - copied verbatim from MARTINI 3 solvents reference.
  - `system.py` - `Bead`, `Bond`, `Constraint`, `Angle`, `Dihedral`, `Exclusion`, `System` dataclasses shared between builder, water packer, and writer.
  - `builder.py` - `build_protein_system(snapshot, sequence_3letter, library)` walks each chain residue-by-residue, places BB beads at BFM-node-interpolated positions, jitters SC beads, applies terminal patches, emits intra- and inter-residue bonded terms, and translates BFM `reactions` into SC4-SC4 dityrosine bonds. Uses minimum-image vectors across periodic boundaries.
  - `water.py` - voxel-grid W bead packer with cell-list-accelerated exclusion check around protein beads (`pack_water(system, library, density_w_per_nm3, exclusion_radius_ang, seed)`).
  - `lammps_writer.py` - emits four files: `<base>.data` (atom_style full, Bonds/Angles/Dihedrals/Impropers + Masses), `<base>.in.settings` (explicit pair_coeff for every bead-type pair, plus bond/angle/dihedral/improper coeffs), `<base>.in.groups` (protein vs. water), `<base>.in` (units real, lj/cut/coul/cut, dielectric=15, minimize -> NVT -> NPT recipe). Constraints emitted as ULTRA-stiff harmonic bonds matching the reference's `#ifdef FLEXIBLE` block.
  - `workflow.py` - `run_protein_network(...)` end-to-end orchestrator.
  - `cli.py` + `__main__.py` - argparse CLI (`python -m topon.protein_network {generate, sweep, topology}`).
- **`tools/extract_residues_from_itp.py`** - one-shot data-file generator. Reads a polyply-generated protein ITP (default: `tests/_martini_extracted/Martini_Ahmet/itp_files/nat_pro.itp`) and emits both `topon/protein_network/residues.py` and the pruned protein FF ITP. Honours `#ifdef FLEXIBLE` / `#ifndef FLEXIBLE`. Re-running it after dropping a new polyply ITP into the source directory regenerates the residue table without code edits.
- **`docs/protein_network.md`** - user-facing reference (CLI, Python API, output layout, known approximations, citations).
- **`tests/output/v33_protein_network_resilin_dry/`** - frozen LAMMPS golden for the regression test.
- **`tests/regression/test_protein_network.py`** - byte-equivalence regression on a 4-chain x 2-repeat resilin system, seed=42.
- **`tests/unit/test_protein_network_*.py`** - 6 new unit-test files (residues, FF, BFM/sequence, builder, water, writer) totalling 65 assertions.

### Approximations
- Reaction-field electrostatics approximated as `dielectric 15.0` + `pair_style lj/cut/coul/cut` (RF correction term omitted; documented in `protein_network.md`).
- GROMACS angle funct=10 (restricted bending) emitted as `angle_style cosine/squared` (omits the `1/sin^2 t` factor).
- `[ exclusions ]` collapsed into `special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0` (1-4 inclusive, more aggressive than the reference's `nrexcl=1`).

### Notes
- Existing DREIDING + KG regression goldens (`tests/output/v21_*`, `tests/output/simbox_crosslink/v4/*`) are unchanged; the protein-network module is purely additive and shares no writer code with the existing pipeline.
- `topon/cli.py`, `topon/pipeline.py`, `topon/forcefield/{dreiding.py,kremer_grest.py}`, and all of `topon/writers/` remain byte-identical (verified by Phase-0 SHA-256 baseline at `tests/.baseline_hashes.txt`).
- `topon.singlechain` and `topon.simbox` are not affected. The MARTINI work does NOT add a `martini` option to either of them.
- This bumps the conceptual project version from V35 to V36; the changelog version numbers are independent of any formal release.

---

## [V35] - 2026-05-01 (peroxide fix)

### Fixed
- **CRITICAL** (`chemistry/builder.py::_create_chain_from_smiles`): For O-terminal monomer SMILES (PDMS, PTFPMS), the chain SMILES had a stray trailing `[O]` placeholder appended after `(unit + linker) * dp`. Combined with the auto-bridge in `_build_chain_atomistic`, this produced a spurious peroxide `Si-O-O-Si` bond at every chain tail. Symptom in regression outputs: exactly one bond of type `O_3-O_3` per chain in `tests/output/solvent_effects/{PDMS,PTFPMS}/*/dp_*/system.data`, regardless of DP. Resolves §P0-1 of `docs/topon_issues_from_solubility_v2.md`.

### Added
- **`chemistry/builder.py::_PEROXIDE_FIX_APPLIED`** module-level marker, `True`. Downstream projects (e.g. the `solubility` package's `_topon_local/` fork) can feature-detect upstream availability of the fix and retire their local patch.
- **`tests/verify_peroxide_fix.py`** standalone check that builds a PDMS DP=10 chain via `ChemistryBuilder` and asserts no O-O bond in the resulting Mol.

### Notes
- Reference LAMMPS data files under `tests/output/solvent_effects/PDMS/*` and `tests/output/solvent_effects/PTFPMS/*` will diff after this change: −1 O atom and −1 bond per chain. Regenerate when convenient.

---

## [V34] - 2026-03 (solvent-effects study infrastructure)

### Added
- **Multi-solvent mixtures** (`singlechain/workflow.py`, `cli.py`): `solvent_mixture` parameter accepts a list of `{"smiles": "...", "weight_fraction": float}` dicts. Molecule counts auto-calculated from target box size and weight fractions.
- **Chain centering** (`singlechain/workflow.py`): Polymer chain is now centered at the box midpoint after packing for cleaner analysis.
- **Solubility parameter module** (`singlechain/solubility.py`): Hoftyzer-Van Krevelen group-contribution HSP estimation using RDKit SMARTS. Supports homopolymers, copolymers (weighted average), end-group correction, and Hansen distance R_a calculation.
- **Improved chain end-group detection** (`singlechain/workflow.py`): `_guess_head` and `_guess_tail` now handle hydrocarbon monomers (EPDM, PIB), fluoropolymers (FKM), and nitrile-containing monomers (NBR).
- **`--solvent-mixture` CLI** (`cli.py`): JSON-based multi-solvent specification for `topon chain`.

### Changed
- `--solvent-smiles` defaults to `None` (internally defaults to toluene). `--n-solvent` auto-calculated when omitted.
- Chain workflow now returns `solvent_species` metadata in result dict.

---


## [V33] - 2026-03 (package hardening + single-chain)

### Added
- **BCC and FCC lattices** (`topology/generator_python.py`): Full Python port of C `create_bcc_lattice` (8-neighbor) and `create_fcc_lattice` (12-neighbor), using high-resolution coordinate map with periodic BC. Both match C output exactly.
- **Graft assignment** (`assignment/manager.py`): `_assign_grafts()` implemented — randomly selects backbone bead positions per edge using `graft_density`, writes `graft_positions`, `graft_dp`, `graft_monomer` edge attributes.
- **Copolymer assignment** (`assignment/manager.py`): `_assign_copolymers()` implemented — generates per-bead monomer sequences via `generate_monomer_sequence()`, writes `monomer_sequence` edge attribute.
- **Gradient copolymer arrangement** (`chemistry/sequences.py`): Linear weight interpolation from head to tail; unknown arrangements now emit `RuntimeWarning`.
- **ChemistryBuilder CG graft/copolymer consumption** (`chemistry/builder.py`): `_build_chain_cg()` reads `monomer_sequence` for per-bead `bead_type`; reads graft attributes to attach linear side-chain bead lists tracked in `builder.graft_atom_map`.
- **`topon simbox` CLI** (`cli.py`, `simbox/workflow.py`): `run_workflow` moved from test script into package; CLI exposes `--n-epoxy`, `--n-amino`, `--n-poss`, `--density`, `--seed`.
- **`topon analyze` CLI** (`cli.py`): Wired to `analysis/report.py`; accepts `.gpickle`, `.nodes`, `.edges`; supports `--format text|json`.
- **Single-chain in solvent** (`singlechain/workflow.py`, `cli.py`): New `topon chain` command; builds a single atomistic DREIDING polymer chain with optional grafts/copolymers and packs it with arbitrary solvent molecules.
- **`CLAUDE.md`**: Root-level architecture rules.
- **`docs/cli.md`**, **`docs/config_reference.md`**: Complete reference documentation.

### Fixed
- **CRITICAL** (`chemistry/builder.py`): SMILES chain fallback silently returned 1-unit chain on failure; now raises `ValueError` with caller `RuntimeWarning`.
- **HIGH** (`forcefield/dreiding.py`): Removed debug `print`; parameter file now resolved from package `__file__` (no CWD dependency).
- **HIGH** (`forcefield/dreiding.py`): Wildcard entries were tuples, not dicts — caused KeyError at write time.
- **MEDIUM** (`assignment/manager.py`): `max_entanglements` was `0` (falsy/misleading); changed to `None`.

### Changed
- **Pipeline** (`pipeline.py`): All 6 stages wired — `topon generate` CLI functional.
- **Docs**: Retired 5 stale docs; `docs/planning/` removed.
- **Tests**: Added POSS 50%/100% regression tests. Full suite: 116/116 pass.

---

## [V32] - 2026-03 (simbox crosslink)

### Added
- **simbox v3**: Full LAMMPS data + input script generation for crosslinking studies.
    - `generate_all.py` and `setup_crosslink.py` workflow scripts.
    - POSS fraction sweeps with 7 compositions (poss_0 through poss_6).

---

## [V31] - 2026-03

### Fixed
- **POSS hydrogens**: Correct H atom placement and valence in AM0270 cage (`v31_fixed`).
- Debug iterations for POSS H conformer generation (`v31_debug`, `v31_debug2`).

---

## [V30] - 2026-03

### Added
- **POSS with explicit H**: AM0270 POSS cage conformer now includes explicit hydrogen atoms.

---

## [V29] - 2026-03

### Fixed
- **POSS conformation**: Improved overlap resolution for POSS-containing systems.

---

## [V28] - 2026-03

### Added
- **POSS atomistic chemistry**: AM0270-POSS junction support in `ChemistryBuilder`.
    - `generate_atomistic_poss.py` workflow.
    - `run_poss_sweep.py` for parametric POSS fraction sweeps.
    - `generate_poss_dataset.py` for dataset-scale POSS generation.
    - `_place_poss_am0270()` builder method (Si8O12 cage + aminopropyl + 7 isooctyl arms).
- **simbox v2**: General-purpose molecule packing sub-system.
    - `Molecule`, `BoxPacker`, `AssembledSystem`, `MoleculeLibrary` classes.
    - `generate_simbox_crosslink.py` workflow for Epoxy-PDMS / Amino-PDMS / POSS systems.
    - 15 iterative sub-versions (v2 through v2.15) refining crosslink setup.

---

## [V27] - 2026-02

### Added
- **Custom Topology Input**: Load arbitrary user-defined graphs directly.
    - `generate_v27_custom_topology.py` workflow.
    - Supports any `.nodes`/`.edges` or `.graphml` graph as topology source.

---

## [V26] - 2026-02

### Added
- **Multi-Entanglement**: Multiple simultaneous entanglement pairs per system.
    - `generate_v26_multi_entanglement.py` — 3× entanglement density.
    - Fast mode variant (`v26_multi_entanglement_fast`) with repulsive pair style.
- **ABC Block Copolymer Entanglements**: Three-component copolymer with entanglements.
    - `generate_v26_abc_entangled.py` workflow.

---

## [V25] - 2026-02

### Added
- **CG Baseline Refinements**: Improved CG generation with repulsive pair style option.
    - `generate_v25_cg.py` workflow.
- **GraphML Workflow**: Full round-trip topology storage via GraphML.
    - `generate_graphml_test.py`, `verify_graphml_format.py` workflows.
    - `generate_cg_from_graphml.py`, `generate_atomistic_from_graphml.py` — load topology from GraphML.
    - `GraphMLWriter.write_graphml()` — exports dual-graph representation with entanglement edges.
- **Ensemble & Study Workflows**:
    - `run_ensemble_study.py` — parametric study runner.
    - `run_study_v2.py` and `visualize_study_v2.py` — enhanced study with metrics.
    - `verify_ensemble_metrics.py` — validates topology statistics across ensemble.
    - `generate_cg_ensemble.py` — multiprocessing ensemble runner (1000 systems, 8 workers).

---

## [V24.2] - 2026-01-23

### Changed
- **Martini Implementation Rolled Back**:
    - Experimental Martini 3 code moved to `older_versions/martini_experiment_2026/`.
    - Cleaned up partial implementation of `topon.chemistry.martini`.
    - Removed `Auto_Martini` and `Polyply` integration hooks.
    - Verified that core Atomistic and Kremer-Grest workflows are fully functional and stable.

---

## [V24.1] - 2026-01-22

### Experimental (Archived)
- Attempted integration of Martini 3 Coarse-Grained model.
- Created `polyply_lib.py`, `generator.py`, `manager.py` in `chemistry/martini`.
- Verified heteropolymer parameter generation (Benzene/PDMS/FPDMS).
- Shelved due to complexity in parameter mixing (Polyply limitations) and prioritization of existing stable models.

---

## [V22.0] - 2026-01-20

### Added
- **Python Topology Generator Port**:
    - Replicated C generator logic in pure Python (`topon.topology.generator_python`).
    - Implemented "Strict Sculpting" and "Systematic Search" algorithms.
    - Achieved 60% success rate on "Hard Case" benchmark (0.3s median time).
    - Validated against original C output.

---

## [V21.2] - 2026-01-21

### Added
- **Defect Injection**: Support for "primary loops" (parallel edges between same node pair).
    - `inject_primary_loops()` algorithm in `topon.assignment.defects`.
    - Integrated with `AssignmentManager`.
    - Verification workflows for CG and Atomistic.
- **Valence Protection**: `max_degree` filter prevents hypervalence (d>4) in atomistic Si networks.

### Verified
- **Baseline CG**: 125 nodes, 210 edges (Energy ~214).
- **Defected CG**: 215 edges, 5 loops (Energy ~219).
- **Defected Atomistic**: 11k atoms, 5 loops. 
    - Resolved RDKit valence errors by skipping fully-coordinated nodes.
    - Confirmed stable energy minimization (-38k kcal/mol).
    - Verified exact atom conservation (+215 net atoms for 5 loops).

---

## [V21.1] - 2026-01-20

### Fixed
- **Entanglement Geometry**:
    - Fixed zero-length bond issue by ensuring bead spacing N+2.
    - Corrected visualization of kinked chains.

---

## [V21.0] - 2026-01-20

### Added
- **Combined Features**: 
    - Simultaneous Entanglements + Grafts.
    - Validated geometry transformation (Linear -> Kinked -> Grafted).
    - Tested with 5 entanglements + 0.2 density grafts.

---

## [V20.0] - 2026-01-20

### Changed
- **Dynamic Graft Scaling**:
    - Graft geometry now scales ratio of Graft DP to Backbone DP.
    - `Effective_Factor = min(Extension_Factor, Graft_DP / Backbone_DP)`.
    - Prevents "spindly" grafts on short backbones.

---

## [V19.0] - 2026-01-20

### Changed
- **Scaled Grafts**:
    - Introduced `graft_extension_factor` (default 0.5).
    - Graft length scales with edge length instead of absolute bond units.

---

## [V18.0] - 2026-01-20

### Added
- **Grafts (Side Chains)**:
    - Support for grafted side chains (e.g. Siloxane on PDMS).
    - CG: Probabilistic attachment (Type 'G').
    - Atomistic: Methyl substitution.

---

## [V17.0] - 2026-01-20

### Added
- **Copolymers**:
    - Support for Block, Random, Alternating, and Gradient sequences.
    - `topon.chemistry.sequences` module.
    - Verified A-B block copolymers in both CG and Atomistic.

---

## [V16.0] - 2026-01-15

### Added
- **Atomistic Entanglements**:
    - Ported Gaussian Kink geometry to DREIDING model.
    - Verified stable overlap resolution for entangled atomistic chains.

---

## [V15.0] - 2026-01-15

### Added
- **Entanglements (CG)**:
    - Gaussian Kink implementation implementation for topological knots.
    - Selects spatially adjacent but topologically distant edges (~10 candidates).
    - Verified kink geometry generation.

---

## [V14.0] - 2026-01-15

### Changed
- **Externalized Parameters**:
     - `run_steps` and dynamics settings moved to config/JSON.
     - Created `experimental.json` for production settings.

---

## [V13.0] - 2026-01-15

### Added
- **Pair Style Configuration**:
    - Added `pair_style` option: `repulsive` (WCA, 1.12 cutoff) vs `attractive` (LJ, 2.5 cutoff).
    - Repulsive mode ~4x faster for equilibration.

---

## [V12] - 2026-01-15

### Added
- `SimulationRunner` module for automated LAMMPS execution
- `--run` flag for `generate_cg.py`
- `execution` config section with `auto_run`, `executable`, `n_procs`
- `GraphAttributor` module for graph re-attribution workflow
- Default assignment files: `node_degree.json`, `edge_uniform.json`
- `assign_graph.py` example workflow
- `verify_workflows.py` comprehensive verification script

### Changed
- Separated `config_cg.json` and `config_atomistic.json`
- Overlap resolution now configurable: CG (0.01, 10 iters), Atomistic (0.1 Å, 20 iters)
- Updated node degree mapping: 1=end, 2-8=A

---

## [V11] - 2026-01-15

### Added
- Global angle toggle (`include_angles`) in config
- If `false`, no angles defined anywhere (data file, input scripts)

### Changed
- `delete_bonds` only generated if `include_angles=true` AND `remove_cg_angles=true`

---

## [V10] - 2026-01-15

### Added
- Configurable angle removal via `remove_cg_angles` in config

---

## [V9] - 2026-01-15

### Fixed
- `delete_bonds` command syntax: `delete_bonds all angle 1-1 remove`

---

## [V8] - 2026-01-15

### Fixed
- Attempted `delete_bonds all angle * remove` (failed)

---

## [V7] - 2026-01-15

### Fixed
- Added `remove` keyword to `delete_bonds` command

---

## [V6] - 2026-01-15

### Fixed
- Explicit `bond_coeff` and `angle_coeff` in `minimize_3_parallel.in` Phase A

---

## [V5] - 2026-01-15

### Added
- Gradual FENE switch workflow:
  - Phase A: Harmonic pre-minimization
  - Phase B: `delete_bonds` to remove angles, switch to FENE
  - Phase C: Dynamics

---

## [V4] - 2026-01-15

### Removed
- Annealing/equilibration script generation for CG (min-only workflow)

---

## [V3] - 2026-01-15

### Changed
- CG minimization scripts renamed to `_parallel.in` convention
- Added explicit `theta0=180` to angle coefficients
- Assigned 'J' bead type to junctions, 'A' to chain atoms

---

## [V2] - 2026-01-14

### Fixed
- CG angle coefficients now include theta0
- Atomistic DP reduced to 5 (~10k atoms)

---

## [V1] - 2026-01-14

### Added
- Initial topon package structure
- Core modules: topology, assignment, chemistry, conformation, writers
- CLI with generate, validate, init commands
- Pydantic config schema
- CG and Atomistic workflow examples
