# Topon Package Development

## Current Objective
Create a clean, modular Python package `topon` for polymer network generation that:
- Preserves all existing functionality from older versions
- Is easy-to-use without requiring source code changes
- Supports both local and HPC/supercomputer workflows
- Uses JSON configuration with web-based GUI (Streamlit + notebook)
- Is modular for future extensions (polyply, ML integration)

---

## User Decisions ✓
- **C Generator**: Preserve, invoke via subprocess (no modifications)
- **GUI**: Streamlit web app + equivalent Jupyter notebook
- **Testing**: Both unit and integration tests
- **Force Fields**: Bundle DREIDING + KG parameters with package
- **Config**: Main JSON + experimental settings JSON

---

## Phase 1: Planning & Design
- [x] Analyze `SystemGenerator_portable` structure
- [x] Analyze `polymer_world` modular architecture
- [x] Review C-code topology generator (`generator_serial_debug11.c`)
- [x] Identify force field implementations (DREIDING, Kremer-Grest)
- [x] Review LAMMPS writers and minimization pipelines
- [x] Write implementation plan v2
- [x] User approval of implementation plan

## Phase 2: Core Package Structure
- [x] Create package directory structure
- [x] Set up `pyproject.toml` with dependencies
- [x] Create main config schema (Pydantic)
- [x] Create experimental settings schema
- [x] Implement config validation

## Phase 3: Topology & Assignment Module
- [x] Create C generator subprocess wrapper
- [x] Create SLURM script generator for HPC
- [x] Implement graph loader (.nodes/.edges/.gpickle)
- [x] Node type assignment (degree/positional/random/explicit)
- [x] Edge type assignment (uniform/random/composite)
- [x] DP distribution with PDI (Schulz-Zimm)
- [x] Entanglement selection (nearest disjoint neighbor)
- [x] Analysis module (max defects/entanglements)
- [x] Defect injection (primary loops enabled) - **VERIFIED v21.2**

## Phase 4: Chemistry Module
- [x] Port `ChemistryBuilder` from polymer_world
- [x] Implement SMILES library with validation
- [x] Support atomistic (DREIDING) mode
- [x] Support CG (Kremer-Grest) mode with options
- [x] Support copolymers, grafts, edge-type mapping
- [ ] Add polyply integration placeholder

## Phase 5: Conformation Module
- [x] Port displacement file generation
- [x] Entanglement kink geometry
- [x] Implement coordinate wrapping

## Phase 6: Writers Module
- [x] Port LAMMPS data writers (DreidingWriter, CGWriter from polymer_world)
- [x] Port network_helpers and network_config to utils/
- [x] Port ConformationManager from polymer_world
- [x] Implement staged minimization pipeline (LammpsInputGenerator)
- [ ] Externalize LAMMPS parameters to experimental.json

## Phase 7: GUI & CLI
- [x] Create CLI entry points
- [ ] Create Streamlit web app
- [ ] Create equivalent Jupyter notebook
- [ ] Implement SMILES validation UI
- [ ] Implement copolymer/edge-type mutual exclusion

## Phase 8: Testing & Documentation
- [x] Copy sample graphs for test data
- [ ] Create unit tests
- [x] Create integration tests (Verified 5x5x5 workflow)
- [x] Refine output formats (3 minimization scripts, valid settings file)
- [x] Fix CG Angle Coeffs (Added Theta0)
- [x] Resize Atomistic Test (DP=5 ~10k atoms)
- [x] CG: Assign 'J' types to junctions, 'A' to chains
- [x] CG: Rename output files to `_parallel.in` convention
- [x] CG: Remove Annealing/Equilibration script generation (Min-only)
- [x] CG: Start with Harmonic Bond Style (Stage 1 & 2)
- [x] CG: Stage 3 Gradual Switch: Harmonic Min -> Remove Angles -> FENE Min
- [x] Verified LAMMPS execution (Min 1, 2, 3 sequence) -> Success
- [x] Made Angle Removal in Stage 3 configurable via JSON (`remove_cg_angles`)
- [x] Refactored Config: Distinct `config_cg.json` and `config_atomistic.json`
- [x] Configurable Overlap Resolution (cutoff, iters)
- [x] Implemented Global Angle Toggle (`include_angles`)
- [x] Created `SimulationRunner` for automated execution
- [x] Integrated execution into `generate_cg.py` (`--run` flag)
- [x] Tested Full Execution Pipeline (V12)
- [x] Updated DEVLOG, README with config documentation
- [x] Created Graph Attributor module (`topon/assignment/attributor.py`)
- [x] Created default assignment files (`node_degree.json`, `edge_uniform.json`)
- [x] Created `assign_graph.py` example workflow
- [x] End-to-End Verification: CG + Atomistic workflows PASSED
- [x] Created `docs/development/` folder structure
- [x] Created `legacy_comparison.md` with parameter tables
- [x] Created `changelog.md` with version history (V1-V12)
- [x] Created `tasks.md` with accomplished/planned items
- [x] Updated atomistic workflow to use `config_atomistic.json`
- [x] Added configurable pair cutoff (repulsive/attractive) to CG workflow
- [x] V13 Test: CG + Atomistic PASSED with new pair_style config
- [x] Write README with quickstart
- [ ] Create usage examples

## Phase 6: Advanced Features Verification (v15-v21)
- [x] **Entanglements (v15-v16)**:
    - [x] Gaussian Kink implementation (Corrected Orientation & Phase)
    - [x] Verification: CG (v15.2) & Atomistic (v16.0)
    - [x] Fix: Bead spacing logic (N+2 points) (v21.1)
- [x] **Copolymers (v17)**:
    - [x] Random, Block, Alternating Sequence Generation
    - [x] Verification: CG & Atomistic Workflows
- [x] **Grafts (v18-v20)**:
    - [x] Replaced Methyl with Siloxane Side Chains
    - [x] Scaled Geometry (v19) and Dynamic DP Scaling (v20)
    - [x] Verification: CG & Atomistic Workflows
- [x] **Combined Features (v21)**:
    - [x] Simultaneous Entanglements + Grafts
    - [x] Reduced Graft Density (0.05) to avoid crowding
    - [x] Verification: CG & Atomistic Workflows (Fixed Geometry)
- [x] **Defect Injection (v21.2)**:
    - [x] Primary Loops Logic (`defects.py`)
    - [x] Valence Protection (Max Degree Filter)
    - [x] Verification: CG (Stable) & Atomistic (Stable, -38k kcal/mol)

## Phase 7: Topology Generation (C Integration)
- [x] **Setup**:
    - [x] Locate/Compile C Executable (`generator.exe`)
    - [x] Verify functionality (Argument signature match)
- [x] **Verification**:
    - [x] Create `test_topology_generation.py` workflow
    - [x] Verify graph generation (Nodes/Edges output)
    - [x] Verify generated graph validity (Load in Python)

## Phase 8: Python Benchmark & Port (Complete)
- [x] **Implementation**:
    - [x] Port `generator_serial_debug11.c` logic to `topon/topology/generator_python.py`
    - [x] Implement `Strict Sculpting` stages (1, 2, 3)
    - [x] Implement `Systematic Search` stage (4)
- [x] **Benchmarking**:
    - [x] Create `benchmark_topology.py`
    - [x] Compare C vs Python execution time for identical inputs
    - [x] **Result**: Python SOLVED "Hard Case" in 0.23s. Python is fully viable.
    - [x] **Batch Verification**: Tested 20 cases from `network_candidates_SC_6x6x6_v2.txt`.
        - Success Rate: 60% (12/20) with 15s timeout.
        - Median Success Time: ~0.3s.

## Phase 9: Martini 3 Coarse-Grained Integration (Archived)
- [x] **Experiment (V24.1)**:
    - [x] Implement `MartiniManager` and basic generators.
    - [x] Verify heteropolymer parameter generation.
- [x] **Decision (V24.2)**:
    - [x] Work archived to `older_versions/martini_experiment_2026/`.
    - [x] Rolled back to focus on Atomistic/KG stability.
    - [x] Reason: Complexity of automated parameter mixing (Polyply limitations).

## Phase 10: Testing & Documentation (Pending)
- [x] Copy sample graphs for test data
- [x] Create unit tests (Config, Topology, Assignment, Chemistry)
- [x] Create integration tests (Verified 5x5x5 workflow)
- [x] Refine output formats (3 minimization scripts, valid settings file)
- [ ] Create usage examples/tutorials
- [ ] Write API documentation

## Phase 11: GUI & Usability (Pending)
- [x] Create CLI entry points
- [x] Wire `Pipeline` stages to real modules (`topon/pipeline.py` — 2026-03-23)
- [ ] Create Streamlit web app (`topon/gui/app.py` missing)
- [ ] Create equivalent Jupyter notebook
- [ ] Implement SMILES validation UI

## Phase 12: Advanced Features (Completed)
- [x] **GraphML Workflow (v25)**:
    - [x] `GraphMLWriter` — dual-graph export with entanglement edges
    - [x] `generate_cg_from_graphml.py`, `generate_atomistic_from_graphml.py` workflows
- [x] **Multi-Entanglement (v26)**:
    - [x] Multiple simultaneous entanglement pairs per system
    - [x] ABC block copolymer + entanglement combined
- [x] **Custom Topology (v27)**:
    - [x] Load arbitrary user graphs as topology source
- [x] **POSS Junctions (v28–v32)**:
    - [x] `_place_poss_am0270()` in `ChemistryBuilder` (Si8O12 cage + arms)
    - [x] `generate_atomistic_poss.py` workflow
    - [x] `run_poss_sweep.py` parametric POSS fraction sweeps
    - [x] `generate_poss_dataset.py` dataset-scale generation
    - [x] Explicit H atom placement verified (v30–v31)
- [x] **Ensemble & Study Workflows (v25)**:
    - [x] `run_cg_ensemble.py` — multiprocessing batch runner (1000 systems, 8 workers)
    - [x] `run_ensemble_study.py`, `run_study_v2.py` — parametric study runners
    - [x] `verify_ensemble_metrics.py` — ensemble statistics validation
    - [x] `docs/cg_ensemble_execution.md` written

## Phase 13: simbox Sub-System (v28+)
- [x] `simbox/molecule.py` — `Molecule` class (SMILES/PDB/RDKit, reactive-site detection)
- [x] `simbox/packer.py` — `BoxPacker` (density-based box sizing, spatial hash overlap)
- [x] `simbox/system.py` — `AssembledSystem` (merged mol, reactive-site registry)
- [x] `simbox/library.py` — `MoleculeLibrary` (Epoxy-PDMS, Amino-PDMS, AM0270-POSS)
- [x] `simbox/writer.py` — LAMMPS data file writer; fixed header bug (`max` not `len`); added `ff_coeffs.in`
- [x] `simbox/inputs.py` — 4-stage LAMMPS scripts (push-off → minimize → NVT → NPT → crosslink template)
- [x] `generate_simbox_crosslink.py` workflow (canonical, v4)
- [x] `docs/simbox.md` — full API documentation with output structure and version history
- [x] `tests/unit/test_simbox.py` — 31 unit tests (Molecule, BoxPacker, assemble)
- [x] `tests/regression/test_simbox_crosslink.py` — 10 regression tests (poss_0, v4 reference)
- [x] `tests/regression/test_simbox_poss.py` — 26 regression tests (poss_50 + poss_100)
- [x] `tests/output/simbox_crosslink/v4/` — reference outputs: poss_0, poss_50, poss_100
- [x] CLI entry point for simbox (`topon simbox` — `topon/simbox/workflow.py` + `topon/cli.py`)

## Known Code Issues (from 2026-03-23 audit)
- [x] **CRITICAL**: `chemistry/builder.py` ~line 509–534 — SMILES fallback now raises `ValueError` instead of silently returning a single monomer; caller emits `RuntimeWarning` and skips the edge
- [x] **HIGH**: `forcefield/dreiding.py` lines 23–24 — removed debug `print`; default param path now resolved from `__file__` via `_BUNDLED_PARAM_FILE`
- [x] **HIGH**: `forcefield/dreiding.py` lines 265–282 — all wildcard entries now use `dict` format matching the parser (bond/angle/improper were tuples)
- [x] **MEDIUM**: `analysis/report.py` line 159 — `count//2` is the correct topological maximum (each entanglement consumes 2 chains); not an estimate
- [x] **MEDIUM**: `assignment/manager.py` — `max_entanglements` changed from `0` to `None` to correctly signal "not computed" (0 was falsy but misleading)
- [x] **MEDIUM**: `chemistry/sequences.py` line 72 — gradient copolymer arrangement implemented (linear weight interpolation); unknown arrangements now emit `RuntimeWarning`
- [x] **MEDIUM**: `topology/generator_python.py` — BCC (8-neighbor) and FCC (12-neighbor) lattices implemented, matching C `create_bcc/fcc_lattice` exactly
- [x] `assignment/manager.py` — `_assign_grafts()` and `_assign_copolymers()` implemented; `ChemistryBuilder._build_chain_cg()` consumes `monomer_sequence`, `graft_positions`, `graft_dp`, `graft_monomer`

---

## Notes
- Do NOT modify files in `older_versions/` folder
- C generator uses `generator_serial_debug11.c` (final version)
- KG options: angle potentials, repulsive/attractive LJ
- GUI shows validation, mutual exclusivity for copolymer/edge-typing
- **CRITICAL**: Coordinates go directly into data file (Python-side), NOT via LAMMPS include commands
- Default density: 0.9 g/cm³
- Charges: Optional toggle, default ON for atomistic
- `tests/protein_network/topro/` is a separate BFM protein network project — not part of topon
- Root-level `output/` directory contains pre-generated topology files (should move to `tests/data/networks/`)
