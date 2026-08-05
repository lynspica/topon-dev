# topon — Architecture

A graph-first Python toolkit for building polymer networks — and, via **topro**, protein networks — ready to drive LAMMPS molecular dynamics. Topology is decided as a graph first, then mapped into chemistry and coordinates; the same graph can produce a coarse-grained Kremer-Grest system or a fully atomistic DREIDING system without changing the topology code.

This document is the **primary onboarding doc for new contributors and AI agents**. Read it first; then go to [USAGE.md](USAGE.md) for how to run things, [DEVELOPMENT.md](DEVELOPMENT.md) for the history of how the package got here, and [`CLAUDE.md`](../CLAUDE.md) at the repo root for the rules every change must follow.

---

## 1. The topon family

Three peer sub-systems share the same Python package. They do related but distinct jobs.

| | **core topon** (polymer networks) | **topro** (protein networks) | **simbox** (molecule packing) |
|---|---|---|---|
| Implementation path | `topon/{topology,assignment,chemistry,conformation,writers}/` | `topon/protein_network/` | `topon/simbox/` |
| Entry point | `topon` CLI / `topon.pipeline.Pipeline` | `python -m topon.protein_network` | `topon simbox` CLI / `topon.simbox` API |
| Resolution | atomistic (DREIDING) or coarse-grained (Kremer-Grest) | coarse-grained (MARTINI 3) | atomistic (DREIDING) |
| Force field | DREIDING / Kremer-Grest | MARTINI 3 (vendored from `Martini_Ahmet/nat_pro.itp`) | DREIDING |
| LAMMPS data | `atom_style full`, wrap-only (7-column atom rows, no image flags) | `atom_style full`, wrap-only (same convention as core topon) | `atom_style full`, wrap-only |
| Pipeline | the six stages below | BFM lattice → JSON snapshot → 3-stage relax (parallel pipeline) | independent packing flow |
| Crosslinks | Y-merge (CG) / chemistry-defined (atomistic) | dityrosine SC4–SC4 (TN6 bead) | reaction templates (epoxy/amine, etc.) |
| Topology shape | lattice graph (SC / BCC / FCC, configurable functionality) | BFM self-avoiding-walk lattice | molecule library + grid packing |

A fourth, smaller utility — `topon/singlechain/` — handles single-chain solubility computations and is used by neither the main pipeline nor topro directly.

> **Naming note.** The integrated module is `topon.protein_network`. The user-facing name for that capability is **topro** (topological protein network — parallel to "topon", topological polymer network). Use "topro" in docs and discussion; use `topon.protein_network` in import statements.

---

## 2. The six-stage pipeline

The `Pipeline` class in [`topon/pipeline.py`](../topon/pipeline.py) is the single orchestrator for **core topon**. Its `run()` method calls six stages in order. Topro does not use this orchestrator — it has its own pipeline (see §6).

```
config.json
    │
    ▼
ConfigLoader.load() ──► ToponConfig (Pydantic) ──► Pipeline.run()
                                                       │
   Stage 1: Topology ──────────────────────────────────┤
   Stage 2: Analysis ──────────────────────────────────┤
   Stage 3: Assignment ────────────────────────────────┤
   Stage 4: Chemistry ─────────────────────────────────┤
   Stage 5: Conformation ──────────────────────────────┤
   Stage 6: Output ────────────────────────────────────┘
                                                       │
                                                       ▼
                                          <study>/<name>/
                                            ├── topology/
                                            ├── 02_Chemistry/
                                            ├── 03_Conformation/
                                            └── 04_Simulation/
```

### Stage 1 — Topology
**Module:** `topon/topology/` &nbsp;**Code:** `pipeline.py:91-140`

Generates or loads a NetworkX `MultiGraph`. Nodes are network junctions; edges are polymer chains. Two paths:
- `source="generate"` — calls the C generator (`generator.exe`) via `topon.topology.generator.run_generator`, then loads `.nodes`/`.edges` files.
- `source="load"` — reads existing `.gpickle` or `.nodes`+`.edges` files via `topon.topology.loader.load_graph`.

Produces: `self.graph` (annotated `MultiGraph`) and `self.dims` (box size as `np.ndarray`).

**The periodic cell is recorded, not inferred.** Generators write the exact
repeat distance into `G.graph["box"]`, and `.nodes` files carry it in a
`# BOX Lx Ly Lz` header. `infer_dims_from_graph` returns that value when it
is present and only falls back to estimating the cell from the coordinate
extent (`max - min + 1`) for graphs written before this existed. The estimate
is exact for SC, whose sites are integer-spaced, but overshoots any lattice
with fractional basis sites: BCC and FCC body/face sites sit at +0.5 and
Diamond sites at quarter-cell offsets, so a 4x4x4 BCC or FCC reported 4.5.
Because `self.dims` is the box every minimum-image calculation uses, that
overshoot violated the `bond < box/2` invariant in Design Principle 3 and
sent roughly a third of BCC edges (a quarter of FCC) to the wrong periodic
replica, where they were built at twice their true bond length. Any new
lattice with non-integer sites must record its cell for the same reason.

### Stage 2 — Analysis
**Module:** `topon/assignment/manager.py:analyze` &nbsp;**Code:** `pipeline.py:146-154`

Computes graph statistics — degree distribution, connectivity, defect/entanglement *capacity* of the topology — before any annotations are written. The analysis result is held in `self.analysis_report`. Read-only; does not modify the graph.

(Distinct from the `topon/analysis/` package, which exposes `analyze_graph()` for the standalone `topon analyze` CLI sub-command. `Pipeline` does not import from `topon.analysis`.)

### Stage 3 — Assignment
**Module:** `topon/assignment/` &nbsp;**Code:** `pipeline.py:160-163`

Annotates the graph in place. `AssignmentManager` (`assignment/manager.py`) orchestrates several module-level assigners:

| Concern | API | Lives in |
|---|---|---|
| Orchestration | `AssignmentManager.run()` | `assignment/manager.py` |
| Node types (degree-based, positional, random, explicit) | `assign_node_types()` | `assignment/node_types.py` |
| Edge types (uniform, random, composite) | `assign_edge_types()` | `assignment/edge_types.py` |
| Degree-of-polymerization distribution (Schulz-Zimm PDI) | `assign_dp()` | `assignment/dp_distribution.py` |
| Copolymer sequences (Block, Random, Alternating, Gradient) | `generate_monomer_sequence()` | `assignment/sequences.py` |
| Defects (primary loops / parallel edges) | `inject_primary_loops()` | `assignment/defects.py` |
| Entanglements (Gaussian Kink geometry) | `select_entanglements()` | `assignment/entanglements.py` |
| Re-attribute existing graph (round-trip workflows) | `GraphAttributor` (class) | `assignment/attributor.py` |

After this stage the graph carries every annotation downstream stages need.

### Stage 4 — Chemistry
**Module:** `topon/chemistry/` &nbsp;**Code:** `pipeline.py:169-245`

`ChemistryBuilder` (in `chemistry/builder.py`, one unified entry point — *not* split into CG/atomistic classes) builds an RDKit `Mol` with 3D coordinates derived from the graph. Supporting modules:

| Module | Role |
|---|---|
| `chemistry/builder.py` | The `ChemistryBuilder` class; bulk of stage-4 logic |
| `chemistry/sequences.py` | Shared monomer sequence helpers used by the builder |
| `chemistry/embed.py` | ETKDGv3 per-chain conformer embedding |
| `chemistry/{dreiding,kg,charmm}/` | **Stub sub-namespaces** reserved for future force-field-specific chemistry; not yet implemented (`kg/` has only a docstring; the others are 1-line stubs). Today's CG vs atomistic switch happens inside `builder.py` via `config.chemistry.model_type`. |

This stage also writes the first set of LAMMPS-relevant outputs to `<output_dir>/02_Chemistry/`:

- `system.data` — atom positions and connectivity (via `CGWriter` or `DreidingWriter` in `topon/writers/`)
- `system.in.settings` — force-field stub (filled by stage 6 for atomistic)
- `system.groups` — `nodes` and `beads` LAMMPS group definitions
- `system_nodes.displace`, `system_beads.displace` — displacement files for stage 5

### Stage 5 — Conformation
**Module:** `topon/conformation/` &nbsp;**Code:** `pipeline.py:251-271`

`ConformationManager` (`conformation/manager.py`) is the entire stage-5 implementation: it reads the chemistry-stage data file, applies the displacement files, adds a small Gaussian noise to break degeneracy, and resolves overlaps iteratively. The sub-namespaces `conformation/{entanglement,packing,placement}/` are 1-line stubs reserved for a future refactor of `manager.py` into smaller pieces.

Conformation defaults (`pipeline.py:47`):

```python
{"overlap_cutoff": 0.01, "overlap_max_iters": 10, "noise_magnitude": 1e-4}
```

Output: `<output_dir>/03_Conformation/system_relaxed.data`.

**The simulation box comes from stage 1, not from the coordinates.** Callers
pass `lattice_box=dims` into `apply_displacements`, giving a box of
`dims * scale`; only when it is omitted does the manager fall back to
estimating `(max node coord + 1) * scale` from the `.displace` files. This
has to match the cell stage 4 used to route chains across the periodic
boundary. When the two disagree, a chain that wraps under one period lands
in a box of another and its closing bond is left stretched across the
system. The two estimates happen to coincide for SC, which is why the
fallback survived so long, and passing the box also makes the written box
exactly `volume^(1/3)`, so the target density is hit on every lattice.

### Stage 6 — Output
**Module:** `topon/writers/` (`LammpsInputGenerator`) &nbsp;**Code:** `pipeline.py:277-300`

Writes the LAMMPS input scripts that drive the simulation:
- a serial-soft-minimization script (relaxes overlaps with a soft pair potential)
- a parallel-production script (runs the actual KG/DREIDING dynamics)

Output: `<output_dir>/04_Simulation/*.in`.

The simulation itself is *not* run by `Pipeline.run()`; that is the job of `topon/simulation/SimulationRunner` (or a manual `lmp` invocation).

---

## 3. Module map

Concise tour of every directory under `topon/`. Every module above the dashed line is part of the main pipeline; sub-systems below are independent.

| Directory | Files / Subdirs | Role |
|---|---|---|
| `topology/` | 4 files + `network/`, `sequence/`, `simple/` (sub-dirs are stubs) | Graph generation and loading. Stage 1. |
| `analysis/` | 2 files | Read-only graph statistics for the `topon analyze` CLI. **Not used by `Pipeline` stage 2** — that calls `AssignmentManager.analyze()`. |
| `assignment/` | 8 files | Graph annotation: node/edge types, DP, defects, entanglements, copolymers. Stage 3. |
| `chemistry/` | 4 files (`builder.py` is most of stage 4); `dreiding/`/`kg/`/`charmm/` are stubs | RDKit Mol construction with 3D coords. Stage 4. |
| `conformation/` | 2 files (`manager.py` is most of stage 5); `entanglement/`/`packing/`/`placement/` are stubs | Overlap resolution and conformation. Stage 5. |
| `writers/` | 6 files | LAMMPS data and input-script writers. Stage 4 + Stage 6. |
| `forcefield/` | 4 files | DREIDING parameter parser; Kremer-Grest parameters. Read by chemistry/writers. |
| `config/` | 4 files | Pydantic `ToponConfig` schema; `ConfigLoader.load()`. |
| `core/` | 2 files | Shared types and protocols. |
| `utils/` | 4 files | Shared helpers (e.g. `write_lammps_displacement_file`). |
| `pipeline.py` | — | The `Pipeline` orchestrator class. |
| `cli.py`, `__main__.py` | — | CLI dispatch. |
| `workflows/` | 4 files | High-level workflow helpers used by `tests/workflows/` scripts. |
| — | — | — |
| `protein_network/` | 16 files + `data/` | **topro**: BFM-based CG MARTINI 3 protein networks. Parallel pipeline. |
| `simbox/` | 8 files | Independent molecule packer + crosslink-template emitter. DREIDING-only. |
| `singlechain/` | 3 files | Single-chain solubility utility. |
| `simulation/` | 2 files + `protocols/` | LAMMPS subprocess runner. |

---

## 4. Design principles

These are the rules every change is expected to follow. Several of them were learned the hard way — see the `legacy/` archive (gitignored) and [DEVELOPMENT.md](DEVELOPMENT.md) for the history.

1. **Graph-first separation.** Topology is decided as a NetworkX graph before chemistry is built. The same graph can produce a CG system or an atomistic system by switching the chemistry builder. *Topology never sees atom types, coordinates, or force-field details.*

2. **Six-stage pipeline, one-way data flow.** Each stage has a single responsibility. Downstream stages do not reach back into upstream state. The stage table in §2 and the module-boundary table in [`CLAUDE.md`](../CLAUDE.md) are authoritative.

3. **LAMMPS data convention is module-specific** as of 2026-05. The core topon and simbox writers (`topon/writers/lammps_data.py`, `topon/simbox/writer.py`) emit 7-column Atoms rows (no `ix iy iz`) and rely on LAMMPS's neighbor / ghost-atom system to handle PBC via min-image, which is correct as long as every bond is shorter than `box/2` (true for KG / DREIDING coarse-grained networks where chains rarely wrap differently from each other). The protein-network writer (`topon/protein_network/lammps_writer.py`) emits **10-column Atoms rows** with `ix iy iz` computed by a **priority-weighted MST** (Kruskal, sort key `(priority, length)` with priority 0 for non-crosslink bonds and priority 1 for crosslinks) over the molecular bond graph — this is necessary because MARTINI 3 protein networks on a BFM lattice contain many chains that wrap multiple times around the box and are crosslinked across those wraps, producing bond pairs whose image-flag-implied unwrapped distance exceeds the ghost-shell cutoff in parallel MPI runs unless image flags are written explicitly. The priority key is required because BFM merges two dityrosine SC4 beads onto the same lattice node, giving every crosslink a wrapped-MIC distance ≈ 0.05 Å (much shorter than the ≈ 6.7 Å BB-BB bonds at the projected scale); a pure length-sorted MST would add crosslinks to the spanning tree first and then drop BB-BB bonds as winding-cycle back-edges. The priority key inverts that so crosslinks are processed last and the cycle-closing back-edges are always the crosslinks themselves — matching the design intent that crosslinks are the redundant, sacrificial elements. A hard assertion guarantees no real funct=integer non-crosslink bond may ever drop. See `topon/protein_network/lammps_writer.py:_kruskal_image_flags_and_drop` and `docs/JOURNAL.md` 2026-05-19 entries.

4. **xyz perturbation handles BFM degeneracy in topro.** The BFM lattice can place two crosslink endpoints from different SAW walks at the same lattice node; in the wrapped frame those bead positions are coincident, and the bond/angle gradient hits division-by-zero in stage 1. The fix is a tiny Gaussian jitter (`coord_perturbation_ang`, σ ≈ 0.05–0.1 Å) at write time — see `topon/protein_network/lammps_writer.py:188-217`. This mirrors core topon's hierarchical-relaxation precondition. **Do not** modify the calibrated stage-1 LAMMPS scripts to work around degeneracies — fix the input geometry instead.

5. **Calibrated LAMMPS scripts are not to be modified.** The `_stage1_soft`, `_stage1_hierarchical`, `_stage2_ljramp`, `_stage3_min_nvt_npt` writers in `topon/protein_network/lammps_writer.py` are tuned by hand. If a run shows force overflows, hot ramps, or NaNs, the fix is in the **input geometry** — xyz perturbation (`coord_perturbation_ang`), hierarchical stage 1 (`hierarchical_stage1=True`), or BFM topology tightening — not in the integrator, thermostat, timestep, or pair_style.

6. **Configuration via Pydantic, no globals in stage code.** Every stage module reads from a `ToponConfig` (or its raw-dict supplements). No `os.environ` lookups, no module-level constants that change behaviour, no hard-coded paths inside `topon/`. Config is loaded in `Pipeline` or in workflow scripts under `tests/workflows/`.

7. **Module boundaries** (from [`CLAUDE.md`](../CLAUDE.md)):

| Module | Does | Does NOT |
|---|---|---|
| `topology/` | Generate connectivity graph | Assign types, coordinates, force field |
| `assignment/` | Assign node/edge types, DP, defects, entanglements | Generate chemistry or coordinates |
| `chemistry/` | Build RDKit molecular structure from graph | Generate placement coordinates |
| `conformation/` | Place atoms in 3D, resolve overlaps | Assign force-field types |
| `writers/` | Format and write LAMMPS files | Computation of any kind |
| `analysis/` | Compute graph statistics | Modify the graph |
| `simbox/` | Independent molecule packing sub-system | Interact with the main pipeline |

8. **Output convention.** End-to-end run artifacts live in `tests/output/<vNN>/<cell>/` (gitignored as of 2026-05-08). Sweep drivers go in `tests/workflows/run_<name>.py`. Do not invent `runs/` or other parallel folders.

9. **Regression coverage on writer changes.** Per `CLAUDE.md`, before modifying anything in `topon/writers/` or `topon/simbox/writer.py`, run `pytest tests/regression/` first to confirm the baseline; make the change; re-run. For `topon/protein_network/lammps_writer.py` the minimum is `pytest tests/unit/`.

---

## 5. Configuration model

The package uses **Pydantic** (`topon/config/`) to validate JSON configuration. Some sections that have not yet been migrated to schema (`conformation`, `simulation`, `experimental`) are passed through as a raw dict alongside the validated config — see `Pipeline.__init__` (`pipeline.py:53`).

Top-level `ToponConfig` sections, in order of pipeline use:

| Section | Used by | Notes |
|---|---|---|
| `study` | all stages | `study.name`, `study.output_dir` |
| `topology` | Stage 1 | `source` (`"generate"` / `"load"`), `lattice_size`, `degree_distribution`, etc. |
| `assignment` | Stage 3 | sub-objects for entanglements, defects, grafts, copolymer sequences |
| `chemistry` | Stage 4 | `model_type` (`"coarse_grained"` / `"atomistic"`), `target_density` |
| `conformation` (raw) | Stage 5 | `overlap_cutoff`, `overlap_max_iters`, `noise_magnitude` |
| `simulation` (raw) | Stage 6 | LAMMPS pair_style, angle handling |
| `execution` (raw) | post-Stage 6 | LAMMPS subprocess settings |
| `experimental` (raw) | Stage 6 | feature-flagged extras |

For the full key-by-key reference and example configs, see [USAGE.md](USAGE.md).

---

## 6. Topro — the parallel protein-network pipeline

`topon/protein_network/` does not run through `Pipeline`. It has its own three-stage flow, designed for MARTINI 3 IDP protein networks.

```
CLI args (block_seq, n_repeats, n_chains, water_density, ...)
    │
    ▼
BFM SAW lattice topology  ──►  JSON topology snapshot
    │
    ▼
Chemistry build (residue lookup from Martini_Ahmet/nat_pro.itp)
    │
    ▼
LAMMPS data + 3-stage relaxation scripts
    │
    ▼
(optional) 14-stage MARTINI annealing port
```

Key topro design decisions, captured for reference:

- **Library lookup, not auto_martini.** Residue topologies are extracted from `tests/_martini_extracted/Martini_Ahmet/itp_files/{nat_pro,high_pro,no_pro}.itp` by `tools/extract_residues_from_itp.py`. Eight amino acids are currently covered (GLY, ALA, ARG, PRO, SER, ASP, TYR, ASN — the resilin reference set). To add residues, re-run the extractor against a polyply ITP that contains them.
- **Wrap-only data + xyz perturbation.** topro uses the same wrap-only convention as core topon (no image flags). The xyz-perturbation hack handles BFM lattice degeneracy at stage-1 write time — see Design Principle 4.
- **Oversized crosslinks are dropped.** BFM merge-site semantics let two chains' SC4 atoms land at the same lattice node from different walks (different image cells). Crosslinks where the two atoms are >`box/4` apart unwrapped are dropped with a warning. Intra-chain bonds always survive.
- **Reaction-field electrostatics approximated.** GROMACS RF (`coulombtype = reaction-field, ε_r=15`) has no exact LAMMPS equivalent. Approximated with `pair_style lj/cut/coul/cut 12.0` + `dielectric 15.0`; loses the RF correction term but uses stock LAMMPS only.
- **Restricted-bending angles approximated.** MARTINI 3 IDP backbone uses GROMACS `funct=10` (`½ K (cos θ − cos θ₀)² / sin²θ`). Stock LAMMPS lacks the `1/sin²θ` factor; we use `cosine/squared`, which matches the leading term and is good enough for IDP relaxation. Folded proteins should use a different angle style.

Known topro limitations (no virtual sites, no elastic network, no NPT for dry boxes, NaCl ions not packed, residue coverage limited to 8 AA) are tracked in [`internal/DEVELOPMENT_INTERNAL.md`](../internal/DEVELOPMENT_INTERNAL.md).

---

## 7. CLI surface

The `topon` CLI (`topon/cli.py`) dispatches to:

| Sub-command | Backend | Status |
|---|---|---|
| `topon generate` | `Pipeline.run()` | Working |
| `topon validate` | `ConfigLoader.load()` | Working |
| `topon init` | scaffolds a new study folder | Working |
| `topon simbox` | `topon.simbox` API | Working |
| `topon chain` | `topon.singlechain` | Working |
| `topon analyze` | `topon.analysis` | **Placeholder** (prints stub) |
| `topon gui` | — | **Not implemented** |
| `python -m topon.protein_network` | topro pipeline | Working |

Full flag-by-flag reference: [USAGE.md](USAGE.md).

---

## 8. Testing layout

| Test class | Location | Run with | Wall-clock |
|---|---|---|---|
| Unit | `tests/unit/` | `pytest tests/unit/` | ~5 s |
| Regression (golden byte-equivalence on writers) | `tests/regression/` | `pytest tests/regression/` | ~1.5 h |
| End-to-end workflows (generate runnable systems) | `tests/workflows/run_<name>.py` | `python <script>` | minutes–hours |

Run outputs land in `tests/output/<vNN>/<cell>/` (gitignored). Sweep drivers (`run_v41_matrix.py`, `run_v42_matrix.py`, `run_v43_core_topon.py`) parameterize the workflow scripts and write into versioned subdirectories.

Regression references are frozen with `seed=42` — see the regression-test surface notes in agent memory before any refactor.

---

## 9. Where to go next

| Question | Doc |
|---|---|
| How do I run topon end-to-end? | [USAGE.md](USAGE.md) |
| Which CLI flag does X? | [USAGE.md](USAGE.md) |
| What's the JSON config schema? | [USAGE.md](USAGE.md) (config-reference appendix) |
| What's changed across versions? | [DEVELOPMENT.md](DEVELOPMENT.md) |
| What's the dev methodology? | [DEVELOPMENT.md](DEVELOPMENT.md) (header) |
| What's planned next, what's broken? | [`internal/DEVELOPMENT_INTERNAL.md`](../internal/DEVELOPMENT_INTERNAL.md) (gitignored, local only) |
| What rules apply to AI-agent edits? | [`CLAUDE.md`](../CLAUDE.md) |
| Code review for a change? | invoke the `topon-reviewer` agent |
| Independent audit of a change? | invoke the `investigator` agent |
