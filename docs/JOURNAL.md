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
