# AGENTS.md — Read this first

> **For any AI agent (Claude Code, ChatGPT, Cursor, Copilot, etc.) starting work on the topon project.** Read this file end-to-end, *then* the four canonical docs in §3 below. After that, you have enough context to take instructions from the user.

---

## 1. What topon is

A modular Python package for generating polymer-network and protein-network structures for LAMMPS molecular dynamics. Topology is decided as a graph first, then mapped to chemistry and 3D coordinates. Three sub-systems share the same package:

- **core topon** — polymer networks, atomistic (DREIDING) or coarse-grained (Kremer-Grest)
- **topro** — `topon.protein_network`, MARTINI 3 protein networks (V36)
- **simbox** — independent molecule packer for crosslinking studies

History, methodology, and design rationale live in the docs in §3.

---

## 2. File-tree map (committed-public unless noted)

```
topon/
├── AGENTS.md                    ← this file
├── README.md                    project front page
├── CLAUDE.md                    rules every change must follow (read second)
├── pyproject.toml               package metadata
├── .gitignore                   gitignored locations
├── topon/                       the package itself
│   ├── pipeline.py              the 6-stage Pipeline orchestrator
│   ├── cli.py                   `topon` CLI dispatch
│   ├── topology/                stage 1
│   ├── analysis/                stage 2 helper (also CLI `topon analyze`)
│   ├── assignment/              stage 3
│   ├── chemistry/               stage 4 (builder.py is the real work)
│   ├── conformation/            stage 5 (manager.py is the real work)
│   ├── writers/                 stage 6
│   ├── protein_network/         topro (parallel pipeline, V36)
│   ├── simbox/                  independent molecule packer
│   ├── singlechain/             single-chain solubility utility
│   ├── config/                  Pydantic schemas + load_config
│   ├── forcefield/              DREIDING, Kremer-Grest parameters
│   ├── workflows/               higher-level orchestrators used by tests/workflows/
│   ├── simulation/              SimulationRunner (LAMMPS subprocess)
│   ├── core/                    shared types
│   └── utils/                   shared helpers
├── docs/
│   ├── USAGE.md                 CLI + APIs + recipes + JSON-config schema
│   ├── ARCHITECTURE.md          6-stage pipeline + module map + design principles
│   ├── DEVELOPMENT.md           V1–V36 changelog + Q→R→R→R methodology
│   ├── manuscript.pdf           npj submission preprint
│   └── martini_annealing_template.in  LAMMPS template referenced by topro
├── examples/
│   ├── README.md                atlas of all demos
│   ├── run_via_api.py           generic Python-API runner (equivalent to `topon generate`)
│   ├── templates/               minimal.json + full.json starters
│   ├── defaults/                shared assignment-fragment configs
│   ├── showcase/                small reference data files (input format)
│   ├── demos/
│   │   ├── polymer/             atomistic + CG, all architecture knobs
│   │   ├── protein/             charmm (legacy pointer) + martini (current)
│   │   ├── topology/            end_linking (C vs Python) + bfm
│   │   └── poss/                POSS-junction atomistic
│   └── npjcompmat/              paper companion dataset (1001 files)
├── tests/
│   ├── unit/                    pytest tests/unit/ (~5s)
│   ├── regression/              pytest tests/regression/ (~1.5h)
│   ├── workflows/               run with `python <path>` (sweep drivers)
│   ├── data/                    pre-generated input topologies
│   ├── sample_graphs/           small sample networks
│   └── output/                  GITIGNORED — accumulated run artifacts
├── tools/                       maintenance scripts (e.g. extract_residues_from_itp.py)
├── .claude/
│   ├── agents/topon-reviewer.md      code-review agent (has Edit/Write)
│   └── agents/investigator.md        unbiased read-only auditor
├── internal/                    GITIGNORED — owner-local notes (see §5)
└── legacy/                      ARCHIVED OUT OF REPO 2026-05-09 → `C:\Users\ahmet\topon_archive\`
```

---

## 3. Reading order — after you finish AGENTS.md

1. **[CLAUDE.md](CLAUDE.md)** — short, hard rules; module boundaries; regression-test requirements; output conventions; "do not modify" list.
2. **[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)** — six-stage pipeline detail; module-by-module map; design principles (wrap-only LAMMPS data, calibrated scripts, etc.); topro section (§6); CLI surface (§7).
3. **[docs/USAGE.md](docs/USAGE.md)** — CLI flags, sub-system APIs, recipes by demo, full JSON-config schema (Appendix A).
4. **[docs/DEVELOPMENT.md](docs/DEVELOPMENT.md)** — version-by-version history (V1–V36), Q→R→R→R methodology, simbox sub-version timeline.

These four are sufficient for ~90 % of tasks. For local-only context (open issues, planned next steps, source-side drift), see [`internal/DEVELOPMENT_INTERNAL.md`](internal/DEVELOPMENT_INTERNAL.md) — present locally on the owner's machine, gitignored, may not exist for new contributors.

---

## 4. Hard conventions (enforced; CLAUDE.md is authoritative)

- **Six-stage pipeline.** Topology → Analysis → Assignment → Chemistry → Conformation → Output. Each stage has one job; downstream never reaches into upstream state.
- **Wrap-only LAMMPS data.** Atom positions in `[0, box)`; 7-column atom rows, no `ix iy iz` image flags. Applies to core topon, simbox, *and* topro. The xyz-perturbation hack handles BFM degeneracy in topro — see ARCHITECTURE §6.
- **Don't modify the calibrated LAMMPS scripts** in `topon/protein_network/lammps_writer.py` (`_stage1_soft`, `_stage1_hierarchical`, `_stage2_ljramp`, `_stage3_min_nvt_npt`). If a run misbehaves, fix the *input geometry*.
- **Configuration via Pydantic** (`topon/config/schema.py`). No globals, no `os.environ`, no hard-coded paths inside stage code.
- **Output convention.** End-to-end run artifacts live in `tests/output/<vNN>/<cell>/` (gitignored). Sweep drivers go in `tests/workflows/run_<name>.py`. Don't invent `runs/` or other parallel folders.
- **Regression coverage on writer changes.** Run `pytest tests/regression/` before touching `topon/writers/` or `topon/simbox/writer.py`.

---

## 5. Local-only state

Several locations are gitignored — they exist on the project owner's machine but won't be present in fresh clones:

- `internal/` — owner's working notes, planned next steps, P0/P1 issue tracker
- `internal/specs/` — feature specs from collaborators (e.g. NPZ output format)
- `legacy/` — moved out of the repo entirely on 2026-05-09 to `C:\Users\ahmet\topon_archive\` (8.3 GB of frozen earlier attempts; reference for fork-provenance comments only)
- `tests/output/` — accumulated LAMMPS run artifacts
- `.vscode/`, `.pytest_cache/`, build artefacts

If a tracked file references a `legacy/...` path or `internal/...` path, treat that as documentation/provenance only — do not assume the path resolves on a fresh clone.

---

## 6. Available agents (in `.claude/agents/`)

- **`topon-reviewer`** — strict project-aware code reviewer. Has `Edit`/`Write`. Use for code-level fixes, LAMMPS-output debugging, regression triage.
- **`investigator`** — unbiased read-only auditor. No `Edit`/`Write` by design. Use for: scientific-correctness audits, doc-consolidation reviews, cross-module change verification, pre-commit sanity checks. Pattern: *draft → investigator review → fix → commit.*

The investigator is the recommended last step before any non-trivial commit.

---

## 7. Git topology

- Default branch: `main`.
- Two remotes:
  - `personal` → `https://github.com/lynspica/topon-dev` — frequent dev pushes, public.
  - `stable` → `https://github.com/keten-group/topon` — paper-companion (v0.1.0); only push there at deliberate release milestones.
- `git push` (no remote) defaults to `personal`. Push to `stable` only with explicit `git push stable main`.

---

## 8. Known open issues — read `internal/DEVELOPMENT_INTERNAL.md`

Several latent issues are tracked in `internal/DEVELOPMENT_INTERNAL.md` (gitignored). The biggest ones a new agent should know:

- **P0-A** — schema gap: `topon generate <config>` strict-validates against `ToponConfig` (extra-forbid), but most demo configs have `conformation`/`simulation`/`execution` sections that aren't in the schema. Those configs only run via the workflow scripts under `tests/workflows/` or via direct `Pipeline(ToponConfig, raw_config={...})` Python wiring. The `examples/demos/topology/` and `examples/demos/poss/` configs are schema-clean and work via `topon generate` today.
- **P0-2** — silent `node_type` fallthrough to `Si` for unknown node types in `topon/chemistry/builder.py::_build_nodes`.

If the user asks you to do something that touches one of these, surface the trade-off explicitly rather than working around it silently.

---

## 9. Saying hi to the user

Once you've finished reading AGENTS.md and the four canonical docs, ask the user what they want to work on. Don't pre-emptively make changes; the user gives instructions and you execute them with full context.
