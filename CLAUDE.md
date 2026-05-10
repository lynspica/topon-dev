# CLAUDE.md — Development Rules for Topon

These rules apply when making changes to this codebase.

---

## Architecture

Topon has a six-stage pipeline:

```
Topology → Analysis → Assignment → Chemistry → Conformation → Output
```

Each stage is a separate module. Keep responsibilities clean:

| Module | Does | Does NOT |
|---|---|---|
| `topon/topology/` | Generate connectivity graph | Assign types, coordinates, force field |
| `topon/assignment/` | Assign node/edge types, DP, defects, entanglements | Generate chemistry or coordinates |
| `topon/chemistry/` | Build RDKit molecular structure from graph | Generate coordinates |
| `topon/conformation/` | Place atoms in 3D, resolve overlaps | Assign force field types |
| `topon/writers/` | Format and write LAMMPS files | Computation of any kind |
| `topon/analysis/` | Compute graph statistics | Modify the graph |
| `topon/simbox/` | Independent molecule packing sub-system | Interact with the main pipeline |

---

## Where New Code Goes

| What you're adding | Where it goes |
|---|---|
| New lattice type or topology algorithm | `topon/topology/` |
| New graph attribute assignment | `topon/assignment/` |
| New monomer / molecule building logic | `topon/chemistry/` |
| New coordinate generation method | `topon/conformation/` |
| New LAMMPS output format | `topon/writers/` |
| New analysis metric | `topon/analysis/` |
| New molecule for simbox packing | `topon/simbox/library.py` |
| New end-to-end workflow script | `tests/workflows/` |

Do not create new top-level modules without discussing first.

---

## LAMMPS Output Is Regression-Tested

The LAMMPS data file format is sensitive. Small changes (whitespace, column order, section ordering) silently break simulations.

Before modifying anything in `topon/writers/` or `topon/simbox/writer.py`:
1. Run `pytest tests/regression/` — confirm passing
2. Make the change
3. Run regression tests again — confirm still passing
4. If no regression test exists for that writer, write one first

---

## Configuration

All behaviour is controlled through Pydantic config objects (`topon/config/schema.py`). Do not use hard-coded paths, global variables, or `os.environ` inside stage modules. Config loading happens in `pipeline.py` or in workflow scripts, not inside stage modules.

---

## Do Not Modify

- `tests/output/` — regression reference outputs (update only by running the generator with `--update-refs`)

> Note: `legacy/` was moved out of the repo on 2026-05-09 to `C:\Users\ahmet\topon_archive\` (with `older_versions/` inside it). Don't re-import either tree into the working copy; if you need a historical reference, open the file in place from the archive.

---

## Testing

Tier markers are registered in `pyproject.toml` and auto-applied by `tests/conftest.py` based on each test's parent directory:

```bash
pytest -m fast                       # fast unit tests (~5s)
pytest tests/unit/chemistry/         # focused on a component you just changed
pytest -m "fast or smoke"            # major pre-push check (smoke includes LAMMPS if `lmp` is on PATH)
pytest -m regression                 # full byte-equivalence suite (~1.5h)
```

Unit tests live in `tests/unit/<component>/` (one subdir per `topon/` sub-package). Smoke tests live in `tests/smoke/` and run end-to-end pipeline + LAMMPS at small scale. Regression tests live in `tests/regression/`. Integration / workflow scripts live in `tests/workflows/` (not pytest — run directly with Python).

`@pytest.mark.requires_lammps` is auto-skipped when `lmp` is not on PATH.

---

## Docs

User-facing docs live in `docs/` (committed) and `internal/` (gitignored, owner-local). Keep them updated when changing behaviour:

- [`docs/USAGE.md`](docs/USAGE.md) — update when changing CLI options, config schema keys, or sub-system APIs
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) — update when changing module structure, pipeline stages, or design principles
- [`docs/DEVELOPMENT.md`](docs/DEVELOPMENT.md) — add a changelog entry for every significant change
- [`docs/JOURNAL.md`](docs/JOURNAL.md) — append a dated journal entry (Change / Why / Issue+solution) for any non-trivial work
- [`internal/DEVELOPMENT_INTERNAL.md`](internal/DEVELOPMENT_INTERNAL.md) (local only) — track open issues, planned next steps, source-side drift

When a change touches multiple of these, update them in the same commit. New AI agents in fresh sessions should read `ARCHITECTURE.md` first.

---

## Agents

Two specialized review agents live under `.claude/agents/`:

- [`topon-reviewer`](.claude/agents/topon-reviewer.md) — topon-specific code reviewer (LAMMPS conventions, force-field rules, has `Edit`/`Write` authority for fixes).
- [`investigator`](.claude/agents/investigator.md) — unbiased read-only auditor for scientific correctness, doc consistency, and cross-module changes (no `Edit`/`Write` tools by design).

Pattern: draft → `investigator` review → fix → commit. Use `topon-reviewer` for code-level fixes; use `investigator` before declaring a non-trivial change complete.
