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

- `legacy/` and `older_versions/` — frozen reference code
- `tests/output/` — regression reference outputs (update only by running the generator with `--update-refs`)

---

## Testing

```bash
pytest tests/unit/         # fast, no LAMMPS needed (~5s)
pytest tests/regression/   # slow, generates full LAMMPS systems (~1.5h)
```

Unit tests live in `tests/unit/`. Regression tests live in `tests/regression/`. Integration / workflow scripts live in `tests/workflows/` (not pytest — run directly with Python).

---

## Docs

User-facing docs live in `docs/`. Keep them updated when changing behaviour:

- `docs/cli.md` — update when adding/changing CLI options
- `docs/config_reference.md` — update when adding config keys to `schema.py`
- `docs/simbox.md` — update when changing simbox API
- `docs/development/changelog.md` — add entry for every significant change
- `docs/development/tasks.md` — mark tasks complete as you finish them
