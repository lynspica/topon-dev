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
| **V47** | 2026-08-05 | **Diamond added to C, per-axis periodicity added to Python** — the last two asymmetries. Both are deterministic, so parity is now exact: all 16 lattice x periodicity combinations produce identical edge sets, and Diamond matches node-for-node including ids. `lattice_type: "Diamond"` works through the config on both paths. Also fixed a C seeding bug that made **every run started in the same second produce the same network**. Found that open boundaries make `"0:0,1:0"` unsatisfiable on BCC and Diamond (free-surface sites have degree 1). |
| **V46** | 2026-08-05 | **C/Python parity swept across 24 configurations** (lattices, sizes, mixtures, distribution modes) by the new `tests/workflows/compare_generators.py`; 22 agree and the other 2 are targets both correctly refuse. Two C bugs found: a `strncmp` prefix match let `MIXED`/`MIXTURE` silently build a pure-SC lattice, and the completion check never read targets above `max_func`, so the generator printed "SUCCESS" over networks that did not satisfy the request. Python gained the matching fail-fast guard. |
| **V45** | 2026-08-05 | **Corrected the vendored C source and sped the Python generator up ~7x.** V44 vendored the newest C by timestamp, but that is an experimental variant that sculpts 1/6 standard configs where the 2025-11-03 version does 6/6; re-vendored the latter. Separately, 99.5% of the Python generator's runtime was a NetworkX subgraph *view* re-evaluating its node filter; a direct adjacency walk is 6-8x faster and provably yields identical networks. |
| **V44** | 2026-08-05 | **Mixed SC/BCC/FCC lattices** (`lattice_type: "MIX"`) in both generators, and the C generator vendored into `topon/topology/csrc/`. Mixing takes a lattice from one edge-length shell to four. `MIX` at `{"SC": 1}` reproduces `SC` exactly; the other corners deliberately do not match the canonical builders. 40 new tests, including C-vs-Python parity checks that compile the source on the fly. |
| **V43** | 2026-08-05 | **Generators record their periodic cell.** `infer_dims_from_graph` used to estimate the box as `max - min + 1`, exact for SC but half a cell too large for BCC/FCC/Diamond, which sent 33% of BCC edges (23% of FCC) to the wrong periodic replica at twice their true bond length. Generators now set `G.graph["box"]`; `.nodes` files carry a `# BOX` header. SC output unchanged. 20 unit tests. |
| **V42** | 2026-07-17 | **README rewritten** around the algorithm and arc animations; "Sub-systems" dropped, protein path (topro) given its own section with MARTINI 3 + CHARMM36m worked examples. Six false claims caught and fixed pre-commit. |
| **V41** | 2026-07-17 | **Entanglement + graft animations** — two-panel single-entanglement (`ent_arc`, `make_ent_movie.py`) and a graft showcase (`graft_arc`, new `grafted` system). All arc animations slowed to ~0.66×. |
| **V40** | 2026-07-17 | `topology/generator_python.py` (+ `generator_python_diamond.py`) — fail fast on unreachable `degree_distribution` targets. An `e:N` above the lattice's edge count (e.g. `e:128` on an 81-edge 3×3×3 SC), or a per-degree target above the node count / max degree, now raises a clear `ValueError` instead of churning through doomed trials. 11 unit tests added. |
| **V39** | 2026-07-17 | Arc **animations** for the gallery — real MD trajectories (lattice → minimised → equilibrated) as boomerang loops, one builder `make_arc_movie.py`. |
| **V38** | 2026-07-16 | Showcase **gallery** (`assets/gallery/`, eleven panels on strict-sculpted heterogeneous networks) + README rewrite. |
| **V37** | 2026-05-20 | In-situ crosslinking (`bfm.generate_topology(crosslink_method="none")`) + opt-in physically-correct CHARMM builds (`--physical-backbone`). NPZ node-feature schema v2. Five CHARMM36m parameter-injection fixes. |

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

## 4. Changelog (V1 – V47)

Notable changes are documented in reverse chronological order.

### [V47] — 2026-08-05 — Diamond in C, per-axis periodicity in Python

#### Added
- **`csrc/generator.c`** — `create_diamond_lattice` plus a `Diamond` /
  `DIAMOND` dispatch. Eight sites per cubic cell, 4-coordinated by
  construction, so a `max_func=4` network needs no pruning. Walks the
  basis in the same order as `generator_python_diamond.py`, so the two
  number their nodes identically.
- **`generator_python.py`** — reads `periodicity` (previously ignored
  entirely) and honours it per axis in the SC, BCC, FCC and MIX
  builders. `generator_python_diamond.create_diamond_lattice` takes the
  same argument. Accepts the C's `"110"` digit string, a bool, or any
  3-element iterable; an unparseable value falls back to fully periodic
  with a warning, since guessing "open" would silently introduce free
  surfaces.
- **`lattice_type: "Diamond"`** in the config schema and in
  `_create_lattice`, which dispatches to the standalone module rather
  than inlining it — the Diamond logic stays reviewable in isolation.
  Also offered by `topon init --interactive` (MIX is not: it needs
  fractions, which is more than that prompt should ask).
- **29 tests** in `tests/unit/topology/test_periodicity.py`, plus 18
  exact-parity cases in `test_c_generator.py`.

#### Verified
Both features are deterministic on both sides, so parity is exact rather
than statistical: across **four lattice types x four periodicity
settings, all 16 produce identical edge sets**. Diamond additionally
matches node-for-node including ids.

#### Behaviour change
A config setting a non-default `periodicity` on the Python path used to
be silently ignored and produce a fully periodic lattice; it now builds
the slab it asked for. Configs at the `"111"` default are unaffected,
which `test_default_is_fully_periodic_and_unchanged` pins.

#### Fixed
- **`csrc/generator.c` seeded only from `time(NULL)`**, which advances
  once a second, so **every C run started within the same second
  produced byte-identical output**. A script looping the executable to
  collect N networks was silently collecting N copies of one network.
  Caught by the comparison sweep, where the "independent" reps of a
  mixture were not independent and skewed the site-count statistics. The
  pid is now mixed in, and `TOPON_SEED` overrides for a reproducible run
  (the Python equivalent is seeding `random` before calling `generate`).

#### Found while testing
- Opening an axis puts corner sites on a free surface with very low
  coordination. On 4x4x4 with every axis open the minimum degree is 3 for
  SC and FCC but **1 for BCC (2 sites) and Diamond (22 sites)**, which
  makes the usual `"0:0,1:0"` unsatisfiable there — clearing a degree-1
  site means cutting its last bond, producing the degree-0 node the same
  request forbids. Both generators decline. Documented in `USAGE.md` with
  the per-lattice table.
- **The C searcher has no wall-clock limit**, so a structurally-doomed
  request looks like a hang. Not an infinite loop — each trial does
  terminate — but stage 4's systematic search scales badly on a graph it
  can never satisfy: for `"0:0,1:0"` on an open Diamond, 2x2x2 is
  instant, 3x3x3 exceeds 90 s for three trials, and a **single 4x4x4
  trial did not finish in 300 s**, where Python gives up in ~0.25 s per
  trial via its `time_limit`. Pre-existing, newly reachable. Two fixes
  suggest themselves and neither is applied here: a pre-flight guard
  comparing requested `d0`/`d1` counts against the sites the boundary
  forces below that degree, and a wall-clock cap matching Python's.

### [V46] — 2026-08-05 — C/Python parity sweep; two C bugs

#### Added
- **`tests/workflows/compare_generators.py`** — runs both generators over
  24 configurations (SC/BCC/FCC/MIX, four sizes including a non-cubic
  3x4x5, seven mixtures, per-degree and `e:N` distribution modes) and
  compares site counts, mean degree, edge-length shells and the recorded
  box, then pushes a subset through the pipeline to a LAMMPS stage-1
  minimize. **22 of 24 agree**; the other two are targets neither
  generator can satisfy and both now refuse. All six LAMMPS pathways
  complete: SC, BCC and FCC pruned to `max_func=4`, a 0.2/0.4/0.4
  mixture, a non-cubic 3x4x5, and an `e:200` edge-count target.

#### Fixed
- **`csrc/generator.c` — prefix match on the MIX dispatch.**
  `strncmp(lattice_type, "MIX", 3) == 0` also matched `MIXED`, `MIXTURE`
  and any other typo starting with those letters, and quietly built a
  pure-SC lattice instead of erroring. Now requires `:` or end-of-string
  after `MIX`.
- **`csrc/generator.c` — targets above `max_func` were never checked.**
  The stage-4 completion check loops only to `max_func`, but the parser
  accepts higher degrees, so `7:5` on a `max_func=4` run was stored and
  then ignored: the generator reported *"SUCCESS: Target distribution
  met!"* over a network with no degree-7 nodes at all. Same for `6:100`,
  where the lattice offers degree 6 but sculpting cannot leave it there.
  No node can finish above `max_func`, so such targets are now rejected
  up front.
- **`generator_python.py`** gained the matching guard, ordered after the
  existing lattice-coordination check so the more fundamental reason
  still produces the error message. Previously these requests churned
  through every trial before giving up.

- **`pyproject.toml`** — `package-data` listed only `data/*`, so an
  installed topon would have shipped no C source at all: `topology/csrc/`
  has no `__init__.py` and is therefore not found as a package. Now
  includes `topology/csrc/*.c` and `*.md`.

#### Documented, not changed
- `is_sc_lattice` is `strcmp(lattice_type, "SC") == 0`, so the degree-2
  sculpting guard is off for `MIX:1,0,0` even though it builds the
  identical SC lattice. On 5x5x5 pruning to `max_func=4`, `SC` averages
  221 edges against `MIX:1,0,0`'s 216 — both valid, drawn from slightly
  different distributions.
- Diamond exists only in Python; per-axis periodicity only in C.

### [V45] — 2026-08-05 — Right C source; Python generator ~7x faster

#### Fixed
- **`topon/topology/csrc/generator.c`** re-vendored from md5 `e7631f4b`
  (2025-11-03) in place of `83d7f9d3` (2026-02-27). The latter is the
  newest file by timestamp, confirmed by editor history, but it is an
  experimental variant under
  `experiments/pruning_research/pruning_algorithm_math*`: it replaces the
  per-degree count check in `is_move_safe` with a cumulative one and
  sculpts **1/6** standard SC configurations where the vendored version
  and the Python port both do **6/6**, failing whenever `max_func` is
  below the lattice coordination. Three signals were missed first time
  round: the shipped `generator.exe` compiles from the 2025-11-03 source,
  the Python port implements the per-degree rule, and the directory name
  marks it as research. `test_c_sculpts_the_configs_python_sculpts` now
  fails on that variant.

#### Changed
- **`generator_python.py:_is_subgraph_connected`** walks `g._adj`
  directly instead of building `g.subgraph(...)` and calling
  `nx.is_connected`. A NetworkX subgraph is a view that re-evaluates its
  node filter on every neighbour access: 6.0M `new_node_ok` calls for one
  1000-node lattice. That call was 99.5% of the generator's runtime, so
  the whole generator is now 6-8x faster (12³: 10.95 s to 1.54 s).
  Verified identical: 1200 fuzzed states against the NetworkX oracle with
  zero disagreements, and identical edge sets across four sizes and three
  seeds with the two implementations swapped. Both are now tests.

#### Clarified
- The C generator and the Python generator are **independent programs
  with different jobs**, documented in `csrc/README.md` and
  `ARCHITECTURE.md`. C is the standalone searcher for long exhaustive
  runs; Python is the quick in-process path. The C is not called from
  Python and must not grow a Python binding. Only lattice construction
  and the `.nodes`/`.edges` format are shared surface. Benchmarked: at 6³
  Python wins on wall clock (the C pays process startup), by 12³ the C is
  14x ahead, and at 24³ it finishes in 6.7 s where Python would run for
  hours.

### [V44] — 2026-08-05 — Mixed SC/BCC/FCC lattices; C generator vendored

#### Added
- **`lattice_type: "MIX"`** in both generators. All three cubic lattices
  share the cell corner and each adds sites on top of it, so `MIX` places
  the corner in every cell, the BCC body centre with probability
  `mix_fractions["BCC"]` and each of the three FCC face centres with
  probability `mix_fractions["FCC"]`. The `"SC"` entry is the remainder and
  places no site of its own, which is what makes the three a partition
  summing to 1. Expected site count `Nx*Ny*Nz * (1 + f_bcc + 3*f_fcc)`
  recovers N / 2N / 4N at the pure corners.
- Edges join every pair within `mix_cutoff` (default 1.0 cell, the
  simple-cubic nearest-neighbour distance) under the minimum image, since a
  mixed point set has no single neighbour shell to enumerate.
- **`topon/topology/csrc/`** — the C generator vendored from
  `generator_serial_debug11.c` (md5 `83d7f9d3`, the copy present in five
  archive locations including the most recent). It gained `MIX` and the
  `# BOX` header from V43. The fractions ride inside the existing
  `lattice_type` argument as `MIX:<sc>,<bcc>,<fcc>[,<cutoff>]`, so the
  eight-positional-argument CLI every existing caller and SLURM script
  writes is unchanged.
- **40 tests**: `test_mixed_lattice.py` (26) and `test_c_generator.py` (14).
  The C tests compile the source on the fly and skip without a compiler.

#### Fixed
- **`topon/topology/csrc/generator.c`** — `target_counts` was sized from a
  per-lattice constant (`max(max_func, 12)`) but is indexed by node degree
  in `run_single_trial`. The pure lattices top out at 12 so it exactly
  fit, but a mixture reaches degree 20 and `mix_cutoff` is user-settable,
  so no constant is safe. It is now sized from the maximum degree of the
  graph actually built, which required moving the lattice construction
  ahead of the argument parsing.
- `srand` moved above the lattice build so `MIX` can draw its sites. A
  no-op for SC/BCC/FCC, which consume no randomness while being built.
- **Heap overflow in `create_mixed_lattice`.** The site scratch buffer was
  sized at `4 * Ncells`, the FCC site count, but a cell can hold five
  sites (corner + body + 3 faces). It looks safe on the mean and is not:
  at `f_bcc=0.1, f_fcc=0.9` roughly 0.2% of 4x4x4 draws exceed `4*N`.
  Now sized at `5 * Ncells`, with a repeated-draw regression test on that
  exact mixture.

#### Known, deliberate
- `MIX` at `{"SC": 1}` reproduces `SC` exactly, down to node ids. It does
  **not** at the other corners: the cutoff also admits the corner-corner
  shell, giving 14 neighbours at `{"BCC": 1}` against canonical BCC's 8,
  and 18 at `{"FCC": 1}` against FCC's 12. The site sets do match. Pinned
  by `test_mix_bcc_corner_is_not_canonical_bcc`.
- A body centre and a face centre can sit 0.5 cells apart, half the SC
  spacing, so at fixed DP the bond length spreads by up to 2x.
- Two pre-existing C/Python divergences are documented in
  `topon/topology/csrc/README.md` and left alone: the C honours per-axis
  periodicity while Python always wraps, and the sculpting degree-2 guard
  is SC-gated in C but unconditional in Python.

### [V43] — 2026-08-05 — Generators record their periodic cell

#### Fixed
- **`topon/topology/loader.py`** — `infer_dims_from_graph` now returns the
  cell recorded by the generator in `G.graph["box"]`, falling back to the
  old `max - min + 1` estimate only when no box is present. The estimate is
  exact for SC, whose sites are integer-spaced with unit separation, but
  overshoots every lattice with fractional basis sites: BCC and FCC body and
  face sites sit at +0.5, Diamond sites at quarter cells, so the coordinates
  stop short of the cell edge and a 4x4x4 BCC or FCC was reported as 4.5.
- Because that value is the box for every minimum-image calculation in the
  pipeline, the overshoot broke the `bond < box/2` invariant of Design
  Principle 3 and pushed edges onto the wrong periodic replica. Measured on
  4x4x4: **169 of 512 BCC edges (33%) and 360 of 1536 FCC edges (23%)** were
  built at exactly twice their true bond length (BCC 1.732 vs 0.866, FCC
  1.414 vs 0.707). Diamond was affected identically.
- A box recorded on the graph now also overrides a `dims` stored alongside
  an older gpickle, which may have come from the estimate.
- **`topon/conformation/manager.py`** — the box was estimated in a *second*,
  independent place: stage 5 re-derived it from the `.displace` files as
  `(max node coord + 1) * scale`. `apply_displacements` gained a
  `lattice_box` argument, and `pipeline.py` plus both workflows now pass the
  same `dims` stage 4 routed the chains with. Fixing only the topology side
  left the two disagreeing (routing on period 4.0, box written at 4.5) and
  produced **63 bonds up to 14.5 Å** where the unfixed code had none over
  1.85 Å. As a side effect the box is now `dims * scale = volume^(1/3)`
  whatever `dims` holds, so the target density is hit exactly on every
  lattice instead of only on SC.

#### Added
- **`topon/topology/generator_python.py`**, **`generator_python_diamond.py`** —
  SC, BCC, FCC and Diamond builders set `G.graph["box"]` to the true lattice
  repeat. The attribute survives `Graph.copy()` (sculpting), the
  `nx.MultiGraph` conversion, vacancy removal and pickling.
- **`.nodes` files** gained an optional `# BOX Lx Ly Lz` header, with
  `read_box_header` / `format_box_header` to parse and render it and a
  `save_nodes_edges` writer that emits it. Files without the header load
  exactly as before.
- **`tests/unit/topology/test_lattice_box.py`** — 20 tests covering the
  recorded cell per lattice, uniform nearest-neighbour edge lengths,
  survival through the pipeline transforms, `.nodes` round-trip, malformed
  and absent headers, and the gpickle precedence rule. One test pins the
  old fallback behaviour so a future change to it is visible rather than
  silent.
- **`tests/workflows/verify_lattice_box.py`** — end-to-end script that
  audits all four lattices, then runs the same BCC network through the
  pipeline twice (with and without the header) and LAMMPS stage-1 minimize.

#### Compatibility
SC output is unchanged by construction, since both estimates were already
exact there: on the frozen 5x5x5 sample graph the old `max + 1` and the new
`dims` are both 5.0. Verified by running `pytest tests/regression/` with the
change stashed and again with it applied: the per-test status diff is empty
(43 passed, 8 failed, 14 skipped, 3 errors in both, the failures
pre-existing). The two reference-generating scripts under `tests/workflows/`
were left passing no `lattice_box` on purpose, since they regenerate SC
references where the fallback is exact. The C generator does not yet emit
the header, so its `.nodes` output still uses the fallback for BCC/FCC.

### [V42] — 2026-07-17 — README rewritten around the animations

#### Changed
- **`README.md`** rewritten: plainer prose, led by the algorithm and arc
  animations rather than a feature list. Dropped the "Sub-systems" section
  (simbox, singlechain) — the protein path is presented as topro, in its own
  section, with MARTINI 3 and CHARMM36m worked examples that take a sequence and
  repeat count and produce a solvated, salted network.

#### Fixed (claims an investigator pass caught in the first draft)
- "POSS **junctions**" → POSS **chain caps**. `diagnostics/rules.py:37` records
  POSS at degree ≥2 as known bug P1-H; the README was advertising the bug.
- `--water-density 0.9` → `10`. The flag is beads/nm³ (`water.py:25`, bulk ≈ 10),
  so 0.9 was ~9% of bulk — the command promised a solvated system and delivered a
  nearly dry one.
- "historically called topro" → topro is the *current* name (`topon topro` is a
  live command, and ARCHITECTURE/AGENTS both use it).
- The CHARMM example could not be copy-pasted: `--topology` is required and must
  already exist. Now shown as two commands, with `--output`.
- "gel point at 0.55" → 0.55 is where the *last* chain joins the cluster; the
  giant component passes half the chains at ≈0.14. Both are now stated.
- "about five times too short" scoped to CG (all-atom is 2.5×), and the "(BFM)"
  initialism dropped — `bfm.py` is a 6-neighbour cubic lattice with one monomer
  per site, not the Carmesin–Kremer bond-fluctuation model.
- **`docs/USAGE.md`** — added the `gradient` caveat, which the README states and
  USAGE previously contradicted.

### [V41] — 2026-07-17 — Entanglement + graft animations, slower arcs

#### Added
- **`assets/gallery/anim/graft_arc.{gif,mp4}`** — a graft showcase arc: a new
  `grafted` system (4×4×4 sculpted to 128 edges, 292 side chains of DP 6 ≈ 40% of
  the beads, no entanglements), backbone blue, side chains teal, junctions dark,
  relaxing lattice → melt. Added a `graft` paint mode + config to
  `make_arc_movie.py` (which now covers cg / copoly / atom / graft; the
  single-entanglement animation stays in its own `make_ent_movie.py`).
- **`assets/gallery/anim/ent_arc.{gif,mp4}`** + `make_ent_movie.py`,
  `movie_ent.in` — a **single entanglement** shown two ways at once, both playing
  the lattice → minimised → equilibrated arc: LEFT the full sculpted network with
  that one entanglement highlighted and a locator box around it; RIGHT that box
  zoomed in. **Exactly three colours**: chain A gold, chain B violet — each
  *including its own side chains* (gold/violet branches, not a fourth colour) —
  everything else faint grey. The boomerang **holds** on the lattice and the melt
  so the "before"/"after" register. `movie_ent.in` reads the pristine
  `03_Conformation` lattice (not `1.restart`, whose soft push has already loosened
  the 0.39 σ crossing to ~1.0 σ) so frame 0 is the tight entanglement.

#### Changed
- **All four arc GIFs/MP4s slowed to ~0.66× speed** (GIF 60 → 91 ms/frame, MP4
  20 → 13 fps), rebuilt from cached frames.
- The README entanglement section leads with the animation; the two stills move
  below it as the as-built reference.

#### Notes for the record — five things the obvious version gets wrong
- **Colour discipline.** Colouring every graft teal, or every entangled strand
  gold, made one entanglement read as many. Only the two focus chains and their
  own grafts carry colour; side chains inherit the parent chain's colour so it
  stays *three* colours and the side groups are still findable (`chain_grafts`
  returns per-chain graft sets).
- **The pair sits on a periodic face**, so a wrapped centroid flips between box
  faces frame to frame. Both panels *follow* the pair (recentre its centroid to
  the box centre each frame); the locator's centre is then fixed.
- **The zoom is adaptive.** The pair grows ~2.5 σ → ~8 σ from crossing to melt
  coil, so a fixed crop shrinks the lattice to a dot or clips the melt. The crop
  radius and camera FOV track the pair's own 88th-percentile extent per frame
  (lightly smoothed), so it fills the panel throughout; the locator box on the
  full panel grows to match (`project()` reads the camera back after `zoom_all`
  and gives pixels-per-σ).
- **A dense KG melt buries any single chain**, so the matrix is a faint thin web
  and the pair a bold string of large beads.
- **Holds at both ends** turn the arc from a blur into a legible before → morph →
  after (identical hold frames coalesce under GIF optimisation, so they cost
  duration, not bytes).

### [V40] — 2026-07-17 — Fail fast on unreachable topology targets

#### Fixed
- **`topology/generator_python.py` + `topology/generator_python_diamond.py` — an unreachable `degree_distribution` now raises instead of hanging.** Both generators validate the requested target against the freshly-built lattice *before* the trial loop (new `_validate_targets_reachable`). Sculpting only ever *removes* edges, so the full lattice is a hard ceiling: an over-target request such as `e:128` on a 3×3×3 SC lattice (81 edges) used to churn through hundreds of thousands of doomed trials with no output — indistinguishable from a hang. It now raises a `ValueError` naming both numbers (e.g. `"degree_distribution e:128 exceeds the 81 edges of a 3x3x3 SC lattice; sculpting only removes edges…"`). The same guard rejects per-degree targets that exceed the node count (`d:N` with `N > nodes`) or ask for a degree above the lattice's maximum (`d > max degree`). Bounds are read from the constructed graph, not a `3·nx·ny·nz` formula, so periodic-boundary edge collapse on tiny lattices (a 2×2×2 SC has 12 edges, not 24) is handled correctly. Reachable targets (`e:81`, `e:70`, `0:2`, empty) are unchanged. 11 unit tests added under `tests/unit/topology/` (the diamond generator previously had none). The archived reference C generator (`generator_serial_debug11.c`, `~/topon_archive/`) and the gitignored experimental `internal/…/generator_serial_diamond.c` share the pattern but are out of the tracked tree and untouched.

### [V39] — 2026-07-17 — Arc animations + entanglement colour rule

#### Added
- **`assets/gallery/anim/{cg,copoly,atom,ent}_arc.{gif,mp4}`** — real MD
  trajectories (lattice → minimised → equilibrated) as seamless boomerang loops,
  ~0.66× speed. The README arc section shows the GIFs instead of the frozen
  stills. (An initial `movie_render.py` + `make_gif.py` pair was written and then
  consolidated the same day into the two builders below; those two scripts no
  longer exist.)
- **Smooth, high-quality arc animations for all three resolutions**
  (`movie_cg.in`, `movie_atom.in` + one builder `make_arc_movie.py`). The
  lattice→melt transition completes inside a single `minimize` under the
  production decks, so a fixed-interval dump jump-cuts it. The movie decks inflate
  the lattice under a ramped soft push with bounded dynamics (`nve/limit`), dumped
  every step; `make_arc_movie.py` selects ~35 frames evenly in a monotone progress
  metric — **bond-orientation disorder** for CG/copolymer (rises through the
  inflation *and* the coiling, unlike bond length, which saturates), **mean
  minimum-image displacement** for atomistic (whose methyl-H bonds give the
  disorder metric a nonzero baseline). Largest per-frame visual gap 0.011 (CG) /
  0.018 (copolymer) / 0.14 (atom) vs a 0.83 cliff; 660 px render, 16 AO samples.
  The README arc row shows the CG + atomistic GIFs; the copolymer section gains a
  block-copolymer arc GIF (A/B halves preserved into the melt).

#### Changed
- **Entanglement panels → a 3-colour rule** (per user): one entanglement's two
  chains in two distinct colours, the whole rest of the network in a single third
  colour, in both the full and zoom views. The CG arc is now single-colour too
  (junctions by size). Topological colouring is reserved for the entanglement
  case; chemistry colouring (copolymer A/B, atomistic elements) is kept.

#### Found while animating
- **Atomistic stage 1 is NOT a no-op** (unlike CG): its staged soft-pushes expand
  the bonds 0.44 → 1.09 Å. Dumping only stages 2–3 gave an atomistic movie with no
  lattice; the fix dumps stage 1 and prepends the pristine `03_Conformation`
  frame. `reset_timestep` also aborts with an active dump, so the instrumenter
  now `undump`s before it.

### [V38] — 2026-07-16 — Gallery + README rewrite

#### Added
- **`assets/gallery/`** — eleven panels, all on **strict-sculpted heterogeneous
  networks** (junction functionality 2–6, mean 4.0, via `degree_distribution =
  "e:N"`), not perfect grids, plus vendored `gen_systems.py` + `render_gallery.py`
  and a README recording provenance and every caption number. Three sections: the
  arc at both resolutions (CG `sculpt_250` **0.17 → 0.96 → 0.97 σ**; atomistic
  `atom_sculpt` **0.44 → 1.09 → 1.12 Å**); the copolymer sequences (sculpted
  4×4×4); entanglements + side chains (sculpted 5×5×5: 12 requested → 24 kinked
  strands, 259 grafts, closest pair 0.39 σ). Both arcs are generated fresh and run
  through LAMMPS; the copolymer and entanglement panels are lattice-state only. No
  `tests/output/` golden is used, and no single-chain panel (not network
  generation).
- **`render_gallery.py` helpers**, each earned by a wrong picture:
  `element_palette()` (colours keyed off the file's own mass table — type IDs are
  *not* portable across systems), `node_ids()` (junctions from topon's own
  `system.groups`, because in a copolymer the junctions and the A monomers are
  both type 1), `rebuild_bond_pbc()` (recompute per-bond PBC shifts after any
  wrap), `half_cell_offset()` (see below), and `walk_chains()` (strip grafts,
  then walk — it returns exactly the 375 edges of a 5×5×5 SC lattice).

#### Found by rendering — all three were marked VERIFIED by a code-read audit
- **Every CG `minimize` is a no-op.** `writers/lammps_inputs.py` hardcodes
  `etol=1.0e-4`; LAMMPS tests |ΔE|/(|E₁|+|E₂|), which at the CG system's E ≈ 3.9e6
  already scores ~4e-7, so each `minimize` stops after **1 iteration**. Measured:
  `system_after_soft` has moved **2e-4 σ** max from the as-built lattice — even the
  soft push-off does nothing. All of the 0.198 → 0.963 relaxation is the
  `run 20000` NVE ramp in stage 2 (ending at T = 18.4). The atomistic writer uses
  `etol=1e-8` and genuinely minimises. Found because the CG "after minimisation"
  panel was indistinguishable from the as-built one.
- **`arrangement: "gradient"` ignores the requested composition.**
  `sequences.py`'s `w = max(0, 1 - |t - pivot| * n)` scales the window by `n`, so
  overlap needs n > 2: with two monomers nothing ever blends and *every*
  composition comes out 50:50 (ask A=0.1, get A=0.50 over 200 seeds). At equal
  fractions and even DP that is byte-identical to `block`. No gradient panel; the
  README no longer claims it works.
- **20% of bonds rendered as stubs**, because the as-built lattice puts nodes at
  x=0 — exactly *on* the periodic face — and the conformation jitter then sends
  neighbouring beads to opposite sides on wrap. 2 761 of 7 625 beads sat within
  0.01 of a face; 1 596 of 7 875 bonds spanned the box. `half_cell_offset()` shifts
  the planes into the interior, leaving the **75** genuinely periodic bonds a
  5×5×5 lattice must have (25 crossing edges per axis × 3). This is what made the
  first draft of the panels look "dotted" — no bead radius could have fixed it.

#### Noted, not established
- **The `v21_cg_combined` golden looks stale** — 190 bonds at 7.45 σ (=√3 × its
  4.2996 spacing) and 10 at 10.53 σ, bead-ends pinned to lattice nodes; fresh
  output shows nothing like it. But the golden *loads* a 210-edge sculpted graph
  rather than generating a 375-edge lattice, and its config
  (`examples/config_cg_combined.json`) no longer exists, so it cannot be re-run
  for a controlled comparison. The solid evidence is narrower: the bead-spacing
  formula moved `a/dp` → `a/(dp+1)`. An earlier draft tied this to
  `test_reproducible_with_seed` failing — wrong: that test dies on the missing
  config, and only ever compared atom counts between two fresh runs.

#### Changed
- **`README.md`** rewritten around the gallery: honest feature list, the
  **CHARMM36m atomistic protein builder** (previously unmentioned), BFM lattice +
  gel-point union-find, `fix bond/react` crosslinking, the pure-Python generator
  (no compiler required), GraphML/NPZ export, test markers, and the real CLI
  surface.

#### Fixed
- **Clone URL** — `README.md` pointed at `https://github.com/lynspica/topon.git`,
  which does not exist. The remotes are `lynspica/topon-dev` (personal) and
  `keten-group/topon` (stable); `pyproject.toml`'s Homepage already said
  `topon-dev`.
- **`assets/logo/README.md`** — corrected a claim introduced in V37 that no
  pipeline code realises kink geometry. It does (`pipeline.py:489`, `:589`). The
  real, narrower defect is now recorded there: `params` / `num_entanglements` are
  never passed to `calculate_entangled_kink`, so `KinkParams` is silently ignored
  (schema defaults coincide with the hardcoded ones, masking it) and
  multi-entanglement geometry is never realised.

#### Known-unsupported claims removed from the README
Each is a real gap, not just wording. Two `investigator` passes; the second found
that the first draft of this very rewrite had introduced five new false claims.
- **Dityrosine `fix bond/react` is not in this repo** — no `*dity*` file is
  tracked; the only `bond/react` templates are simbox's epoxy-amine. topon places
  dityrosine as build-time harmonic SC4–SC4 bonds (`protein_network/builder.py:489`).
- **`bfm.py` is not a bond-fluctuation model** — 6-neighbour cubic lattice, one
  monomer per site, fixed bond length 1 (`bfm.py:1`, `:105-115`, `:131-157`).
  The name and the three `examples/demos/topology/bfm/` echoes are misleading.
- **Kinks are applied in the chemistry stage**, not conformation (`pipeline.py:489`,
  `:589` are inside `_run_chemistry_stage`) — which also contradicts CLAUDE.md's
  "chemistry does NOT generate coordinates".
- `--physical-backbone` is **opt-in, default OFF** (`build_systems.py:115`).
- `lattice_type` accepts only `SC`/`BCC`/`FCC`; `generator_python_diamond.py` is
  not reachable from a config. BCC/FCC are untested.
- `.graphml` / `.npz` load only through the Python API (`run_from_graph`), not
  `ExistingFilesConfig`.
- `defects.secondary_loops` is schema-validated and gated but never injected;
  and what `assignment/defects.py` calls a *primary* loop is a *secondary* loop
  in the literature — and is the same object `analysis/report.py` calls
  secondary.
- Entanglement "multi" support is metadata only (see Fixed, above).
- "Dynamic DP scaling" — DP is static config; what scales is graft geometric
  extension. Atomistic grafts are additionally PDMS-only
  (`chemistry/builder.py:556`).
- `pytest -m fast` is ~13–22 s, not the ~5 s claimed in `pyproject.toml:79`,
  `tests/conftest.py:9` and CLAUDE.md.

#### Not fixed here — needs a decision
- **The regression tier is red**, and was already red at `fa86928` before this
  docs-only work: 9 failed / 3 errors / 14 skipped. `poss_100
  test_system_data_identical` mismatches all 10 470 atoms. CLAUDE.md's writer
  rule ("run `pytest tests/regression/` — confirm passing") cannot currently be
  satisfied.
- **`tests/output/` is gitignored** — goldens are local-only, so a clone skips
  the regression tier entirely.

### [V37] — 2026-05-20 — In-situ crosslinking + physically correct CHARMM builds

#### Added
- **`bfm.generate_topology(crosslink_method="none")`** — emits a single conv=0 snapshot labelled `uncrosslinked` (`reactions=[]`), skipping the crosslink loop and candidate search. The starting point for **in-situ** crosslinking: form dityrosine bonds *during* MD (LAMMPS `fix bond/react`) rather than a priori at build time. No CHARMM-builder change was needed — the builder infers a crosslink only where two Y nodes share a lattice site, and BFM's excluded volume guarantees none do on an unreacted snapshot. Existing methods and snapshot labels are untouched.
- **`charmm/build_systems.py --physical-backbone` / `--xpro-cis-fraction`** (opt-in) — build physically correct starting geometry:
  - `charmm_ff.py` parses the RTF internal-coordinate (IC) tables (370 entries) into `residues[..]["ics"]`.
  - `builder._build_physical_positions` lays backbone N/CA/C along the (coiled) CA trace at real bond lengths + ~111° N–CA–C, then NeRF-builds every remaining atom from the residue's IC table → ideal bonds/angles, planar impropers, real sidechain rotamers. Sidechains are reflected across the N–CA–C plane where the lattice's hairpins would flip chirality → **100 % L**.
  - `builder._coil_positions` decompresses interior residues to ~3.8 Å CA–CA (crosslinker Y residues stay on their lattice nodes, so a-priori crosslink geometry is preserved), so stage-1 needs no violent expansion.
  - `lammps_writer.find_omega_dihedrals` / `find_chirality_impropers` + a stage-1..3 `fix restrain` side-include (`*.in.omega`) holding peptide omega trans (minus an `--xpro-cis-fraction` X-Pro subset) and the N–C–CA–CB chirality improper at L, released before each stage's dynamics. CHARMM has no CA-chirality improper, so without this the soft-min inverts ~20 % of centres to D.
- **`charmm/build_systems.py --no-image-flags`** — legacy 7-column Atoms, keeping all crosslinks (exact topology, single-rank only). Default emits 10-column image flags via a priority-weighted MST and drops only winding-cycle crosslinks (MPI-safe).
- **NPZ node-feature schema v2** (`writers/npz_writer.py`, `SCHEMA_VERSION = 2`) — 8 → 10 columns: `node_degree` split into `chem_degree` / `phys_degree`, new `frac_ext`, and the conformation-derived columns (`contour_length`, `rg`, `COMX/Y/Z`, `frac_ext`) are NaN at write time (they are filled in by the downstream LAMMPS run, not the topology generator).

#### Fixed
- **CHARMM36m `(DEFAULT)` parameter injection** (`charmm_ff.py`, `charmm/builder.py`, `charmm/lammps_writer.py`) — five force-field-correctness bugs that silently substituted generic parameters: dihedral wildcard fallback (`X t2 t3 X`), bidirectional improper lookup (central atom first *or* last), N-terminal proline patch (`PROP`, not `NTER`), multi-term proper dihedrals (101/575 keys were truncated to their first Fourier term), and HIS → HSD remap (bare `HIS` is not in the RTF and was silently dropped, fusing its neighbours). After the fixes both resilin sequences regenerate with **0 DEFAULTs and net charge exactly 0.0000**.
- **`writers/npz_writer.py` `edge_index`** — remap from the original sparse/offset ID space to 0-based row positions into `node_features` (PyG requirement; raw IDs exceeded N and triggered out-of-bounds).
- **`topology/loader.py` `load_npz`** — reconstruct the original simple-cubic lattice coordinates from `node_ids` + `box` when the v2 COM columns are NaN (validated: every chain must join two lattice neighbours under PBC; returns NaN + a warning for non-SC graphs rather than a silently invalid build).
- **`conformation/manager.py`** — initialise `moved_count` before the loop (`max_iters=0` no longer raises `UnboundLocalError`).

#### Notes
- The `--physical-backbone` path is **opt-in**; the default CHARMM build is verified byte-for-byte identical to the previous builder (data + stage scripts).
- Stage 1 is serial by design (`pair_style soft 1.0` → ~3 Å comm cutoff < ~5.5 Å longest bond); stages 2/3 are MPI-safe. See [USAGE.md](USAGE.md) §4.1.

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
