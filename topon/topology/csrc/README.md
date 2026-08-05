# C topology generator

`generator.c` is the standalone searcher. It runs on its own, without
Python, and is the tool for long exhaustive searches over many trials.

The pure-Python [`generator_python.py`](../generator_python.py) is a
separate program with a different job: quick in-process generation of
likely networks, no compiler needed, and it is the pipeline default.

**These are two independent tools, not a library and a wrapper.** Nothing
here is called from Python and nothing here should grow a Python binding.
The division is deliberate: long searches belong to C because it is
faster at them, everyday generation belongs to Python because it is
immediate. What they share is the lattice construction and the
`.nodes`/`.edges` format, and only that shared surface has to stay in
step.

## Build

```bash
gcc -O2 -o generator.exe generator.c -lm
```

The binary is gitignored. Point topon at it with
`topology.generator.exe_path`, or leave that `null` to use the Python
generator.

## Command line

```
generator.exe <dims> <periodicity> <max_func> <max_trials> <max_saves> "<degree_dist>" <logging> <lattice_type>
```

| argument | meaning |
|---|---|
| `dims` | `NxxNyxNz`, e.g. `8x6x8` |
| `periodicity` | one digit per axis, `1` periodic and `0` open, e.g. `111` |
| `max_func` | maximum crosslink degree |
| `max_trials` | trials before giving up |
| `max_saves` | networks to write |
| `degree_dist` | `"d:N,..."` per-degree counts, or `"e:N"` for a total edge count |
| `logging` | `0` or `1` |
| `lattice_type` | `SC`, `BCC`, `FCC`, or `MIX:<sc>,<bcc>,<fcc>[,<cutoff>]` |

The mix fractions ride inside the `lattice_type` argument rather than
taking a ninth position, so scripts written against the original
eight-argument CLI keep working unchanged.

`MIX` builds its edges by an all-pairs neighbour search, which is O(N²) in
the site count. That is unnoticeable at the cell counts topon normally
uses (a 6x6x6 mixture is a few thousand sites) but grows quickly: expect
seconds around 20x20x20. The pure lattices are unaffected, since they
enumerate a fixed neighbour pattern instead.

```bash
./generator.exe 6x6x6 111 4 1000 1 "0:0,1:0" 0 SC
./generator.exe 6x6x6 111 4 1000 1 "0:0,1:0" 0 MIX:0.2,0.4,0.4
```

## Output

Writes `output/network_N<dims>_trial<n>.{nodes,edges}`. The `.nodes` file
opens with a `# BOX Lx Ly Lz` header recording the true periodic cell,
which `topon.topology.loader` reads back. That header matters: without it
the loader falls back to estimating the cell from the coordinate extent,
which is exact only for SC and overshoots any lattice with fractional
basis sites.

## What has to stay in step

Only the shared surface: lattice construction, and the `.nodes` /
`.edges` format. Changes there land in both this file and
`generator_python.py`. The sculpting search itself is each program's own
business.

`tests/unit/topology/test_c_generator.py` compiles this source and checks
site counts, coordinates, the `# BOX` header, and that the C sculpts the
same configurations Python does. It skips when no compiler is on PATH.

The two do not produce identical draws. This one seeds from
`srand(time(NULL))` and draws from `rand()`, so parity is over
distributions, not individual networks.

### Known divergences

Both predate the vendoring and are deliberately left alone:

- **Periodicity.** This file honours `p_dims` per axis; the Python
  builders always wrap, ignoring the `periodicity` config value.
- **The degree-2 guard.** In the sculpting stages it is gated on
  `is_sc_lattice` here, but applied unconditionally in Python, so the two
  sculpt BCC and FCC differently.

## Provenance

Vendored 2026-08-05 from `generator_serial_debug11.c`, md5
`e7631f4bbcb963d50c382721de3b3c18`, dated 2025-11-03. That is the version
the shipped `generator.exe` was built from and the one
`generator_python.py` was ported from.

### The variant that was NOT taken

A later `generator_serial_debug11.c` exists (md5 `83d7f9d37c72`,
2026-02-27) under `experiments/pruning_research/pruning_algorithm_math*`,
replicated to five archive locations. Editor history confirms it is the
newest source by timestamp, and it was vendored first on that basis. It
is wrong for this slot.

It replaces the per-degree count check in `is_move_safe` with a
cumulative one (`N_leq_v + increase > T_leq_v`). Measured across six
standard SC configurations it sculpts **1/6** where this version and the
Python port both do **6/6**, failing every case where `max_func` sits
below the lattice coordination, which is the ordinary use. Treat it as an
open experiment rather than a newer release. If the cumulative rule turns
out to be the correct one, the Python port has to move with it and the
sculpting failures need fixing first.

`test_c_sculpts_the_configs_python_sculpts` pins this: it fails on the
`83d7f9d3` variant and passes on this one.
