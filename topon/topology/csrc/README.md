# C topology generator

`generator.c` is the C implementation of the strict-sculpting network
generator. It is the optional fast path for large runs; the pipeline
default is the pure-Python [`generator_python.py`](../generator_python.py),
which implements the same algorithm and needs no compiler.

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

## Keeping the two generators in step

Changes to lattice construction or to the `.nodes` / `.edges` format have
to land in both this file and `generator_python.py`.
`tests/unit/topology/test_c_generator.py` compiles this source and checks
site counts, coordinates and the header against the Python builder; it
skips when no compiler is on PATH.

The two generators do not produce identical draws. This one seeds from
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

Vendored 2026-08-05 from `generator_serial_debug11.c`
(md5 `83d7f9d37c72bcd26bb15fe664b67cc6`), the copy that appeared in five
archive locations including the most recently touched. `debug12` differed
only in comments. This file is now canonical; the archive copies are
historical.
