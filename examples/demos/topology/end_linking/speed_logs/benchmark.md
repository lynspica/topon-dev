# Topology generator benchmark

Comparison between the C generator (`generator.exe`, written by the original
project) and the pure-Python port (`topon.topology.generator_python`).
The Python port reproduces the C output bit-for-bit on success cases.

## Hard-case test

A 6x6x6 SC lattice with constraint `degree_distribution="0:13,1:25"`,
`max_functionality=4`. This is the test the V22 changelog refers to as
the "Hard Case".

| Generator | Solved? | Wall time |
|---|---|---|
| Python (`generator_python.py`) | yes | 0.23 s |
| C (`generator.exe`) | yes | TBD — re-run with the binary in `examples/demos/topology/end_linking/c/run.py` to fill this row |

## Batch benchmark — 20 cases from `network_candidates_SC_6x6x6_v2.txt`

15-second timeout per case.

| Generator | Success rate | Median time | Notes |
|---|---|---|---|
| Python (`generator_python.py`) | **60 % (12 / 20)** | ~0.3 s on solved cases | the 8 timeouts are constraint configurations the sculpting + systematic-search algorithm cannot solve in 15 s |
| C (`generator.exe`) | TBD | TBD | re-run benchmark needed |

## Reproducing

```bash
# Python
python examples/demos/topology/end_linking/python/run.py

# C — first set TOPON_GENERATOR_EXE to your compiled binary path
export TOPON_GENERATOR_EXE=/abs/path/to/generator.exe
python examples/demos/topology/end_linking/c/run.py
```

The two scripts use the same `run_generator(...)` entry point from
`topon.topology.generator`; the difference is the `exe_path` argument
(`None` for the Python port, a binary path for the C generator).

## Provenance

The Python port + benchmark were added in V22 (2026-01-20). See the V22
entry of [`docs/DEVELOPMENT.md`](../../../../../docs/DEVELOPMENT.md).
