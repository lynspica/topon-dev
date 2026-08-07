# Hard-core three-stage protocol (experimental)

An alternative to the scripts `topon/writers/lammps_inputs.py` generates.
**Not** a replacement: the generated scripts are unchanged and remain the
default everywhere. These are used only by `tests/workflows/entangle_steps.py`
under `--protocol hardcore`, for the entanglement work.

## What differs

Both changes are commented-out lines, not rewrites, so the diff against the
generated scripts is readable.

| stage | generated | here |
|---|---|---|
| 1 | `pair_style soft 1.0` ramping the prefactor 0 to 30 | soft push disabled; `lj/cut 1.122462` (WCA) throughout |
| 2 | `pair_style soft 1.0`, then LJ epsilon ramped 0.001 to 1.0 | both disabled; WCA throughout |
| 3 | unchanged | unchanged |

## Why

`pair_style soft` has finite energy at zero separation, so two beads may sit
on top of one another at bounded cost. That is what the stage-1 soft push is
for: it resolves a tangled starting geometry by letting chains pass through
each other. For a network carrying *prescribed* entanglements, that is the
one move that undoes the work.

WCA has a hard core. Nothing crosses it, so topology built into the starting
structure is still there at the end.

Measured on one prescribed pair, one winding, same build, same network:

| stage | generated (soft) | here (hard core) |
|---|---|---|
| as built | −1.13, sep 1.39 | −1.13, sep 1.39 |
| stage 1 | +0.32, sep 1.98 | **−1.13, sep 1.35** |
| stage 2 | −0.05, sep 15.12 | **−0.95, sep 0.94** |
| stage 3 | +0.09, sep 4.29 | **−0.92, sep 2.05** |

The generated protocol destroys the entanglement and leaves the partners
15 sigma apart. This one keeps it, same sign and magnitude, through all
three stages. Final bond statistics are identical either way
(0.968 median, 1.106 max), with no errors and no lost atoms.

## Known caveat

Without the soft push there is nothing to relieve a bad starting overlap
gently. Stage 1 here runs bonds at median 0.498 / max 1.697, so a bond
briefly passes the FENE limit of 1.5 before recovering by stage 2. It did
not error on this system. A denser or larger system is where that would
bite, and it is the first thing to check if one fails.

## Status

Experimental. Before this could become the default it needs to pass the
full suite, including the atomistic side-chain cases, which it has not been
run against.

Origin: hand-edited by the repository owner from the generated scripts of
`tests/output/entangle_steps/step2_band1_b90_wca`, and kept here so runs are
reproducible rather than depending on a local copy.

## Confirmed by primitive-path analysis

Z1+ (Kroger, Comput. Phys. Commun. 283 (2023) 108567) on the prescribed pair
with the rest of the network removed, so anything it reports is between
those two chains and nothing else:

| stage | hard core | generated (soft) |
|---|---|---|
| as built | Z = 1, 1 | Z = 1, 1 |
| after stage 1 | **1, 1** | **0, 0** |
| equilibrated | **1, 1** | **0, 0** |

Both protocols start with exactly the one entanglement that was designed.
The soft push removes it during stage 1. The hard core keeps it, and the
primitive path length barely moves (20.6 to 19.2).

Z1+ is not in this repository: its README says it should not be
re-distributed. Obtain it from the CPC library or mk@mat.ethz.ch. Only a
Linux binary ships for the core module, so on Windows it runs under WSL;
`tests/workflows/run_z1.sh` drives it over a directory of .Z1 files.
