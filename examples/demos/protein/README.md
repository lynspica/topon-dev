# Protein-network demos (topro)

The protein-network capability in topon — user-facing name **topro** ("topological protein network") — currently has two flavours, captured in the sub-folders here:

- [`charmm/`](charmm/) — the **legacy** atomistic CHARMM protein-network workflow that the original topro package shipped. Not invoked from the current `topon` package; pointers to where the code lives.
- [`martini/`](martini/) — the **current** coarse-grained MARTINI 3 protein-network generator (`topon.protein_network`, V36). Working, regression-tested.

For the design rationale of why MARTINI 3 superseded the CHARMM path (and the approximations involved in porting from GROMACS to LAMMPS), see [`docs/ARCHITECTURE.md`](../../../docs/ARCHITECTURE.md) §6 and [`docs/USAGE.md`](../../../docs/USAGE.md) §4.1.
