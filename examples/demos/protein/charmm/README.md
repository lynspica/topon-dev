# Topro — CHARMM atomistic protein networks (legacy)

The original topro package emitted full-atomistic CHARMM (CHARMM36 / 36m + CMAP) protein networks: BFM SAW lattice topology, dityrosine crosslinks via `CE2-CE2` + DITY patch, TIP3P water + SHAKE, flat soft-push relaxation. The topology pipeline (BFM, gel-point detection, JSON snapshot) was forked verbatim into the current `topon.protein_network` (V36); the **chemistry stage** — atomistic CHARMM placement — was not.

## Why this folder is a stub

The CHARMM atomistic chemistry stage was not migrated into the integrated topon. The current protein-network capability under `topon.protein_network` produces MARTINI 3 coarse-grained output instead (see the sibling `martini/` folder). The reasons are documented in [`docs/ARCHITECTURE.md`](../../../../docs/ARCHITECTURE.md) §6 and the topon-family table in §1.

## Where the legacy code lives

The original topro package source lives under `legacy/subprojects/protein_network/topro/` (gitignored — owner-local). Three known issues from that code path informed the design of the integrated MARTINI port:

1. `pos % box` wrap breaks bonds at large `lattice_scale`.
2. Interpolation overwrites anchor positions with chain-walk values.
3. `prefactor = ramp(0, 30)` no-op during minimize.

These issues are why the integrated MARTINI port (`topon.protein_network`) uses a wrap-only convention with an xyz perturbation hack rather than the legacy approach.

## Reviving CHARMM atomistic

If a downstream project needs CHARMM atomistic protein networks, the cleanest path is to write a new `topon.protein_network.chemistry_charmm` module that produces atomistic placements from the same BFM topology snapshot used by the MARTINI builder. The topology, water packer, and ion packer are all reusable.
