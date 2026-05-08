# Topon — Version History

Maps canonical version tags to their reference output directories and the workflow scripts that generated them.

For detailed per-version changelogs see `docs/development/changelog.md`.

---

## Network Pipeline Versions

| Version | Output directory | Workflow script | Notes |
|---------|-----------------|-----------------|-------|
| v21 CG combined | `tests/output/v21_cg_combined/` | `generate_cg_entangled.py` | CG + entanglements + grafts |
| v21 Atomistic combined | `tests/output/v21_atomistic_combined/` | `generate_atomistic_combined.py` | Atomistic + entanglements + grafts |

**Regression reference:** The regression tests (`tests/regression/test_cg_network.py`, `tests/regression/test_atomistic_network.py`) compare fresh workflow output against the reference directories above.

---

## simbox Crosslink Versions

| Version | Output directory | Workflow script | Key changes |
|---------|-----------------|-----------------|-------------|
| v1 (root) | `tests/output/simbox_crosslink/` | `generate_all.py` (archived) | Original prototype; has `.mol` files for fix/bond/react |
| v2.0–v2.9 | `tests/output/simbox_crosslink/v2*/` | `generate_simbox_crosslink.py` variants | Iterative refinements to packing, atom typing, settings |
| v4/poss_0 | `tests/output/simbox_crosslink/v4/poss_0/` | `generate_simbox_crosslink.py` | **Canonical poss_0.** Fixed writer header bug; added `ff_coeffs.in`; N_3=4, H_=5; seed=42, n_epoxy=600, n_amino=300, n_poss=0 |
| v4/poss_50 | `tests/output/simbox_crosslink/v4/poss_50/` | `generate_simbox_crosslink.py` | 50% POSS fraction; seed=42, n_epoxy=60, n_amino=15, n_poss=15 |
| v4/poss_100 | `tests/output/simbox_crosslink/v4/poss_100/` | `generate_simbox_crosslink.py` | 100% POSS (no Amino-PDMS); seed=42, n_epoxy=60, n_amino=0, n_poss=30 |

**Regression references:**
- `tests/regression/test_simbox_crosslink.py` → `v4/poss_0/`
- `tests/regression/test_simbox_poss.py` → `v4/poss_50/` and `v4/poss_100/`

### Compositions

| Reference | n_epoxy | n_amino | n_poss | Total atoms | POSS fraction |
|-----------|---------|---------|--------|-------------|---------------|
| v4/poss_0   | 600 | 300 |   0 | 77,700 | 0% |
| v4/poss_50  |  60 |  15 |  15 |  9,120 | 50% |
| v4/poss_100 |  60 |   0 |  30 | 10,470 | 100% |

All v4 references: seed=42, density=0.85 g/cm³, atom types Si3=1/O_3=2/C_3=3/N_3=4/H_=5, 34 dihedral types.

---

## Updating the Reference

When the generation code changes in a way that intentionally alters output:

1. Run the workflow with a fixed seed to produce new output.
2. Verify the output is correct (atom count, type ordering, header counts).
3. Copy output to a new version subdirectory (e.g. `v5/poss_0/`).
4. Update `REF_DIR` in the relevant regression test file.
5. Add an entry to this file and to `docs/development/changelog.md`.
