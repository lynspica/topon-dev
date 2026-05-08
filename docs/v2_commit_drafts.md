# Draft commit messages for solubility v2 topon upstreams

These are the two commits the solubility v2 work wants topon to carry.
They are separate and independent.  Both sit uncommitted in the working
tree; push them at your discretion after code review.

---

## Commit 1 — "Add ETKDGv3 chain embedder (fixes branched-chain coord collapse)"

```
Add ETKDGv3 chain embedder (fixes branched-chain coord collapse)

`_assign_extended_linear_coords` (singlechain/workflow.py) walks only
the backbone heavy-atom path.  Pendant atoms — notably the three
methyls on each Si of a PDMS repeat, but also the four carbons on the
quaternary C of polyisobutylene — are placed at their parent atom's
position with a random displacement.  At DP=30 PDMS this produces an
observed minimum pair-atom distance of 0.436 Å (Si-methyl C on top of
backbone O), which translates to LJ pair potentials of ~-10²⁷ kcal/mol
and makes LAMMPS abort at step 1 with "Non-numeric box dimensions —
simulation unstable".

Replacement: `topon.chemistry.embed.embed_with_etkdg(mol, seed)` runs
RDKit's ETKDGv3 (`useRandomCoords=True` for robustness on large
molecules) followed by MMFF94 (200 iters) with UFF (150 iters) as a
fallback when MMFF cannot parameterise the molecule (typical for
Si-containing chains).

Backwards compatibility:
  - `_assign_extended_linear_coords` remains in place but now emits
    DeprecationWarning with a pointer at the replacement.
  - `run_workflow` gains an `embedder="linear"|"etkdg"` parameter,
    defaulting to "linear" so every pre-existing caller is unchanged.
    New code should pass `embedder="etkdg"` for any polymer with side
    chains or stereochemistry.

Regression tests in tests/unit/test_embed.py:
  - ETKDGv3 gives min_pairwise_distance > 0.8 Å on DP=10 PDMS.
  - Deterministic across same-seed runs.
  - Legacy placer is explicitly pinned as producing min dist < 0.8 Å,
    so the test is a red flag if the linear placer is ever silently
    "fixed".
  - Embedding an empty RWMol raises.

Related downstream: the solubility package (v0.2.0) imports this
function from topon and removes its local copy of the embedder.

Closes DOW/Studies/solubility v2_roadmap.md priority 1a.
```

Files touched:
* `topon/chemistry/embed.py`           — new
* `topon/chemistry/__init__.py`        — re-export
* `topon/singlechain/workflow.py`      — add `embedder` param + deprecate linear placer
* `tests/unit/test_embed.py`           — new

---

## Commit 2 — "Accept legacy 'type' graph attribute in ChemistryBuilder"

```
Accept legacy 'type' graph attribute in ChemistryBuilder

ChemistryBuilder._build_nodes reads the graph attribute 'node_type'
(with a hard-coded "A" default).  Several legacy callers — most
prominently DOW/solvent_effects/generate_matrix.py::build_chain — set
'type' instead of 'node_type'.  Before this commit, topon silently
ignored the legacy attribute and fell through to the "A" default,
which maps to the Si-atom simple-atom fallback.  For hydrocarbon
polymers (Butyl, EPDM, NBR, FKM, Polyacrylate) this leaked Si atoms
into the built chain.

After this commit, the builder reads 'node_type' preferentially, then
'type' as a back-compat fallback (emitting a DeprecationWarning), then
falls through to "A" only if neither is present.  New callers should
use 'node_type'; the warning makes the migration path obvious.

Regression tests in tests/unit/test_node_attributes.py:
  - Butyl DP=4 built from a `type="end"` graph has no Si atoms.
  - Building with `type` emits DeprecationWarning mentioning 'node_type'.
  - Building with `node_type` emits no warning.

Closes DOW/Studies/solubility v2_roadmap.md priority 1b.
```

Files touched:
* `topon/chemistry/builder.py`           — fallback + deprecation warning
* `tests/unit/test_node_attributes.py`   — new

---

## How to commit

```bash
cd /c/Users/ahmet/OneDrive\ -\ Northwestern\ University/DOW-Ahmet/topon

# Commit 1
git add topon/chemistry/embed.py \
        topon/chemistry/__init__.py \
        topon/singlechain/workflow.py \
        tests/unit/test_embed.py
git commit -F <(sed -n '/^## Commit 1/,/^Files touched/p' docs/v2_commit_drafts.md \
                 | sed '1,2d;/^Files touched/d' | sed -n '/^```$/,/^```$/p' | sed '1d;$d')

# Commit 2
git add topon/chemistry/builder.py \
        tests/unit/test_node_attributes.py
git commit -F <(sed -n '/^## Commit 2/,/^Files touched/p' docs/v2_commit_drafts.md \
                 | sed '1,2d;/^Files touched/d' | sed -n '/^```$/,/^```$/p' | sed '1d;$d')
```

Or simply copy the message bodies into `git commit -m` manually.
