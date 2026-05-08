# topon issues surfaced by the solubility v2 work

Written by the solubility package maintainer after completing v0.2.0
(2026-04-24).  Two issues (ETKDGv3 embedder, `node_type` fallback) were
upstreamed into topon as v2 commits 1a and 1b.  This document collects
the **remaining** issues that the solubility package works around but
should be fixed inside topon so that every downstream project benefits
(notably `DOW/solvent_effects/generate_matrix.py`, which hits all the
same bugs silently).

Each entry records:

1. **What** — the observable symptom.
2. **Where** — the function and (approximate) line.
3. **Impact** — what the downstream user sees if unaware.
4. **Evidence** — the v2 test / build run that surfaced it.
5. **Proposed fix** — a concrete change sketch.

Issues are grouped by severity.  Priority labels map to v3 candidates:

* 🔴 **P0 — correctness bug, silent failure that affects physical results.**
* 🟠 **P1 — robustness issue, user-visible when they look but easy to miss.**
* 🟡 **P2 — developer experience, not a science bug.**

---

## ✅ P0-1. Auto-bridge creates peroxide `-O-O-` on PDMS / PTFPMS chains [RESOLVED]

**Status (2026-05-01):** fixed in `topon/chemistry/builder.py` by removing the
trailing `[O]` placeholder in `_create_chain_from_smiles`. Module now exports
`_PEROXIDE_FIX_APPLIED = True` so downstream projects (e.g.
`solubility/_topon_local/`) can feature-detect and retire their local fork.
Verified: PDMS DP=10 build now yields canonical SMILES with no `OO` substring,
0 O-O bonds, and correct atom counts (Si=12, O=11, C=26 with TMS caps).
Reference LAMMPS data files in `tests/output/solvent_effects/PDMS/*` and
`tests/output/solvent_effects/PTFPMS/*` will need regenerating once — they
will lose 1 O atom and 1 O-O bond per chain.

**What.** For repeat-unit SMILES that already end in a linkable atom
(e.g. PDMS `[Si](C)(C)O`), the chain built by `ChemistryBuilder` gets
an extra bridge O at the head end:

    ...-[Si(CH3)3]-O-O-[Si(CH3)2-O]n-...
                    ^^^ peroxide artifact

**Where.** `topon/chemistry/builder.py::_create_chain_from_smiles`
(appends a trailing `[O]` marker when monomer ends in `O`) combined
with `_build_chain_atomistic` (auto-bridge logic at the head end).

**Impact.**

* In OPLS typing, the extra O is classified as a siloxane O
  (charge −0.30 e), leaving every PDMS chain **net −0.30 e**.  Summed
  over 10 chains that's −3 e excess in the simulation cell — LAMMPS
  warns but runs, and the extra net charge shifts δ_ele_short /
  δ_ele_long by a few tenths of MPa^½ on screening runs.
* Chemically the peroxide is not present in real PDMS.  Any bond/angle
  FF parameters for Si-O-O-Si don't exist in real OPLS-AA; topon's
  builder currently forces downstream code to invent them.
* Affects every caller that uses the default auto-bridge path with
  O-terminal monomer SMILES.  Silent for non-OPLS tiers where the
  extra O's Gasteiger charge happens to balance out.

**Evidence.**

* `solubility/KNOWN_ISSUES.md §A (residual)`: documented after v2.
* Reproducer: build a PDMS DP=5 chain via
  `solubility.chemistry.homopolymer.build_homopolymer_chain(get_polymer('PDMS'), dp=5, ...)`
  and inspect the SMILES of the resulting chain — you'll see
  `C[Si](C)(C)OO[Si](C)(C)...` (double-O).
* Same artifact appears in `E:/PhD/DOW/solvent_effects/generate_matrix.py`
  builds, currently unreported because generate_matrix.py doesn't type
  with OPLS.

**Proposed fix.** Restructure `_create_chain_from_smiles` to emit the
heavy-atom chain without the trailing `[O]` placeholder:

1. Drop the `(unit + linker) * dp + "[O]"` pattern.  Build the chain
   as `(unit + linker) * (dp - 1) + unit` so the last repeat's linker
   O becomes the attachment point naturally, with no separate `[O]`.
2. In `_build_chain_atomistic`, auto-bridge should check whether the
   chain's head atom is already bonded to an atom of the end-cap's
   symbol *through the monomer SMILES*; if so, skip the bridge.
3. Add a regression test: build PDMS DP=5, sanitize, assert the SMILES
   contains no `OO` substring and atoms/bonds match the expected 5-mer
   stoichiometry.

---

## 🔴 P0-2. Default node_type falls through to "Si" simple-atom

**What.** When a graph has no `node_type_map` entry for a node type (or
no `node_type` attribute at all — see v2 priority 1b), topon silently
places a single Si atom at that node.

**Where.** `topon/chemistry/builder.py::_build_nodes` lines ~130–140:

    if not node_config:
        molecule = "Si"
        is_end_cap = degree == 1

**Impact.** Hydrocarbon polymers get silently contaminated with Si
atoms when the caller forgets to supply a `node_type_map`.  v2
priority 1b partially fixed the attribute-name mismatch (`type` vs
`node_type`), but the default-to-Si fallthrough is still live.  Even
with correct attributes, a typo in `node_type_map` keys
(`{"End": ...}` vs `"end"`) silently re-triggers the same bug.

**Evidence.**

* `solubility/SPEC_QUESTIONS.md #1` and
  `solubility/SPEC_QUESTIONS.md #2` — Si leaked into Butyl / EPDM /
  NBR / FKM / Polyacrylate chains during v1, producing
  `"No types for atom ..."` foyer errors on OPLS typing.
* Reproducer in `solubility/tests/test_chemistry.py::test_homopolymer_chain_builds`
  (the assertion `n_si2 == 0` would fail without the v1 workaround).

**Proposed fix.**

1. Change the default from a silent `molecule = "Si"` to one of:
   - Raise `ChemistryBuilderError` with an explicit message naming
     the unmatched node type and listing the keys present in
     `node_type_map`.
   - Or warn loudly (not silent) and fall through to a hydrogen atom
     (neutral, non-reactive) instead of Si.
2. Add a config-validation hook that checks every node's `node_type`
   against `config.node_type_map.keys()` before building.

---

## 🟠 P1-3. `Chem.MolFromSmiles("Si")` fails silently in `_place_simple_atom`

**What.** When the end-cap fallback is hit with atom symbol `"Si"`
(the P0-2 default), topon calls `Chem.MolFromSmiles("Si")` — which
RDKit rejects because `Si` is two letters (S + I) and not a valid
SMILES unless wrapped `[Si]`.  The function catches the parse error
and tries `Chem.Atom("Si")` next, which succeeds.

**Where.** `topon/chemistry/builder.py::_place_simple_atom` lines
~155–170.

**Impact.** Not a correctness bug (the fallback path works), but
during v2 we observed repeated stderr noise:

    SMILES Parse Error: Failed parsing SMILES 'Si' for input: 'Si'

printed many times per build.  Confuses users who think the build is
failing.

**Evidence.** Every hydrocarbon-polymer build in v1 produced a torrent
of these messages.  Documented in
`solubility/SPEC_QUESTIONS.md` (implicit — the build flow discussed
there manifested this).

**Proposed fix.** Inside `_place_simple_atom`, when `atom_symbol` is a
one-or-two-letter element symbol, don't call `MolFromSmiles` — just
construct `Chem.Atom(atom_symbol)` directly (this is already the path
the code takes for `atom_symbol.isalpha()`; the bug is that the
MolFromSmiles path is *also* tried first when called from some code
paths).  Quick audit of all `Chem.MolFromSmiles(<element-like-string>)`
call sites recommended.

---

## 🟠 P1-4. `_guess_head` / `_guess_tail` is a regex that can silently return wrong atom

**What.** `topon/singlechain/workflow.py::_guess_head` scans the SMILES
for `[A-Z][a-z]?` and returns the first matching element, without
checking that the element is actually at a chain-head position.

**Where.** Lines ~667–723 of `singlechain/workflow.py`.

**Impact.** For unusual repeat units (stereo markers, ring
specifications, leading `(`), the guess silently picks the wrong atom.
This then feeds `MonomerConfig(chain_head=wrong_symbol)` and the
builder connects via the wrong atom, producing a malformed chain.
Currently every solubility v2 monomer happens to have a
straightforwardly-parseable SMILES so the bug doesn't bite, but a
general-purpose polymer library would trip over it.

**Evidence.** No failures in the v2 test suite, but the heuristic is
inspected in `solubility/SPEC_QUESTIONS.md #1` and flagged.  Test on
e.g. `smiles="C[C@H](C)C(=O)O"` (a chiral amino-acid-like monomer) —
the regex returns `"C"` but the real chain-head depends on which
carbon carries the extension.

**Proposed fix.** Replace the regex heuristic with an RDKit-based
walk:

1. Parse the SMILES to a Mol.
2. Identify the two "terminal" atoms of the backbone (atoms with
   degree ≥ 1 that sit on a longest-path through heavy atoms).
3. Return their element symbols.

This matches what `_find_backbone_path` already does for extended
linear coords; extract a shared helper.

---

## 🟠 P1-5. `BoxPacker` default `min_dist=0.0` produces overlapping atoms

**What.** `topon.simbox.packer.BoxPacker(min_dist=0.0, ...)` by default
disables overlap detection completely.  Packed molecules can be placed
with atoms at identical coordinates, which causes LAMMPS to abort with
`Non-numeric box dimensions — simulation unstable` at step 1.

**Where.** `topon/simbox/packer.py::BoxPacker.__init__` line 164.

**Impact.** Every caller that uses the default is silently playing
Russian roulette with LAMMPS stability.  In v2 the solubility package
explicitly passes `min_dist=2.0` to work around it
(`solubility/cellbuilder/polymer_bulk.py` line 100).

**Evidence.** During v1 integration, the default min_dist=0.0 caused
LAMMPS to abort on every PDMS / FKM cell.  Fix: bump to 2.0 in
solubility.

**Proposed fix.** Change the default to `min_dist=2.0` (a
chemically-sensible C-C non-bonded lower limit).  Callers that
deliberately want to allow overlaps (e.g. for SoftPacking pre-stages)
can still pass 0.0 explicitly.  Add a deprecation period: warn
when `min_dist=0.0` is used explicitly so downstream callers migrate.

---

## 🟡 P2-6. Verbose `print()` throughout `ChemistryBuilder` and `BoxPacker`

**What.** Every chain build emits:

    Building chemistry...
      Model type: atomistic
      Building nodes...
        Placed 2 node structures
      Building chains...
        Built 1 chains
      Total atoms: XX
      Total bonds: YY

And every packing run emits:

    [BoxPacker] Target density = ...
    [BoxPacker] Initial box = ...
    [BoxPacker] Placed N/M molecules

For a 27-config × 3-replica × 10-chain solubility build, that's
~900 `Building chemistry...` blocks in stdout.  Pipe-buffered through
`tee | grep`, it's the dominant reason stdout looks frozen during
builds.

**Where.** Scattered through `topon/chemistry/builder.py` and
`topon/simbox/packer.py`.

**Impact.** Dev UX only — no science effect.

**Proposed fix.** Replace every `print()` with a module-level
`logger = logging.getLogger(__name__)` and use `logger.info(...)`.
Downstream callers then silence topon with:

    logging.getLogger("topon").setLevel(logging.WARNING)

…which is exactly what `solubility build-all` wants to do.

---

## 🟡 P2-7. Heavy imports inside hot functions

**What.** `_build_nodes`, `_build_chain_atomistic`, and several other
inner-loop functions import RDKit at every call:

    def _build_nodes(self):
        from rdkit import Chem  # ← inside the loop body

For a 900-chain build that's 900 `from rdkit import Chem` calls.
Python caches imports, so the *cost* is negligible, but it makes the
code harder to read and fogs the dependency graph.

**Where.** Multiple places in `topon/chemistry/builder.py` and
`topon/singlechain/workflow.py`.

**Impact.** Cosmetic.

**Proposed fix.** Move all intra-function `from rdkit import Chem` to
the module top.  Use `TYPE_CHECKING` guards only for the genuine
optional-dependency case (topon's cli may work without rdkit for CG
models).

---

## 🟡 P2-8. `ChemistryConfig.target_density` is required but unused for single-chain builds

**What.** Every single-chain caller (including generate_matrix.py and
every solubility polymer builder) has to supply `target_density` to
`ChemistryConfig`, but the builder never reads it when building one
chain — density is used later by `BoxPacker`.

**Where.** `topon/config/schema.py::ChemistryConfig` field declaration.

**Impact.** Cognitive overhead — callers think the density matters for
chain build (it doesn't) and pick arbitrary values that then get
logged as "target density" in the chain-build output.

**Proposed fix.** Make `target_density` optional in `ChemistryConfig`
(default None).  `BoxPacker` already takes density separately.  If
both are set, validate they match.

---

## 🟡 P2-9. Silent `except Exception: pass` in `build()` around SanitizeMol

**What.** In topon builders and in many downstream callers, post-build
sanitization is wrapped in a catch-all:

    try: Chem.SanitizeMol(rw)
    except: pass

This hides real chemistry errors (e.g. valence violations introduced
by auto-bridge bugs).

**Where.** Not in topon itself but in every downstream caller that
copied the pattern, including (previously)
`solubility/chemistry/homopolymer.py` until v2 removed it.

**Impact.** Downstream maintainers propagate the pattern instead of
fixing the underlying sanitization issues.

**Proposed fix.** Inside topon, after `_build_chain_atomistic`, do a
`SanitizeMol` with the set of safe flags (skip kekulization for
aromatic-heavy cases, e.g. `SANITIZE_ALL ^ SANITIZE_KEKULIZE`) and
*raise* on failure.  Document which monomer SMILES forms are expected
to pass.  This surfaces bugs like P0-1 at build time rather than at
LAMMPS time.

---

## 🟡 P2-10. No version or changelog — downstream pinning is fragile

**What.** `topon/pyproject.toml` says `version = "0.1.0"` and hasn't
moved across the last two topon edit rounds.  Downstream
(`solubility/pyproject.toml`) says `topon >= 0.1.0`, so solubility
v0.2.0 *accidentally* depends on topon v2's ETKDGv3 embedder but
nothing enforces it.

**Where.** `topon/pyproject.toml`.

**Impact.** A fresh solubility install against a pre-v2 topon will
silently fall back to the legacy linear placer and produce the
0.4-Å-collapsed coordinates we spent v2 fixing.

**Proposed fix.**

1. Bump `topon` version to `0.2.0` when the v2 commits (ETKDGv3 +
   `node_type` fallback) land, `0.3.0` when the P0-1 / P0-2 fixes
   from this document land.
2. Maintain `topon/CHANGELOG.md` with a bullet per commit referencing
   the downstream project that caught each bug.
3. Pin `solubility/pyproject.toml::dependencies`:
   `"topon>=0.3.0"` once P0 fixes land.

---

## Summary of solubility v3 asks

In priority order of what would most unblock the solubility package:

| P | Issue | v3 action in solubility |
|---|---|---|
| 🔴 P0-1 | Auto-bridge peroxide | Upstream fix lands → delete `KNOWN_ISSUES.md §A (residual)` |
| 🔴 P0-2 | Silent Si default | Upstream fix lands → remove dual-attribute workaround |
| 🟠 P1-5 | `BoxPacker` default | Upstream fix lands → remove `min_dist=2.0` override |
| 🟠 P1-3 | `MolFromSmiles("Si")` noise | Upstream fix lands → remove stderr redirect in solubility CLI |
| 🟠 P1-4 | `_guess_head` regex | Required before adding any new polymer with non-trivial SMILES |
| 🟡 P2-6 | Prints → logger | Quality-of-life; enables silencing in `build-all` |
| 🟡 P2-10 | Version + changelog | Fix the dependency pin hygiene |
| 🟡 P2-7, -8, -9 | Code hygiene | Do as part of any larger topon refactor |

The v3 candidate list (solubility side) explicitly names only the
auto-bridge fix as #2; after reviewing this document the user can
decide whether to also promote P0-2 / P1-5 into v3 scope.
