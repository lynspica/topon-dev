# topon TODOs — Mirrored from Anthology (2026-05-07)

*Mirrored from Anthology TODAY.md (2026-04-24 user content carried forward through 2026-05-07 inbox processing). Captured here so the items live in topon's devlog rather than only in Anthology.*

---

## Open TODOs (from 2026-04-24 user-written TODAY.md, paired with NPJ release)

1. **SELFIES support after SMILES** — extend `topon.chemistry` to accept SELFIES strings as an alternative to SMILES. Useful for downstream ML pipelines (some generative models output SELFIES natively).
2. **Fix GitHub link of topon examples** — the atomistic / CG examples links are broken. Repo housekeeping.
3. **Chain-end problem** — solubility folder fixed; md written for topon describing the issue. Need to confirm the fix and document the canonical handling. *Anthology context: Xu 2025 Sci. Adv. (distilled 2026-05-07) gives a physical mechanism for why chain ends matter — N_EG (chain ends per cooperatively rearranging region) controls Tg via cooperative-motion unjamming. Worth using this physics as the framing for topon's chain-end documentation/handling decisions.*
4. **Martini integration** — add MARTINI force field support to topon's chemistry-mapping subpackage. *Anthology context: Beltukov 2019 (Eur. Phys. J. D, distilled 2026-05-07) demonstrates MARTINI works for PS elastic moduli via internal-pressure protocol — μs-feasible. Salerno 2016 PRL (also distilled 2026-05-07) cautions that 4-CH₂/bead loses some dynamics fidelity but quasi-static moduli are preserved. The `Martini_Ahmet.zip` archive at `Anthology/archive/martini-ahmet/Martini_Ahmet.zip` likely holds Ahmet's prep work.*
5. **npz integration** — add npz file format support (likely for descriptor outputs / dataset packaging compatible with the Ramprasad / camel-vector-1 ML pipelines).
6. **Long HSP calculation (weekend)** — finish the long-running HSP descriptor calculations across the topon polymer set. *Anthology context: Lindvig 2002 (distilled 2026-05-07) provides the predictive HSP→χ recipe; could be used as a calibration baseline against topon's HSP outputs.*
7. **RESP for DREIDING — add to solubility package** — add RESP charge fitting as an option in topon's solubility-package chemistry mapping for DREIDING mode.
8. **Poster** — research-output deliverable. Likely related to the SES 2026 abstract (Ahmet's own conference output, preserved at `Anthology/notes/conference-abstracts/SES_2026.md`) "Topological control of mechanical response in polymer networks" with Arora, Ahn, Keten.

---

## Anthology Sparks Worth Flagging for topon Design Decisions

These are sparks recorded in `Anthology/SPARKS.md` that bear directly on topon design or descriptor-stack choices:

- **N_EG (chain ends per CRR) as new topon descriptor** (May 7, Xu 2025) — graph-countable; should join {F_eff, κ, α/β, CB, SP, Fiedler}. Mechanistically explains POSS Tg shift as topology-intrinsic, not interface-driven.
- **Three independent paths to χ** (May 7) — Lindvig HSP table / Venetsanos KB-MD / ML routes. topon should likely export the descriptors needed for all three (HSP groups, pair-correlation prerequisites, ML-fingerprint inputs).
- **SMiPoly closes the inverse-design loop** (May 7) — SMiPoly's 22 polymerization rules + 1083 monomer catalog is a defensible starting catalog for topon's chemistry-mapping subpackage's monomer set.
- **Sequence distribution at constant composition tunes properties** (May 7) — for the side-chain-crosslinking pattern in topro: end-biased / center-biased / random crosslink placement at constant total density should give different mechanical response. Sequence-design ML on top of this is a natural camel-vector-2 instantiation.
- **Three converging vectors for camel** (Apr 24) — topon outputs (topology, descriptors, networks) feed all three vectors. Vector-3 (LLM orchestration) is the natural lane for topon-MCP integration.

---

## Pointers Back to Anthology

- **INTERESTS.md** — `Anthology/INTERESTS.md`. Active Project #1 (NPJ paper) lists topon cleanup as the companion sprint. Active Project #2 (camel) names topon as the topology-feature provider.
- **Daily digest 2026-05-07** — `Anthology/digests/daily/2026-05-07.md`. Full record of what was processed and how it bears on topon.
- **`Anthology/archive/martini-ahmet/Martini_Ahmet.zip`** — Ahmet's MARTINI prep files. Useful for the Martini integration item.

---

*Update this file as items are completed in topon. Cross-reference back to Anthology for context-rich background on individual items.*
