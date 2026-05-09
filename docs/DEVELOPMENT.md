# topon — Development

> **Stub.** Being populated as part of the doc consolidation (2026-05-08).
> Until this is filled in, the canonical sources are:
> - Version-by-version history: [`docs/development/changelog.md`](development/changelog.md)
> - Closed and open phase work: [`docs/development/tasks.md`](development/tasks.md)
> - simbox version log: tail of [`docs/simbox.md`](simbox.md)

## Methodology (placeholder)

Each non-trivial change should walk through:

**Questions → Research → Requirements → Roadmap**

1. **Questions** — what is unclear, what assumption could be wrong, what would make this concrete?
2. **Research** — read the relevant code and docs; cite file:line; check the regression suite.
3. **Requirements** — write down what the change must do and must NOT do, and the regression boundary.
4. **Roadmap** — sequence the change into reviewable steps; each step a separate commit.

This section will be expanded with concrete examples and links to past changes that followed the methodology.

See [ARCHITECTURE.md](ARCHITECTURE.md) for the package layout.
