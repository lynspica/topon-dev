"""Config-linting rules surfaced by `topon doctor`.

A rule is a function ``(toponConfig, raw_dict) -> Iterable[Issue]``. Each
``Issue`` has a ``level`` (``"ok"``, ``"warn"``, ``"error"``), a short
``rule`` name, a human-readable ``message``, and an optional ``fix`` hint.

Rules are stateless and cheap (no chemistry build, no LAMMPS). They exist
to catch the known footguns in `internal/DEVELOPMENT_INTERNAL.md` before
the user spends 30 minutes on a build that produces nonsense.
"""
from .rules import Issue, run_all_rules, RULE_REGISTRY

__all__ = ["Issue", "run_all_rules", "RULE_REGISTRY"]
