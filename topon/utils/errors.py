"""Friendly error formatting for config loading + Pydantic validation.

Wraps a few common Pydantic v2 error shapes into a one-screen text block
that points at the offending JSON path, explains the constraint in plain
language, and (when possible) offers a "did you mean ...?" suggestion.

Used by `topon validate` and `topon generate` to replace raw stack
traces. Falls back gracefully — if the exception isn't a Pydantic
`ValidationError`, the original `str(exc)` is returned unchanged.
"""
from __future__ import annotations

import difflib
import json
from pathlib import Path
from typing import Iterable


# Common typo hints we've seen in demo configs / forum posts.
_KNOWN_ALIASES = {
    "degree_of_polymerizatin": "degree_of_polymerization",
    "target_densty": "target_density",
    "target_density_gcc": "target_density",
    "bead_density": "target_density",  # legacy alias (already auto-hoisted)
    "lattice": "lattice_size",
    "n_func": "max_functionality",
    "max_funct": "max_functionality",
    "model": "model_type",
}


def _did_you_mean(key: str, candidates: Iterable[str]) -> str | None:
    """Return the closest match from candidates, or None if nothing close."""
    if key in _KNOWN_ALIASES and _KNOWN_ALIASES[key] in candidates:
        return _KNOWN_ALIASES[key]
    matches = difflib.get_close_matches(key, list(candidates), n=1, cutoff=0.7)
    return matches[0] if matches else None


def format_pydantic_error(exc: Exception, config_path: Path | str | None = None) -> str:
    """Pretty-print a Pydantic v2 ValidationError as multi-line text.

    Non-Pydantic exceptions fall through to their default ``str``. Always
    safe to call.
    """
    try:
        from pydantic import ValidationError
    except ImportError:  # very defensive — pydantic is a hard dep
        return str(exc)

    if not isinstance(exc, ValidationError):
        return str(exc)

    lines: list[str] = []
    header = "Configuration error"
    if config_path is not None:
        header += f" in {Path(config_path).resolve()}"
    lines.append(header)
    lines.append("-" * len(header))

    for err in exc.errors():
        loc = ".".join(str(p) for p in err.get("loc", ())) or "<root>"
        msg = err.get("msg", "(no message)")
        etype = err.get("type", "")
        # The pydantic "Literal" message already contains the allowed set.
        lines.append(f"  {loc}: {msg}")
        # Did-you-mean for unknown-key extras.
        if etype == "extra_forbidden" and isinstance(err.get("input"), dict) is False:
            # the path's last segment is the unknown key
            key = err["loc"][-1] if err["loc"] else None
            if key:
                suggestion = _did_you_mean(str(key), [])  # no schema introspection here
                if suggestion:
                    lines.append(f"        did you mean '{suggestion}'?")

    lines.append("")
    lines.append("Hints:")
    lines.append("  - `topon doctor <config.json>` lints common footguns")
    lines.append("  - `docs/USAGE.md` Appendix A lists every schema field")
    return "\n".join(lines)


def load_config_or_die(config_path: str | Path) -> tuple[object, dict]:
    """Load + Pydantic-validate a config, printing a friendly error and
    raising SystemExit(1) on failure. Used by the CLI commands."""
    import sys
    from topon.config import load_config_full

    try:
        return load_config_full(config_path)
    except FileNotFoundError:
        print(f"Configuration file not found: {Path(config_path).resolve()}",
              file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as exc:
        print(f"Configuration is not valid JSON ({Path(config_path).resolve()}):",
              file=sys.stderr)
        print(f"  Line {exc.lineno}, col {exc.colno}: {exc.msg}", file=sys.stderr)
        sys.exit(1)
    except Exception as exc:
        print(format_pydantic_error(exc, config_path), file=sys.stderr)
        sys.exit(1)
