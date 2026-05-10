"""Project-wide pytest configuration.

Auto-applies tier markers (`fast`, `smoke`, `regression`) based on each
test's directory, so individual test files don't need explicit
``pytestmark = ...`` lines. Also auto-skips any test marked
``requires_lammps`` when the ``lmp`` binary is not on ``PATH``.

Run patterns:
    pytest -m fast                          # quick check (~5s)
    pytest tests/unit/chemistry/            # focused on a component
    pytest -m "fast or smoke"               # major pre-push check
    pytest -m regression                    # nightly (~1.5h)

Markers are registered in ``pyproject.toml`` under
``[tool.pytest.ini_options].markers``.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest


_TIER_BY_DIR = {
    "unit": "fast",
    "smoke": "smoke",
    "regression": "regression",
}


def pytest_collection_modifyitems(config, items):
    """Auto-apply tier markers by test path; skip requires_lammps if no lmp."""
    has_lammps = shutil.which("lmp") is not None
    skip_no_lammps = pytest.mark.skip(reason="LAMMPS (lmp) not on PATH")

    for item in items:
        # Find which tier directory the test lives under.
        parts = Path(str(item.path)).parts
        for i, p in enumerate(parts):
            if p == "tests" and i + 1 < len(parts):
                tier = _TIER_BY_DIR.get(parts[i + 1])
                if tier is not None:
                    item.add_marker(getattr(pytest.mark, tier))
                break

        # Skip LAMMPS-running tests when binary is unavailable.
        if not has_lammps and "requires_lammps" in item.keywords:
            item.add_marker(skip_no_lammps)
