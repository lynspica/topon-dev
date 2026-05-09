"""Generic Python-API runner for any topon demo config.

Equivalent to ``topon generate <config_path>`` but exercises the
``topon.pipeline.Pipeline`` class directly — useful as a starting point
for users who want to script the pipeline (e.g. in a Jupyter notebook
or as part of a larger Python workflow) instead of the CLI.

Usage:
    python examples/run_via_api.py demos/polymer/coarse_grained/basic/config.json
    python examples/run_via_api.py demos/polymer/atomistic/combined/config.json
"""
from __future__ import annotations

import sys
from pathlib import Path

from topon.config import load_config
from topon.pipeline import Pipeline


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2

    cfg_path = Path(argv[1])
    if not cfg_path.is_absolute():
        cfg_path = Path(__file__).parent / cfg_path
    cfg_path = cfg_path.resolve()

    if not cfg_path.exists():
        print(f"Error: config not found: {cfg_path}")
        return 1

    print(f"Loading config: {cfg_path}")
    config = load_config(str(cfg_path))
    Pipeline(config).run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
