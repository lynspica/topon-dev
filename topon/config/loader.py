"""
Configuration loader for Topon.

Handles loading configuration from JSON files and merging with defaults.
"""

import json
from pathlib import Path
from typing import Tuple, Union

from topon.config.schema import ToponConfig


# Top-level keys covered by the Pydantic ToponConfig schema. Any other
# top-level keys in a JSON config are treated as "raw" sections and
# returned separately by load_config_full so callers can forward them
# to Pipeline(..., raw_config=...).
_SCHEMA_KEYS = {"study", "topology", "assignment", "chemistry", "output"}


def load_config_full(
    config_path: Union[str, Path],
) -> Tuple[ToponConfig, dict]:
    """
    Load and split a JSON config into (validated schema, raw extras).

    Top-level keys covered by ToponConfig (study, topology, assignment,
    chemistry, output) are validated. Any other top-level keys (e.g.
    conformation, simulation, execution, experimental) are returned in
    the raw dict, suitable for ``Pipeline(config, raw_config=raw)``.

    Args:
        config_path: Path to the JSON configuration file.

    Returns:
        Tuple ``(ToponConfig, raw_dict)``.

    Raises:
        FileNotFoundError: If config file doesn't exist.
        pydantic.ValidationError: If a schema-covered section is invalid.
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_path, "r", encoding="utf-8") as f:
        config_data = json.load(f)

    # Backward-compat hoists for deprecated key locations:
    #   chemistry.degree_of_polymerization -> assignment.dp_distribution.default.mean
    #   chemistry.bead_density            -> chemistry.target_density (CG)
    # See P2-G note in internal/DEVELOPMENT_INTERNAL.md sec.1.
    chem = config_data.get("chemistry", {})
    if "degree_of_polymerization" in chem:
        dp = chem.pop("degree_of_polymerization")
        assignment = config_data.setdefault("assignment", {})
        dp_dist = assignment.setdefault("dp_distribution", {})
        default = dp_dist.setdefault("default", {})
        default["mean"] = float(dp)
    if "bead_density" in chem:
        chem["target_density"] = float(chem.pop("bead_density"))

    schema_data = {k: v for k, v in config_data.items() if k in _SCHEMA_KEYS}
    raw_data = {k: v for k, v in config_data.items() if k not in _SCHEMA_KEYS}

    return ToponConfig(**schema_data), raw_data


def load_config(config_path: Union[str, Path]) -> ToponConfig:
    """
    Load configuration from a JSON file (schema-only view).

    Backward-compatible entry point: returns just the validated
    ToponConfig and silently drops any raw-extras sections (conformation,
    simulation, execution, experimental). For full access to those
    sections, use :func:`load_config_full`.

    Args:
        config_path: Path to the JSON configuration file.

    Returns:
        Validated ToponConfig object.

    Raises:
        FileNotFoundError: If config file doesn't exist.
        pydantic.ValidationError: If a schema-covered section is invalid.
    """
    config, _ = load_config_full(config_path)
    return config


def merge_configs(*configs: dict) -> dict:
    """
    Deep merge multiple configuration dictionaries.
    
    Later configs override earlier ones.
    
    Args:
        *configs: Configuration dictionaries to merge.
        
    Returns:
        Merged configuration dictionary.
    """
    result = {}
    
    for config in configs:
        result = _deep_merge(result, config)
    
    return result


def _deep_merge(base: dict, override: dict) -> dict:
    """
    Deep merge two dictionaries.
    
    Args:
        base: Base dictionary.
        override: Dictionary with values to override.
        
    Returns:
        Merged dictionary.
    """
    result = base.copy()
    
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    
    return result


def save_config(config: ToponConfig, output_path: Union[str, Path]) -> None:
    """
    Save configuration to a JSON file.
    
    Args:
        config: ToponConfig object to save.
        output_path: Path to save the JSON file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(config.model_dump(), f, indent=2)


def create_default_config() -> ToponConfig:
    """
    Create a default configuration.
    
    Returns:
        ToponConfig with all defaults.
    """
    return ToponConfig()
