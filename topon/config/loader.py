"""
Configuration loader for Topon.

Handles loading configuration from JSON files and merging with defaults.
"""

import json
from pathlib import Path
from typing import Union

from topon.config.schema import ToponConfig


def load_config(config_path: Union[str, Path]) -> ToponConfig:
    """
    Load configuration from a JSON file.
    
    Args:
        config_path: Path to the JSON configuration file.
        
    Returns:
        Validated ToponConfig object.
        
    Raises:
        FileNotFoundError: If config file doesn't exist.
        ValidationError: If config is invalid.
    """
    config_path = Path(config_path)
    
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, "r", encoding="utf-8") as f:
        config_data = json.load(f)
    
    return ToponConfig(**config_data)


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
