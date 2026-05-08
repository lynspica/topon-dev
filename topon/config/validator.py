"""
Configuration validator for Topon.

Handles validation of configuration beyond Pydantic schema validation,
including conflict detection and constraint checking.
"""

from typing import Optional
from topon.config.schema import ToponConfig


class ConfigValidationError(Exception):
    """Raised when configuration validation fails."""
    
    def __init__(self, errors: list[str]):
        self.errors = errors
        message = "Configuration validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
        super().__init__(message)


def validate_config(config: ToponConfig, max_possible: Optional[dict] = None) -> list[str]:
    """
    Validate configuration for conflicts and constraint violations.
    
    Args:
        config: ToponConfig object to validate.
        max_possible: Optional dict with max possible values from graph analysis.
                      Keys: 'primary_loops', 'secondary_loops', 'entanglements'
                      
    Returns:
        List of error messages (empty if valid).
    """
    errors = []
    
    # Check for copolymer + heterogeneous edge types conflict
    errors.extend(_check_copolymer_edge_conflict(config))
    
    # Check for missing type mappings
    errors.extend(_check_type_mappings(config))
    
    # Check for missing monomers
    errors.extend(_check_monomer_references(config))
    
    # Check target constraints if max_possible provided
    if max_possible:
        errors.extend(_check_target_constraints(config, max_possible))
    
    return errors


def validate_config_strict(config: ToponConfig, max_possible: Optional[dict] = None) -> None:
    """
    Validate configuration and raise exception if invalid.
    
    Args:
        config: ToponConfig object to validate.
        max_possible: Optional dict with max possible values.
        
    Raises:
        ConfigValidationError: If validation fails.
    """
    errors = validate_config(config, max_possible)
    if errors:
        raise ConfigValidationError(errors)


def _check_copolymer_edge_conflict(config: ToponConfig) -> list[str]:
    """Check for copolymer + heterogeneous edge types conflict."""
    errors = []
    
    copolymer_enabled = config.assignment.copolymer.enabled
    edge_method = config.assignment.edge_types.method
    
    if copolymer_enabled and edge_method in ["random", "composite"]:
        errors.append(
            "CONFLICT: Cannot use copolymer with heterogeneous edge types. "
            "Choose one: either copolymer (same edge type, different monomers within chain) "
            "OR heterogeneous edges (different edge types with different chemistry)."
        )
    
    return errors


def _check_type_mappings(config: ToponConfig) -> list[str]:
    """Check that all assigned types have corresponding chemistry mappings."""
    errors = []
    
    # Get all possible node types from assignment config
    node_types_used = set()
    
    # From degree mapping
    for node_type in config.assignment.node_types.degree.mapping.values():
        node_types_used.add(node_type)
    
    # From positional
    for node_type in config.assignment.node_types.positional.layer_types:
        node_types_used.add(node_type)
    
    # From random
    for node_type in config.assignment.node_types.random.type_ratios.keys():
        node_types_used.add(node_type)
    
    # Check all are in node_type_map
    for node_type in node_types_used:
        if node_type not in config.chemistry.node_type_map:
            errors.append(
                f"Node type '{node_type}' is used in assignment but not defined in chemistry.node_type_map"
            )
    
    # Get all possible edge types from assignment config
    edge_types_used = set()
    
    # From uniform
    edge_types_used.add(config.assignment.edge_types.uniform.type)
    
    # From random
    for edge_type in config.assignment.edge_types.random.type_ratios.keys():
        edge_types_used.add(edge_type)
    
    # From composite
    for edge_type in config.assignment.edge_types.composite.layer_types:
        edge_types_used.add(edge_type)
    
    # Check all are in edge_type_map
    for edge_type in edge_types_used:
        if edge_type not in config.chemistry.edge_type_map:
            errors.append(
                f"Edge type '{edge_type}' is used in assignment but not defined in chemistry.edge_type_map"
            )
    
    return errors


def _check_monomer_references(config: ToponConfig) -> list[str]:
    """Check that all referenced monomers exist in the library."""
    errors = []
    
    available_monomers = set(config.chemistry.monomers.keys())
    
    # Check edge_type_map references
    for edge_type, edge_chem in config.chemistry.edge_type_map.items():
        if edge_chem.monomer not in available_monomers:
            errors.append(
                f"Monomer '{edge_chem.monomer}' referenced by edge type '{edge_type}' "
                f"is not defined in chemistry.monomers"
            )
    
    # Check graft references
    if config.assignment.grafts.enabled:
        for edge_type, graft_config in config.assignment.grafts.per_edge_type.items():
            if graft_config.side_chain_monomer not in available_monomers:
                errors.append(
                    f"Graft monomer '{graft_config.side_chain_monomer}' for edge type '{edge_type}' "
                    f"is not defined in chemistry.monomers"
                )
    
    # Check copolymer references
    if config.assignment.copolymer.enabled:
        for edge_type, copoly_config in config.assignment.copolymer.per_edge_type.items():
            for comp in copoly_config.composition:
                if comp.monomer not in available_monomers:
                    errors.append(
                        f"Copolymer monomer '{comp.monomer}' for edge type '{edge_type}' "
                        f"is not defined in chemistry.monomers"
                    )
    
    return errors


def _check_target_constraints(config: ToponConfig, max_possible: dict) -> list[str]:
    """Check that target values don't exceed max possible."""
    errors = []
    
    # Primary loops
    if config.assignment.defects.primary_loops.enabled:
        target = config.assignment.defects.primary_loops.target
        target_type = config.assignment.defects.primary_loops.target_type
        max_val = max_possible.get("primary_loops", 0)
        
        if target_type == "count" and target > max_val:
            errors.append(
                f"Requested {target} primary loops but only {max_val} possible in this graph"
            )
        elif target_type == "percentage" and target > 100:
            errors.append(f"Primary loops percentage cannot exceed 100 (got {target})")
    
    # Secondary loops
    if config.assignment.defects.secondary_loops.enabled:
        target = config.assignment.defects.secondary_loops.target
        target_type = config.assignment.defects.secondary_loops.target_type
        max_val = max_possible.get("secondary_loops", 0)
        
        if target_type == "count" and target > max_val:
            errors.append(
                f"Requested {target} secondary loops but only {max_val} possible in this graph"
            )
        elif target_type == "percentage" and target > 100:
            errors.append(f"Secondary loops percentage cannot exceed 100 (got {target})")
    
    # Entanglements
    if config.assignment.entanglements.enabled:
        target = config.assignment.entanglements.target
        target_type = config.assignment.entanglements.target_type
        max_val = max_possible.get("entanglements", 0)
        
        if target_type == "count" and target > max_val:
            errors.append(
                f"Requested {target} entanglements but only {max_val} possible in this graph"
            )
        elif target_type == "percentage" and target > 100:
            errors.append(f"Entanglements percentage cannot exceed 100 (got {target})")
    
    return errors
