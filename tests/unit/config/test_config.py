"""
Unit tests for Configuration module.
Validation of schema constraints and conflict detection.
"""

import pytest
from pydantic import ValidationError
from topon.config.schema import (
    ToponConfig, 
    StudyConfig, 
    TopologyConfig,
    AssignmentConfig,
    ChemistryConfig,
    OutputConfig
)
from topon.config.validator import validate_config

# =============================================================================
# Schema Tests (Pydantic)
# =============================================================================

def test_default_config_validity():
    """Test that default configuration is valid."""
    config = ToponConfig()
    assert config.study.name == "my_network"
    assert config.topology.generator.lattice_size == "6x6x6"
    assert config.chemistry.model_type == "coarse_grained"

def test_config_invalid_lattice_type():
    """Test validation of lattice_type enum."""
    with pytest.raises(ValidationError):
        ToponConfig(
            topology=TopologyConfig(
                generator={"lattice_type": "INVALID"}
            )
        )

def test_config_range_validation():
    """Test numeric range validation."""
    # Test PDI >= 1.0
    with pytest.raises(ValidationError):
        ToponConfig(
            assignment=AssignmentConfig(
                dp_distribution={"default": {"pdi": 0.5}}
            )
        )
    
    # Test Fraction <= 1.0
    with pytest.raises(ValidationError):
        ToponConfig(
            assignment=AssignmentConfig(
                copolymer={
                    "per_edge_type": {
                        "A": {
                            "composition": [{"monomer": "M1", "fraction": 1.5}]
                        }
                    }
                }
            )
        )

# =============================================================================
# Logical Validation Tests (validator.py)
# =============================================================================

@pytest.fixture
def base_config():
    """Return a valid base configuration."""
    return ToponConfig()

def test_copolymer_edge_conflict(base_config):
    """Test conflict between enabled copolymer and heterogeneous edge method."""
    # Enable copolymer
    base_config.assignment.copolymer.enabled = True
    # Set incompatible edge method
    base_config.assignment.edge_types.method = "random"
    
    errors = validate_config(base_config)
    assert any("CONFLICT" in e for e in errors)
    assert any("heterogeneous edge types" in e for e in errors)

def test_missing_node_type_mapping(base_config):
    """Test detection of used node types missing from chemistry map."""
    # Switch to the random method and use a type 'C' not in the chemical map
    base_config.assignment.node_types.method = "random"
    base_config.assignment.node_types.random.type_ratios = {"A": 50, "C": 50}

    # 'C' is not in default chemical mapping
    errors = validate_config(base_config)
    assert any("Node type 'C'" in e for e in errors)

def test_missing_monomer_reference(base_config):
    """Test detection of missing monomers."""
    # Reference unknown monomer in edge type
    base_config.chemistry.edge_type_map["A"].monomer = "UNKNOWN_MONOMER"
    
    errors = validate_config(base_config)
    assert any("Monomer 'UNKNOWN_MONOMER'" in e for e in errors)

def test_target_constraints(base_config):
    """Test validation of target values against max possible."""
    base_config.assignment.entanglements.enabled = True
    base_config.assignment.entanglements.target = 100
    
    max_possible = {"entanglements": 50}
    errors = validate_config(base_config, max_possible=max_possible)
    
    assert any("Requested 100 entanglements" in e for e in errors)
    assert any("only 50 possible" in e for e in errors)
