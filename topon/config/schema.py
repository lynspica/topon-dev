"""
Pydantic schema definitions for Topon configuration.

This module defines all configuration models using Pydantic for validation.
The design follows the principle: Assignment = Abstract Types, Chemistry = Concrete Molecules.
"""

from typing import Optional, Literal
from pydantic import BaseModel, Field, field_validator


# =============================================================================
# STUDY CONFIG
# =============================================================================

class StudyConfig(BaseModel):
    """Study-level configuration."""
    name: str = Field(default="my_network", description="Name of the study")
    output_dir: str = Field(default="./output", description="Output directory path")


# =============================================================================
# TOPOLOGY CONFIG
# =============================================================================

class GeneratorConfig(BaseModel):
    """C generator configuration."""
    exe_path: Optional[str] = Field(default=None, description="Path to generator executable")
    lattice_size: str = Field(default="6x6x6", description="Lattice dimensions (e.g., '6x6x6')")
    lattice_type: Literal["SC", "BCC", "FCC"] = Field(default="SC", description="Lattice type")
    periodicity: str = Field(default="111", description="Periodicity string (e.g., '111' for all periodic)")
    max_functionality: int = Field(default=6, ge=1, description="Maximum node functionality")
    max_trials: int = Field(default=1000000, ge=1, description="Maximum trials for generator")
    max_saves: int = Field(default=1, ge=1, description="Number of networks to save")
    degree_distribution: str = Field(
        default="0:0,1:0",
        description="Target degree distribution (e.g., '0:15,1:30,3:43')"
    )


class ExistingFilesConfig(BaseModel):
    """Configuration for loading existing topology files."""
    nodes_file: Optional[str] = Field(default=None, description="Path to .nodes file")
    edges_file: Optional[str] = Field(default=None, description="Path to .edges file")
    gpickle_file: Optional[str] = Field(default=None, description="Path to .gpickle file")


class TopologyConfig(BaseModel):
    """Topology generation/loading configuration."""
    source: Literal["generate", "load"] = Field(
        default="load",
        description="Whether to generate new topology or load existing"
    )
    generator: GeneratorConfig = Field(default_factory=GeneratorConfig)
    existing_files: ExistingFilesConfig = Field(default_factory=ExistingFilesConfig)


# =============================================================================
# ASSIGNMENT CONFIG
# =============================================================================

class DegreeNodeTypeConfig(BaseModel):
    """Degree-based node type assignment."""
    mapping: dict[str, str] = Field(
        default={"1": "end", "2": "A", "3": "A", "4": "A", "5": "A", "6": "A"},
        description="Map degree -> node type"
    )


class PositionalConfig(BaseModel):
    """Positional (layer-based) assignment config."""
    dimension: Literal["x", "y", "z"] = Field(default="z")
    num_layers: int = Field(default=2, ge=1)
    layer_types: list[str] = Field(default=["A", "B"])


class RandomTypeConfig(BaseModel):
    """Random type assignment config."""
    type_ratios: dict[str, float] = Field(
        default={"A": 100},
        description="Type ratios (will be normalized)"
    )


class NodeTypesConfig(BaseModel):
    """Node type assignment configuration."""
    method: Literal["degree", "positional", "random", "explicit"] = Field(default="degree")
    degree: DegreeNodeTypeConfig = Field(default_factory=DegreeNodeTypeConfig)
    positional: PositionalConfig = Field(default_factory=PositionalConfig)
    random: RandomTypeConfig = Field(default_factory=RandomTypeConfig)
    explicit: dict[int, str] = Field(default_factory=dict, description="Per-node ID type assignment")


class UniformEdgeConfig(BaseModel):
    """Uniform edge type config."""
    type: str = Field(default="A")


class CompositeEdgeConfig(BaseModel):
    """Composite/lamellar edge type config."""
    dimension: Literal["x", "y", "z"] = Field(default="z")
    num_layers: int = Field(default=2, ge=1)
    layer_types: list[str] = Field(default=["A", "B"])


class EdgeTypesConfig(BaseModel):
    """Edge type assignment configuration."""
    method: Literal["uniform", "random", "composite"] = Field(default="uniform")
    uniform: UniformEdgeConfig = Field(default_factory=UniformEdgeConfig)
    random: RandomTypeConfig = Field(default_factory=RandomTypeConfig)
    composite: CompositeEdgeConfig = Field(default_factory=CompositeEdgeConfig)


class DPConfig(BaseModel):
    """DP distribution for a type."""
    mean: float = Field(default=25.0, gt=0)
    pdi: float = Field(default=1.0, ge=1.0, description="Polydispersity index (1.0 = monodisperse)")


class DPDistributionConfig(BaseModel):
    """DP distribution configuration."""
    default: DPConfig = Field(default_factory=DPConfig)
    per_edge_type: dict[str, DPConfig] = Field(default_factory=dict)


class TargetConfig(BaseModel):
    """Target specification for defects/entanglements."""
    enabled: bool = Field(default=False)
    target: int = Field(default=0, ge=0, description="Target count or percentage value")
    target_type: Literal["count", "percentage"] = Field(default="count")


class DefectsConfig(BaseModel):
    """Defect injection configuration."""
    primary_loops: TargetConfig = Field(default_factory=TargetConfig)
    secondary_loops: TargetConfig = Field(default_factory=TargetConfig)


class KinkParams(BaseModel):
    """Parameters for entanglement kink geometry."""
    overshoot: float = Field(default=0.2, ge=0, le=1)
    z_amp: float = Field(default=0.5, ge=0)
    sigma: float = Field(default=0.15, gt=0)


class EntanglementsConfig(BaseModel):
    """Entanglement configuration."""
    enabled: bool = Field(default=False)
    target: int = Field(default=0, ge=0)
    target_type: Literal["count", "percentage"] = Field(default="count")

    # Distribution mode: specify average crosslinks per chain
    # Formula: total_draws = avg_crosslinks_per_chain * 0.5 * num_chains
    avg_crosslinks_per_chain: Optional[float] = Field(
        default=None, ge=0,
        description="Average crosslinks per chain. If set, uses distribution mode with replacement."
    )

    kink_params: KinkParams = Field(default_factory=KinkParams)

    # Spatial placement bias. Default "uniform" gives the legacy
    # behaviour (homogeneous random selection from crossing
    # candidates). Other kinds reweight the draw pool by a spatial
    # function of each candidate's midpoint center.
    #
    # See :func:`topon.assignment.entanglements.compute_bias_weights`
    # for the supported kinds and their params keys:
    #   * region      — uniform inside a sphere, low outside.
    #                    params: center (fractional [0..1] x/y/z),
    #                            radius (fraction of min(dims)),
    #                            strength (in/out density ratio).
    #   * anti_region — depleted inside a sphere, normal outside.
    #                    params: center, radius, strength.
    #   * gradient    — power-law gradient along an axis.
    #                    params: axis ("x"/"y"/"z"), strength (exponent).
    #   * clusters    — gaussian peaks at multiple centers.
    #                    params: centers (list of fractional xyz),
    #                            sigma (fraction of min(dims)),
    #                            strength.
    placement_bias_kind: Literal[
        "uniform", "region", "anti_region", "gradient", "clusters"
    ] = Field(
        default="uniform",
        description="Spatial bias applied to the entanglement-candidate draw.",
    )
    placement_bias_params: dict = Field(
        default_factory=dict,
        description="Parameters consumed by the placement-bias function.",
    )


class GraftConfig(BaseModel):
    """Graft configuration for a single edge type."""
    graft_density: float = Field(default=0.5, ge=0, le=1)
    side_chain_monomer: str = Field(default="PDMS")
    side_chain_dp: int = Field(default=5, ge=1)
    cap_atom: str | None = Field(
        default=None,
        description=(
            "Atom symbol for the unreacted side-chain end. None (default) "
            "means: cap with the same chemistry as the backbone monomer "
            "(e.g. trimethylsilyl for PDMS = add one extra methyl to the "
            "terminal Si). Set to 'H' to leave the terminal atom hydrogen-"
            "capped (RDKit default valence fill). Currently only used by "
            "the atomistic generator; CG ignores this field."
        ),
    )


class GraftsConfig(BaseModel):
    """Grafts configuration."""
    enabled: bool = Field(default=False)
    per_edge_type: dict[str, GraftConfig] = Field(default_factory=dict)


class CopolymerComposition(BaseModel):
    """Single monomer in copolymer composition."""
    monomer: str
    fraction: float = Field(ge=0, le=1)


class CopolymerTypeConfig(BaseModel):
    """Copolymer config for a single edge type."""
    arrangement: Literal["block", "alternating", "random", "gradient"] = Field(default="block")
    composition: list[CopolymerComposition] = Field(default_factory=list)


class CopolymerConfig(BaseModel):
    """Copolymer configuration."""
    enabled: bool = Field(default=False)
    per_edge_type: dict[str, CopolymerTypeConfig] = Field(default_factory=dict)


class AssignmentConfig(BaseModel):
    """Complete assignment configuration."""
    node_types: NodeTypesConfig = Field(default_factory=NodeTypesConfig)
    edge_types: EdgeTypesConfig = Field(default_factory=EdgeTypesConfig)
    dp_distribution: DPDistributionConfig = Field(default_factory=DPDistributionConfig)
    defects: DefectsConfig = Field(default_factory=DefectsConfig)
    entanglements: EntanglementsConfig = Field(default_factory=EntanglementsConfig)
    grafts: GraftsConfig = Field(default_factory=GraftsConfig)
    copolymer: CopolymerConfig = Field(default_factory=CopolymerConfig)


# =============================================================================
# CHEMISTRY CONFIG
# =============================================================================

class NodeMoleculeConfig(BaseModel):
    """Configuration for a node type's molecule."""
    molecule: str = Field(description="SMILES or molecule name (e.g., 'Si', 'POSS')")
    is_end_cap: bool = Field(default=False, description="Whether this is an end-cap molecule")


class MonomerConfig(BaseModel):
    """Monomer definition."""
    smiles: str = Field(description="SMILES string for the repeating unit")
    chain_head: str = Field(default="Si", description="Atom type at chain head")
    chain_tail: str = Field(default="O", description="Atom type at chain tail")


class EdgeChemistryConfig(BaseModel):
    """Chemistry configuration for an edge type."""
    monomer: str = Field(description="Monomer name from monomers library")


class ConnectionConfig(BaseModel):
    """Chain-node connection configuration."""
    auto_bridge: bool = Field(
        default=True,
        description="Automatically insert bridge atom when needed"
    )
    default_bridge_atom: str = Field(
        default="O",
        description="Default bridge atom when auto_bridge is True"
    )


class ChemistryConfig(BaseModel):
    """Chemistry configuration."""
    model_type: Literal["atomistic", "coarse_grained"] = Field(default="coarse_grained")
    target_density: float = Field(default=0.9, gt=0, description="Target density in g/cm³")
    
    node_type_map: dict[str, NodeMoleculeConfig] = Field(
        default={
            "end": NodeMoleculeConfig(molecule="[Si](C)(C)C", is_end_cap=True),
            "A": NodeMoleculeConfig(molecule="Si"),
        },
        description="Map node types to molecules"
    )
    
    edge_type_map: dict[str, EdgeChemistryConfig] = Field(
        default={"A": EdgeChemistryConfig(monomer="PDMS")},
        description="Map edge types to chemistry"
    )
    
    monomers: dict[str, MonomerConfig] = Field(
        default={
            "PDMS": MonomerConfig(smiles="[Si](C)(C)O", chain_head="Si", chain_tail="O"),
            "FPDMS": MonomerConfig(smiles="[Si](C)(CCC(F)(F)F)O", chain_head="Si", chain_tail="O"),
            "Phenyl": MonomerConfig(smiles="[Si](C)(c1ccccc1)O", chain_head="Si", chain_tail="O"),
        },
        description="Monomer library"
    )
    
    connection: ConnectionConfig = Field(default_factory=ConnectionConfig)


# =============================================================================
# OUTPUT CONFIG
# =============================================================================

class OutputConfig(BaseModel):
    """Output configuration."""
    lammps_data: bool = Field(default=True, description="Generate LAMMPS data file")
    lammps_inputs: bool = Field(default=True, description="Generate LAMMPS input scripts")
    visualization: bool = Field(default=True, description="Generate visualization HTML")
    analysis_report: bool = Field(default=True, description="Generate analysis report")
    save_attributed_graph: bool = Field(default=True, description="Save attributed graph as gpickle")
    export_graphml: bool = Field(
        default=False,
        description="Export the dual-graph (chains + entanglement edges) as a GraphML file",
    )
    export_npz: bool = Field(
        default=False,
        description="Export the graph as a compressed NPZ dataset for downstream GNN pipelines",
    )


# =============================================================================
# MAIN CONFIG
# =============================================================================

class ToponConfig(BaseModel):
    """Complete Topon configuration."""
    study: StudyConfig = Field(default_factory=StudyConfig)
    topology: TopologyConfig = Field(default_factory=TopologyConfig)
    assignment: AssignmentConfig = Field(default_factory=AssignmentConfig)
    chemistry: ChemistryConfig = Field(default_factory=ChemistryConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)
    
    model_config = {"extra": "forbid"}
