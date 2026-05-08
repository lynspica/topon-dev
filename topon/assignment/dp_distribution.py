"""
DP (Degree of Polymerization) distribution assignment for Topon.

Assigns DP values to edges based on:
- Default mean and PDI
- Per-edge-type specific distributions
"""

import random
import math
from typing import Optional
import networkx as nx

from topon.config.schema import DPDistributionConfig, DPConfig


def assign_dp(G: nx.MultiGraph, config: DPDistributionConfig) -> None:
    """
    Assign DP values to all edges in the graph.
    
    Args:
        G: Graph to modify (in place).
        config: DP distribution configuration.
    """
    # Count edges by type
    edge_types = {}
    for u, v, key, data in G.edges(keys=True, data=True):
        edge_type = data.get("edge_type", "default")
        if edge_type not in edge_types:
            edge_types[edge_type] = []
        edge_types[edge_type].append((u, v, key))
    
    # Assign DP for each type
    for edge_type, edges in edge_types.items():
        # Get config for this type (or use default)
        if edge_type in config.per_edge_type:
            dp_config = config.per_edge_type[edge_type]
        else:
            dp_config = config.default
        
        # Generate DP values
        mean_dp = dp_config.mean
        pdi = dp_config.pdi
        
        dp_values = generate_dp_distribution(len(edges), mean_dp, pdi)
        
        # Assign to edges
        for (u, v, key), dp in zip(edges, dp_values):
            G.edges[u, v, key]["dp"] = dp
    
    # Report stats
    all_dps = [data.get("dp", 0) for u, v, key, data in G.edges(keys=True, data=True)]
    if all_dps:
        print(f"    DP range: {min(all_dps)} - {max(all_dps)}")
        print(f"    DP mean: {sum(all_dps) / len(all_dps):.1f}")


def generate_dp_distribution(n: int, mean_dp: float, pdi: float) -> list[int]:
    """
    Generate n DP values with target mean and PDI.
    
    For PDI = 1 (monodisperse): all values are the same
    For PDI > 1: uses Schulz-Zimm distribution approximation
    
    Args:
        n: Number of values to generate.
        mean_dp: Target mean DP.
        pdi: Polydispersity index (Mw/Mn), must be >= 1.
        
    Returns:
        List of integer DP values.
    """
    if pdi < 1.0:
        pdi = 1.0
    
    if pdi == 1.0:
        # Monodisperse - all same value
        return [max(1, int(round(mean_dp)))] * n
    
    # Schulz-Zimm distribution
    # Shape parameter k = 1 / (PDI - 1)
    # Scale parameter θ = mean_dp * (PDI - 1)
    k = 1.0 / (pdi - 1.0)
    theta = mean_dp / k
    
    dp_values = []
    for _ in range(n):
        # Gamma distribution
        value = random.gammavariate(k, theta)
        # Ensure minimum DP of 1
        dp = max(1, int(round(value)))
        dp_values.append(dp)
    
    # Adjust to match target mean more closely
    current_mean = sum(dp_values) / n if n > 0 else mean_dp
    if abs(current_mean - mean_dp) > 1:
        adjustment = mean_dp - current_mean
        for i in range(n):
            dp_values[i] = max(1, int(round(dp_values[i] + adjustment)))
    
    return dp_values
