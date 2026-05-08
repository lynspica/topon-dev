"""
Monomer Sequence Generation Logic
"""
import random
import math
import warnings
from typing import List

def generate_monomer_sequence(dp: int, config: dict, default_monomer: str = 'A') -> List[str]:
    """
    Generate a sequence of monomers based on copolymer configuration.
    
    Args:
        dp: Degree of Polymerization (total length)
        config: CopolymerTypeConfig dict (arrangement, composition)
        default_monomer: Default monomer type if config is invalid/empty
        
    Returns:
        List of monomer/bead type strings
    """
    if not config or not config.get('composition'):
        return [default_monomer] * dp

    arrangement = config.get('arrangement', 'random')
    comps = config.get('composition', [])
    
    # Normalize fractions if needed? 
    # Schema says floats. Let's assume they sum to approx 1.0 or weights.
    # If using random.choices, weights are fine.
    
    monomers = [c['monomer'] for c in comps]
    fractions = [c['fraction'] for c in comps]
    
    if not monomers:
        return [default_monomer] * dp

    if arrangement == 'random':
        return random.choices(monomers, weights=fractions, k=dp)
        
    elif arrangement == 'block':
        # Split DP according to fractions
        # e.g. DP=20, A=0.5, B=0.5 -> 10 A, 10 B
        seq = []
        current_count = 0
        
        # Calculate precise counts to sum to DP
        counts = []
        rem = dp
        total_frac = sum(fractions)
        
        for i, f in enumerate(fractions):
            if i == len(fractions) - 1:
                count = rem # Last one takes remainder to ensure sum=DP
            else:
                count = int(round(dp * (f / total_frac)))
                rem -= count
            counts.append(count)
            
        for m, c in zip(monomers, counts):
            seq.extend([m] * c)
            
        return seq[:dp] # Safety clamp
        
    elif arrangement == 'alternating':
        # A, B, C, A, B, C... irrelevant of fractions usually, but maybe weighted alternating?
        # Standard alternating ignores fractions usually.
        seq = []
        n_types = len(monomers)
        for i in range(dp):
            seq.append(monomers[i % n_types])
        return seq

    elif arrangement == 'gradient':
        # Linear gradient: monomer probability shifts smoothly from the first
        # component to the last along the chain length.
        # Position 0 → weight biased toward monomers[0]
        # Position dp-1 → weight biased toward monomers[-1]
        seq = []
        n = len(monomers)
        for i in range(dp):
            t = i / max(dp - 1, 1)  # 0.0 … 1.0
            # Linear interpolation weights: starts at fractions[0], ends at fractions[-1]
            weights = []
            for j, f in enumerate(fractions):
                pivot = j / max(n - 1, 1)
                # Weight peaks at the monomer's pivot position along the gradient
                w = max(0.0, 1.0 - abs(t - pivot) * n)
                weights.append(w * f)
            if sum(weights) == 0:
                weights = fractions  # fallback if all zero
            seq.append(random.choices(monomers, weights=weights, k=1)[0])
        return seq

    warnings.warn(
        f"generate_monomer_sequence: unknown arrangement {arrangement!r}; "
        f"returning homopolymer of default_monomer={default_monomer!r}.",
        RuntimeWarning,
        stacklevel=2,
    )
    return [default_monomer] * dp
