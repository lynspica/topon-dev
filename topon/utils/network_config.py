# polymer_world/utils/network_config.py

# Dictionary of Monomer Repeating Units
# The SMILES string represents the -[Si(R1)(R2)-O]- repeating unit.
PENDANT_GROUP_SMILES = {
    # Base Polymers
    'PDMS': '[Si](C)(C)O',
    'FPDMS': '[Si](C)(CCC(F)(F)F)O',  # Trifluoropropyl
    
    # High Strength / Aromatic
    'Phenyl': '[Si](C)(c1ccccc1)O',
    'Diphenyl': '[Si](c1ccccc1)(c1ccccc1)O',
    'Methylphenyl': '[Si](C)(c1ccc(C)cc1)O',
    'Chlorophenyl': '[Si](C)(c1ccc(Cl)cc1)O',
    'Trifluoromethylphenyl': '[Si](C)(c1ccc(C(F)(F)F)cc1)O',
    
    # New Additions (for Bottlebrushes)
    'Polyaramid': 'Nc1ccc(cc1)NC(=O)c1ccc(cc1)C(=O)',
    'PEO': 'CCO',
    'PS': 'CC(c1ccccc1)',

    # Fluorinated
    'Trifluoropropyl': '[Si](C)(CCC(F)(F)F)O',
    'Tridecafluorooctyl': '[Si](C)(CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)O',
    
    # Reactive / Functional
    'Vinyl': '[Si](C)(C=C)O',
    'Methacryloxypropyl': '[Si](C)(CCCOC(=O)C(C)=C)O',
    'Epoxy': '[Si](C)(CCCOC1OC1)O',
    'Aminopropyl': '[Si](C)(CCCN)O',
    'Carboxypropyl': '[Si](C)(CCC(=O)O)O',
    'Hydroxypropyl': '[Si](C)(CCCO)O',
    'Mercaptopropyl': '[Si](C)(CCCS)O',
    'Cyanopropyl': '[Si](C)(CCC#N)O',
    
    # Alkyl
    'Cyclohexyl': '[Si](C)(C1CCCCC1)O',
    'Hexyl': '[Si](C)(CCCCCC)O',
    'Hexadecyl': '[Si](C)(CCCCCCCCCCCCCCCC)O',
    'Isobutyl': '[Si](C)(CC(C)C)O'
}

# SMILES for end-capping units in uncrosslinked systems
END_CAP_START_SMILES = "C[Si](C)(C)O"
END_CAP_END_SMILES = "[Si](C)(C)C"

# Simulation Parameters
DREIDING_PARAM_FILE = "DreidingX6parameters.txt"
TEMP_DATA_FILE = "temp_datafile.data"
GROUP_DEFS_FILE = "set.groups"

DEFAULT_TARGET_DENSITY = 0.9
BOX_SCALING_FACTORS = [1.001, 1.001, 1.001]