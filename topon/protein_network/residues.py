"""MARTINI 3 protein residue table.

Auto-generated from nat_pro.itp by tools/extract_residues_from_itp.py.
Edit the extractor and re-run; do not hand-edit this file.

Detected repeating block: 'GGRPSDSYGAPGGGN' x 18.
Bead types referenced (14): P2, Q5, SC3, SP1, SP2, SP2a, SP5, SQ3p, SQ5n, TC3, TC4, TC5, TN6, TP1.
"""
from __future__ import annotations

SOURCE_ITP = 'nat_pro.itp'
REFERENCE_BLOCK = 'GGRPSDSYGAPGGGN'
REFERENCE_N_REPEATS = 18

THREE_TO_ONE = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLN': 'Q',
    'GLU': 'E',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

# Canonical mid-chain bead pattern per residue: [(atom_name, bead_type, charge), ...]
RESIDUES: dict[str, dict] = {
    'ALA': {
        "beads": [('BB', 'SP2', 0.0), ('SC1', 'TC3', 0.0)],
        "intra_bonds": [],
        "intra_constraints": [('BB', 'SC1', 1, 0.27)],
        "intra_angles": [],
        "intra_dihedrals_proper": [],
        "intra_dihedrals_improper": [],
        "intra_exclusions": [],
        "occurrences_in_reference": [10, 25, 40, 55, 70, 85, 100, 115, 130, 145, 160, 175, 190, 205, 220, 235, 250, 265],
    },
    'ARG': {
        "beads": [('BB', 'P2', 0.0), ('SC1', 'SC3', 0.0), ('SC2', 'SQ3p', 1.0)],
        "intra_bonds": [('BB', 'SC1', 1, 0.33, 5000.0), ('SC1', 'SC2', 1, 0.38, 5000.0)],
        "intra_constraints": [],
        "intra_angles": [('BB', 'SC1', 'SC2', 2, 180.0, 25.0)],
        "intra_dihedrals_proper": [],
        "intra_dihedrals_improper": [],
        "intra_exclusions": [],
        "occurrences_in_reference": [3, 18, 33, 48, 63, 78, 93, 108, 123, 138, 153, 168, 183, 198, 213, 228, 243, 258],
    },
    'ASN': {
        "beads": [('BB', 'P2', 0.0), ('SC1', 'SP5', 0.0)],
        "intra_bonds": [('BB', 'SC1', 1, 0.352, 5000.0)],
        "intra_constraints": [],
        "intra_angles": [],
        "intra_dihedrals_proper": [],
        "intra_dihedrals_improper": [],
        "intra_exclusions": [],
        "occurrences_in_reference": [15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270],
    },
    'ASP': {
        "beads": [('BB', 'P2', 0.0), ('SC1', 'SQ5n', -1.0)],
        "intra_bonds": [('BB', 'SC1', 1, 0.352, 7500.0)],
        "intra_constraints": [],
        "intra_angles": [],
        "intra_dihedrals_proper": [],
        "intra_dihedrals_improper": [],
        "intra_exclusions": [],
        "occurrences_in_reference": [6, 21, 36, 51, 66, 81, 96, 111, 126, 141, 156, 171, 186, 201, 216, 231, 246, 261],
    },
    'GLY': {
        "beads": [('BB', 'SP1', 0.0)],
        "intra_bonds": [],
        "intra_constraints": [],
        "intra_angles": [],
        "intra_dihedrals_proper": [],
        "intra_dihedrals_improper": [],
        "intra_exclusions": [],
        "occurrences_in_reference": [1, 2, 9, 12, 13, 14, 16, 17, 24, 27, 28, 29, 31, 32, 39, 42, 43, 44, 46, 47, 54, 57, 58, 59, 61, 62, 69, 72, 73, 74, 76, 77, 84, 87, 88, 89, 91, 92, 99, 102, 103, 104, 106, 107, 114, 117, 118, 119, 121, 122, 129, 132, 133, 134, 136, 137, 144, 147, 148, 149, 151, 152, 159, 162, 163, 164, 166, 167, 174, 177, 178, 179, 181, 182, 189, 192, 193, 194, 196, 197, 204, 207, 208, 209, 211, 212, 219, 222, 223, 224, 226, 227, 234, 237, 238, 239, 241, 242, 249, 252, 253, 254, 256, 257, 264, 267, 268, 269],
    },
    'PRO': {
        "beads": [('BB', 'SP2a', 0.0), ('SC1', 'SC3', 0.0)],
        "intra_bonds": [('BB', 'SC1', 1, 0.33, 7500.0)],
        "intra_constraints": [],
        "intra_angles": [],
        "intra_dihedrals_proper": [],
        "intra_dihedrals_improper": [],
        "intra_exclusions": [],
        "occurrences_in_reference": [4, 11, 19, 26, 34, 41, 49, 56, 64, 71, 79, 86, 94, 101, 109, 116, 124, 131, 139, 146, 154, 161, 169, 176, 184, 191, 199, 206, 214, 221, 229, 236, 244, 251, 259, 266],
    },
    'SER': {
        "beads": [('BB', 'P2', 0.0), ('SC1', 'TP1', 0.0)],
        "intra_bonds": [('BB', 'SC1', 1, 0.287, 7500.0)],
        "intra_constraints": [],
        "intra_angles": [],
        "intra_dihedrals_proper": [],
        "intra_dihedrals_improper": [],
        "intra_exclusions": [],
        "occurrences_in_reference": [5, 7, 20, 22, 35, 37, 50, 52, 65, 67, 80, 82, 95, 97, 110, 112, 125, 127, 140, 142, 155, 157, 170, 172, 185, 187, 200, 202, 215, 217, 230, 232, 245, 247, 260, 262],
    },
    'TYR': {
        "beads": [('BB', 'P2', 0.0), ('SC1', 'TC4', 0.0), ('SC2', 'TC5', 0.0), ('SC3', 'TC5', 0.0), ('SC4', 'TN6', 0.0)],
        "intra_bonds": [('BB', 'SC1', 1, 0.325, 5000.0)],
        "intra_constraints": [('SC1', 'SC2', 1, 0.3), ('SC1', 'SC3', 1, 0.3), ('SC2', 'SC3', 1, 0.3), ('SC2', 'SC4', 1, 0.285), ('SC3', 'SC4', 1, 0.285)],
        "intra_angles": [('BB', 'SC1', 'SC2', 2, 120.0, 60.0), ('BB', 'SC1', 'SC3', 2, 120.0, 60.0)],
        "intra_dihedrals_proper": [],
        "intra_dihedrals_improper": [('SC4', 'SC2', 'SC3', 'SC1', 2, 180.0, 50.0)],
        "intra_exclusions": [('BB', 'SC1', 'SC2', 'SC3', 'SC4'), ('SC1', 'SC2', 'SC3', 'SC4'), ('SC2', 'SC3', 'SC4'), ('SC3', 'SC4')],
        "occurrences_in_reference": [8, 23, 38, 53, 68, 83, 98, 113, 128, 143, 158, 173, 188, 203, 218, 233, 248, 263],
    },
}

TERMINAL_PATCHES = {
    "N_term": [('BB', 'Q5', 1.0)],
    "C_term": [('BB', 'Q5', -1.0), ('SC1', 'SP5', 0.0)],
}

# Backbone BB-BB bonds extracted as (resname_a, resname_b, funct, length_nm, k_kJ_mol_nm2).
BACKBONE_BB_BONDS = [
    ('ALA', 'PRO', 1, 0.36, 10000.0),
    ('ARG', 'PRO', 1, 0.36, 10000.0),
    ('ASN', 'GLY', 1, 0.36, 8000.0),
    ('ASP', 'SER', 1, 0.36, 8000.0),
    ('GLY', 'ALA', 1, 0.36, 8000.0),
    ('GLY', 'ARG', 1, 0.36, 8000.0),
    ('GLY', 'ASN', 1, 0.36, 8000.0),
    ('GLY', 'GLY', 1, 0.36, 8000.0),
    ('PRO', 'GLY', 1, 0.305, 10000.0),
    ('PRO', 'SER', 1, 0.305, 10000.0),
    ('SER', 'ASP', 1, 0.36, 8000.0),
    ('SER', 'TYR', 1, 0.36, 8000.0),
    ('TYR', 'GLY', 1, 0.36, 8000.0),
]

BACKBONE_BBB_ANGLE = (10, 137.0, 25.0)  # (funct, angle_deg, k_kJ_mol_rad2)

# Backbone BBBB dihedrals bucketed by polyply-comment label (e.g. 'BBBB', 'GGGX').
BACKBONE_BBBB_DIHEDRALS: dict[str, list[tuple]] = {
    'BBBB': [
        (9, 60.0, 2.8, 1),
        (9, 150.0, -0.6, 1),
        (9, 130.0, -1.2, 2),
    ],
    'GGGX': [
        (9, 160.0, 0.8, 1),
        (9, -160.0, 0.8, 1),
        (9, 0.0, 1.2, 2),
    ],
    'GGXG': [
        (9, -160.0, 0.8, 1),
        (9, -80.0, 0.8, 2),
        (9, 0.0, 1.2, 2),
    ],
    'GXGG': [
        (9, 60.0, 2.0, 1),
        (9, -100.0, 0.8, 1),
    ],
    'XGGG': [
        (9, 160.0, 0.8, 1),
        (9, -160.0, 0.8, 1),
        (9, 0.0, 1.2, 2),
    ],
}

