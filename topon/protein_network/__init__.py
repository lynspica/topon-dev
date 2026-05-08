"""MARTINI 3 protein-network generator for topon.

Sequence-driven coarse-grained polymer-network builder, peer to topon.simbox
and topon.singlechain. Mirrors the topro CHARMM atomistic generator's two-stage
shape: BFM lattice topology -> JSON -> sequence + MARTINI FF + water -> LAMMPS.
"""
