"""
topon.chemistry.kg
==================
Stage 2 backend: Kremer-Grest (coarse-grained) chemistry.

Responsible for:
  - Assigning bead types (J = junction, A = backbone, G = graft)
  - Building the RDKit mol from a NetworkGraph
  - Writing LAMMPS data file via CGWriter
  - Computing box scale from bead density
  - Writing displacement files for Stage 3

Status: logic lives inline in topon/workflows/cg_network.py (Stage 2 section).
A standalone assign(network_graph, config) → CGSystem function will be extracted
here once the workflow regression tests are stable.

Parameters (from config JSON)
------------------------------
  chemistry.degree_of_polymerization : int   — beads per chain
  chemistry.bead_density             : float — LJ units
  assignment.grafts                  : dict  — optional graft config
  assignment.entanglements           : dict  — optional entanglement config
  simulation.include_angles          : bool  — whether to include angle terms
"""
