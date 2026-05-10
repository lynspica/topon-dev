"""
Main pipeline orchestrator for Topon.

Coordinates all stages of the polymer network generation process:
  1. Topology  — generate via C executable or load from .nodes/.edges/.gpickle
  2. Analysis  — graph statistics, defect/entanglement capacity
  3. Assignment — node/edge types, DP, defects, entanglements
  4. Chemistry  — build RDKit molecular structure
  5. Conformation — place atoms, resolve overlaps
  6. Output    — write LAMMPS data + input scripts

Usage::

    from topon.config.loader import ConfigLoader
    from topon.pipeline import Pipeline

    config = ConfigLoader.load("examples/config_cg.json")
    pipe = Pipeline(config)
    pipe.run()

NOTE: The topology stage with ``source="generate"`` uses the C subprocess
generator when ``topology.generator.exe_path`` is set (faster), and the
pure-Python ``PythonTopologyGenerator`` otherwise (no compiler required).
"""

from pathlib import Path
from typing import Optional

import numpy as np

from topon.config.schema import ToponConfig


class Pipeline:
    """
    Main pipeline for polymer network generation.

    Parameters
    ----------
    config : ToponConfig
        Validated config object from :func:`topon.config.loader.ConfigLoader.load`.
    raw_config : dict, optional
        Raw JSON dict for sections not yet covered by the Pydantic schema
        (keys: ``conformation``, ``simulation``, ``experimental``).
    """

    # Conformation defaults (not yet in Pydantic schema)
    _DEFAULT_CONFORMATION = {
        "overlap_cutoff": 0.01,
        "overlap_max_iters": 10,
        "noise_magnitude": 1e-4,
    }

    def __init__(self, config: ToponConfig, raw_config: Optional[dict] = None):
        self.config = config
        self.raw_config = raw_config or {}

        self.graph = None
        self.dims: Optional[np.ndarray] = None
        self.analysis_report: Optional[dict] = None
        self.chemical_space = None
        self._builder = None
        self._assignment_manager = None

        self.output_dir = Path(config.study.output_dir) / config.study.name
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def run(self) -> None:
        """Run the complete pipeline end-to-end."""
        print(f"=== Topon Pipeline: {self.config.study.name} ===")
        print(f"Output directory: {self.output_dir}")
        print()

        self._run_topology_stage()
        self._run_analysis_stage()
        self._run_assignment_stage()
        self._run_chemistry_stage()
        self._run_conformation_stage()
        self._run_output_stage()

        print()
        print("=== Pipeline Complete ===")

    # ------------------------------------------------------------------
    # Stage 1: Topology
    # ------------------------------------------------------------------

    def _run_topology_stage(self) -> None:
        print("--- Stage 1: Topology ---")
        if self.config.topology.source == "generate":
            self._generate_topology()
        else:
            self._load_existing_topology()
        print(f"  Nodes: {self.graph.number_of_nodes()}")
        print(f"  Edges: {self.graph.number_of_edges()}")
        print()

    def _generate_topology(self) -> None:
        import networkx as nx

        from topon.topology.generator import run_generator
        from topon.topology.generator_python import PythonTopologyGenerator
        from topon.topology.loader import (
            infer_dims_from_graph,
            load_graph,
            remove_vacancies,
        )

        gen_cfg = self.config.topology.generator
        topology_dir = self.output_dir / "topology"
        topology_dir.mkdir(parents=True, exist_ok=True)

        if gen_cfg.exe_path:
            # C subprocess path: writes <topology_dir>/output/*.nodes + *.edges,
            # then re-load via the standard loader.
            nodes_path, edges_path = run_generator(
                gen_cfg, topology_dir, exe_path=gen_cfg.exe_path
            )
            self.graph, self.dims = load_graph(
                nodes_path=str(nodes_path),
                edges_path=str(edges_path),
            )
        else:
            # Pure-Python path: in-memory graph, no file round-trip.
            gen = PythonTopologyGenerator(gen_cfg)
            graphs = gen.generate(
                trials=gen_cfg.max_trials,
                max_saves=gen_cfg.max_saves,
            )
            if not graphs:
                raise RuntimeError(
                    f"PythonTopologyGenerator produced no graphs after "
                    f"{gen_cfg.max_trials} trials (constraints may be too "
                    f"strict for lattice {gen_cfg.lattice_size} "
                    f"with degree_distribution={gen_cfg.degree_distribution!r})."
                )
            G = graphs[0]
            if not isinstance(G, nx.MultiGraph):
                G = nx.MultiGraph(G)
            remove_vacancies(G)
            self.graph = G
            self.dims = infer_dims_from_graph(G)

    def _load_existing_topology(self) -> None:
        from topon.topology.loader import load_graph

        files = self.config.topology.existing_files
        if files.gpickle_file:
            self.graph, self.dims = load_graph(gpickle_path=files.gpickle_file)
        elif files.nodes_file and files.edges_file:
            self.graph, self.dims = load_graph(
                nodes_path=files.nodes_file,
                edges_path=files.edges_file,
            )
        else:
            raise ValueError(
                "No topology files specified. Provide gpickle_file or "
                "both nodes_file and edges_file."
            )

    # ------------------------------------------------------------------
    # Stage 2: Analysis
    # ------------------------------------------------------------------

    def _run_analysis_stage(self) -> None:
        print("--- Stage 2: Analysis ---")
        from topon.assignment.manager import AssignmentManager

        self._assignment_manager = AssignmentManager(
            self.graph, self.dims, self.config.assignment
        )
        self.analysis_report = self._assignment_manager.analyze()
        print()

    # ------------------------------------------------------------------
    # Stage 3: Assignment
    # ------------------------------------------------------------------

    def _run_assignment_stage(self) -> None:
        print("--- Stage 3: Assignment ---")
        self._assignment_manager.run()
        print()

    # ------------------------------------------------------------------
    # Stage 4: Chemistry
    # ------------------------------------------------------------------

    def _run_chemistry_stage(self) -> None:
        print("--- Stage 4: Chemistry ---")
        from topon.chemistry.builder import ChemistryBuilder
        from topon.writers import CGWriter, DreidingWriter
        from topon.utils import write_lammps_displacement_file

        self._builder = ChemistryBuilder(
            self.graph, self.dims, self.config.chemistry
        )
        self.chemical_space = self._builder.build()

        chem_dir = self.output_dir / "02_Chemistry"
        chem_dir.mkdir(parents=True, exist_ok=True)
        data_path = str(chem_dir / "system.data")

        model = self.config.chemistry.model_type
        if model == "coarse_grained":
            writer = CGWriter(self.chemical_space, data_path)
            writer.write()
        else:
            writer = DreidingWriter(self.chemical_space, data_path)
            writer.write()

        # Displacement files for conformation stage
        n_atoms = self.chemical_space.GetNumAtoms()
        density = self.config.chemistry.target_density
        vol = n_atoms / density
        scale = (vol / float(np.prod(self.dims))) ** (1.0 / 3.0)

        node_coords: dict[int, tuple] = {}
        for node, atom_ref in self._builder.node_map.items():
            pos = self.graph.nodes[node].get("pos", (0.0, 0.0, 0.0))
            primary_idx = atom_ref[0] if isinstance(atom_ref, (list, tuple)) else atom_ref
            node_coords[primary_idx] = tuple(pos)

        write_lammps_displacement_file(
            node_coords, scale, scale, scale,
            str(chem_dir / "system_nodes.displace"), "nodes"
        )

        # edge_atom_map keys are (u, v, key) tuples from the MultiGraph
        # (see chemistry/builder.py::_build_chains). Iterate via dict keys
        # rather than treating them as sequential int indices into the
        # graph.edges() list.
        chain_coords: dict[int, tuple] = {}
        for (u, v, _key), atom_indices in self._builder.edge_atom_map.items():
            pos_u = np.array(self.graph.nodes[u].get("pos", (0.0, 0.0, 0.0)))
            pos_v = np.array(self.graph.nodes[v].get("pos", (0.0, 0.0, 0.0)))
            vec = pos_v - pos_u
            mic = vec - self.dims * np.round(vec / self.dims)
            for j, atom_idx in enumerate(atom_indices):
                frac = (j + 1) / (len(atom_indices) + 1)
                chain_coords[atom_idx] = tuple(pos_u + frac * mic)

        write_lammps_displacement_file(
            chain_coords, scale, scale, scale,
            str(chem_dir / "system_beads.displace"), "beads"
        )

        node_atom_ids = []
        for atom_ref in self._builder.node_map.values():
            if isinstance(atom_ref, (list, tuple)):
                node_atom_ids.extend(int(i) + 1 for i in atom_ref)
            else:
                node_atom_ids.append(int(atom_ref) + 1)
        node_atom_ids.sort()

        with open(chem_dir / "system.groups", "w") as fh:
            fh.write("# LAMMPS group definitions\n")
            fh.write(f"group nodes id {' '.join(str(x) for x in node_atom_ids)}\n")
            fh.write("group beads subtract all nodes\n")

        settings_path = chem_dir / "system.in.settings"
        if not settings_path.exists():
            with open(settings_path, "w") as fh:
                fh.write("# Force field settings (auto-generated stub)\n")
        print()

    # ------------------------------------------------------------------
    # Stage 5: Conformation
    # ------------------------------------------------------------------

    def _run_conformation_stage(self) -> None:
        print("--- Stage 5: Conformation ---")
        from topon.conformation import ConformationManager

        conf_params = {
            **self._DEFAULT_CONFORMATION,
            **self.raw_config.get("conformation", {}),
        }
        cm = ConformationManager(
            str(self.config.study.output_dir),
            self.config.study.name,
        )
        conformed, roles = cm.apply_displacements("system.data")
        noisy = cm.apply_noise(conformed, magnitude=conf_params["noise_magnitude"])
        cm.resolve_overlaps(
            noisy,
            roles,
            cutoff=conf_params["overlap_cutoff"],
            max_iters=conf_params["overlap_max_iters"],
        )
        print()

    # ------------------------------------------------------------------
    # Stage 6: Output
    # ------------------------------------------------------------------

    def _run_output_stage(self) -> None:
        print("--- Stage 6: Output ---")
        from topon.writers import LammpsInputGenerator

        sim_cfg = self.raw_config.get("simulation", {})
        experimental = self.raw_config.get("experimental", {})
        # LammpsInputGenerator branches on "cg" vs "atomistic" literals
        # (see topon/writers/lammps_inputs.py); the schema's chemistry.model_type
        # uses "coarse_grained" / "atomistic". Map at the call site rather than
        # touching every comparison in the writer.
        model = "cg" if self.config.chemistry.model_type == "coarse_grained" else "atomistic"

        # Pass the BASE output_dir (not self.output_dir which already includes
        # study.name) — LammpsInputGenerator re-appends study.name internally.
        # Matches the ConformationManager call pattern at line 259-263.
        gen = LammpsInputGenerator(
            str(self.config.study.output_dir),
            self.config.study.name,
            config=sim_cfg,
            experimental=experimental,
        )
        gen.write_serial_soft_minimization(
            settings_file="system.in.settings",
            model_type=model,
        )
        gen.write_parallel_production(
            settings_file="system.in.settings",
            model_type=model,
        )
        print(f"  LAMMPS scripts written to: {self.output_dir / '04_Simulation'}")
        print()
