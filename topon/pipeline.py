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


def _local_perp_unit(backbone_xyz, k: int, fallback_unit, rand_vec) -> np.ndarray:
    """Return a unit vector for a graft to stick out perpendicular to the
    backbone at index ``k``, biased *outward* on kinked sections.

    Algorithm:
      1. Tangent T = central difference (P[k+1] - P[k-1]); endpoint cases
         use one-sided differences. Falls back to the per-edge chord unit
         when the backbone has fewer than 2 atoms.
      2. Curvature direction K = P[k+1] - 2*P[k] + P[k-1] (points *into*
         the bend, i.e. toward the centre of curvature).
      3. If K is large enough (the chain is genuinely curving — entangled
         kinks always trip this), the graft direction is the **outward**
         normal -K, projected to be perpendicular to T. This places grafts
         on the convex side of the bend so they never dive into the chain.
      4. If K is tiny (a straight backbone), fall back to the per-edge
         random ``rand_vec`` projected perpendicular to T (same answer as
         the pre-2026-05-11 chord-perpendicular behaviour on straight
         chains).
    """
    n = len(backbone_xyz)
    if n < 2:
        local_unit = fallback_unit
        curv = None
    else:
        if k <= 0:
            tangent = np.asarray(backbone_xyz[1]) - np.asarray(backbone_xyz[0])
            curv = None  # no second-derivative at endpoint
        elif k >= n - 1:
            tangent = np.asarray(backbone_xyz[-1]) - np.asarray(backbone_xyz[-2])
            curv = None
        else:
            tangent = np.asarray(backbone_xyz[k + 1]) - np.asarray(backbone_xyz[k - 1])
            curv = (
                np.asarray(backbone_xyz[k + 1])
                - 2.0 * np.asarray(backbone_xyz[k])
                + np.asarray(backbone_xyz[k - 1])
            )
        tn = float(np.linalg.norm(tangent))
        local_unit = tangent / tn if tn > 1e-9 else fallback_unit

    # Outward normal from curvature, if the chain bends enough at this point.
    if curv is not None:
        cn = float(np.linalg.norm(curv))
        if cn > 1e-3:
            outward = -curv / cn
            # Project out the tangent component so the graft is strictly perp.
            outward = outward - np.dot(outward, local_unit) * local_unit
            on = float(np.linalg.norm(outward))
            if on > 1e-6:
                return outward / on

    # Straight (or near-straight) backbone: random perp from rand_vec.
    perp = np.cross(local_unit, rand_vec)
    pn = float(np.linalg.norm(perp))
    if pn < 1e-6:
        perp = np.cross(local_unit, np.array([1.0, 0.0, 0.0]))
        pn = float(np.linalg.norm(perp))
        if pn < 1e-6:
            perp = np.cross(local_unit, np.array([0.0, 1.0, 0.0]))
            pn = float(np.linalg.norm(perp))
    return perp / (pn + 1e-12)


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

    def run_from_graph(
        self,
        graph,
        dims,
    ) -> None:
        """Run only chemistry + conformation + output stages.

        Use when the graph is already fully prepared (e.g. loaded from a
        graphml / npz dual-graph file via
        :func:`topon.topology.loader.load_graphml` /
        :func:`topon.topology.loader.load_npz`) and stages 1-3 (topology,
        analysis, assignment) should be skipped.

        For deterministic output across different sources of the same
        topology (e.g. graphml-load vs npz-load), seed both
        :mod:`random` and :mod:`numpy.random` before calling this. The
        chemistry stage uses ``np.random.randn`` for graft-perp
        directions, and the conformation stage uses noise from
        ``numpy.random`` -- without seeding, byte-equivalence cannot be
        guaranteed.

        Args:
            graph: NetworkX MultiGraph with crosslink nodes (``pos``) and
                chain edges (``dp``, ``entangled_with``,
                ``entanglement_count``). Edge / node ``type`` attributes
                are optional and default to ``"A"``.
            dims: Box dimensions array ``[Lx, Ly, Lz]``.
        """
        self.graph = graph
        self.dims = (
            np.asarray(dims, dtype=float)
            if not isinstance(dims, np.ndarray) else dims.astype(float)
        )
        print(f"=== Topon Pipeline (rebuild from graph): "
              f"{self.config.study.name} ===")
        print(f"Output directory: {self.output_dir}")
        print(f"  Skipping stages 1-3 (graph supplied directly).")
        print(f"  Nodes: {self.graph.number_of_nodes()}, "
              f"Edges: {self.graph.number_of_edges()}")
        print()

        self._run_chemistry_stage()
        self._run_conformation_stage()
        self._run_output_stage()

        print()
        print("=== Pipeline Complete (rebuild) ===")

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
        from topon.utils.network_helpers import (
            generate_approximate_side_chain_coords,
        )
        from topon.conformation.entanglement.realize import (
            entangled_backbone_paths,
        )

        self._builder = ChemistryBuilder(
            self.graph, self.dims, self.config.chemistry
        )
        self.chemical_space = self._builder.build()

        chem_dir = self.output_dir / "02_Chemistry"
        chem_dir.mkdir(parents=True, exist_ok=True)
        data_path = str(chem_dir / "system.data")

        model = self.config.chemistry.model_type
        density = self.config.chemistry.target_density

        if model == "coarse_grained":
            # Honour raw_config's simulation.include_angles flag (default True
            # to preserve historic behaviour). Mirrors the topon.workflows.
            # cg_network call pattern and the schema doc-string in
            # topon/chemistry/kg/__init__.py.
            sim_cfg_for_writer = self.raw_config.get("simulation", {})
            include_angles = sim_cfg_for_writer.get("include_angles", True)
            writer = CGWriter(
                self.chemical_space, data_path,
                include_angles=include_angles,
            )
            writer.write()
            # Count-based volume for CG (matches v21 cg_network reference).
            n_atoms = self.chemical_space.GetNumAtoms()
            vol = n_atoms / density
        else:
            # Atomistic: mirror topon.workflows.atomistic_network's tail —
            # the canonical path that produces healthy LAMMPS stage 2/3
            # output (v21/v43 reference). ChemistryBuilder.build() returns
            # a heavy-atom-only RWMol; we Sanitize -> AddHs -> Gasteiger so
            # the data file (a) has H_ atom-type rows DREIDING needs, and
            # (b) is charge-neutral so PPPM auto-gewald doesn't crash. AddHs
            # preserves heavy-atom indices, so _builder.node_map and
            # edge_atom_map remain valid. Sanitize can fail on demos that
            # produce over-valent atoms (e.g. defect demo's degree-6 Si);
            # fall back to writing the heavy-atom mol uncharged in that
            # case — system.data is then DREIDING-incomplete (no H_) but
            # the chemistry stage still completes for inspection.
            from rdkit import Chem
            from rdkit.Chem import AllChem
            try:
                try:
                    Chem.SanitizeMol(self.chemical_space)
                except Chem.AtomValenceException:
                    # Defect demos can produce degree-6 Si junctions that
                    # exceed RDKit's permitted valence (max 6 by table, but
                    # the strict check trips on Si@6 with no charge). Skip
                    # just the valence-property check; the rest of sanitize
                    # (kekulize, ring-find, etc.) still runs. AddHs then
                    # assigns 0 implicit H to those Si atoms.
                    Chem.SanitizeMol(
                        self.chemical_space,
                        sanitizeOps=(
                            Chem.SanitizeFlags.SANITIZE_ALL
                            ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                        ),
                    )
                mol_h = Chem.AddHs(self.chemical_space)
                AllChem.ComputeGasteigerCharges(mol_h)
                # Over-valent atoms (defect's degree-6 Si) make Gasteiger
                # emit NaN for those atoms and their neighbours. Scrub to
                # 0 so the LAMMPS data file is numeric; the residual net
                # charge logged below is usually within PPPM's tolerance.
                import math
                nan_count = 0
                for atom in mol_h.GetAtoms():
                    if atom.HasProp("_GasteigerCharge"):
                        q = atom.GetDoubleProp("_GasteigerCharge")
                        if math.isnan(q) or math.isinf(q):
                            atom.SetDoubleProp("_GasteigerCharge", 0.0)
                            nan_count += 1
                if nan_count:
                    print(f"  [WARN] Gasteiger NaN/Inf on {nan_count} atoms "
                          f"(over-valent neighbours); zeroed.")

                # Background charge neutralization: redistribute residual net
                # charge uniformly across all atoms with a valid Gasteiger
                # value. Two things make this load-bearing:
                #   (1) The defect demo has over-valent Si (degree 5-6 from
                #       parallel-edge defects). NaN-zeroing those atoms can
                #       leave several e residual; PPPM then prints "System
                #       is not charge neutral" and runs slowly.
                #   (2) Gasteiger output itself has a tiny non-zero residual
                #       (~1e-12 e) from finite-precision iteration on every
                #       molecule. Cheap to scrub here.
                valid_atoms = [
                    a for a in mol_h.GetAtoms() if a.HasProp("_GasteigerCharge")
                ]
                total_q = sum(
                    a.GetDoubleProp("_GasteigerCharge") for a in valid_atoms
                )
                if valid_atoms and abs(total_q) > 1e-6:
                    delta = -total_q / len(valid_atoms)
                    for a in valid_atoms:
                        a.SetDoubleProp(
                            "_GasteigerCharge",
                            a.GetDoubleProp("_GasteigerCharge") + delta,
                        )
                    print(f"  Charge-neutralized: spread {-total_q:+.4f} e "
                          f"across {len(valid_atoms)} atoms "
                          f"(delta = {delta:+.2e} e/atom)")

                self.chemical_space = mol_h
                # Mass-based volume (matches canonical workflow at v43).
                mass = sum(a.GetMass() for a in mol_h.GetAtoms())
                vol = (mass / density) * 1.66054  # A^3 / Da at g/cm^3
                writer = DreidingWriter(mol_h, data_path, use_charges=True)
            except Exception as exc:
                print(f"  [WARN] AddHs/Gasteiger skipped ({exc}); "
                      f"writing heavy-atom data file uncharged.")
                mol_h = self.chemical_space
                n_atoms = mol_h.GetNumAtoms()
                vol = n_atoms / density
                writer = DreidingWriter(mol_h, data_path, use_charges=False)
            writer.write()

        scale = (vol / float(np.prod(self.dims))) ** (1.0 / 3.0)
        sx = sy = sz = scale

        # Nodes displacement (both CG and atomistic). POSS node_map values
        # can be lists (cage); preserve the isinstance branch.
        node_coords: dict[int, tuple] = {}
        for node, atom_ref in self._builder.node_map.items():
            pos = self.graph.nodes[node].get("pos", (0.0, 0.0, 0.0))
            primary_idx = atom_ref[0] if isinstance(atom_ref, (list, tuple)) else atom_ref
            node_coords[primary_idx] = tuple(pos)

        write_lammps_displacement_file(
            node_coords, sx, sy, sz,
            str(chem_dir / "system_nodes.displace"), "nodes"
        )

        # Entangled edges' paths, both branches, computed once. The method
        # is the config's: "waypoint" (default) draws each pair together
        # with a prescribed winding count; "kink" is the legacy Gaussian
        # bump. Edges not in this dict are linear, as they always were.
        ent_cfg = self.config.assignment.entanglements
        ent_paths = entangled_backbone_paths(
            self.graph, self.dims, self._builder.edge_atom_map,
            method=ent_cfg.method,
            kink_params=ent_cfg.kink_params.model_dump(),
        )

        if model == "coarse_grained":
            # CG: same backbone + grafts loop as atomistic (entanglement-aware
            # winding for backbone, perpendicular placement with 3-way length
            # cap for grafts). CG doesn't have pendant/H passes — graft beads
            # come directly from `_builder.graft_atom_map`.
            backbone_coords: dict[int, tuple] = {}
            graft_coords: dict[int, tuple] = {}
            graft_atom_map = getattr(self._builder, "graft_atom_map", {}) or {}
            ext_factor = 0.5

            for (u, v, key), atoms in self._builder.edge_atom_map.items():
                data = self.graph[u][v][key]
                pos_u = np.array(self.graph.nodes[u].get("pos", (0.0, 0.0, 0.0)))
                pos_v = np.array(self.graph.nodes[v].get("pos", (0.0, 0.0, 0.0)))
                vec = pos_v - pos_u
                mic = vec - self.dims * np.round(vec / self.dims)
                edge_len = float(np.linalg.norm(mic))
                unit_vec = mic / (edge_len + 1e-9)
                rand_vec = np.random.randn(3)
                perp = np.cross(unit_vec, rand_vec)
                if np.linalg.norm(perp) < 1e-6:
                    perp = np.cross(unit_vec, np.array([1.0, 0.0, 0.0]))
                perp_unit = perp / (np.linalg.norm(perp) + 1e-9)

                backbone_xyz = ent_paths.get((u, v, key))
                if backbone_xyz is None:
                    backbone_xyz = []
                    for j in range(len(atoms)):
                        frac = (j + 1) / (len(atoms) + 1)
                        backbone_xyz.append(pos_u + frac * mic)

                for j, a_idx in enumerate(atoms):
                    backbone_coords[a_idx] = tuple(backbone_xyz[j])

                # Grafts (CG): perp to the local backbone *tangent* at the
                # anchor (kinked backbones twist; chord-perp grafts would
                # otherwise dive back into the chain on entangled edges).
                edge_grafts = graft_atom_map.get((u, v, key))
                if edge_grafts:
                    lattice_spacing = (
                        float(min(self.dims)) if self.dims is not None else edge_len
                    )
                    for frac, g_atoms in edge_grafts:
                        k_float = frac * (len(atoms) + 1) - 1
                        k = max(0, min(int(round(k_float)), len(backbone_xyz) - 1))
                        anchor_pos = backbone_xyz[k]
                        local_perp = _local_perp_unit(
                            backbone_xyz, k, unit_vec, rand_vec
                        )
                        graft_dp_eff = max(len(g_atoms), 1)
                        backbone_dp = max(len(atoms), 1)
                        eff_factor = min(
                            ext_factor,
                            graft_dp_eff / backbone_dp,
                            0.5 * lattice_spacing / max(edge_len, 1e-9),
                        )
                        graft_vec = local_perp * (edge_len * eff_factor)
                        for m, g_idx in enumerate(g_atoms):
                            g_frac = (m + 1) / len(g_atoms)
                            graft_coords[g_idx] = tuple(anchor_pos + g_frac * graft_vec)

            write_lammps_displacement_file(
                backbone_coords, sx, sy, sz,
                str(chem_dir / "system_backbone.displace"), "backbone"
            )
            write_lammps_displacement_file(
                graft_coords, sx, sy, sz,
                str(chem_dir / "system_grafts.displace"), "grafts"
            )
        else:
            # Atomistic: backbone + grafts + pendant + hydrogens. Backbone
            # path consults `entangled_with` on the (u, v, key) edge data
            # for kinked-chain placement (v21.1 N+2 fix). graft_atom_map is
            # currently only populated by ChemistryBuilder._build_chain_cg
            # — for atomistic it's empty, so system_grafts.displace will be
            # empty and graft side-chain atoms instead get coords from the
            # pendant pass (neighbor propagation through mol_h).
            backbone_coords: dict[int, tuple] = {}
            graft_coords: dict[int, tuple] = {}
            graft_atom_map = getattr(self._builder, "graft_atom_map", {}) or {}
            ext_factor = 0.5  # canonical workflow default

            for (u, v, key), atoms in self._builder.edge_atom_map.items():
                data = self.graph[u][v][key]
                pos_u = np.array(self.graph.nodes[u].get("pos", (0.0, 0.0, 0.0)))
                pos_v = np.array(self.graph.nodes[v].get("pos", (0.0, 0.0, 0.0)))
                vec = pos_v - pos_u
                mic = vec - self.dims * np.round(vec / self.dims)
                edge_len = float(np.linalg.norm(mic))
                unit_vec = mic / (edge_len + 1e-9)
                rand_vec = np.random.randn(3)
                perp = np.cross(unit_vec, rand_vec)
                if np.linalg.norm(perp) < 1e-6:
                    perp = np.cross(unit_vec, np.array([1.0, 0.0, 0.0]))
                perp_unit = perp / (np.linalg.norm(perp) + 1e-9)

                backbone_xyz = ent_paths.get((u, v, key))
                if backbone_xyz is None:
                    backbone_xyz = []
                    for j in range(len(atoms)):
                        frac = (j + 1) / (len(atoms) + 1)
                        backbone_xyz.append(pos_u + frac * mic)

                for j, a_idx in enumerate(atoms):
                    backbone_coords[a_idx] = tuple(backbone_xyz[j])

                # Grafts: place perpendicular to the local backbone *tangent*
                # at the anchor (per-anchor finite difference, not per-edge
                # chord). On a kinked entangled chain the local tangent
                # curves away from the chord, so chord-perp grafts would
                # otherwise dive back into the chain. Length is the minimum
                # of three competing constraints:
                #   (a) extension_factor (default 0.5) — half the edge len
                #   (b) graft_dp / backbone_dp           — chain-length scaling
                #   (c) 0.5 * lattice_spacing / edge_len — never past half
                #       a lattice cell into neighbouring cells
                edge_grafts = graft_atom_map.get((u, v, key))
                if edge_grafts:
                    lattice_spacing = float(min(self.dims)) if self.dims is not None else edge_len
                    for frac, g_atoms in edge_grafts:
                        k_float = frac * (len(atoms) + 1) - 1
                        k = max(0, min(int(round(k_float)), len(backbone_xyz) - 1))
                        anchor_pos = backbone_xyz[k]
                        local_perp = _local_perp_unit(
                            backbone_xyz, k, unit_vec, rand_vec
                        )
                        graft_dp_eff = max(len(g_atoms), 1)
                        backbone_dp = max(len(atoms), 1)
                        eff_factor = min(
                            ext_factor,
                            graft_dp_eff / backbone_dp,
                            0.5 * lattice_spacing / max(edge_len, 1e-9),
                        )
                        graft_vec = local_perp * (edge_len * eff_factor)
                        for m, g_idx in enumerate(g_atoms):
                            g_frac = (m + 1) / len(g_atoms)
                            graft_coords[g_idx] = tuple(anchor_pos + g_frac * graft_vec)

            write_lammps_displacement_file(
                backbone_coords, sx, sy, sz,
                str(chem_dir / "system_backbone.displace"), "backbone"
            )
            write_lammps_displacement_file(
                graft_coords, sx, sy, sz,
                str(chem_dir / "system_grafts.displace"), "grafts"
            )

            # Pendant heavy + hydrogens via neighbor propagation through mol_h.
            known = {**node_coords, **backbone_coords, **graft_coords}
            side_coords = generate_approximate_side_chain_coords(mol_h, known)
            h_coords = {k: v for k, v in side_coords.items()
                        if mol_h.GetAtomWithIdx(k).GetSymbol() == "H"}
            p_coords = {k: v for k, v in side_coords.items()
                        if mol_h.GetAtomWithIdx(k).GetSymbol() != "H"}

            write_lammps_displacement_file(
                p_coords, sx, sy, sz,
                str(chem_dir / "system_pendant.displace"), "pendant"
            )
            write_lammps_displacement_file(
                h_coords, sx, sy, sz,
                str(chem_dir / "system_hydrogens.displace"), "hydrogens"
            )

        # Groups: nodes (junction atoms) vs beads (everything else).
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

    def _graph_periodicity(self):
        """Per-axis boundaries the topology was built with.

        Recorded by the generators and carried in ``.nodes`` files by a
        ``# PERIODICITY`` header. Returns None when unknown, which every
        consumer reads as fully periodic -- the behaviour before open
        boundaries were supported.
        """
        from topon.topology.loader import graph_periodicity

        return graph_periodicity(self.graph)

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
        # Hand down the same cell stage 4 routed the chains with, so a
        # chain that wraps the boundary lands in a box of the same period,
        # and the boundary conditions so open axes are not wrapped at all.
        periodicity = self._graph_periodicity()
        conformed, roles = cm.apply_displacements(
            "system.data",
            lattice_box=None if self.dims is None else tuple(self.dims),
            periodicity=periodicity,
        )
        noisy = cm.apply_noise(conformed, magnitude=conf_params["noise_magnitude"])
        cm.resolve_overlaps(
            noisy,
            roles,
            cutoff=conf_params["overlap_cutoff"],
            max_iters=conf_params["overlap_max_iters"],
            periodicity=periodicity,
        )
        print()

    # ------------------------------------------------------------------
    # Stage 6: Output
    # ------------------------------------------------------------------

    def _run_output_stage(self) -> None:
        print("--- Stage 6: Output ---")
        from topon.writers import LammpsInputGenerator

        # Optional graph-format exports (GraphML, NPZ).
        if self.config.output.export_graphml:
            from topon.writers.graphml_writer import write_graphml
            graphml_path = self.output_dir / f"{self.config.study.name}.graphml"
            mean_dp = int(self.config.assignment.dp_distribution.default.mean)
            write_graphml(
                self.graph,
                str(graphml_path),
                dp=mean_dp,
                dims=self.dims,
            )
            print(f"  GraphML written to: {graphml_path}")
        if self.config.output.export_npz:
            try:
                from topon.writers.npz_writer import write_npz
            except ImportError:
                print(
                    "  [skip] NPZ export requested but topon.writers.npz_writer "
                    "is not yet implemented (see internal/specs/npz_format.md)."
                )
            else:
                npz_path = self.output_dir / f"{self.config.study.name}.npz"
                write_npz(self.graph, str(npz_path), dims=self.dims)
                print(f"  NPZ written to: {npz_path}")

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
