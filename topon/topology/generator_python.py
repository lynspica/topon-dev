
import random
import time
import math
import networkx as nx
from collections import defaultdict, deque

class PythonTopologyGenerator:
    """
    A Python implementation of the 'Strict Sculpting' algorithm for polymer network generation.
    Designed to exactly match the logic of the C-based generator (generator_serial_debug11.c).
    """

    def __init__(self, config):
        """
        Initialize with a topology configuration.
        Expected config attributes:
        - lattice_source: "SC" (Simple Cubic), "BCC", "FCC" (Currently only SC implemented for benchmark)
        - dimension: tuple (nx, ny, nz)
        - periodicity: bool or tuple
        - degree_distribution: str (e.g., "0:13,1:25,..." or "e:371")
        - functionality: int (max_func)
        """
        self.config = config
        self.dims = getattr(config, 'dimension', getattr(config, 'lattice_size', (6, 6, 6)))
        if isinstance(self.dims, str):
             # Parse "6x6x6" string
             try:
                 parts = self.dims.lower().split('x')
                 self.dims = (int(parts[0]), int(parts[1]), int(parts[2]))
             except:
                 print(f"Warning: Could not parse dimension string '{self.dims}', using default (6,6,6)")
                 self.dims = (6, 6, 6)
        
        # Accept either the schema's `lattice_type` (current) or the legacy
        # `lattice_source` attribute name (older callers / namedtuples).
        self.lattice_type = getattr(
            config, 'lattice_type',
            getattr(config, 'lattice_source', 'SC'),
        )
        self.max_func = getattr(config, 'max_functionality', getattr(config, 'functionality', 4))

        # Mixed-lattice knobs. Only read when lattice_type == "MIX"; the
        # defaults reproduce a plain simple-cubic lattice.
        self.mix_fractions = dict(
            getattr(config, 'mix_fractions', None)
            or {"SC": 1.0, "BCC": 0.0, "FCC": 0.0}
        )
        self.mix_cutoff = float(getattr(config, 'mix_cutoff', 1.0) or 1.0)

        # Per-axis periodic boundaries, matching the C searcher's p_dims.
        # Defaults to fully periodic, which is what every builder did
        # unconditionally before this was read at all.
        self.periodicity = self._parse_periodicity(
            getattr(config, 'periodicity', '111')
        )
        
        # Parse degree distribution string
        self.target_counts = defaultdict(lambda: -2)  # -2 means not specified
        self.target_edge_count = -1
        self._parse_degree_distribution(getattr(config, 'degree_distribution', ""))

    @staticmethod
    def _wrap_hr(coord, hr_dims, periodic):
        """Wrap a high-res neighbour address, or None if it left an open face.

        Mirrors the C searcher's neighbour step: wrap the axis when it is
        periodic, otherwise keep the raw index and let the bounds check
        drop it.
        """
        out = []
        for value, extent, is_periodic in zip(coord, hr_dims, periodic):
            if is_periodic:
                out.append(value % extent)
            elif 0 <= value < extent:
                out.append(value)
            else:
                return None
        return tuple(out)

    @staticmethod
    def _parse_periodicity(value):
        """Normalise a periodicity spec to ``(px, py, pz)`` booleans.

        Accepts the C searcher's ``"111"`` / ``"110"`` digit string, a
        single bool applied to all three axes, or any 3-element iterable
        of truthy values. Anything unrecognised falls back to fully
        periodic with a warning rather than silently producing an open
        lattice, since a surprise free surface changes the physics.

        A non-periodic axis simply omits the wrap-around bonds, so the
        lattice grows a free surface there and the sites on it have
        reduced coordination.
        """
        if value is None:
            return (True, True, True)
        if isinstance(value, bool):
            return (value, value, value)
        if isinstance(value, str):
            digits = [c for c in value.strip() if c in "01"]
            if len(digits) == 3:
                return tuple(c == "1" for c in digits)
            print(f"Warning: could not parse periodicity {value!r}; "
                  f"using fully periodic '111'.")
            return (True, True, True)
        try:
            axes = tuple(bool(v) for v in value)
        except TypeError:
            print(f"Warning: could not parse periodicity {value!r}; "
                  f"using fully periodic '111'.")
            return (True, True, True)
        if len(axes) != 3:
            print(f"Warning: periodicity {value!r} does not have 3 axes; "
                  f"using fully periodic '111'.")
            return (True, True, True)
        return axes

    def _parse_degree_distribution(self, dist_str):
        if not dist_str:
            return
        
        parts = dist_str.split(',')
        for part in parts:
            part = part.strip()
            if part.startswith('e:'):
                self.target_edge_count = int(part.split(':')[1])
            elif ':' in part:
                d_str, n_str = part.split(':')
                d = int(d_str.replace('d', '')) # Handle "d3" or "3"
                self.target_counts[d] = int(n_str)

    def _lattice_label(self):
        """Human-readable "<nx>x<ny>x<nz> <TYPE>" tag for error messages."""
        label = f"{self.dims[0]}x{self.dims[1]}x{self.dims[2]} {self.lattice_type}"
        if self.lattice_type == "MIX":
            frac = ",".join(
                f"{k}:{self.mix_fractions.get(k, 0.0):g}"
                for k in ("SC", "BCC", "FCC")
            )
            label += f" ({frac})"
        return label

    def _validate_targets_reachable(self, base_graph):
        """Fail fast when the requested ``degree_distribution`` can never be met.

        Sculpting only ever REMOVES edges from the freshly-built lattice, so
        that full lattice is a hard ceiling: its edge count bounds any ``e:N``
        target and its maximum node degree bounds any per-degree ``d:N`` target.
        Without this guard an over-target request makes :meth:`generate` grind
        through every trial (each one doomed) before giving up — for a large
        ``trials`` count that looks like an indefinite hang. Bounds are read
        from the actual constructed graph, not a ``3*nx*ny*nz`` formula, so
        periodic-boundary edge collapse on tiny lattices (e.g. a 2x2x2 SC has
        12 edges, not 24) is accounted for.

        Raises
        ------
        ValueError
            If the target edge count exceeds the lattice's edges, or a
            per-degree target asks for more nodes than exist or for a degree
            higher than any node in the base lattice.
        """
        base_edges = base_graph.number_of_edges()
        base_nodes = base_graph.number_of_nodes()
        label = self._lattice_label()

        # --- e:N  (total edge-count target) ---
        if self.target_edge_count != -1 and self.target_edge_count > base_edges:
            raise ValueError(
                f"degree_distribution e:{self.target_edge_count} exceeds the "
                f"{base_edges} edges of a {label} lattice; sculpting only "
                f"removes edges, so this target is unreachable."
            )

        # --- d:N  (per-degree targets) ---
        # target_counts is a defaultdict; iterate a snapshot of only the
        # explicitly-parsed entries. count == 0 means "forbidden" (reachable)
        # and the -2 sentinel means "unspecified"; both are skipped.
        if self.target_counts:
            max_base_degree = max((d for _, d in base_graph.degree()), default=0)
            for degree, count in list(self.target_counts.items()):
                if count <= 0:
                    continue
                if count > base_nodes:
                    raise ValueError(
                        f"degree_distribution {degree}:{count} exceeds the "
                        f"{base_nodes} nodes of a {label} lattice; a lattice "
                        f"cannot hold more nodes of degree {degree} than it has "
                        f"nodes, so this target is unreachable."
                    )
                if degree > max_base_degree:
                    raise ValueError(
                        f"degree_distribution {degree}:{count} requires degree-"
                        f"{degree} nodes, but the maximum degree in a {label} "
                        f"lattice is {max_base_degree}; sculpting only removes "
                        f"edges, so this target is unreachable."
                    )
                # Checked after the lattice bound, which is the more
                # fundamental reason when both apply. Stage 3 prunes every
                # node to max_func and stage 4 refuses to finish while any
                # sits above it, so a target above the ceiling can never be
                # met however rich the lattice was. Without this the run
                # burns through every trial before giving up, and the C
                # searcher used to report success on such a request
                # outright -- see topon/topology/csrc/README.md.
                if degree > self.max_func:
                    raise ValueError(
                        f"degree_distribution {degree}:{count} requires degree-"
                        f"{degree} nodes on a {label} lattice, but "
                        f"max_functionality is {self.max_func}; sculpting "
                        f"enforces that ceiling, so no node can finish with "
                        f"degree {degree}."
                    )

    def generate(self, trials=1, max_saves=1, time_limit=None):
        """
        Run multiple trials to generate a valid network.
        Returns a list of successful graphs (networkx.Graph objects).
        """
        successful_graphs = []

        base_graph = self._create_lattice(self.dims, self.lattice_type)
        # Reject structurally-unreachable targets before churning through
        # trials (sculpting only removes edges — see _validate_targets_reachable).
        self._validate_targets_reachable(base_graph)
        print(f"DEBUG: Entering generate loop (trials={trials})")
        
        start_time = time.time()
        
        for trial in range(trials):
            if time_limit and (time.time() - start_time > time_limit):
                print(f"  [Python] Time limit reached ({time_limit}s). Stopping.")
                break
                
            if trial % 100 == 0:
                 print(f"  [Python] Trial {trial}/{trials}...")
            g = self.run_single_trial(base_graph, trial)
            if g is not None:
                print(f"  [Python] Success on trial {trial}!")
                successful_graphs.append(g)
                if len(successful_graphs) >= max_saves:
                    break
                    
        return successful_graphs

    def _create_lattice(self, dims, lattice_type):
        """
        Creates the initial full lattice.
        Matches C generator's create_sc/bcc/fcc_lattice functions.
        """
        nx_val, ny_val, nz_val = dims

        if lattice_type == "SC":
            return self._create_sc_lattice(nx_val, ny_val, nz_val)
        elif lattice_type == "BCC":
            return self._create_bcc_lattice(nx_val, ny_val, nz_val)
        elif lattice_type == "FCC":
            return self._create_fcc_lattice(nx_val, ny_val, nz_val)
        elif lattice_type == "MIX":
            return self._create_mixed_lattice(nx_val, ny_val, nz_val)
        elif lattice_type in ("Diamond", "DIAMOND"):
            # Delegates to the standalone module rather than inlining the
            # basis here: the Diamond logic deliberately lives in its own
            # file so it can be reviewed in isolation. This is a dispatch
            # bridge so a config can name Diamond on either generator, not
            # a merge of the two.
            from topon.topology.generator_python_diamond import (
                create_diamond_lattice,
            )
            return create_diamond_lattice(
                nx_val, ny_val, nz_val, self.periodicity
            )
        else:
            raise NotImplementedError(
                f"Lattice type {lattice_type} not supported. "
                f"Use SC, BCC, FCC, Diamond, or MIX."
            )

    def _create_sc_lattice(self, nx_val, ny_val, nz_val):
        """Simple Cubic: N nodes, 6 neighbours each when fully periodic.

        Honours ``self.periodicity`` per axis, matching the C searcher's
        ``if (p_dims[0] || x < Nx - 1)`` guard: an open axis simply omits
        the wrap-around bond, leaving a free surface whose sites have
        reduced coordination.
        """
        px, py, pz = self.periodicity
        g = nx.Graph()
        g.graph["box"] = (float(nx_val), float(ny_val), float(nz_val))
        g.graph["periodicity"] = (px, py, pz)

        total_nodes = nx_val * ny_val * nz_val
        for i in range(total_nodes):
            z = i // (nx_val * ny_val)
            rem = i % (nx_val * ny_val)
            y = rem // nx_val
            x = rem % nx_val
            g.add_node(i, pos=(float(x), float(y), float(z)))

        for z in range(nz_val):
            for y in range(ny_val):
                for x in range(nx_val):
                    u = z * (nx_val * ny_val) + y * nx_val + x

                    if px or x < nx_val - 1:
                        v_x = z * (nx_val * ny_val) + y * nx_val + (x + 1) % nx_val
                        if u != v_x and not g.has_edge(u, v_x):
                            g.add_edge(u, v_x)

                    if py or y < ny_val - 1:
                        v_y = z * (nx_val * ny_val) + ((y + 1) % ny_val) * nx_val + x
                        if u != v_y and not g.has_edge(u, v_y):
                            g.add_edge(u, v_y)

                    if pz or z < nz_val - 1:
                        v_z = ((z + 1) % nz_val) * (nx_val * ny_val) + y * nx_val + x
                        if u != v_z and not g.has_edge(u, v_z):
                            g.add_edge(u, v_z)

        return g

    def _create_bcc_lattice(self, nx_val, ny_val, nz_val):
        """Body-Centered Cubic: 2*N nodes, 8 neighbors each (periodic BC).

        Mirrors C create_bcc_lattice:
        - Corner atoms at (i, j, k), high-res coords (2i, 2j, 2k)
        - Body atoms at (i+0.5, j+0.5, k+0.5), high-res coords (2i+1, 2j+1, 2k+1)
        - Neighbors via all 8 (±1, ±1, ±1) offsets in high-res space

        The periodic cell is (nx, ny, nz), not the extent of the site
        coordinates: body-centre sites sit at +0.5, so the coordinates
        only reach nx-0.5 and a max-min+1 estimate would overshoot by
        half a cell. Recording the cell explicitly keeps every downstream
        minimum-image calculation on the right periodic replica.
        """
        px, py, pz = self.periodicity
        g = nx.Graph()
        g.graph["box"] = (float(nx_val), float(ny_val), float(nz_val))
        g.graph["periodicity"] = (px, py, pz)

        hr_nx = 2 * nx_val
        hr_ny = 2 * ny_val
        hr_nz = 2 * nz_val

        # Map from high-res (hx, hy, hz) index -> node id
        coord_map = {}
        node_idx = 0

        for k in range(nz_val):
            for j in range(ny_val):
                for i in range(nx_val):
                    # Corner node
                    cx, cy, cz = 2 * i, 2 * j, 2 * k
                    g.add_node(node_idx, pos=(float(i), float(j), float(k)))
                    coord_map[(cx, cy, cz)] = node_idx
                    node_idx += 1

                    # Body-center node
                    bx, by, bz = 2 * i + 1, 2 * j + 1, 2 * k + 1
                    g.add_node(node_idx, pos=(i + 0.5, j + 0.5, k + 0.5))
                    coord_map[(bx, by, bz)] = node_idx
                    node_idx += 1

        # Connect: each node links to 8 diagonal neighbors in high-res space.
        # An open axis does not wrap, so neighbours off that face simply
        # do not exist and the surface sites lose coordination.
        for (hx, hy, hz), uid in coord_map.items():
            for dx in (-1, 1):
                for dy in (-1, 1):
                    for dz in (-1, 1):
                        nbr = self._wrap_hr(
                            (hx + dx, hy + dy, hz + dz),
                            (hr_nx, hr_ny, hr_nz), (px, py, pz),
                        )
                        if nbr is None:
                            continue
                        vid = coord_map.get(nbr)
                        if vid is not None and uid < vid:
                            g.add_edge(uid, vid)

        return g

    def _create_fcc_lattice(self, nx_val, ny_val, nz_val):
        """Face-Centered Cubic: 4*N nodes, 12 neighbors each (periodic BC).

        Mirrors C create_fcc_lattice:
        - Corner at (2i, 2j, 2k)
        - Face-XY at (2i+1, 2j+1, 2k)
        - Face-XZ at (2i+1, 2j, 2k+1)
        - Face-YZ at (2i, 2j+1, 2k+1)
        - Neighbors via 12 face-diagonal offsets: XY(±1,±1,0), XZ(±1,0,±1), YZ(0,±1,±1)

        As for BCC, the periodic cell is (nx, ny, nz) while the face-site
        coordinates only reach nx-0.5, so the cell is recorded rather than
        inferred from the coordinate extent.
        """
        px, py, pz = self.periodicity
        g = nx.Graph()
        g.graph["box"] = (float(nx_val), float(ny_val), float(nz_val))
        g.graph["periodicity"] = (px, py, pz)

        hr_nx = 2 * nx_val
        hr_ny = 2 * ny_val
        hr_nz = 2 * nz_val

        coord_map = {}
        node_idx = 0

        for k in range(nz_val):
            for j in range(ny_val):
                for i in range(nx_val):
                    # Corner
                    coord_map[(2 * i, 2 * j, 2 * k)] = node_idx
                    g.add_node(node_idx, pos=(float(i), float(j), float(k)))
                    node_idx += 1
                    # Face XY (z shared)
                    coord_map[(2 * i + 1, 2 * j + 1, 2 * k)] = node_idx
                    g.add_node(node_idx, pos=(i + 0.5, j + 0.5, float(k)))
                    node_idx += 1
                    # Face XZ (y shared)
                    coord_map[(2 * i + 1, 2 * j, 2 * k + 1)] = node_idx
                    g.add_node(node_idx, pos=(i + 0.5, float(j), k + 0.5))
                    node_idx += 1
                    # Face YZ (x shared)
                    coord_map[(2 * i, 2 * j + 1, 2 * k + 1)] = node_idx
                    g.add_node(node_idx, pos=(float(i), j + 0.5, k + 0.5))
                    node_idx += 1

        # 12 nearest-neighbor offsets in high-res space (face diagonals)
        fcc_offsets = [
            (1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0),  # XY plane
            (1, 0, 1), (1, 0, -1), (-1, 0, 1), (-1, 0, -1),  # XZ plane
            (0, 1, 1), (0, 1, -1), (0, -1, 1), (0, -1, -1),  # YZ plane
        ]

        for (hx, hy, hz), uid in coord_map.items():
            for dx, dy, dz in fcc_offsets:
                nbr = self._wrap_hr(
                    (hx + dx, hy + dy, hz + dz),
                    (hr_nx, hr_ny, hr_nz), (px, py, pz),
                )
                if nbr is None:
                    continue
                vid = coord_map.get(nbr)
                if vid is not None and uid < vid:
                    g.add_edge(uid, vid)

        return g

    def _create_mixed_lattice(self, nx_val, ny_val, nz_val):
        """Overlay of SC / BCC / FCC basis sites in one cubic cell.

        All three lattices share the cell corner, and each adds its own
        sites on top of it: BCC one body centre, FCC three face centres.
        So the corner is placed in every cell, the body centre with
        probability ``mix_fractions["BCC"]`` and each face centre with
        probability ``mix_fractions["FCC"]``. The ``"SC"`` fraction is the
        remainder, contributing no site of its own, which is what makes
        the three fractions a partition summing to 1.

        Expected site count is ``Nx*Ny*Nz * (1 + f_bcc + 3*f_fcc)``, which
        recovers the exact counts of the pure lattices: ``N`` for SC,
        ``2N`` for BCC, ``4N`` for FCC.

        Edges join every pair within ``mix_cutoff`` under the minimum
        image, rather than the fixed offset patterns the pure builders
        use, because on a mixed point set there is no single neighbour
        shell. The 1.0 default is the simple-cubic nearest-neighbour
        distance, which keeps the always-present corner sublattice
        connected however few body and face sites are drawn.

        Two consequences worth knowing:

        * ``MIX`` at fractions ``(1, 0, 0)`` reproduces ``SC`` exactly,
          same node ids, positions and edges. It is **not** true of the
          other two corners: at ``(0, 1, 0)`` the cutoff also admits the
          corner-corner shell at 1.0, so nodes carry 14 neighbours rather
          than BCC's 8, and at ``(0, 0, 1)`` 18 rather than FCC's 12. Use
          ``lattice_type`` SC / BCC / FCC when the canonical coordination
          is what you want; ``MIX`` is for genuine mixtures.
        * Body and face sites can land 0.5 cells apart, closer than any
          pure lattice's nearest-neighbour distance (SC 1.0, BCC 0.866,
          FCC 0.707). Since DP is assigned independently of edge length,
          that widens the spread of bond lengths a strand of given DP is
          built at.
        """
        import numpy as np

        f_bcc = float(self.mix_fractions.get("BCC", 0.0))
        f_fcc = float(self.mix_fractions.get("FCC", 0.0))

        px, py, pz = self.periodicity
        g = nx.Graph()
        g.graph["box"] = (float(nx_val), float(ny_val), float(nz_val))
        g.graph["periodicity"] = (px, py, pz)
        g.graph["mix_fractions"] = dict(self.mix_fractions)

        # Face-centre offsets, in the same XY / XZ / YZ order the pure FCC
        # builder uses so a fraction of 1.0 gives the same site set.
        face_offsets = ((0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))

        positions = []
        # Cell order matches the SC builder (z outer, then y, then x) so
        # that fractions (1, 0, 0) yields identical node ids.
        for k in range(nz_val):
            for j in range(ny_val):
                for i in range(nx_val):
                    positions.append((float(i), float(j), float(k)))
                    if f_bcc > 0.0 and random.random() < f_bcc:
                        positions.append((i + 0.5, j + 0.5, k + 0.5))
                    if f_fcc > 0.0:
                        for fx, fy, fz in face_offsets:
                            if random.random() < f_fcc:
                                positions.append((i + fx, j + fy, k + fz))

        for idx, pos in enumerate(positions):
            g.add_node(idx, pos=pos)

        # Neighbour search under the minimum image. Blocked rather than
        # one big (N, N, 3) array so a large lattice does not blow up
        # memory; N stays in the low thousands for realistic cell counts.
        pts = np.asarray(positions, dtype=float)
        box = np.array([nx_val, ny_val, nz_val], dtype=float)
        # Only periodic axes take the minimum image; an open axis keeps
        # the raw separation, so nothing bonds across that face.
        wrap = np.array([px, py, pz], dtype=bool)
        cutoff_sq = self.mix_cutoff ** 2
        n = len(pts)
        block = 512
        for start in range(0, n, block):
            chunk = pts[start:start + block]
            delta = chunk[:, None, :] - pts[None, :, :]
            delta -= np.where(wrap, box * np.round(delta / box), 0.0)
            dist_sq = (delta * delta).sum(axis=-1)
            rows, cols = np.nonzero(
                (dist_sq <= cutoff_sq + 1e-12) & (dist_sq > 1e-12)
            )
            for r, c in zip(rows + start, cols):
                if r < c:
                    g.add_edge(int(r), int(c))

        return g

    def run_single_trial(self, base_graph, trial_num):
        """
        Runs the Strict Sculpting algorithm stages on a copy of the base graph.
        """
        g = base_graph.copy()
        total_nodes = g.number_of_nodes()
        
        # Track edge removal history for visualization
        move_history = []
        
        # Node Status: 0=ACTIVE, 1=IS_DEGREE_0, 2=IS_DEGREE_1
        # In Python we can use a dict or node attribute
        # Default active
        node_status = {n: "ACTIVE" for n in g.nodes()}
        
        # Shuffle node indices
        node_indices = list(g.nodes())
        random.shuffle(node_indices)
        
        current_node_offset = 0
        
        # Targets
        n0_target = max(0, self.target_counts[0]) if self.target_counts[0] != -2 else 0
        n1_target = max(0, self.target_counts[1]) if self.target_counts[1] != -2 else 0
        
        target_degree_sum = self.target_edge_count * 2
        
        # --- Stage 1: Set d0 (Strict) ---
        for _ in range(n0_target):
            if current_node_offset >= total_nodes: break
            node_idx = node_indices[current_node_offset]
            current_node_offset += 1
            
            while g.degree[node_idx] > 0:
                neighbors = list(g.neighbors(node_idx))
                removed = False
                for neighbor in neighbors:
                    # SC Optimization in C: if neighbor degree <= 2, skip to avoid breaking chains too much
                    if g.degree[neighbor] <= 2:
                        continue
                        
                    if not self._is_move_safe(g, node_idx, neighbor, stage=1, 
                                              target_degree_sum=target_degree_sum, 
                                              current_total_degree_sum=-1): # sum not needed for stg 1
                        continue
                        
                    g.remove_edge(node_idx, neighbor)
                    move_history.append({'stage': 1, 'edge': (node_idx, neighbor), 'reason': 'd0'})
                    removed = True
                    break
                
                if not removed:
                    return None # Failed to isolate node
            
            node_status[node_idx] = "IS_DEGREE_0"
            
        # --- Stage 2: Set d1 (Strict) ---
        for _ in range(n1_target):
            if current_node_offset >= total_nodes: break
            node_idx = node_indices[current_node_offset]
            if node_status[node_idx] != "ACTIVE": 
                # Should not happen as we iterate sequential offset, but good check
                pass
            current_node_offset += 1
            
            while g.degree[node_idx] > 1:
                neighbors = list(g.neighbors(node_idx))
                random.shuffle(neighbors)
                removed = False
                for neighbor in neighbors:
                    if g.degree[neighbor] <= 2:
                        continue

                    if not self._is_move_safe(g, node_idx, neighbor, stage=2, 
                                              target_degree_sum=target_degree_sum, 
                                              current_total_degree_sum=-1):
                        continue
                        
                    g.remove_edge(node_idx, neighbor)
                    
                    if self._is_subgraph_connected(g, node_status):
                        move_history.append({'stage': 2, 'edge': (node_idx, neighbor), 'reason': 'd1'})
                        removed = True
                        break
                    else:
                        g.add_edge(node_idx, neighbor) # Backtrack
                        
                if not removed:
                    return None # Failed to reduce to d1
            
            node_status[node_idx] = "IS_DEGREE_1"

        # --- Stage 3: Enforce Max Functionality (Strict) ---
        for i in range(total_nodes):
            node_idx = node_indices[i]
            if node_status[node_idx] != "ACTIVE": continue
            
            while g.degree[node_idx] > self.max_func:
                neighbors = list(g.neighbors(node_idx))
                random.shuffle(neighbors)
                removed = False
                for neighbor in neighbors:
                    if g.degree[neighbor] <= 2:
                        continue
                        
                    if not self._is_move_safe(g, node_idx, neighbor, stage=3, 
                                              target_degree_sum=target_degree_sum, 
                                              current_total_degree_sum=-1):
                        continue
                        
                    g.remove_edge(node_idx, neighbor)
                    if self._is_subgraph_connected(g, node_status):
                        move_history.append({'stage': 3, 'edge': (node_idx, neighbor), 'reason': 'max_func'})
                        removed = True
                        break
                    else:
                        g.add_edge(node_idx, neighbor) # Backtrack
                
                if not removed:
                    return None # Failed to enforce max func

        # --- Stage 4: Systematic Search Loop ---
        while True:
            # Check current distribution
            current_degree_sum = sum(d for n, d in g.degree())
            current_counts = defaultdict(int)
            has_high_degree = False
            for n, d in g.degree():
                current_counts[d] += 1
                if node_status[n] == "ACTIVE" and d > self.max_func:
                    has_high_degree = True
            
            is_done = True
            
            if has_high_degree:
                is_done = False
            else:
                # 1. Check explicit targets
                for d, count in self.target_counts.items():
                    if count >= 0 and current_counts[d] != count:
                        is_done = False
                        break
                
                # 2. Check total edge count / connectivity depending on mode
                if is_done:
                    if self.target_edge_count != -1:
                         # e:N mode
                        if current_degree_sum != target_degree_sum:
                            is_done = False
                        if not self._is_subgraph_connected(g, node_status):
                            is_done = False
                    else:
                        # Legacy mode (d0 already met, just need connectivity)
                        if self._is_subgraph_connected(g, node_status):
                             is_done = True
                        else:
                             return None # Failed connectivity check at end
            
            if is_done:
                # Attach move history to graph
                g.graph['move_history'] = move_history
                return g
            
            # Not done, perform systematic edge removal
            edges = list(g.edges())
            random.shuffle(edges)
            
            move_made = False
            
            for u, v in edges:
                u_deg = g.degree[u]
                v_deg = g.degree[v]
                
                if u_deg <= 1 or v_deg <= 1: continue
                
                # Legacy d2 check
                if u_deg == 2 or v_deg == 2:
                    if self.target_counts[1] != -1: # if d1 count is tracked
                        d1_count = sum(1 for n in g.nodes() if node_status[n] == "ACTIVE" and g.degree[n] == 1)
                        if d1_count >= self.target_counts[1]:
                            continue
                            
                if not self._is_move_safe(g, u, v, stage=4, 
                                          target_degree_sum=target_degree_sum, 
                                          current_total_degree_sum=current_degree_sum):
                    continue
                
                g.remove_edge(u, v)
                
                if self._is_subgraph_connected(g, node_status):
                    move_history.append({'stage': 4, 'edge': (u, v), 'reason': 'systematic'})
                    move_made = True
                    break # Restart loop
                else:
                    g.add_edge(u, v) # Backtrack
            
            if not move_made:
                return None # Stuck

    def _is_subgraph_connected(self, g, node_status):
        """
        Checks if the subgraph of ACTIVE nodes is connected.
        Ignores IS_DEGREE_0 nodes (they are isolated by definition).
        IS_DEGREE_1 nodes are part of the active graph usually?
        Wait, C code: `if (node_status[i] == ACTIVE) ...`
        Wait, in C `IS_DEGREE_1` nodes are *excluded* from the connectivity check loop?
        
        Let's look at C code `is_subgraph_connected`:
        `if (node_status[i] == ACTIVE) ...`
        Yes! C code ONLY checks connectivity among "ACTIVE" nodes.
        Nodes marked IS_DEGREE_0 or IS_DEGREE_1 are NOT part of the connectivity check.
        They are considered "done" and "removed" from the main component logic?
        
        Wait, `IS_DEGREE_1` nodes (dangling ends) *should* be connected to the main component.
        If we exclude them from the check, we only ensure the core is connected.
        Let's double check C code line 304: `if (node_status[node_status[i] == ACTIVE])`.
        
        In Stage 2 (Set d1), we mark nodes as `IS_DEGREE_1`.
        If they are excluded from connectivity check, that means we only care if the *remaining* network is connected.
        Dangling ends are by definition connected to *something* (degree 1), so as long as that something is in the main component, they are fine.

        Implementation note: this walks the adjacency mapping directly
        rather than building ``g.subgraph(active)`` and calling
        ``nx.is_connected``. The two give identical answers, but a
        NetworkX subgraph is a *view* that re-evaluates its node filter on
        every neighbour access, which costs about six million predicate
        calls per check on a 1000-node lattice. Since this function is
        roughly 99% of the generator's runtime, that made the whole
        generator about eight times slower than it needed to be.
        """
        # Edges count only when BOTH ends are ACTIVE, matching the C
        # searcher: `if (node_status[pCrawl->dest] == ACTIVE) unite_sets(...)`.
        # Raw {node: {neighbour: attrs}}. `_adj` is NetworkX-internal but
        # stable and ~1.6x faster than the public `adj` view, so take it
        # when present and fall back if a future release renames it.
        adj = getattr(g, "_adj", None)
        if adj is None:
            adj = g.adj
        start = None
        n_active = 0
        for n in adj:
            if node_status[n] == "ACTIVE":
                n_active += 1
                if start is None:
                    start = n
        if start is None:
            return True

        seen = {start}
        stack = [start]
        while stack:
            for nbr in adj[stack.pop()]:
                if nbr not in seen and node_status[nbr] == "ACTIVE":
                    seen.add(nbr)
                    stack.append(nbr)
        return len(seen) == n_active


    def _is_move_safe(self, g, u, v, stage, target_degree_sum, current_total_degree_sum):
        """
        Equivalent to C `is_move_safe`.
        """
        
        # --- Target Edge Count Check (Stage 4 only) ---
        if stage == 4 and target_degree_sum != -2: # -2 is check for "not set"
             # If removing this edge (degree sum - 2) drops us below target
             if current_total_degree_sum <= target_degree_sum:
                 return False

        u_new_degree = g.degree[u] - 1
        v_new_degree = g.degree[v] - 1
        
        # --- Check 1: Victim 'v' ---
        
        # 1a. Forbidden Degree (target=0)
        # In Python target_counts[d] returns -2 if not set.
        # If set to 0, it means forbidden.
        if v_new_degree >= 0 and self.target_counts[v_new_degree] == 0:
            return False
            
        # 1b. Overshooting
        if v_new_degree >= 0 and self.target_counts[v_new_degree] > 0:
            # Case A: d0/d1 (Sacred)
            if v_new_degree <= 1:
                # Count current
                current_count = sum(1 for n in g.nodes() if g.degree[n] == v_new_degree)
                if current_count >= self.target_counts[v_new_degree]:
                    return False
            # Case B: d2+ (Only Stage 4)
            elif stage == 4:
                current_count = sum(1 for n in g.nodes() if g.degree[n] == v_new_degree)
                if current_count >= self.target_counts[v_new_degree]:
                    return False
                    
        # --- Check 2: Actor 'u' (Only Stage 4) ---
        if stage == 4:
            # 2a. Forbidden
            if u_new_degree >= 0 and self.target_counts[u_new_degree] == 0:
                return False
                
            # 2b. Overshooting
            if u_new_degree >= 0 and self.target_counts[u_new_degree] > 0:
                 current_count = sum(1 for n in g.nodes() if g.degree[n] == u_new_degree)
                 if current_count >= self.target_counts[u_new_degree]:
                     return False
                     
        return True
