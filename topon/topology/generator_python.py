
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
        
        self.lattice_type = getattr(config, 'lattice_source', 'SC')
        self.max_func = getattr(config, 'max_functionality', getattr(config, 'functionality', 4))
        
        # Parse degree distribution string
        self.target_counts = defaultdict(lambda: -2)  # -2 means not specified
        self.target_edge_count = -1
        self._parse_degree_distribution(getattr(config, 'degree_distribution', ""))

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
        
    def generate(self, trials=1, max_saves=1, time_limit=None):
        """
        Run multiple trials to generate a valid network.
        Returns a list of successful graphs (networkx.Graph objects).
        """
        successful_graphs = []
        
        base_graph = self._create_lattice(self.dims, self.lattice_type)
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
        else:
            raise NotImplementedError(f"Lattice type {lattice_type} not supported. Use SC, BCC, or FCC.")

    def _create_sc_lattice(self, nx_val, ny_val, nz_val):
        """Simple Cubic: N nodes, 6 neighbors each (periodic BC)."""
        g = nx.Graph()

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

                    v_x = z * (nx_val * ny_val) + y * nx_val + (x + 1) % nx_val
                    if not g.has_edge(u, v_x):
                        g.add_edge(u, v_x)

                    v_y = z * (nx_val * ny_val) + ((y + 1) % ny_val) * nx_val + x
                    if not g.has_edge(u, v_y):
                        g.add_edge(u, v_y)

                    v_z = ((z + 1) % nz_val) * (nx_val * ny_val) + y * nx_val + x
                    if not g.has_edge(u, v_z):
                        g.add_edge(u, v_z)

        return g

    def _create_bcc_lattice(self, nx_val, ny_val, nz_val):
        """Body-Centered Cubic: 2*N nodes, 8 neighbors each (periodic BC).

        Mirrors C create_bcc_lattice:
        - Corner atoms at (i, j, k), high-res coords (2i, 2j, 2k)
        - Body atoms at (i+0.5, j+0.5, k+0.5), high-res coords (2i+1, 2j+1, 2k+1)
        - Neighbors via all 8 (±1, ±1, ±1) offsets in high-res space
        """
        g = nx.Graph()

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

        # Connect: each node links to 8 diagonal neighbors in high-res space
        for (hx, hy, hz), uid in coord_map.items():
            for dx in (-1, 1):
                for dy in (-1, 1):
                    for dz in (-1, 1):
                        nx_ = (hx + dx) % hr_nx
                        ny_ = (hy + dy) % hr_ny
                        nz_ = (hz + dz) % hr_nz
                        vid = coord_map.get((nx_, ny_, nz_))
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
        """
        g = nx.Graph()

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
                nx_ = (hx + dx) % hr_nx
                ny_ = (hy + dy) % hr_ny
                nz_ = (hz + dz) % hr_nz
                vid = coord_map.get((nx_, ny_, nz_))
                if vid is not None and uid < vid:
                    g.add_edge(uid, vid)

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
        """
        active_nodes = [n for n in g.nodes() if node_status[n] == "ACTIVE"]
        if not active_nodes: return True
        
        # Determine subgraph of active nodes
        # Note: This means we only traverse edges where BOTH ends are ACTIVE?
        # C code:
        # `if (node_status[pCrawl->dest] == ACTIVE) unite_sets(...)`
        # Yes, edges are only considered if both nodes are ACTIVE.
        
        subg = g.subgraph(active_nodes)
        return nx.is_connected(subg)


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
