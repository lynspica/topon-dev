"""Diamond-lattice topology generator (Python).

A self-contained, standalone module that builds the **Diamond** lattice
and runs the same strict-sculpting algorithm as
``generator_python.py``. Lives in its own file (NOT a modification of
the existing generator) so the user can review Diamond logic in
isolation before deciding whether to merge it into the main generator.

The Diamond lattice is **the only common cubic lattice where every
node is exactly 4-coordinated by construction**. Two interpenetrating
FCC sublattices, offset by (¼, ¼, ¼) along the body diagonal:

   sublattice A:  (0,0,0), (½,½,0), (½,0,½), (0,½,½)
   sublattice B:  (¼,¼,¼), (¾,¾,¼), (¾,¼,¾), (¼,¾,¾)

   → 8 atoms per conventional cubic unit cell
   → every atom has exactly 4 nearest neighbors (tetrahedral coordination)
   → unconstrained network is exactly 4-regular: avg_degree = 4.00

This means for the polymer-network use case, Diamond is the cleanest
backbone for a max_func=4 network: there is no need for the strict-
sculpting algorithm to throw away high-degree neighbors, so the
"native" graph already satisfies max_func=4. With ``degree_distribution
= ""`` the algorithm short-circuits and the raw lattice is returned.

Usage (pure-Python; mirrors the SC/BCC/FCC API)::

    from topon.topology.generator_python_diamond import (
        DiamondTopologyGenerator,
    )
    from types import SimpleNamespace

    cfg = SimpleNamespace(
        lattice_size="6x6x6",         # 8*216 = 1728 atoms
        lattice_type="Diamond",
        max_functionality=4,
        degree_distribution="",        # no constraints -> native 4-regular
    )
    gen = DiamondTopologyGenerator(cfg)
    graphs = gen.generate(trials=1, max_saves=1)
    G = graphs[0]
    assert all(d == 4 for _, d in G.degree())

A one-liner sanity check is provided at the bottom of this file
(``verify_diamond_lattice``) — call it from a Python REPL to confirm
edge count, degree distribution, and connectivity on a fresh
``create_diamond_lattice(N, N, N)``.
"""
from __future__ import annotations

import random
import time
from collections import defaultdict

import networkx as nx


# Tetrahedral neighbor offsets in the high-resolution (×4) integer grid.
# Each diamond atom has 4 neighbors at offsets of the form (±1, ±1, ±1)
# where the *product of signs* selects which sublattice the neighbor
# belongs to. A-type atoms (hr-coord sum ≡ 0 mod 4) reach B-type
# neighbors with sign-product = +1; B-type atoms reach A-type with
# sign-product = -1.
_DIAMOND_OFFSETS_A_TO_B = (   # sign product = +1, jumps A -> B
    (+1, +1, +1),
    (+1, -1, -1),
    (-1, +1, -1),
    (-1, -1, +1),
)
_DIAMOND_OFFSETS_B_TO_A = (   # sign product = -1, jumps B -> A
    (-1, -1, -1),
    (-1, +1, +1),
    (+1, -1, +1),
    (+1, +1, -1),
)


def _diamond_sublattice_offsets(hx: int, hy: int, hz: int):
    """Return the 4 high-res neighbor offsets for the atom at (hx,hy,hz).

    Sublattice is determined by (hx+hy+hz) mod 4:
      0 -> A (corner / face-center positions)        -> A_TO_B offsets
      3 -> B (¼-shifted positions)                    -> B_TO_A offsets

    Diamond places atoms ONLY at sites where the coord-sum is 0 or 3
    mod 4; sites with sum 1 or 2 are empty.
    """
    s = (hx + hy + hz) % 4
    if s == 0:
        return _DIAMOND_OFFSETS_A_TO_B
    if s == 3:
        return _DIAMOND_OFFSETS_B_TO_A
    raise ValueError(
        f"({hx},{hy},{hz}) is not a Diamond lattice site "
        f"(sum mod 4 = {s}, expected 0 or 3)"
    )


def create_diamond_lattice(Nx: int, Ny: int, Nz: int) -> nx.Graph:
    """Build the full Diamond lattice as a NetworkX Graph (periodic BC).

    Returns a graph with ``8*Nx*Ny*Nz`` nodes, each with attribute
    ``pos = (x, y, z)`` in units of the conventional cubic unit cell
    (so the box spans ``[0, Nx] x [0, Ny] x [0, Nz]``).

    Every node has exactly 4 neighbors → ``avg_degree = 4.0``.

    The periodic cell (Nx, Ny, Nz) is recorded on the graph. Diamond
    basis sites sit at quarter-cell offsets, so the coordinates only
    reach Nx-0.25 and estimating the box from their extent would
    overshoot it.
    """
    g = nx.Graph()
    g.graph["box"] = (float(Nx), float(Ny), float(Nz))

    # 8 atoms per conventional unit cell, in (low-res, low-res, low-res)
    # × 4 -> integer high-res coords. The 8 positions are the 4 FCC sites
    # (×4 = (0,0,0), (2,2,0), (2,0,2), (0,2,2)) plus the 4 ¼-offset sites
    # (×4 = (1,1,1), (3,3,1), (3,1,3), (1,3,3)).
    basis_hr = (
        (0, 0, 0), (2, 2, 0), (2, 0, 2), (0, 2, 2),     # A sublattice
        (1, 1, 1), (3, 3, 1), (3, 1, 3), (1, 3, 3),     # B sublattice
    )
    basis_pos = tuple(
        (b[0] / 4.0, b[1] / 4.0, b[2] / 4.0) for b in basis_hr
    )

    hr_Nx, hr_Ny, hr_Nz = 4 * Nx, 4 * Ny, 4 * Nz

    # Step 1: place atoms.
    # We use a hash from (hx, hy, hz) -> node id so the neighbor lookup
    # below can resolve periodic-image addresses by modding.
    coord_to_id: dict[tuple[int, int, int], int] = {}
    node_idx = 0
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                for b_idx, (bx, by, bz) in enumerate(basis_hr):
                    hx = 4 * i + bx
                    hy = 4 * j + by
                    hz = 4 * k + bz
                    pos = (
                        float(i) + basis_pos[b_idx][0],
                        float(j) + basis_pos[b_idx][1],
                        float(k) + basis_pos[b_idx][2],
                    )
                    g.add_node(node_idx, pos=pos)
                    coord_to_id[(hx, hy, hz)] = node_idx
                    node_idx += 1

    # Step 2: connect each atom to its 4 nearest neighbors (PBC via
    # modulo on the high-res grid). The (uid < vid) guard prevents
    # double-adding undirected edges.
    for (hx, hy, hz), uid in coord_to_id.items():
        for dx, dy, dz in _diamond_sublattice_offsets(hx, hy, hz):
            nhx = (hx + dx) % hr_Nx
            nhy = (hy + dy) % hr_Ny
            nhz = (hz + dz) % hr_Nz
            vid = coord_to_id.get((nhx, nhy, nhz))
            if vid is None:
                # Should not happen on a fully populated, periodic
                # Diamond lattice -- but guard anyway.
                continue
            if uid < vid:
                g.add_edge(uid, vid)

    return g


def verify_diamond_lattice(Nx: int = 4, Ny: int = 4, Nz: int = 4) -> dict:
    """Build a fresh Diamond lattice and assert its invariants.

    Returns a small summary dict (for printing in tests / REPL); raises
    AssertionError if any invariant is broken.
    """
    G = create_diamond_lattice(Nx, Ny, Nz)
    n = G.number_of_nodes()
    e = G.number_of_edges()
    expected_n = 8 * Nx * Ny * Nz
    # Each atom has 4 neighbors → total degree sum = 4n → edges = 2n
    expected_e = 2 * n
    degrees = [d for _, d in G.degree()]
    assert n == expected_n, f"node count {n} != expected {expected_n}"
    assert e == expected_e, f"edge count {e} != expected {expected_e}"
    assert all(d == 4 for d in degrees), (
        f"degree distribution not uniformly 4: "
        f"{dict((d, degrees.count(d)) for d in set(degrees))}"
    )
    assert nx.is_connected(G), "Diamond lattice is not connected!"
    return {
        "Nx": Nx, "Ny": Ny, "Nz": Nz,
        "nodes": n, "edges": e,
        "avg_degree": 2 * e / n,
        "max_degree": max(degrees),
        "min_degree": min(degrees),
        "connected": True,
    }


# ----------------------------------------------------------------------
# Strict-sculpting algorithm (cloned from generator_python.py, lightly
# adapted -- no SC/BCC/FCC code paths, just Diamond). This is here so a
# user can pass *constraints* (e.g. some d2 / d3 defects) on top of the
# 4-regular base. With an empty ``degree_distribution``, the raw lattice
# is returned untouched.
# ----------------------------------------------------------------------


class DiamondTopologyGenerator:
    """Diamond-lattice analogue of ``PythonTopologyGenerator``.

    Expected ``config`` attributes (any duck-typed object with these):
      - ``lattice_size``: str "NxNxN" or tuple (Nx, Ny, Nz)
      - ``max_functionality``: int (defaults to 4; cannot exceed 4 on
        Diamond since no node can have more than 4 neighbors)
      - ``degree_distribution``: str like "0:13,2:5" or "" for none
    """

    def __init__(self, config):
        self.config = config

        dims_raw = getattr(
            config, "lattice_size",
            getattr(config, "dimension", (4, 4, 4)),
        )
        if isinstance(dims_raw, str):
            parts = dims_raw.lower().split("x")
            self.dims = (int(parts[0]), int(parts[1]), int(parts[2]))
        else:
            self.dims = tuple(dims_raw)

        # Diamond's lattice_type field is fixed; we still accept it
        # being passed for API symmetry but reject anything else.
        lt = getattr(
            config, "lattice_type",
            getattr(config, "lattice_source", "Diamond"),
        )
        if lt != "Diamond":
            raise ValueError(
                f"DiamondTopologyGenerator only handles Diamond, got {lt!r}"
            )
        self.lattice_type = "Diamond"

        self.max_func = getattr(
            config, "max_functionality",
            getattr(config, "functionality", 4),
        )
        if self.max_func > 4:
            # Allowed but useless: Diamond can never exceed 4-coordination.
            print(
                f"  [Diamond] WARN: max_functionality={self.max_func} > 4; "
                f"clamping effective limit to 4 (lattice is 4-regular)."
            )
            self.max_func = 4

        self.target_counts = defaultdict(lambda: -2)   # -2 = unset
        self.target_edge_count = -1
        self._parse_degree_distribution(
            getattr(config, "degree_distribution", "") or ""
        )

    def _parse_degree_distribution(self, dist_str: str) -> None:
        for part in (p.strip() for p in dist_str.split(",") if p.strip()):
            if part.startswith("e:"):
                self.target_edge_count = int(part.split(":", 1)[1])
            elif ":" in part:
                d_str, n_str = part.split(":", 1)
                d = int(d_str.replace("d", ""))
                self.target_counts[d] = int(n_str)

    def _lattice_label(self) -> str:
        """Human-readable "<nx>x<ny>x<nz> Diamond" tag for error messages."""
        return f"{self.dims[0]}x{self.dims[1]}x{self.dims[2]} {self.lattice_type}"

    def _validate_targets_reachable(self, base_graph: nx.Graph) -> None:
        """Fail fast when the requested ``degree_distribution`` can't be met.

        Sculpting only ever REMOVES edges, so the freshly-built lattice is a
        hard ceiling: its edge count bounds any ``e:N`` target and its maximum
        node degree bounds any per-degree ``d:N`` target. Without this guard an
        over-target request grinds through every trial (each one doomed) before
        giving up, which looks like a hang. Bounds are read from the actual
        constructed graph, not a formula, so periodic-boundary edge collapse on
        tiny lattices is handled. Mirrors
        ``PythonTopologyGenerator._validate_targets_reachable``.
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

        # --- d:N  (per-degree targets; 0 = forbidden, -2 = unspecified) ---
        if self.target_counts:
            max_base_degree = max((d for _, d in base_graph.degree()), default=0)
            for degree, count in list(self.target_counts.items()):
                if count <= 0:
                    continue
                if count > base_nodes:
                    raise ValueError(
                        f"degree_distribution {degree}:{count} exceeds the "
                        f"{base_nodes} nodes of a {label} lattice; a lattice "
                        f"cannot hold more nodes of degree {degree} than it "
                        f"has nodes, so this target is unreachable."
                    )
                if degree > max_base_degree:
                    raise ValueError(
                        f"degree_distribution {degree}:{count} requires degree-"
                        f"{degree} nodes, but the maximum degree in a {label} "
                        f"lattice is {max_base_degree}; sculpting only removes "
                        f"edges, so this target is unreachable."
                    )

    def generate(self, trials: int = 1, max_saves: int = 1,
                 time_limit: float | None = None) -> list[nx.Graph]:
        """Build the base Diamond lattice and run ``run_single_trial``.

        Short-circuits and returns the raw lattice when no constraints
        were specified -- the common case for "give me a clean
        4-regular network."
        """
        base = create_diamond_lattice(*self.dims)
        # Reject structurally-unreachable targets before churning through
        # trials (sculpting only removes edges).
        self._validate_targets_reachable(base)

        # Short-circuit: empty constraint set + no edge target.
        if (
            not any(v >= 0 for v in self.target_counts.values())
            and self.target_edge_count == -1
        ):
            return [base.copy()]

        out: list[nx.Graph] = []
        t0 = time.time()
        for trial in range(trials):
            if time_limit and time.time() - t0 > time_limit:
                print(f"  [Diamond] Time limit ({time_limit}s) reached.")
                break
            g = self._run_single_trial(base, trial)
            if g is not None:
                out.append(g)
                if len(out) >= max_saves:
                    break
        return out

    # -- Strict-sculpting trial (compact clone of generator_python.py) --

    def _run_single_trial(self, base_graph: nx.Graph,
                          trial_num: int) -> nx.Graph | None:
        g = base_graph.copy()
        total_nodes = g.number_of_nodes()
        node_status = {n: "ACTIVE" for n in g.nodes()}
        node_indices = list(g.nodes())
        random.shuffle(node_indices)

        n0_target = max(0, self.target_counts[0]) if self.target_counts[0] != -2 else 0
        n1_target = max(0, self.target_counts[1]) if self.target_counts[1] != -2 else 0
        target_degree_sum = self.target_edge_count * 2 if self.target_edge_count > 0 else -1

        cursor = 0
        # Stage 1: enforce d0 (vacancies).
        for _ in range(n0_target):
            if cursor >= total_nodes:
                break
            u = node_indices[cursor]
            cursor += 1
            while g.degree[u] > 0:
                done = False
                for v in list(g.neighbors(u)):
                    if g.degree[v] <= 2:
                        continue
                    if not self._is_move_safe(g, u, v, 1, target_degree_sum, -1):
                        continue
                    g.remove_edge(u, v)
                    done = True
                    break
                if not done:
                    return None
            node_status[u] = "IS_DEGREE_0"

        # Stage 2: enforce d1 (dangling).
        for _ in range(n1_target):
            if cursor >= total_nodes:
                break
            u = node_indices[cursor]
            cursor += 1
            while g.degree[u] > 1:
                neighbors = list(g.neighbors(u))
                random.shuffle(neighbors)
                done = False
                for v in neighbors:
                    if g.degree[v] <= 2:
                        continue
                    if not self._is_move_safe(g, u, v, 2, target_degree_sum, -1):
                        continue
                    g.remove_edge(u, v)
                    if self._is_subgraph_connected(g, node_status):
                        done = True
                        break
                    g.add_edge(u, v)
                if not done:
                    return None
            node_status[u] = "IS_DEGREE_1"

        # Stage 3: enforce max functionality.
        # Diamond starts 4-regular, so this is a no-op when max_func >= 4.
        for u in node_indices:
            if node_status[u] != "ACTIVE":
                continue
            while g.degree[u] > self.max_func:
                neighbors = list(g.neighbors(u))
                random.shuffle(neighbors)
                done = False
                for v in neighbors:
                    if g.degree[v] <= 2:
                        continue
                    if not self._is_move_safe(g, u, v, 3, target_degree_sum, -1):
                        continue
                    g.remove_edge(u, v)
                    if self._is_subgraph_connected(g, node_status):
                        done = True
                        break
                    g.add_edge(u, v)
                if not done:
                    return None

        # Stage 4: systematic search for explicit d2+/edge-count targets.
        while True:
            current_sum = sum(d for _, d in g.degree())
            counts = defaultdict(int)
            for _, d in g.degree():
                counts[d] += 1

            is_done = True
            for d, cnt in self.target_counts.items():
                if cnt >= 0 and counts[d] != cnt:
                    is_done = False
                    break
            if is_done and target_degree_sum != -1:
                if current_sum != target_degree_sum:
                    is_done = False
                if not self._is_subgraph_connected(g, node_status):
                    is_done = False
            if is_done and not self._is_subgraph_connected(g, node_status):
                return None

            if is_done:
                return g

            edges = list(g.edges())
            random.shuffle(edges)
            made_move = False
            for u, v in edges:
                if g.degree[u] <= 1 or g.degree[v] <= 1:
                    continue
                if not self._is_move_safe(g, u, v, 4, target_degree_sum, current_sum):
                    continue
                g.remove_edge(u, v)
                if self._is_subgraph_connected(g, node_status):
                    made_move = True
                    break
                g.add_edge(u, v)
            if not made_move:
                return None

    def _is_subgraph_connected(self, g: nx.Graph, node_status: dict) -> bool:
        active = [n for n in g.nodes() if node_status[n] == "ACTIVE"]
        if not active:
            return True
        return nx.is_connected(g.subgraph(active))

    def _is_move_safe(self, g, u, v, stage, target_degree_sum,
                      current_total_degree_sum) -> bool:
        if stage == 4 and target_degree_sum != -1:
            if current_total_degree_sum <= target_degree_sum:
                return False
        u_new = g.degree[u] - 1
        v_new = g.degree[v] - 1
        # Victim 'v'.
        if v_new >= 0 and self.target_counts[v_new] == 0:
            return False
        if v_new >= 0 and self.target_counts[v_new] > 0:
            if v_new <= 1 or stage == 4:
                cur = sum(1 for n in g.nodes() if g.degree[n] == v_new)
                if cur >= self.target_counts[v_new]:
                    return False
        # Actor 'u' (stage 4 only).
        if stage == 4:
            if u_new >= 0 and self.target_counts[u_new] == 0:
                return False
            if u_new >= 0 and self.target_counts[u_new] > 0:
                cur = sum(1 for n in g.nodes() if g.degree[n] == u_new)
                if cur >= self.target_counts[u_new]:
                    return False
        return True


if __name__ == "__main__":
    # Self-test on a small lattice
    import json
    info = verify_diamond_lattice(4, 4, 4)
    print(json.dumps(info, indent=2))
