"""
Graph loader for Topon.

Handles loading topology from various file formats:

* gpickle / .nodes+.edges -- basic topology (lattice + connectivity).
* graphml / npz           -- post-assignment dual graph: chains + crosslinks
                              with DP, entanglements (multiplicity preserved),
                              and crosslink positions. Suitable for skipping
                              the topology + analysis + assignment stages and
                              going straight into chemistry/conformation/output.
"""

import pickle
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Optional, Union

import networkx as nx
import numpy as np
import pandas as pd


def load_graph(
    gpickle_path: Optional[Union[str, Path]] = None,
    nodes_path: Optional[Union[str, Path]] = None,
    edges_path: Optional[Union[str, Path]] = None,
) -> tuple[nx.MultiGraph, Optional[np.ndarray]]:
    """
    Load a topology graph from file(s).
    
    Args:
        gpickle_path: Path to a .gpickle file (takes precedence).
        nodes_path: Path to a .nodes file.
        edges_path: Path to a .edges file.
        
    Returns:
        Tuple of (NetworkX MultiGraph, box dimensions array or None).
        
    Raises:
        ValueError: If no valid file paths provided.
        FileNotFoundError: If specified files don't exist.
    """
    if gpickle_path:
        return _load_from_gpickle(gpickle_path)
    elif nodes_path and edges_path:
        return _load_from_nodes_edges(nodes_path, edges_path)
    else:
        raise ValueError(
            "Must provide either gpickle_path or both nodes_path and edges_path"
        )


def _load_from_gpickle(path: Union[str, Path]) -> tuple[nx.MultiGraph, Optional[np.ndarray]]:
    """
    Load graph from a .gpickle file.
    
    Args:
        path: Path to .gpickle file.
        
    Returns:
        Tuple of (graph, dims).
    """
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"Gpickle file not found: {path}")
    
    with open(path, "rb") as f:
        data = pickle.load(f)
    
    # Handle different gpickle formats
    if isinstance(data, tuple) and len(data) == 2:
        # Format: (graph, dims)
        G, dims = data
    elif isinstance(data, nx.Graph):
        # Format: just the graph
        G = data
        dims = infer_dims_from_graph(G)
    elif isinstance(data, dict) and "graph" in data:
        # Format: dict with 'graph' and optionally 'dims'
        G = data["graph"]
        dims = data.get("dims")
    else:
        raise ValueError(f"Unrecognized gpickle format in {path}")
    
    # Ensure it's a MultiGraph
    if not isinstance(G, nx.MultiGraph):
        G = nx.MultiGraph(G)

    # A recorded box wins over a stored dims: it is the generator's exact
    # cell, whereas a dims saved alongside an older graph may have come
    # from the positional fallback (and be half a cell too large on any
    # lattice with fractional basis sites).
    if G.graph.get("box") is not None:
        dims = infer_dims_from_graph(G)

    # Remove vacancies (degree-0 nodes)
    n_removed = remove_vacancies(G)
    
    print(f"Loaded graph from {path.name}")
    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    if n_removed > 0:
        print(f"  Removed {n_removed} vacancies (degree-0 nodes)")
    
    return G, dims


def _load_from_nodes_edges(
    nodes_path: Union[str, Path],
    edges_path: Union[str, Path]
) -> tuple[nx.MultiGraph, Optional[np.ndarray]]:
    """
    Load graph from .nodes and .edges files.
    
    File formats:
    - .nodes: NodeID X Y Z Degree (whitespace-separated, # comments)
    - .edges: Node1 Node2 (whitespace-separated, # comments)
    
    Args:
        nodes_path: Path to .nodes file.
        edges_path: Path to .edges file.
        
    Returns:
        Tuple of (graph, dims).
    """
    nodes_path = Path(nodes_path)
    edges_path = Path(edges_path)
    
    if not nodes_path.exists():
        raise FileNotFoundError(f"Nodes file not found: {nodes_path}")
    if not edges_path.exists():
        raise FileNotFoundError(f"Edges file not found: {edges_path}")
    
    # Load nodes
    nodes_df = pd.read_csv(
        nodes_path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["id", "x", "y", "z", "degree"]
    )
    
    # Load edges
    edges_df = pd.read_csv(
        edges_path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["u", "v"]
    )
    
    # Build graph
    G = nx.MultiGraph()

    for _, row in nodes_df.iterrows():
        G.add_node(
            int(row["id"]),
            pos=(float(row["x"]), float(row["y"]), float(row["z"]))
        )

    for _, row in edges_df.iterrows():
        u, v = int(row["u"]), int(row["v"])
        if G.has_node(u) and G.has_node(v):
            G.add_edge(u, v)

    # An optional "# BOX Lx Ly Lz" header carries the true periodic cell.
    # Files written without one fall back to the positional heuristic.
    box = read_box_header(nodes_path)
    if box is not None:
        G.graph["box"] = box

    # "# PERIODICITY 110" records which axes are open. Absent means fully
    # periodic, which is what every file predating the header represents.
    axes = read_periodicity_header(nodes_path)
    if axes is not None:
        G.graph["periodicity"] = axes

    # Infer dimensions from positions
    dims = infer_dims_from_graph(G)
    
    # Remove vacancies (degree-0 nodes)
    n_removed = remove_vacancies(G)
    
    print(f"Loaded graph from {nodes_path.name} + {edges_path.name}")
    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    if n_removed > 0:
        print(f"  Removed {n_removed} vacancies (degree-0 nodes)")
    
    return G, dims


def remove_vacancies(G: nx.Graph) -> int:
    """
    Remove degree-0 nodes (vacancies) from graph.
    
    Vacancies are lattice positions with no edges - they should not
    become atoms in the simulation.
    
    Args:
        G: Graph to modify in-place.
        
    Returns:
        Number of nodes removed.
    """
    vacancies = [n for n in G.nodes() if G.degree(n) == 0]
    if vacancies:
        G.remove_nodes_from(vacancies)
    return len(vacancies)


# Header line the topology generators write into ``.nodes`` files to record
# the exact periodic cell, e.g. "# BOX 6 6 6". Held as a module constant so
# the Python reader/writer and the C generator agree on one spelling.
BOX_HEADER_KEY = "BOX"


PERIODICITY_HEADER_KEY = "PERIODICITY"


def format_periodicity_header(periodicity) -> str:
    """Render per-axis boundaries as a ``.nodes`` header line.

    Args:
        periodicity: Iterable of three truthy values, one per axis.

    Returns:
        The header line, e.g. ``"# PERIODICITY 110"``, without a newline.
    """
    digits = "".join("1" if p else "0" for p in periodicity)
    return f"# {PERIODICITY_HEADER_KEY} {digits}"


def read_periodicity_header(path: Union[str, Path]):
    """Read the ``# PERIODICITY 110`` header from a .nodes file.

    Returns ``(px, py, pz)`` booleans, or None when absent or malformed.
    Absent means fully periodic, which is what every file written before
    this header existed represents.
    """
    try:
        with open(path) as f:
            for line in f:
                if not line.startswith("#"):
                    break
                parts = line.lstrip("#").split()
                if len(parts) == 2 and parts[0].upper() == PERIODICITY_HEADER_KEY:
                    digits = parts[1].strip()
                    if len(digits) == 3 and set(digits) <= {"0", "1"}:
                        return tuple(c == "1" for c in digits)
                    return None
    except OSError:
        return None
    return None


def format_box_header(box) -> str:
    """Render a 3-component box as its ``.nodes`` header line.

    Args:
        box: Iterable of three box lengths in lattice units.

    Returns:
        The header line, without a trailing newline.
    """
    lx, ly, lz = (float(v) for v in box)
    return f"# {BOX_HEADER_KEY} {lx:g} {ly:g} {lz:g}"


def read_box_header(path: Union[str, Path]) -> Optional[tuple[float, float, float]]:
    """Read the ``# BOX Lx Ly Lz`` header from a .nodes file.

    Only the leading comment block is scanned, so this stops after a
    handful of lines on any file. Returns None when the header is absent
    or malformed, which is the expected case for files written before
    generators recorded their box.

    Args:
        path: Path to a ``.nodes`` file.

    Returns:
        ``(Lx, Ly, Lz)`` in lattice units, or None.
    """
    try:
        with open(path) as f:
            for line in f:
                if not line.startswith("#"):
                    break
                parts = line.lstrip("#").split()
                if len(parts) == 4 and parts[0].upper() == BOX_HEADER_KEY:
                    try:
                        lx, ly, lz = (float(v) for v in parts[1:4])
                    except ValueError:
                        return None
                    if min(lx, ly, lz) <= 0:
                        return None
                    return (lx, ly, lz)
    except OSError:
        return None
    return None


def save_nodes_edges(
    G: nx.Graph,
    nodes_path: Union[str, Path],
    edges_path: Union[str, Path],
    box=None,
    periodicity=None,
) -> None:
    """Write a graph in the ``.nodes`` / ``.edges`` format.

    Matches what the C generator emits, plus a ``# BOX`` header carrying
    the exact periodic cell so a reload does not have to guess it, and a
    ``# PERIODICITY`` header when any axis is open.

    Args:
        G: Graph whose nodes carry ``pos``.
        nodes_path: Destination ``.nodes`` path.
        edges_path: Destination ``.edges`` path.
        box: Periodic cell to record. Defaults to ``G.graph["box"]``;
             omitted from the header when neither is available.
        periodicity: Per-axis boundaries. Defaults to
             ``G.graph["periodicity"]``. Written only when an axis is
             actually open, so fully periodic files keep the exact
             format they had before this existed.
    """
    nodes_path = Path(nodes_path)
    edges_path = Path(edges_path)
    nodes_path.parent.mkdir(parents=True, exist_ok=True)
    edges_path.parent.mkdir(parents=True, exist_ok=True)

    if box is None:
        box = G.graph.get("box")
    if periodicity is None:
        periodicity = G.graph.get("periodicity")

    with open(nodes_path, "w") as f:
        if box is not None:
            f.write(format_box_header(box) + "\n")
        if periodicity is not None and not all(periodicity):
            f.write(format_periodicity_header(periodicity) + "\n")
        f.write("# NodeID X Y Z Degree\n")
        for node in sorted(G.nodes()):
            x, y, z = G.nodes[node].get("pos", (0.0, 0.0, 0.0))
            f.write(f"{node} {x:f} {y:f} {z:f} {G.degree(node)}\n")

    with open(edges_path, "w") as f:
        f.write("# Node1 Node2\n")
        for u, v in sorted(G.edges()):
            f.write(f"{u} {v}\n")


def infer_dims_from_graph(G: nx.Graph) -> Optional[np.ndarray]:
    """
    Infer box dimensions from node positions.
    
    Args:
        G: Graph with 'pos' node attributes.
        
    Returns:
        Box dimensions as numpy array, or None if no positions.

    Notes:
        Prefers the exact cell recorded by the generator in
        ``G.graph["box"]``. The ``max - min + 1`` fallback below is only
        correct when every site sits on an integer coordinate with unit
        spacing, i.e. simple cubic. BCC, FCC and Diamond place basis
        sites at fractional offsets, so the fallback overshoots the true
        cell by that offset: a 4x4x4 BCC or FCC lattice reports 4.5
        instead of 4.0. Since this value feeds every minimum-image
        calculation downstream, that overshoot makes roughly a third of
        BCC edges (and a quarter of FCC edges) resolve to the wrong
        periodic replica and be built at twice their true bond length.
        Generators therefore record the true cell explicitly; the
        fallback survives only for graphs written before they did
        (old gpickles, ``.nodes`` files with no ``# BOX`` header).
    """
    box = G.graph.get("box")
    if box is not None:
        arr = np.asarray(box, dtype=float).ravel()
        if arr.size == 3 and np.all(np.isfinite(arr)) and np.all(arr > 0):
            return arr

    positions = []
    for node, data in G.nodes(data=True):
        if "pos" in data:
            positions.append(data["pos"])
    
    if not positions:
        return None
    
    positions = np.array(positions)
    
    # Assume box starts at 0, dimensions are max + 1 (for lattice spacing)
    # This is a heuristic - actual dims should be stored in gpickle
    max_pos = positions.max(axis=0)
    min_pos = positions.min(axis=0)
    
    # For integer lattice positions, dims = max - min + 1
    dims = max_pos - min_pos + 1
    
    return dims


def save_graph(
    G: nx.Graph,
    output_path: Union[str, Path],
    dims: Optional[np.ndarray] = None
) -> None:
    """
    Save graph to a .gpickle file.
    
    Args:
        G: NetworkX graph to save.
        output_path: Path for output file.
        dims: Optional box dimensions to save with graph.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    data = (G, dims) if dims is not None else G
    
    with open(output_path, "wb") as f:
        pickle.dump(data, f)
    
    print(f"Saved graph to {output_path}")


def get_node_positions(G: nx.Graph) -> dict[int, np.ndarray]:
    """
    Extract node positions as a dictionary.

    Args:
        G: Graph with 'pos' node attributes.

    Returns:
        Dict mapping node ID to position array.
    """
    positions = {}
    for node, data in G.nodes(data=True):
        if "pos" in data:
            positions[node] = np.array(data["pos"])
    return positions


# ======================================================================
# Dual-graph loaders (graphml + npz)
#
# Both formats are written by ``topon.writers.{graphml_writer,npz_writer}``
# in the dual representation:
#   * polymer nodes   (one per topon chain) carry DP via "length"
#   * crosslinker nodes (one per junction) carry pos via COMX/Y/Z
#   * chemical edges connect chain <-> its two end-crosslinks
#   * entanglement edges connect chain <-> chain, replicated `count` times
#     to preserve multiplicity
#
# The loaders here invert that transformation to rebuild the topon
# MultiGraph (junctions = nodes, chains = edges with dp / entangled_with /
# entanglement_count attributes). With that graph plus dims, the chemistry
# + conformation + output stages can be re-run without going through
# topology generation or assignment.
# ======================================================================


def load_graphml(
    path: Union[str, Path],
) -> tuple[nx.MultiGraph, Optional[np.ndarray]]:
    """Load a topon graphml dual graph back into a MultiGraph + dims.

    Reverses :func:`topon.writers.graphml_writer.write_graphml`.

    Args:
        path: Path to the ``.graphml`` file.

    Returns:
        ``(G, dims)`` where ``G`` is a ``nx.MultiGraph`` with crosslink
        nodes (carrying ``pos``) and chain edges (carrying ``dp``,
        ``entangled_with``, ``entanglement_count``). ``dims`` is the
        inferred box dimensions (from graph attributes if present, else
        from node positions).
    """
    path = Path(path)
    tree = ET.parse(path)
    root = tree.getroot()
    ns = {"g": "http://graphml.graphdrawing.org/xmlns"}

    # Parse <key> declarations so we can map data-element key IDs back to
    # their semantic names (e.g. "length", "COMX", "edge_type", "xhi").
    key_meta: dict[str, tuple[str, str, str]] = {}
    for key in root.findall("g:key", ns):
        key_meta[key.get("id")] = (
            key.get("for"),
            key.get("attr.name"),
            key.get("attr.type"),
        )

    graph_el = root.find("g:graph", ns)
    if graph_el is None:
        raise ValueError(f"No <graph> element in {path}")

    polymer_attrs: dict[int, dict] = {}     # chain dual-id -> attrs
    crosslink_attrs: dict[int, dict] = {}   # junction id   -> attrs

    for node in graph_el.findall("g:node", ns):
        nid = int(node.get("id"))
        attrs: dict = {}
        for d in node.findall("g:data", ns):
            kid = d.get("key")
            if kid not in key_meta:
                continue
            _, name, type_ = key_meta[kid]
            txt = d.text
            if txt is None or txt == "NaN":
                val = float("nan")
            elif type_ == "long":
                val = int(txt)
            elif type_ in ("double", "float"):
                val = float(txt)
            else:
                val = txt
            attrs[name] = val
        ntype = attrs.get("type", "")
        if ntype == "polymer":
            polymer_attrs[nid] = attrs
        elif ntype == "crosslinker":
            crosslink_attrs[nid] = attrs
        else:
            # Fall back on whether the node has COMX/Y/Z that aren't NaN
            # (crosslinks carry positions; chains do not).
            comx = attrs.get("COMX")
            if isinstance(comx, float) and not np.isnan(comx):
                crosslink_attrs[nid] = attrs
            else:
                polymer_attrs[nid] = attrs

    # Pull edges, split by edge_type. Chemical edges link a chain to a
    # crosslink; entanglement edges link chain to chain (replicated for
    # multiplicity).
    chemical_pairs: list[tuple[int, int]] = []
    entanglement_pairs: list[tuple[int, int]] = []
    for e in graph_el.findall("g:edge", ns):
        src = int(e.get("source"))
        tgt = int(e.get("target"))
        etype = "chemical"
        for d in e.findall("g:data", ns):
            kid = d.get("key")
            if kid in key_meta and key_meta[kid][1] == "edge_type":
                etype = d.text or "chemical"
        if etype == "chemical":
            chemical_pairs.append((src, tgt))
        elif etype == "entanglement":
            entanglement_pairs.append((src, tgt))

    G, cid_to_edge = _build_multigraph_from_dual(
        polymer_attrs, crosslink_attrs, chemical_pairs, entanglement_pairs
    )

    # Graph-level box bounds (xlo/xhi/ylo/yhi/zlo/zhi). All NaN by default
    # in the writer, so fall back to inferring from positions.
    box: dict[str, float] = {}
    for d in graph_el.findall("g:data", ns):
        kid = d.get("key")
        if kid in key_meta and key_meta[kid][0] == "graph":
            try:
                box[key_meta[kid][1]] = float(d.text)
            except (TypeError, ValueError):
                pass
    dims: Optional[np.ndarray] = None
    needed = {"xlo", "xhi", "ylo", "yhi", "zlo", "zhi"}
    if needed <= set(box) and not any(np.isnan(box[k]) for k in needed):
        dims = np.array(
            [box["xhi"] - box["xlo"], box["yhi"] - box["ylo"], box["zhi"] - box["zlo"]],
            dtype=float,
        )
    else:
        dims = infer_dims_from_graph(G)

    n_ent_pairs = sum(
        1 for _, _, _, data in G.edges(keys=True, data=True)
        if data.get("entangled_with")
    ) // 2  # symmetric
    print(f"Loaded graphml from {path.name}")
    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}, "
          f"Entangled pairs: {n_ent_pairs}")
    return G, dims


def _sc_positions_from_ids(
    xl_ids: list[int],
    box: np.ndarray,
    node_features: np.ndarray,
    node_ids: np.ndarray,
    edge_index: np.ndarray,
    edge_type: np.ndarray,
    n_polymer: int,
) -> Optional[dict[int, tuple[float, float, float]]]:
    """Reconstruct simple-cubic lattice coords from node ids + box.

    The C generator lays SC sites out as ``id = x + y*Nx + z*Nx*Ny`` with
    ``coords = (x, y, z)``; ``write_npz`` stores ``box = [0, Nx, 0, Ny,
    0, Nz]``. Inverting is therefore exact -- but only if the graph
    really came from an SC lattice, so we verify: every chain must join
    two sites that are nearest neighbours under periodic boundaries.

    Returns ``{node_id: (x, y, z)}``, or ``None`` when the box is
    unusable, an id falls outside the lattice, or the neighbour check
    fails (e.g. a BCC/FCC/Diamond graph, or a non-lattice topology).
    """
    if box is None or box.size != 6 or bool(np.isnan(box).any()):
        return None
    Nx, Ny, Nz = int(round(float(box[1]))), int(round(float(box[3]))), int(round(float(box[5])))
    if min(Nx, Ny, Nz) <= 0:
        return None
    n_sites = Nx * Ny * Nz
    if not xl_ids or max(xl_ids) >= n_sites or min(xl_ids) < 0:
        return None

    pos = {
        i: (float(i % Nx), float((i // Nx) % Ny), float(i // (Nx * Ny)))
        for i in xl_ids
    }

    # Validate: each chain's two junctions must be lattice neighbours.
    from collections import defaultdict
    chain_to_xl: dict[int, set[int]] = defaultdict(set)
    chem = edge_index[:, edge_type == 0]
    for k in range(chem.shape[1]):
        a, b = int(chem[0, k]), int(chem[1, k])
        if a < n_polymer <= b:
            chain_to_xl[a].add(int(node_ids[b]))
        elif b < n_polymer <= a:
            chain_to_xl[b].add(int(node_ids[a]))
    checked = 0
    for xls in chain_to_xl.values():
        if len(xls) != 2:
            continue
        p, q = (pos[i] for i in xls)
        dims_ = (Nx, Ny, Nz)
        dist = 0.0
        for ax in range(3):
            d = abs(p[ax] - q[ax])
            dist += min(d, dims_[ax] - d)
        if abs(dist - 1.0) > 1e-6:
            return None                      # not an SC nearest-neighbour graph
        checked += 1
    if checked == 0:
        return None
    return pos


def load_npz(
    path: Union[str, Path],
) -> tuple[nx.MultiGraph, Optional[np.ndarray]]:
    """Load a topon npz dual graph back into a MultiGraph + dims.

    Reverses :func:`topon.writers.npz_writer.write_npz`.

    Args:
        path: Path to the ``.npz`` file.

    Returns:
        ``(G, dims)`` -- same shape as :func:`load_graphml`.
    """
    path = Path(path)
    data = np.load(path)
    node_features = data["node_features"]
    node_ids = data["node_ids"]
    edge_index = data["edge_index"]
    edge_type = data["edge_type"]
    box = data["box"]
    n_polymer = int(data["n_polymer"])
    n_crosslinker = int(data["n_crosslinker"])

    # Feature columns (must match npz_writer._FEATURE_COLUMNS):
    #   v1: [type, length, contour_length, rg, COMX, COMY, COMZ, node_degree]
    #   v2: [type, length, contour_length, rg, COMX, COMY, COMZ,
    #        chem_degree, phys_degree, frac_ext]
    # In v2 the COM columns are deliberately NaN (they are conformation
    # outputs, filled in after a LAMMPS run -- not known at generation
    # time). The downstream chemistry/conformation stages nevertheless
    # need a 3-D embedding of the junctions, so when COM is NaN we
    # reconstruct the ORIGINAL LATTICE COORDINATES from node_ids + box.
    #
    # The C generator numbers simple-cubic sites as
    #     id = x + y*Nx + z*Nx*Ny        coords = (x, y, z)
    # and write_npz stores box = [0, Nx, 0, Ny, 0, Nz], so the mapping
    # inverts exactly. ``_sc_positions_from_ids`` validates the result
    # (every chain must join two lattice neighbours) and returns None if
    # the graph is not an SC lattice, in which case COM/NaN is kept.
    polymer_attrs: dict[int, dict] = {}
    crosslink_attrs: dict[int, dict] = {}

    for i in range(n_polymer):
        nid = int(node_ids[i])
        polymer_attrs[nid] = {"length": int(node_features[i, 1])}

    xl_ids = [int(node_ids[i]) for i in range(n_polymer, n_polymer + n_crosslinker)]
    com_is_nan = (
        n_crosslinker > 0
        and bool(np.isnan(node_features[n_polymer:, 4:7]).all())
    )
    recovered = (
        _sc_positions_from_ids(xl_ids, box, node_features, node_ids,
                               edge_index, edge_type, n_polymer)
        if com_is_nan else None
    )
    if com_is_nan and recovered is None:
        print("  WARNING: COM columns are NaN and lattice positions could "
              "not be reconstructed; junction coordinates will be NaN and "
              "any downstream conformation/LAMMPS build will be invalid.")

    for i in range(n_polymer, n_polymer + n_crosslinker):
        nid = int(node_ids[i])
        if recovered is not None:
            cx, cy, cz = recovered[nid]
        else:
            cx = float(node_features[i, 4])
            cy = float(node_features[i, 5])
            cz = float(node_features[i, 6])
        crosslink_attrs[nid] = {"COMX": cx, "COMY": cy, "COMZ": cz}

    # Edges are stored bi-directionally; reduce to unordered pairs.
    #
    # edge_index holds 0-based ROW POSITIONS into node_features (PyG
    # convention -- see the npz_writer remap bugfix). The rest of this
    # loader keys nodes by their original IDs, so map each edge endpoint
    # back through node_ids: node_ids[row_position] -> original ID.
    chemical_pairs_set: set[tuple[int, int]] = set()
    entanglement_pairs_list: list[tuple[int, int]] = []
    n_edges = edge_index.shape[1]
    for k in range(n_edges):
        src = int(node_ids[int(edge_index[0, k])])
        tgt = int(node_ids[int(edge_index[1, k])])
        et = int(edge_type[k])
        if et == 0:                              # chemical (bidirectional copy)
            chemical_pairs_set.add(tuple(sorted((src, tgt))))
        elif et == 1 and src < tgt:              # entanglement: keep one dir
            entanglement_pairs_list.append((src, tgt))

    G, cid_to_edge = _build_multigraph_from_dual(
        polymer_attrs,
        crosslink_attrs,
        list(chemical_pairs_set),
        entanglement_pairs_list,
    )

    # Box: [xlo, xhi, ylo, yhi, zlo, zhi]
    if box.size == 6 and not np.isnan(box[1]):
        dims = np.array(
            [box[1] - box[0], box[3] - box[2], box[5] - box[4]], dtype=float
        )
    else:
        dims = infer_dims_from_graph(G)

    n_ent_pairs = sum(
        1 for _, _, _, data in G.edges(keys=True, data=True)
        if data.get("entangled_with")
    ) // 2
    print(f"Loaded npz from {path.name}")
    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}, "
          f"Entangled pairs: {n_ent_pairs}")
    return G, dims


def _build_multigraph_from_dual(
    polymer_attrs: dict[int, dict],
    crosslink_attrs: dict[int, dict],
    chemical_pairs: list[tuple[int, int]],
    entanglement_pairs: list[tuple[int, int]],
) -> tuple[nx.MultiGraph, dict[int, tuple[int, int, int]]]:
    """Shared dual-graph -> MultiGraph reconstruction (graphml & npz).

    * Each polymer (chain) dual-node becomes an edge in the topon graph;
      its two chemical-edge neighbours are the chain's u/v crosslinks.
    * Crosslink positions come from the COMX/Y/Z fields.
    * Entanglement multiplicity is recovered by counting how many times
      a chain-chain pair appears among ``entanglement_pairs``.

    Returns the rebuilt graph plus a ``{chain_dual_id: (u, v, key)}`` map
    so the caller can wire `entangled_with` after edges are added.
    """
    G = nx.MultiGraph()
    for xid, attrs in crosslink_attrs.items():
        pos = (
            float(attrs.get("COMX", 0.0)),
            float(attrs.get("COMY", 0.0)),
            float(attrs.get("COMZ", 0.0)),
        )
        G.add_node(xid, pos=pos)

    # For each chain, find its two crosslink endpoints from chemical pairs.
    chain_endpoints: dict[int, list[int]] = {}
    for a, b in chemical_pairs:
        chain, junction = (a, b) if a in polymer_attrs else (b, a)
        if chain not in polymer_attrs:
            # Neither end is a known polymer dual-node; skip
            continue
        chain_endpoints.setdefault(chain, []).append(junction)

    # Add edges in chain-id sorted order so the resulting (u, v, key)
    # assignment is deterministic across graphml and npz loads of the
    # same network.
    cid_to_edge: dict[int, tuple[int, int, int]] = {}
    for cid in sorted(chain_endpoints):
        endpoints = chain_endpoints[cid]
        if len(endpoints) != 2:
            print(f"  [warn] chain {cid} has {len(endpoints)} crosslink "
                  f"endpoints (expected 2); skipping")
            continue
        # Canonicalise the (u, v) order so the same chain always yields
        # the same edge key regardless of how chemical_pairs were ordered
        # in the source file.
        u, v = sorted(endpoints)
        dp = int(polymer_attrs[cid].get("length", 1))
        key = G.add_edge(u, v, dp=dp)
        cid_to_edge[cid] = (u, v, key)

    # Recover entanglement multiplicities by counting chain pairs.
    ent_counts: Counter = Counter()
    for cid1, cid2 in entanglement_pairs:
        pair = tuple(sorted((cid1, cid2)))
        ent_counts[pair] += 1

    for (cid1, cid2), count in ent_counts.items():
        if cid1 not in cid_to_edge or cid2 not in cid_to_edge:
            continue
        e1 = cid_to_edge[cid1]
        e2 = cid_to_edge[cid2]
        G[e1[0]][e1[1]][e1[2]]["entangled_with"] = e2
        G[e1[0]][e1[1]][e1[2]]["entanglement_count"] = count
        G[e2[0]][e2[1]][e2[2]]["entangled_with"] = e1
        G[e2[0]][e2[1]][e2[2]]["entanglement_count"] = count

    return G, cid_to_edge
