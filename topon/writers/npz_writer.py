"""
NPZ Writer for Topon
====================
Exports the dual-graph (chains + crosslinks) of a topon polymer network as a
single compressed `.npz` file. Schema and rationale: see
``internal/specs/npz_format.md``. Designed to be loaded directly into
PyG / PyTorch Geometric pipelines.

Mirror of :func:`topon.writers.graphml_writer.write_graphml`'s dual-graph
transformation:
  - Chain nodes (one per edge of the topon MultiGraph)
  - Crosslink nodes (one per junction node in the topon MultiGraph)
  - Chemical edges: each chain connects its two end-crosslinks
  - Entanglement edges: chain–chain links from ``edge.entangled_with``

Stress / strain arrays are returned empty by default — those typically come
from the LAMMPS run, not the generation step. Pass them in if you have them.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import networkx as nx
import numpy as np


# Node-feature column order (must match the spec). Index in node_features[:, k].
_FEATURE_COLUMNS = (
    "type",            # 0 = polymer (chain), 1 = crosslinker
    "length",          # DP for chains, 1 for crosslinkers
    "contour_length",  # DP * bond_length (NaN if unknown)
    "rg",              # radius of gyration (NaN if unknown)
    "COMX",
    "COMY",
    "COMZ",
    "node_degree",
)


def write_npz(
    G: nx.MultiGraph,
    output_path: str,
    *,
    dp: int = 50,
    bond_length: float = 1.0,
    dims: Optional[np.ndarray] = None,
    strain: Optional[np.ndarray] = None,
    stress: Optional[np.ndarray] = None,
) -> Path:
    """Write the dual graph of ``G`` as a compressed NPZ.

    Args:
        G: Topon ``MultiGraph`` (nodes = crosslinks, edges = chains).
        output_path: Path to the output ``.npz`` file.
        dp: Default degree of polymerization if not stored per-edge.
        bond_length: Used to compute ``contour_length = length * bond_length``
            for chain nodes (NaN remains for crosslinker nodes).
        dims: Box dimensions ``[Lx, Ly, Lz]`` (writes box as
            ``[0, Lx, 0, Ly, 0, Lz]``); if ``None``, box is left as NaN.
        strain: Optional 1-D float32 array (e.g. from a LAMMPS dump).
        stress: Optional 1-D float32 array (e.g. from a LAMMPS dump).

    Returns:
        Path to the written ``.npz`` file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    edges = list(G.edges(keys=True, data=True))

    # Chain nodes: one per MultiGraph edge.
    edge_to_id: dict[tuple[int, int, int], int] = {}
    chain_features: list[list[float]] = []   # rows of node_features
    chain_ids: list[int] = []                # original mol_id (per-chain)

    # Offset chain IDs above the max crosslink ID so they don't collide.
    max_xlink_id = max(max(u, v) for u, v, _ in G.edges) if G.edges else 0
    next_id = max_xlink_id + 1

    for u, v, key, data in edges:
        edge_to_id[(u, v, key)] = next_id
        chain_ids.append(next_id)
        length = float(data.get("dp", dp))
        contour = length * float(bond_length)
        # Center of mass for a chain: midpoint of its end-crosslink positions.
        pu = G.nodes.get(u, {}).get("pos")
        pv = G.nodes.get(v, {}).get("pos")
        if pu is not None and pv is not None and len(pu) == 3 and len(pv) == 3:
            comx = 0.5 * (float(pu[0]) + float(pv[0]))
            comy = 0.5 * (float(pu[1]) + float(pv[1]))
            comz = 0.5 * (float(pu[2]) + float(pv[2]))
        else:
            comx = comy = comz = float("nan")
        chain_features.append([
            0.0,            # type = polymer
            length,
            contour,
            float("nan"),   # rg unknown at generation time
            comx,
            comy,
            comz,
            2.0,            # chain nodes have degree 2 in the dual graph
        ])
        next_id += 1

    # Crosslink nodes: one per topon-graph node that participates in an edge.
    crosslink_ids = sorted({u for u, v, _ in G.edges} | {v for u, v, _ in G.edges})
    crosslink_features: list[list[float]] = []
    for xid in crosslink_ids:
        attrs = G.nodes.get(xid, {})
        pos = attrs.get("pos")
        if pos is not None and len(pos) == 3:
            cx, cy, cz = float(pos[0]), float(pos[1]), float(pos[2])
        else:
            cx = cy = cz = float("nan")
        # Crosslink "node_degree" in the dual = number of chains attached
        # = original graph degree. Same as G.degree(xid).
        crosslink_features.append([
            1.0,            # type = crosslinker
            1.0,            # length placeholder
            0.0,            # contour_length
            0.0,            # rg
            cx,
            cy,
            cz,
            float(G.degree(xid)),
        ])

    n_polymer = len(chain_features)
    n_crosslinker = len(crosslink_features)

    # Stack node_features in [polymer rows, crosslinker rows] order.
    if n_polymer + n_crosslinker == 0:
        node_features = np.zeros((0, len(_FEATURE_COLUMNS)), dtype=np.float32)
        node_ids = np.zeros(0, dtype=np.int32)
    else:
        node_features = np.asarray(
            chain_features + crosslink_features, dtype=np.float32
        )
        node_ids = np.asarray(chain_ids + crosslink_ids, dtype=np.int32)

    # Build chemical edges: each chain <-> its two end-crosslinks.
    chemical_edges: list[tuple[int, int]] = []
    for (u, v, key), cid in edge_to_id.items():
        chemical_edges.append((cid, u))
        chemical_edges.append((cid, v))

    # Build entanglement edges (chain <-> chain), without double-counting.
    entanglement_edges: list[tuple[int, int]] = []
    seen: set[frozenset] = set()
    for u, v, key, data in edges:
        ew = data.get("entangled_with")
        if not ew:
            continue
        ew = tuple(ew)
        e1 = (u, v, key)
        e2 = ew
        pair = frozenset([frozenset(e1[:2]), frozenset(e2[:2])])
        if pair in seen:
            continue
        seen.add(pair)
        cid1 = edge_to_id.get(e1) or edge_to_id.get((e1[1], e1[0], e1[2]))
        cid2 = edge_to_id.get(e2) or edge_to_id.get((e2[1], e2[0], e2[2]))
        if cid1 is None or cid2 is None:
            continue
        count = data.get("entanglement_count", 1)
        for _ in range(count):
            entanglement_edges.append((cid1, cid2))

    # Build edge_index and edge_type arrays. PyG expects undirected edges
    # represented as both directions, so we emit (i, j) and (j, i).
    #
    # NOTE: the (src, tgt) values collected here are in the *original
    # ID space* -- ``cid`` chain dual-ids (offset above max_xlink_id) and
    # ``u``/``v`` crosslink node-ids (sparse: vacancy removal leaves
    # gaps). They are NOT yet 0-based row positions into node_features.
    # The remap step below converts them; see the bugfix block.
    edge_pairs: list[tuple[int, int, int]] = []  # (src, tgt, type_int)
    for src, tgt in chemical_edges:
        edge_pairs.append((src, tgt, 0))
        edge_pairs.append((tgt, src, 0))
    for src, tgt in entanglement_edges:
        edge_pairs.append((src, tgt, 1))
        edge_pairs.append((tgt, src, 1))

    if edge_pairs:
        ei = np.asarray([(p[0], p[1]) for p in edge_pairs], dtype=np.int32).T
        et = np.asarray([p[2] for p in edge_pairs], dtype=np.int32)
    else:
        ei = np.zeros((2, 0), dtype=np.int32)
        et = np.zeros((0,), dtype=np.int32)

    # --- BUGFIX: remap edge_index from original-ID space -> 0-based row
    # positions into node_features.
    #
    # node_features rows are [chain rows ... crosslink rows]; node_ids[i]
    # is the original simulation/mol ID of row i. edge_index built above
    # holds those original IDs, but PyG (and the spec) require edge_index
    # to be 0-based positions into node_features. Crosslink IDs are
    # sparse and chain IDs are offset above max_xlink_id, so the raw-ID
    # range exceeds N -> PyTorch Geometric throws CUDA "index out of
    # bounds". Remap with a dense ID->position lookup table.
    if ei.shape[1] > 0:
        max_id = int(node_ids.max()) + 1
        id_to_pos = np.full(max_id, -1, dtype=np.int64)
        id_to_pos[node_ids] = np.arange(len(node_ids), dtype=np.int64)
        ei = id_to_pos[ei].astype(np.int32)
        if (ei < 0).any():
            raise ValueError(
                "npz_writer: edge_index references an ID absent from "
                "node_ids -- cannot remap to a row position. This "
                "indicates an upstream graph-construction bug."
            )

    # Box bounds: [xlo, xhi, ylo, yhi, zlo, zhi].
    if dims is not None and np.asarray(dims).size >= 3:
        d = np.asarray(dims, dtype=np.float32).flatten()
        box = np.asarray(
            [0.0, d[0], 0.0, d[1], 0.0, d[2]], dtype=np.float32
        )
    else:
        box = np.full(6, np.nan, dtype=np.float32)

    strain_arr = (
        np.asarray(strain, dtype=np.float32)
        if strain is not None
        else np.zeros(0, dtype=np.float32)
    )
    stress_arr = (
        np.asarray(stress, dtype=np.float32)
        if stress is not None
        else np.zeros(0, dtype=np.float32)
    )

    np.savez_compressed(
        output_path,
        node_features=node_features,
        node_ids=node_ids,
        edge_index=ei,
        edge_type=et,
        box=box,
        strain=strain_arr,
        stress=stress_arr,
        n_polymer=np.int32(n_polymer),
        n_crosslinker=np.int32(n_crosslinker),
    )

    print(f"  NPZ written to {output_path}")
    print(f"    Polymer nodes : {n_polymer}")
    print(f"    Crosslink nodes : {n_crosslinker}")
    print(f"    Chemical edges (directed) : {len(chemical_edges) * 2}")
    print(f"    Entanglement edges (directed): {len(entanglement_edges) * 2}")

    return output_path
