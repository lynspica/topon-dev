"""
GraphML Writer for Topon
========================
Exports a topon polymer network graph as a GraphML file.

In topon's internal representation:
  - Nodes = crosslink sites
  - Edges = polymer chains

The output GraphML uses the dual (line-graph) representation:
  - Nodes = polymer chains  (with type, length, rg, COMX/Y/Z)
  - Edges = chemical bonds   (two chains sharing a crosslink)
          + entanglement      (entangled chain pairs, with multiplicity)
"""

from pathlib import Path
from typing import Optional
from xml.etree.ElementTree import Element, SubElement, ElementTree
from collections import defaultdict

import networkx as nx
import numpy as np


def write_graphml(
    G: nx.MultiGraph,
    output_path: str,
    dp: int = 50,
    dims: Optional[np.ndarray] = None,
) -> Path:
    """
    Write a topon graph to GraphML with the dual-graph transformation.

    Args:
        G: Topon MultiGraph (nodes=crosslinks, edges=chains).
        output_path: Path for the output .graphml file.
        dp: Default degree of polymerization if not stored per-edge.
        dims: Box dimensions [x, y, z] (written as graph attributes).

    Returns:
        Path to the written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Build edge → chain-node mapping
    # ------------------------------------------------------------------
    edges = list(G.edges(keys=True, data=True))
    edge_to_id = {}          # (u, v, key) → chain_node_id (1-based)
    chain_dp = {}             # chain_node_id → DP

    for idx, (u, v, key, data) in enumerate(edges, start=1):
        edge_key = (u, v, key)
        edge_to_id[edge_key] = idx
        chain_dp[idx] = data.get("dp", dp)

    # ------------------------------------------------------------------
    # 2. Build crosslinker nodes and chemical edges
    #    For every edge in the Topon graph (a chain), it connects two
    #    crosslinks (u and v). In GraphML, the chain is a node (cid),
    #    and the crosslinks (u, v) are ALSO nodes.
    #    We add chemical edges: (cid <-> u) and (cid <-> v).
    # ------------------------------------------------------------------
    crosslink_nodes = set()
    chemical_edges = set() # (src, tgt) where one is cid, one is crosslink id
    
    # We must offset chain IDs from crosslink IDs so they don't collide.
    # Let's say chain IDs start after max crosslink ID.
    max_xlink_id = max(max(u, v) for u, v, _ in G.edges) if G.edges else 0
    
    # Re-map edge_to_id so chain IDs don't collide with crosslinks.
    edge_to_id = {}
    chain_dp = {}
    for idx, (u, v, key, data) in enumerate(edges, start=max_xlink_id + 1):
        edge_key = (u, v, key)
        edge_to_id[edge_key] = idx
        chain_dp[idx] = data.get("dp", dp)
        
        # Track crosslinks
        crosslink_nodes.add(u)
        crosslink_nodes.add(v)
        
        # Chemical edges connect the chain to its two endpoints
        chemical_edges.add((idx, u))
        chemical_edges.add((idx, v))

    # ------------------------------------------------------------------
    # 3. Build entanglement edges  (with multiplicity)
    #    Entanglements connect a chain (cid1) to a chain (cid2).
    # ------------------------------------------------------------------
    entanglement_edges = []   # list of (src, tgt)  — may repeat
    seen_pairs = set()        # avoid double-counting symmetric marking
    for (u, v, key, data) in edges:
        ew = data.get("entangled_with")
        if not ew:
            continue
        ew = tuple(ew)
            
        e1 = (u, v, key)
        e2 = ew
        
        pair = frozenset([frozenset(e1[:2]), frozenset(e2[:2])])
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)

        cid1 = edge_to_id.get(e1) or edge_to_id.get((e1[1], e1[0], e1[2]))
        cid2 = edge_to_id.get(e2) or edge_to_id.get((e2[1], e2[0], e2[2]))
        
        if cid1 is None or cid2 is None:
            continue

        count = data.get("entanglement_count", 1)
        for _ in range(count):
            entanglement_edges.append((cid1, cid2))

    # ------------------------------------------------------------------
    # 4. Build XML
    # ------------------------------------------------------------------
    root = Element("graphml")
    root.set("xmlns", "http://graphml.graphdrawing.org/xmlns")
    root.set("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
    root.set(
        "xsi:schemaLocation",
        "http://graphml.graphdrawing.org/xmlns "
        "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd",
    )

    _add_key(root, "d13", "edge", "edge_type", "string")
    _add_key(root, "d12", "node", "COMZ", "double")
    _add_key(root, "d11", "node", "COMY", "double")
    _add_key(root, "d10", "node", "COMX", "double")
    _add_key(root, "d9",  "node", "rg", "double")
    _add_key(root, "d8",  "node", "contour_length", "double")
    _add_key(root, "d7",  "node", "length", "long")
    _add_key(root, "d6",  "node", "type", "string")
    _add_key(root, "d5",  "graph", "zhi", "float")
    _add_key(root, "d4",  "graph", "zlo", "float")
    _add_key(root, "d3",  "graph", "yhi", "float")
    _add_key(root, "d2",  "graph", "ylo", "float")
    _add_key(root, "d1",  "graph", "xhi", "float")
    _add_key(root, "d0",  "graph", "xlo", "float")

    graph_el = SubElement(root, "graph", edgedefault="undirected")

    # --- Polymer Nodes (chains) ---
    for cid in sorted(edge_to_id.values()):
        node_el = SubElement(graph_el, "node", id=str(cid))
        _add_data(node_el, "d6", "polymer")
        _add_data(node_el, "d7", str(chain_dp[cid]))
        _add_data(node_el, "d8", "NaN")   # contour_length
        _add_data(node_el, "d9", "NaN")   # rg
        _add_data(node_el, "d10", "NaN")  # COMX
        _add_data(node_el, "d11", "NaN")  # COMY
        _add_data(node_el, "d12", "NaN")  # COMZ
        
    # --- Crosslinker Nodes ---
    for xid in sorted(crosslink_nodes):
        node_el = SubElement(graph_el, "node", id=str(xid))
        _add_data(node_el, "d6", "crosslinker")
        _add_data(node_el, "d7", "1")     # Crosslinkers have length 1 in example
        _add_data(node_el, "d8", "0.0")
        _add_data(node_el, "d9", "0.0")
        
        # If the Topon graph has pos attribute, we can add it, but fallback to NaN
        node_data = G.nodes.get(xid, {})
        pos = node_data.get("pos")
        if pos is not None and len(pos) == 3:
            _add_data(node_el, "d10", str(pos[0]))
            _add_data(node_el, "d11", str(pos[1]))
            _add_data(node_el, "d12", str(pos[2]))
        else:
            _add_data(node_el, "d10", "NaN")
            _add_data(node_el, "d11", "NaN")
            _add_data(node_el, "d12", "NaN")

    # --- Chemical edges (chain <-> crosslink) ---
    for src, tgt in sorted(chemical_edges):
        edge_el = SubElement(graph_el, "edge", source=str(src), target=str(tgt))
        _add_data(edge_el, "d13", "chemical")

    # --- Entanglement edges (chain <-> chain) ---
    for src, tgt in entanglement_edges:
        edge_el = SubElement(graph_el, "edge", source=str(src), target=str(tgt))
        _add_data(edge_el, "d13", "entanglement")

    # --- Graph-level data (box bounds) ---
    _add_data(graph_el, "d0", "NaN")  # xlo
    _add_data(graph_el, "d1", "NaN")  # xhi
    _add_data(graph_el, "d2", "NaN")  # ylo
    _add_data(graph_el, "d3", "NaN")  # yhi
    _add_data(graph_el, "d4", "NaN")  # zlo
    _add_data(graph_el, "d5", "NaN")  # zhi

    # ------------------------------------------------------------------
    # 5. Write to disk
    # ------------------------------------------------------------------
    tree = ElementTree(root)
    _indent(root)  # pretty-print
    with open(output_path, "wb") as f:
        tree.write(f, encoding="utf-8", xml_declaration=True)

    # Summary
    n_chains = len(edge_to_id)
    n_xlink = len(crosslink_nodes)
    n_chem = len(chemical_edges)
    n_ent = len(entanglement_edges)
    print(f"  GraphML written to {output_path}")
    print(f"    Polymer nodes : {n_chains}")
    print(f"    Crosslink nodes : {n_xlink}")
    print(f"    Chemical edges (chain-link) : {n_chem}")
    print(f"    Entanglement edges (chain-chain): {n_ent}")

    return output_path


# ======================================================================
# Helpers
# ======================================================================

def _add_key(parent: Element, kid: str, for_: str, name: str, type_: str):
    """Add a <key> element."""
    SubElement(
        parent, "key",
        id=kid, **{"for": for_},
        **{"attr.name": name, "attr.type": type_},
    )


def _add_data(parent: Element, key: str, value: str):
    """Add a <data> element."""
    d = SubElement(parent, "data", key=key)
    d.text = value


def _indent(elem: Element, level: int = 0):
    """Add indentation to XML for pretty-printing."""
    indent = "\n" + "  " * level
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = indent + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = indent
        for child in elem:
            _indent(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = indent
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = indent
    if level == 0:
        elem.tail = "\n"
