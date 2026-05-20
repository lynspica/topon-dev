"""Regression guard: graphml_writer must emit d14=strain + d15=stress
graph-level keys (matching the GNN-pipeline `translate.py` convention),
and must populate the box (d0..d5) from ``dims`` rather than always-NaN.

Without these, the file is "format-incomplete" per the GNN pipeline:
the user's friend's translate.py expects strain/stress as JSON-array
strings on the <graph> element (empty -> "[]"). The bug pre-fix was
that graphml_writer declared only d0..d13 and emitted NaN for box
regardless of dims.
"""
from __future__ import annotations

import json
import xml.etree.ElementTree as ET

import networkx as nx
import numpy as np

from topon.writers.graphml_writer import write_graphml


GRAPHML_NS = "{http://graphml.graphdrawing.org/xmlns}"


def _simple_graph() -> nx.MultiGraph:
    G = nx.MultiGraph()
    G.add_node(0, pos=(0.0, 0.0, 0.0))
    G.add_node(1, pos=(1.0, 0.0, 0.0))
    G.add_node(2, pos=(0.0, 1.0, 0.0))
    G.add_edge(0, 1, dp=10)
    G.add_edge(1, 2, dp=12)
    return G


def _parse_keys(root) -> dict:
    """{attr_name: (id, for, type)}."""
    out = {}
    for k in root.findall(f"{GRAPHML_NS}key"):
        out[k.get("attr.name")] = (k.get("id"), k.get("for"),
                                   k.get("attr.type"))
    return out


def _parse_graph_data(graph_el, keys) -> dict:
    """{attr_name: text} for direct <data> children of <graph>."""
    by_id = {v[0]: name for name, v in keys.items()}
    out = {}
    for d in graph_el.findall(f"{GRAPHML_NS}data"):
        name = by_id.get(d.get("key"))
        if name is not None:
            out[name] = d.text
    return out


def test_graphml_declares_strain_stress_keys(tmp_path):
    """The graphml file must declare d14=strain and d15=stress as
    for='graph' string keys."""
    out = tmp_path / "g.graphml"
    write_graphml(_simple_graph(), str(out),
                  dp=10, dims=np.array([2.0, 2.0, 2.0]))
    root = ET.parse(out).getroot()
    keys = _parse_keys(root)

    assert "strain" in keys, "graph-level strain key missing"
    assert "stress" in keys, "graph-level stress key missing"
    assert keys["strain"] == ("d14", "graph", "string")
    assert keys["stress"] == ("d15", "graph", "string")


def test_graphml_emits_empty_strain_stress_data(tmp_path):
    """Default (no strain/stress passed) must emit value '[]' on the
    <graph> element, matching translate.py's empty convention."""
    out = tmp_path / "g.graphml"
    write_graphml(_simple_graph(), str(out),
                  dp=10, dims=np.array([2.0, 2.0, 2.0]))
    root = ET.parse(out).getroot()
    keys = _parse_keys(root)
    graph_el = root.find(f"{GRAPHML_NS}graph")
    data = _parse_graph_data(graph_el, keys)

    assert data.get("strain") == "[]"
    assert data.get("stress") == "[]"


def test_graphml_emits_populated_strain_stress_data(tmp_path):
    """When arrays are passed, they must serialise as JSON strings."""
    out = tmp_path / "g.graphml"
    strain = np.linspace(0, 1.0, 5)
    stress = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    write_graphml(_simple_graph(), str(out),
                  dp=10, dims=np.array([2.0, 2.0, 2.0]),
                  strain=strain, stress=stress)
    root = ET.parse(out).getroot()
    keys = _parse_keys(root)
    data = _parse_graph_data(root.find(f"{GRAPHML_NS}graph"), keys)

    s_in = json.loads(data["strain"])
    p_in = json.loads(data["stress"])
    np.testing.assert_allclose(s_in, strain, rtol=1e-6)
    np.testing.assert_allclose(p_in, stress, rtol=1e-6)


def test_graphml_box_populated_from_dims(tmp_path):
    """When dims is passed, box bounds (d0..d5) must be 0..Lx, 0..Ly,
    0..Lz -- not 'NaN' (the pre-fix bug)."""
    out = tmp_path / "g.graphml"
    write_graphml(_simple_graph(), str(out),
                  dp=10, dims=np.array([2.5, 3.5, 4.5]))
    root = ET.parse(out).getroot()
    keys = _parse_keys(root)
    data = _parse_graph_data(root.find(f"{GRAPHML_NS}graph"), keys)

    assert data["xlo"] == "0.0"
    assert float(data["xhi"]) == 2.5
    assert data["ylo"] == "0.0"
    assert float(data["yhi"]) == 3.5
    assert data["zlo"] == "0.0"
    assert float(data["zhi"]) == 4.5


def test_graphml_box_nan_when_dims_missing(tmp_path):
    """No dims -> box stays NaN (backwards-compatible)."""
    out = tmp_path / "g.graphml"
    write_graphml(_simple_graph(), str(out), dp=10, dims=None)
    root = ET.parse(out).getroot()
    keys = _parse_keys(root)
    data = _parse_graph_data(root.find(f"{GRAPHML_NS}graph"), keys)

    for k in ("xlo", "xhi", "ylo", "yhi", "zlo", "zhi"):
        assert data[k] == "NaN"
