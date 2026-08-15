"""Backbone paths for entangled edges, by either construction.

One function serves the pipeline and both canonical workflow modules, which
until now each carried their own copy of the kink loop:

    paths = entangled_backbone_paths(graph, dims, edge_atoms,
                                     method="waypoint", kink_params={...})
    # paths[(u, v, key)] -> [xyz, ...], one per chain atom, junctions excluded

Two constructions:

``waypoint`` (the default) is the prescribed winding of
:mod:`topon.conformation.entanglement.waypoints`: the two chains of a pair
are drawn together as splines that spiral about their contact in antiphase,
so the pair carries exactly ``entanglement_count`` windings by construction.
Verified with primitive-path analysis (V49): the requested count is
delivered through the full protocol, and nothing appears on pairs that were
not asked.

``kink`` is the legacy Gaussian bump aimed at the partner's midpoint
(:func:`topon.utils.network_helpers.calculate_entangled_kink`). Each chain
is drawn alone, so what the pair carries after relaxation is statistical
rather than prescribed. Kept for comparison with pre-V49 systems.
"""
from __future__ import annotations

import numpy as np

from topon.conformation.entanglement.waypoints import (
    Site,
    entangled_pair,
    resample_path,
)
from topon.utils.network_helpers import calculate_entangled_kink


def _mic(vec, dims):
    if dims is None:
        return np.asarray(vec, float)
    d = np.asarray(dims, float)
    v = np.asarray(vec, float)
    return v - d * np.round(v / d)


def _perp_of(mic):
    """A deterministic vector perpendicular to the chord, for the rare case
    of a partner whose midpoint coincides with the edge's own."""
    axis = np.zeros(3)
    axis[int(np.argmin(np.abs(mic)))] = 1.0
    p = np.cross(mic, axis)
    n = np.linalg.norm(p)
    return p / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])


def _kink_path(pos_u, mic, n_atoms, orient_vec, count, kink_params):
    """The legacy per-edge kink, exactly as the pipeline drew it."""
    kink_dict = calculate_entangled_kink(
        start_pos=np.zeros(3),
        end_pos=mic,
        num_atoms=n_atoms + 2,          # N+2 fix (v21.1)
        params=kink_params,
        orientation_vec=orient_vec,
        z_phase=1.0,
        num_entanglements=count,
    )
    full = [kink_dict[k] for k in sorted(kink_dict.keys())]
    return [pos_u + np.array(pt) for pt in full[1:-1]]


def entangled_backbone_paths(graph, dims, edge_atoms, method="waypoint",
                             kink_params=None):
    """Interior bead positions for every edge carrying ``entangled_with``.

    ``edge_atoms`` maps ``(u, v, key)`` to that edge's chain atoms (only the
    length is used). Returns ``{edge_key: [xyz, ...]}`` with one position per
    atom, junction ends excluded, for exactly the edges that are entangled;
    the caller places every other chain however it already does.

    The partner's chord is taken in the image nearest the edge's own
    midpoint, which is the same convention the kink always used and what
    makes a pair across the periodic boundary wind rather than reach across
    the box.
    """
    kink_params = kink_params or {}
    out = {}
    seen_pairs = set()

    multi = graph.is_multigraph()
    for edge_key, atoms in edge_atoms.items():
        u, v, key = edge_key
        data = graph[u][v][key] if multi else graph[u][v]
        partner = data.get("entangled_with")
        if partner is None or edge_key in out:
            continue

        pos_u = np.asarray(graph.nodes[u].get("pos", (0.0,) * 3), float)
        pos_v = np.asarray(graph.nodes[v].get("pos", (0.0,) * 3), float)
        mic = _mic(pos_v - pos_u, dims)
        count = int(data.get("entanglement_count", 1))

        p_u, p_v = partner[0], partner[1]
        p_key = partner[2] if len(partner) > 2 else 0
        p_pos_u = np.asarray(graph.nodes[p_u].get("pos", (0.0,) * 3), float)
        p_pos_v = np.asarray(graph.nodes[p_v].get("pos", (0.0,) * 3), float)
        p_mic = _mic(p_pos_v - p_pos_u, dims)

        my_mid = pos_u + 0.5 * mic
        delta = _mic((p_pos_u + 0.5 * p_mic) - my_mid, dims)
        b0 = my_mid + delta - 0.5 * p_mic       # partner chord, my image
        b1 = b0 + p_mic

        if method == "kink":
            orient = (delta if np.linalg.norm(delta) >= 0.01
                      else _perp_of(mic))
            out[edge_key] = _kink_path(pos_u, mic, len(atoms), orient,
                                       count, kink_params)
            continue

        # Waypoint: the pair is drawn together, once. Both edges get their
        # paths here; the partner's entry is filled so the caller's loop
        # finds it whichever edge it reaches first.
        pair_id = frozenset([edge_key, (p_u, p_v, p_key)])
        if pair_id in seen_pairs:
            continue
        seen_pairs.add(pair_id)

        p_atoms = edge_atoms.get((p_u, p_v, p_key))
        if p_atoms is None:
            # Partner edge not built (filtered dangling end): nothing to
            # wind around, fall back to the single-chain kink.
            orient = (delta if np.linalg.norm(delta) >= 0.01
                      else _perp_of(mic))
            out[edge_key] = _kink_path(pos_u, mic, len(atoms), orient,
                                       count, kink_params)
            continue

        na, nb = len(atoms) + 2, len(p_atoms) + 2
        pa, pb, _info = entangled_pair(
            pos_u, pos_u + mic, b0, b1,
            [Site(at=0.5, turns=count)],
            n_beads=max(na, nb))
        # Each chain at its own bead count: the pair is drawn at one density
        # and re-placed by arc length, so unequal DP costs nothing.
        pa = resample_path(pa, na)
        pb = resample_path(pb, nb)
        out[edge_key] = [np.array(q) for q in pa[1:-1]]
        out[(p_u, p_v, p_key)] = [np.array(q) for q in pb[1:-1]]

    return out
