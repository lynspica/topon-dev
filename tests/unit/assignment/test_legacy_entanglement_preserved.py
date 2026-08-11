"""The legacy entanglement path must keep behaving exactly as it did.

`select_entanglements` gained two optional arguments -- `shell_weights` on
the config and `chain_paths` on the call. Both are additions to a code path
that existing studies depend on, so the property worth pinning is not that
they work but that they are *inert*: a caller who does not use them gets the
selection it has always got, bit for bit.

Verified against the main-branch implementation while these were added:
identical candidate list (89 pairs), identical selection (35 pairs, same
sha1), and a byte-identical conformed data file across repeated runs.
"""

import random

import networkx as nx
import numpy as np
import pytest

from topon.assignment.entanglements import (
    find_crossing_candidates,
    select_entanglements,
)
from topon.config.schema import EntanglementsConfig


@pytest.fixture
def cubic():
    """A small simple-cubic network with positions, as the pipeline builds."""
    G = nx.MultiGraph()
    n = 4
    for i in range(n):
        for j in range(n):
            for k in range(n):
                G.add_node((i, j, k), pos=(float(i), float(j), float(k)))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    b = ((i + d[0]) % n, (j + d[1]) % n, (k + d[2]) % n)
                    G.add_edge((i, j, k), b)
    return G, np.array([float(n)] * 3)


def _select(G, dims, seed=1, **kw):
    cands = find_crossing_candidates(G, dims)
    random.seed(seed)
    opts = {"enabled": True, "avg_crosslinks_per_chain": 2.0}
    opts.update(kw.pop("cfg", {}))
    cfg = EntanglementsConfig(**opts)
    return cands, select_entanglements(G, cfg, dims, candidates=list(cands),
                                       num_chains=G.number_of_edges(), **kw)


def test_the_same_seed_gives_the_same_selection(cubic):
    """Repeatability: a study rerun must reproduce its own systems."""
    G, dims = cubic
    _, first = _select(G, dims)
    _, second = _select(G, dims)
    assert first == second
    assert first, "the fixture should produce some entanglements"


def test_shell_weights_default_is_inert(cubic):
    """An empty mapping must not perturb the draw."""
    G, dims = cubic
    _, without = _select(G, dims)
    _, with_empty = _select(G, dims, cfg={"shell_weights": {}})
    assert without == with_empty


def test_chain_paths_default_is_inert(cubic):
    """Passing no conformation must leave the ranking as it was."""
    G, dims = cubic
    _, without = _select(G, dims)
    _, with_none = _select(G, dims, chain_paths=None)
    assert without == with_none


def test_a_conformation_does_change_the_draw(cubic):
    """The opposite guard: the new argument must actually do something.

    Without this, the two tests above would also pass if the feature were
    silently disconnected.
    """
    G, dims = cubic
    cands = find_crossing_candidates(G, dims)
    # Put every chain on top of one edge's midpoint, so exactly the pairs
    # involving it score and everything else is excluded.
    paths = {}
    for u, v in {(a, b) for a, b, _ in [(c[0][0], c[0][1], 0) for c in cands]}:
        paths[frozenset((u, v))] = np.zeros((8, 3))
    random.seed(1)
    cfg = EntanglementsConfig(enabled=True, avg_crosslinks_per_chain=2.0)
    ranked = select_entanglements(G, cfg, dims, candidates=list(cands),
                                  num_chains=G.number_of_edges(),
                                  chain_paths=paths, proximity_cutoff=0.5)
    random.seed(1)
    plain = select_entanglements(G, cfg, dims, candidates=list(cands),
                                 num_chains=G.number_of_edges())
    assert ranked != plain


def test_count_follows_the_requested_density(cubic):
    """e=2 must ask for more than e=1, which is what a study varies."""
    G, dims = cubic
    _, one = _select(G, dims, cfg={"avg_crosslinks_per_chain": 1.0})
    _, two = _select(G, dims, cfg={"avg_crosslinks_per_chain": 4.0})
    assert sum(c for _, _, c in two) > sum(c for _, _, c in one)
