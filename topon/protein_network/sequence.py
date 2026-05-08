"""Sequence-to-BFM-node mapping for arbitrary residue blocks.

Forked from `legacy/subprojects/protein_network/topro/topro/protein/sequence.py`,
generalized to any one-letter block by inferring the crosslinker position from
the block string itself (uses the first 'Y' by default, configurable via
`crosslinker_letter`). The 3-letter mapping comes from this package's residues
module rather than being duplicated.

The mapping function returns ``{node_idx: residue_idx}`` so callers can ask
"which residue does this BFM lattice node represent?". Anchors are fixed first
(chain ends + crosslinker positions); intermediate NC nodes are linearly
interpolated between adjacent anchor pairs.
"""
from __future__ import annotations

from .residues import ONE_TO_THREE, REFERENCE_BLOCK


DEFAULT_BLOCK_SEQ: str = REFERENCE_BLOCK
DEFAULT_CROSSLINKER_LETTER: str = "Y"
DEFAULT_Y_IN_BLOCK: int = REFERENCE_BLOCK.index("Y")  # 7 for resilin GGRPSDSYGAPGGGN


def get_node_residue_mapping(
    n_repeats: int,
    segs_per_block: int,
    y_offset_in_block: int = 0,
    block_seq: str | None = None,
    block_size: int | None = None,
    y_in_block: int | None = None,
    crosslinker_letter: str = DEFAULT_CROSSLINKER_LETTER,
) -> dict[int, int]:
    """Map every BFM chain node index to a residue index.

    Anchors:
      * End nodes  -> residues 0 and (n_repeats * block_size - 1)
      * Y nodes    -> the crosslinker residue within each repeat block

    Intermediate (NC) nodes are linearly interpolated between adjacent anchors.
    """
    if block_seq is not None:
        block_size = len(block_seq)
        if crosslinker_letter in block_seq:
            y_in_block = block_seq.index(crosslinker_letter)
        else:
            y_in_block = block_size // 2
    else:
        if block_size is None:
            block_seq = DEFAULT_BLOCK_SEQ
            block_size = len(DEFAULT_BLOCK_SEQ)
        if y_in_block is None:
            y_in_block = DEFAULT_Y_IN_BLOCK

    n_nodes = n_repeats * segs_per_block + 1
    n_residues = n_repeats * block_size

    mapping: dict[int, int] = {}
    mapping[0] = 0
    mapping[n_nodes - 1] = n_residues - 1
    for k in range(n_repeats):
        y_node = k * segs_per_block + y_offset_in_block + 1
        y_res = k * block_size + y_in_block
        mapping[y_node] = y_res

    sorted_nodes = sorted(mapping.keys())
    for i in range(len(sorted_nodes) - 1):
        n_start = sorted_nodes[i]
        n_end = sorted_nodes[i + 1]
        r_start = mapping[n_start]
        r_end = mapping[n_end]
        span = n_end - n_start
        for j in range(1, span):
            nc_node = n_start + j
            if nc_node not in mapping:
                t = j / span
                mapping[nc_node] = round(r_start + t * (r_end - r_start))
    return mapping


def build_full_sequence(block_seq: str, n_repeats: int) -> list[str]:
    """Tile a one-letter repeat block n_repeats times and convert to 3-letter codes."""
    full_1l = block_seq * n_repeats
    out: list[str] = []
    for aa in full_1l:
        code = ONE_TO_THREE.get(aa.upper())
        if code is None:
            raise ValueError(f"Unknown one-letter amino acid code: {aa!r}")
        out.append(code)
    return out


def get_tyr_node_indices(n_repeats: int, segs_per_block: int, y_offset_in_block: int = 0) -> set[int]:
    """Return the set of chain node indices that are crosslinker positions."""
    return {k * segs_per_block + y_offset_in_block + 1 for k in range(n_repeats)}
