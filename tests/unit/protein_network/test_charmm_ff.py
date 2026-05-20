"""Unit tests for the CHARMM36m force-field lookups + terminal patching.

Guards the 2026-05-19 fixes to `topon/protein_network/charmm/`:
  1. lookup_dihedral wildcard fallback   (proline-ring + sidechain torsions)
  2. lookup_improper bidirectional match (ARG guanidinium, carboxylate)
  3. PROP patch for N-terminal proline   (correct CP1/NP typing + integer charge)

Before these fixes the writer injected K=0 / generic "(DEFAULT)" parameters
for any term the lookup missed, silently zeroing real CHARMM torsions and
producing non-integer protein charge for proline-initial chains.
"""
from __future__ import annotations

from importlib import resources

import pytest

from topon.protein_network import sequence
from topon.protein_network.bfm import generate_topology
from topon.protein_network.charmm import Atom, build_protein_system, build_type_maps
from topon.protein_network.charmm.charmm_ff import CHARMMForceField


@pytest.fixture(scope="module")
def ff() -> CHARMMForceField:
    with resources.as_file(resources.files("topon.protein_network.charmm.data")) as d:
        return CHARMMForceField(
            str(d / "par_all36m_prot_C2L.prm"),
            str(d / "top_all36m_prot_C2L.rtf"),
        )


# ── 1. Dihedral wildcard fallback ────────────────────────────────────────────

@pytest.mark.parametrize("quad", [
    ("C", "CP1", "CP2", "CP2"),     # proline ring, covered by X CP1 CP2 X
    ("CP1", "CP2", "CP2", "CP3"),   #               covered by X CP2 CP2 X
    ("CP2", "CP2", "CP1", "N"),     #               covered by X CP1 CP2 X (rev)
    ("HB1", "CP1", "CP2", "CP2"),
    ("HA2", "CP2", "CP2", "CP3"),
])
def test_dihedral_wildcard_resolves_proline_ring(ff, quad):
    """Proline-ring torsions are defined in the .prm only as wildcard
    terms (X _ _ X). The lookup must find them, not fall through to K=0."""
    res = ff.lookup_dihedral(*quad)
    assert res is not None, f"{quad} should resolve via X-wildcard, got None"
    # CHARMM proline ring torsions are small but non-zero (0.14/0.16).
    k = res[0][0]
    assert k > 0.0, f"{quad} resolved but K={k} (expected non-zero barrier)"


def test_dihedral_exact_takes_priority_over_wildcard(ff):
    """A fully specified dihedral must win over any wildcard term."""
    # C-N-CP1-C is explicitly listed (0.8, 3) in the proline PEP block.
    exact = ff.lookup_dihedral("C", "CP1", "N", "C")
    assert exact is not None
    # The explicit value (0.8) differs from a generic wildcard, confirming
    # the exact branch fired.
    assert any(abs(term[0] - 0.8) < 1e-6 for term in exact)


def test_dihedral_reverse_order_matches(ff):
    """Reversed atom order must resolve to the same parameter."""
    fwd = ff.lookup_dihedral("C", "CP1", "CP2", "CP2")
    rev = ff.lookup_dihedral("CP2", "CP2", "CP1", "C")
    assert fwd == rev and fwd is not None


# ── 2. Improper bidirectional matching ───────────────────────────────────────

def test_improper_guanidinium_resolves(ff):
    """ARG guanidinium improper (central C, three NC2) is stored in the
    .prm as `NC2 X X C` (central atom LAST). The bidirectional lookup
    must find it with K=45, not the K=20 default."""
    res = ff.lookup_improper("C", "NC2", "NC2", "NC2")
    assert res is not None, "guanidinium improper should resolve via reversed match"
    assert abs(res[0] - 45.0) < 1e-6, f"expected K=45 (got {res[0]})"


def test_improper_carboxylate_resolves(ff):
    """ASP/GLU/C-term carboxylate improper (central CC) is stored as
    `OC X X CC`. Bidirectional lookup must find it with K=96."""
    res = ff.lookup_improper("CC", "CT2A", "OC", "OC")
    assert res is not None, "carboxylate improper should resolve via reversed match"
    assert abs(res[0] - 96.0) < 1e-6, f"expected K=96 (got {res[0]})"


def test_improper_forward_central_still_works(ff):
    """A normal peptide improper (central atom first) must still resolve."""
    # O C(=O) backbone improper: C  X  X  O style exists in CHARMM.
    res = ff.lookup_improper("C", "CT1", "NH1", "O")
    assert res is not None


# ── 3. N-terminal proline uses PROP (typing + integer charge) ────────────────

def _build(ff, block_seq, n_repeats=4):
    topo = generate_topology(n_chains=2, n_repeats=n_repeats, segs_per_block=2,
                             equil_steps=1000, seed=7, verbose=False)
    snap = topo["snapshots"][0]
    full_seq = sequence.build_full_sequence(block_seq, n_repeats)
    node_to_res = sequence.get_node_residue_mapping(
        n_repeats=n_repeats, segs_per_block=2, block_seq=block_seq,
    )
    atoms, bonds, *_ = build_protein_system(ff, snap, full_seq, node_to_res)
    return atoms


def test_nterminal_proline_typed_cp1_not_ct1(ff):
    """A chain starting with proline must keep its N-terminal CA as CP1
    (proline ring), not the generic NTER CT1. The CA-CT1 mis-typing was
    the source of the CHARMM-nonexistent angles/dihedrals + charge bug."""
    atoms = _build(ff, "PGRPSDSYPAPGPPN")
    # First residue (res_idx 0) CA of chain 0.
    ca0 = [a for a in atoms if a.res_name == "PRO" and a.name == "CA"
           and a.chain_id == 0]
    assert ca0, "expected a proline CA in chain 0"
    # The very first one is the N-terminal proline.
    first_ca = min(ca0, key=lambda a: a.idx)
    assert first_ca.atype == "CP1", (
        f"N-terminal proline CA typed {first_ca.atype!r}, expected 'CP1' "
        f"(NTER mis-types it as CT1)"
    )


def test_proline_initial_chain_has_integer_charge(ff):
    """The whole proline-initial protein must sum to an integer charge.
    NTER-on-proline produced a fractional total (root of the non-neutral
    system bug); PROP restores integrality."""
    atoms = _build(ff, "PGRPSDSYPAPGPPN")
    total = sum(a.charge for a in atoms)
    assert abs(total - round(total)) < 1e-6, (
        f"proline-initial protein net charge {total:+.4f} is not integer"
    )


def test_glycine_initial_chain_unaffected(ff):
    """Sanity: a glycine-initial chain (GLYP path) is still integer and
    its first CA is the generic CT2 (glycine has no sidechain CB)."""
    atoms = _build(ff, "GGQPSDSYGAPGGGN")
    total = sum(a.charge for a in atoms)
    assert abs(total - round(total)) < 1e-6


# ── 4. Multi-term dihedral expansion (one LAMMPS type per Fourier term) ───────

def _atom(i, atype):
    return Atom(i, atype, atype, 0.0, "X", 1, 0, (float(i), 0.0, 0.0), 1)


def test_build_type_maps_expands_multiterm_dihedral(ff):
    """The omega peptide torsion CT1-C-NH1-CT1 has TWO CHARMM Fourier
    terms (n=1 K=1.6 and n=2 K=2.5). build_type_maps must allocate one
    LAMMPS dihedral type per term, not collapse to a single (truncated)
    type that would silently drop the dominant n=2 planarity term."""
    atoms = [_atom(1, "CT1"), _atom(2, "C"), _atom(3, "NH1"), _atom(4, "CT1")]
    bonds = [(1, 2), (2, 3), (3, 4)]
    maps = build_type_maps(atoms, bonds, [], [(1, 2, 3, 4)], [], ff)
    dmap = maps[3]
    tids = dmap[("CT1", "C", "NH1", "CT1")]
    assert len(tids) == 2, f"omega should expand to 2 types, got {tids}"
    # reverse-key alias points at the same list
    assert dmap[("CT1", "NH1", "C", "CT1")] == tids


def test_build_type_maps_singleterm_dihedral_one_type(ff):
    """A single-term dihedral (C-CP1-N-C, proline PEP 0.8/n=3) maps to
    exactly one type — no spurious expansion."""
    atoms = [_atom(1, "C"), _atom(2, "CP1"), _atom(3, "N"), _atom(4, "C")]
    maps = build_type_maps(atoms, [(1, 2), (2, 3), (3, 4)], [],
                           [(1, 2, 3, 4)], [], ff)
    assert len(maps[3][("C", "CP1", "N", "C")]) == 1


def test_build_type_maps_without_ff_collapses_to_single(ff):
    """Legacy guard: with ff=None every quad collapses to one type
    (the old truncating behaviour) — confirms ff is what enables
    expansion, so the default can never silently change."""
    atoms = [_atom(1, "CT1"), _atom(2, "C"), _atom(3, "NH1"), _atom(4, "CT1")]
    maps = build_type_maps(atoms, [(1, 2), (2, 3), (3, 4)], [],
                           [(1, 2, 3, 4)], [], ff=None)
    assert len(maps[3][("CT1", "C", "NH1", "CT1")]) == 1


# ── 5. Histidine remap (no silent residue drop) ──────────────────────────────

def test_histidine_built_as_hsd_not_dropped(ff):
    """A sequence containing H must build (HIS->HSD default tautomer),
    not silently skip the residue (which would fuse its neighbours)."""
    atoms = _build(ff, "GHGPSDSYGAPGGGN")  # H at block position 2
    assert any(a.res_name == "HSD" for a in atoms), (
        "histidine should be instantiated as HSD, not dropped"
    )
    # net charge still integer (HSD is neutral)
    total = sum(a.charge for a in atoms)
    assert abs(total - round(total)) < 1e-6
