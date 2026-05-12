"""Diagnostics rules for `topon doctor`.

Each rule receives the validated Pydantic `ToponConfig` plus the raw dict
(so we can also inspect schema-gap fields like `conformation`, `simulation`,
`execution`). Rules emit zero or more `Issue` records.

Adding a new rule:
    1. Write a function `check_<name>(cfg, raw) -> list[Issue]`
    2. Append it to RULE_REGISTRY at the bottom of this file
    3. The rule should be self-contained — no chemistry build, no LAMMPS.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, List


@dataclass
class Issue:
    rule: str
    level: str          # "ok" | "warn" | "error"
    message: str
    fix: str | None = None


def _node_type_mapping(cfg) -> dict:
    """Pull the degree->node_type mapping from the active assignment method."""
    assign = cfg.assignment
    if assign.node_types.method == "degree":
        return dict(assign.node_types.degree.mapping or {})
    return {}


# ---------- rules ---------------------------------------------------------

def check_poss_at_internal_junction(cfg, raw) -> List[Issue]:
    """P1-H: POSS at degree-≥2 junctions crashes LAMMPS stage 1."""
    mapping = _node_type_mapping(cfg)
    if not mapping:
        return []
    poss_keys = {
        tname for tname, mol_cfg in (cfg.chemistry.node_type_map or {}).items()
        if mol_cfg.molecule and mol_cfg.molecule.upper().startswith("POSS")
    }
    if not poss_keys:
        return []
    out: List[Issue] = []
    for degree_str, tname in mapping.items():
        if tname in poss_keys:
            try:
                deg = int(degree_str)
            except (TypeError, ValueError):
                continue
            if deg >= 2:
                out.append(Issue(
                    rule="poss_at_internal_junction",
                    level="warn",
                    message=(
                        f"POSS molecule '{tname}' is mapped to degree-{deg} "
                        f"nodes. POSS at internal junctions (degree >= 2) "
                        f"hits known bug P1-H: bond extent > half periodic "
                        f"box at LAMMPS stage 1."
                    ),
                    fix=(
                        f"Map POSS to degree-1 chain caps only "
                        f'(e.g. `"mapping": {{"1": "{tname}"}}`) and use a '
                        f"plain Si node-type for higher-degree junctions."
                    ),
                ))
    return out


def check_unknown_node_type(cfg, raw) -> List[Issue]:
    """P0-2: a `node_type` in assignment.degree.mapping must appear in
    `chemistry.node_type_map` — otherwise the chemistry stage silently
    falls through to a single Si atom, contaminating hydrocarbon polymers.
    """
    mapping = _node_type_mapping(cfg)
    if not mapping:
        return []
    chem_types = set((cfg.chemistry.node_type_map or {}).keys())
    out: List[Issue] = []
    for degree_str, tname in mapping.items():
        if tname not in chem_types:
            out.append(Issue(
                rule="unknown_node_type",
                level="warn",
                message=(
                    f"Assignment maps degree-{degree_str} to '{tname}', but "
                    f"'{tname}' isn't a key in chemistry.node_type_map "
                    f"(present: {sorted(chem_types) or '[]'}). The chemistry "
                    f"stage will silently fall through to a single Si atom."
                ),
                fix=(
                    f'Add a "{tname}" entry under chemistry.node_type_map '
                    f'with a `molecule` SMILES or element symbol.'
                ),
            ))
    return out


def check_atomistic_graft_non_pdms(cfg, raw) -> List[Issue]:
    """P1-L follow-up: atomistic graft is hard-coded to PDMS structure; a
    non-PDMS monomer SMILES with graft_density > 0 silently skips grafts.
    """
    if cfg.chemistry.model_type != "atomistic":
        return []
    grafts_cfg = cfg.assignment.grafts.per_edge_type or {}
    if not grafts_cfg:
        return []
    monomers = cfg.chemistry.monomers or {}
    edge_types = cfg.chemistry.edge_type_map or {}
    out: List[Issue] = []
    for etype, gcfg in grafts_cfg.items():
        if gcfg.graft_density <= 0:
            continue
        edge_chem = edge_types.get(etype)
        if edge_chem is None:
            continue
        mon_cfg = monomers.get(edge_chem.monomer)
        if mon_cfg is None or mon_cfg.smiles == "[Si](C)(C)O":
            continue
        out.append(Issue(
            rule="atomistic_graft_non_pdms",
            level="warn",
            message=(
                f"Edge type '{etype}' has graft_density={gcfg.graft_density} "
                f"but its monomer '{edge_chem.monomer}' has SMILES "
                f"{mon_cfg.smiles!r} (not PDMS '[Si](C)(C)O'). Atomistic "
                f"graft is hard-coded to PDMS structure; non-PDMS monomers "
                f"emit a RuntimeWarning at build time and skip grafts."
            ),
            fix="Use the CG model_type, or change the monomer to PDMS, or "
                "set graft_density=0 for this edge type.",
        ))
    return out


def check_schema_gap_extras(cfg, raw) -> List[Issue]:
    """P0-B (residual): top-level `conformation`/`simulation`/`execution`
    sections aren't Pydantic-validated. CLI path handles them via
    `load_config_full`; direct `Pipeline(ToponConfig(...))` construction
    without `raw_config=...` ignores them silently.
    """
    extras = [k for k in ("conformation", "simulation", "execution") if k in (raw or {})]
    if not extras:
        return []
    return [Issue(
        rule="schema_gap_extras",
        level="ok",
        message=(
            f"Config has unvalidated top-level section(s): {extras}. The "
            f"CLI handles these via `load_config_full`; direct API users "
            f"must pass them as `raw_config={{...}}` to `Pipeline(...)`."
        ),
        fix="If using the API directly, do "
            "`config, raw = load_config_full(path); Pipeline(config, raw_config=raw)`.",
    )]


def check_lattice_size_format(cfg, raw) -> List[Issue]:
    """Catch the common 'lattice_size: 5' (int) vs '5x5x5' (string) confusion."""
    if cfg.topology.source != "generate":
        return []
    gen = cfg.topology.generator
    if gen is None:
        return []
    ls = gen.lattice_size
    if not isinstance(ls, str) or "x" not in ls.lower():
        return [Issue(
            rule="lattice_size_format",
            level="error",
            message=(
                f"topology.generator.lattice_size = {ls!r}; expected a "
                f"string like '5x5x5' or '6x6x6'."
            ),
            fix='Set `"lattice_size": "5x5x5"` (string with two `x` separators).',
        )]
    return []


def check_dp_below_kuhn(cfg, raw) -> List[Issue]:
    """DP < 5 produces chains shorter than a Kuhn length — geometrically valid
    but the entanglement/conformation stages have edge cases."""
    dp_cfg = cfg.assignment.dp_distribution
    if dp_cfg is None or dp_cfg.default is None:
        return []
    mean = dp_cfg.default.mean
    if mean is None or mean >= 5:
        return []
    return [Issue(
        rule="dp_below_kuhn",
        level="warn",
        message=(
            f"DP mean = {mean} is shorter than a typical Kuhn length (~5). "
            f"Conformation overlap-resolution and entanglement detection "
            f"may produce poor geometries at this length."
        ),
        fix="Use DP >= 5 unless intentionally testing the short-chain limit.",
    )]


def check_defects_without_endcap_safety(cfg, raw) -> List[Issue]:
    """P1-K reminder: primary loops were over-valencing end-cap nodes
    pre-2026-05-10. The fix is in `assignment/defects.py` and is always
    applied; this rule just informs the user that defect+endcap is now safe.
    """
    if not cfg.assignment.defects.primary_loops.enabled:
        return []
    return [Issue(
        rule="defects_endcap_safe",
        level="ok",
        message=(
            "Primary-loop defects are enabled. Eligible candidates skip "
            "degree-1 chain caps (max_degree=4 + exclude_node_types=('end',)), "
            "so the chemistry build stays chemically valid."
        ),
    )]


# ---------- registry + runner ---------------------------------------------

RuleFn = Callable[[object, dict], Iterable[Issue]]

RULE_REGISTRY: list[RuleFn] = [
    check_lattice_size_format,
    check_unknown_node_type,
    check_poss_at_internal_junction,
    check_atomistic_graft_non_pdms,
    check_dp_below_kuhn,
    check_defects_without_endcap_safety,
    check_schema_gap_extras,
]


def run_all_rules(cfg, raw: dict | None = None) -> List[Issue]:
    """Run every rule against ``cfg`` (Pydantic ``ToponConfig``)."""
    raw = raw or {}
    issues: List[Issue] = []
    for rule in RULE_REGISTRY:
        try:
            issues.extend(rule(cfg, raw) or [])
        except Exception as exc:  # rules must never crash doctor
            issues.append(Issue(
                rule=rule.__name__,
                level="warn",
                message=f"Rule crashed: {type(exc).__name__}: {exc}",
            ))
    return issues
