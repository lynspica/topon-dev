"""Unit tests for topon.diagnostics rules.

Tests fire each rule against a synthetic config that's tailored to trigger
(or not trigger) the rule. No chemistry build, no LAMMPS — these are
millisecond tests.
"""
from __future__ import annotations

import pytest

from topon.config.schema import (
    AssignmentConfig,
    ChemistryConfig,
    DegreeNodeTypeConfig,
    DPConfig,
    DPDistributionConfig,
    DefectsConfig,
    GeneratorConfig,
    NodeMoleculeConfig,
    NodeTypesConfig,
    OutputConfig,
    StudyConfig,
    TargetConfig,
    ToponConfig,
    TopologyConfig,
)
from topon.diagnostics import run_all_rules


def _make_cfg(**overrides) -> ToponConfig:
    """Build a minimal-but-valid ToponConfig; override fields piecewise."""
    base = dict(
        study=StudyConfig(name="t", output_dir="out"),
        topology=TopologyConfig(
            source="generate",
            generator=GeneratorConfig(
                lattice_size="5x5x5",
                lattice_type="SC",
                max_functionality=4,
                degree_distribution="0:0,1:0",
            ),
        ),
        chemistry=ChemistryConfig(
            model_type="atomistic",
            target_density=1.0,
            node_type_map={"A": NodeMoleculeConfig(molecule="Si")},
        ),
        assignment=AssignmentConfig(
            node_types=NodeTypesConfig(
                method="degree",
                degree=DegreeNodeTypeConfig(mapping={"1": "A", "4": "A"}),
            ),
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=10.0)),
        ),
        output=OutputConfig(),
    )
    base.update(overrides)
    return ToponConfig(**base)


def test_clean_config_emits_no_warns_or_errors():
    issues = run_all_rules(_make_cfg(), raw={})
    assert all(i.level == "ok" for i in issues), [
        f"{i.rule}: {i.level}" for i in issues if i.level != "ok"
    ]


def test_poss_at_degree_4_warns():
    cfg = _make_cfg(
        chemistry=ChemistryConfig(
            model_type="atomistic",
            target_density=1.0,
            node_type_map={
                "A": NodeMoleculeConfig(molecule="Si"),
                "POSS": NodeMoleculeConfig(molecule="POSS_AM0270", is_end_cap=True),
            },
        ),
        assignment=AssignmentConfig(
            node_types=NodeTypesConfig(
                method="degree",
                degree=DegreeNodeTypeConfig(mapping={"1": "A", "4": "POSS"}),
            ),
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=10.0)),
        ),
    )
    issues = run_all_rules(cfg, raw={})
    poss_warns = [i for i in issues if i.rule == "poss_at_internal_junction"]
    assert poss_warns, "Expected a poss_at_internal_junction warn"
    assert poss_warns[0].level == "warn"


def test_unknown_node_type_warns():
    cfg = _make_cfg(
        assignment=AssignmentConfig(
            node_types=NodeTypesConfig(
                method="degree",
                degree=DegreeNodeTypeConfig(mapping={"1": "A", "4": "MISSING"}),
            ),
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=10.0)),
        ),
    )
    issues = run_all_rules(cfg, raw={})
    typed = [i for i in issues if i.rule == "unknown_node_type"]
    assert typed and typed[0].level == "warn"
    assert "MISSING" in typed[0].message


def test_dp_below_kuhn_warns():
    cfg = _make_cfg(
        assignment=AssignmentConfig(
            node_types=NodeTypesConfig(
                method="degree",
                degree=DegreeNodeTypeConfig(mapping={"1": "A", "4": "A"}),
            ),
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=3.0)),
        ),
    )
    issues = run_all_rules(cfg, raw={})
    dp_warns = [i for i in issues if i.rule == "dp_below_kuhn"]
    assert dp_warns and dp_warns[0].level == "warn"


def test_defects_endcap_safe_ok():
    cfg = _make_cfg(
        assignment=AssignmentConfig(
            node_types=NodeTypesConfig(
                method="degree",
                degree=DegreeNodeTypeConfig(mapping={"1": "A", "4": "A"}),
            ),
            dp_distribution=DPDistributionConfig(default=DPConfig(mean=10.0)),
            defects=DefectsConfig(primary_loops=TargetConfig(enabled=True, target=5)),
        ),
    )
    issues = run_all_rules(cfg, raw={})
    ds = [i for i in issues if i.rule == "defects_endcap_safe"]
    assert ds and ds[0].level == "ok"


def test_schema_gap_extras_reports_only_when_present():
    cfg = _make_cfg()
    no_extras = run_all_rules(cfg, raw={})
    assert not [i for i in no_extras if i.rule == "schema_gap_extras"]

    with_extras = run_all_rules(cfg, raw={"conformation": {}, "simulation": {}})
    extras = [i for i in with_extras if i.rule == "schema_gap_extras"]
    assert extras and extras[0].level == "ok"
    assert "conformation" in extras[0].message
