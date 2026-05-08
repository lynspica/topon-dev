"""
Regression tests: CG (Kremer-Grest) polymer network.

Reference: tests/output/v21_cg_combined/cg_combined/
Workflow:  topon.workflows.cg_network (new canonical) +
           tests/workflows/generate_cg_combined.py (legacy, for structure tests)

The v21 reference was generated without a fixed seed, so the atom count
(grafts are stochastic) differs between runs. The regression tests verify:
  - TestCGHeader:     structure of the frozen reference (no new run needed)
  - TestCGNewWorkflow: new canonical workflow produces correct file structure
                       and LAMMPS format (with seed=42 for reproducibility)

Run with:
    pytest tests/regression/test_cg_network.py -v
"""

import sys
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

REF_DIR = ROOT / "tests/output/v21_cg_combined/cg_combined"
GRAPH_NODES = ROOT / "tests/sample_graphs/network_N5x5x5_trial3.nodes"
GRAPH_EDGES = ROOT / "tests/sample_graphs/network_N5x5x5_trial3.edges"
CONFIG = ROOT / "examples/config_cg_combined.json"
EXPERIMENTAL = ROOT / "examples/experimental_test.json"

from tests.regression.lammps_compare import assert_lammps_identical, parse_lammps_data


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ref_chemistry_data():
    """Parsed reference data file from Stage 2 (Chemistry)."""
    path = REF_DIR / "02_Chemistry/system.data"
    if not path.exists():
        pytest.skip(f"Reference file not found: {path}")
    return parse_lammps_data(path)


@pytest.fixture(scope="module")
def new_workflow_output(tmp_path_factory):
    """Run the new canonical CG workflow with seed=42."""
    if not GRAPH_NODES.exists():
        pytest.skip(f"Sample graph not found: {GRAPH_NODES}")
    if not CONFIG.exists():
        pytest.skip(f"Config not found: {CONFIG}")

    from topon.workflows.cg_network import run
    output_dir = tmp_path_factory.mktemp("cg_new")
    root = run(
        nodes_path=GRAPH_NODES,
        edges_path=GRAPH_EDGES,
        config_path=CONFIG,
        experimental_path=EXPERIMENTAL,
        output_dir=output_dir,
        seed=42,
    )
    return root


# ---------------------------------------------------------------------------
# TestCGHeader: reference file introspection (no new run)
# ---------------------------------------------------------------------------

class TestCGHeader:
    def test_atom_count(self, ref_chemistry_data):
        assert ref_chemistry_data.header.get("atoms", 0) > 0

    def test_bond_count(self, ref_chemistry_data):
        assert "bonds" in ref_chemistry_data.header

    def test_atom_types(self, ref_chemistry_data):
        # CG combined has 3 bead types: A (chain), G (graft), J (junction)
        assert ref_chemistry_data.header.get("atom types") == 3

    def test_bond_types(self, ref_chemistry_data):
        assert ref_chemistry_data.header.get("bond types") == 1

    def test_bead_type_names(self, ref_chemistry_data):
        names = {v[1] for v in ref_chemistry_data.masses.values()}
        assert names == {"A", "G", "J"}


# ---------------------------------------------------------------------------
# TestCGNewWorkflow: new workflow produces correct LAMMPS format
# ---------------------------------------------------------------------------

class TestCGNewWorkflow:
    def test_output_files_exist(self, new_workflow_output):
        assert (new_workflow_output / "02_Chemistry/system.data").exists()
        assert (new_workflow_output / "02_Chemistry/system.groups").exists()
        assert (new_workflow_output / "02_Chemistry/system_nodes.displace").exists()
        assert (new_workflow_output / "02_Chemistry/system_beads.displace").exists()
        assert (new_workflow_output / "03_Conformation/system_conformed.data").exists()
        assert (new_workflow_output / "03_Conformation/system_relaxed.data").exists()
        assert (new_workflow_output / "04_Simulation/minimize_1_serial.in").exists()
        assert (new_workflow_output / "04_Simulation/minimize_2_parallel.in").exists()
        assert (new_workflow_output / "04_Simulation/minimize_3_parallel.in").exists()

    def test_bead_types(self, new_workflow_output):
        data = parse_lammps_data(new_workflow_output / "02_Chemistry/system.data")
        assert data.header.get("atom types") == 3
        names = {v[1] for v in data.masses.values()}
        assert names == {"A", "G", "J"}

    def test_bond_types(self, new_workflow_output):
        data = parse_lammps_data(new_workflow_output / "02_Chemistry/system.data")
        assert data.header.get("bond types") == 1
        assert data.header.get("angle types") == 1

    def test_lammps_input_script_format(self, new_workflow_output):
        script = (new_workflow_output / "04_Simulation/minimize_1_serial.in").read_text()
        assert "units           lj" in script
        assert "bond_style      harmonic" in script
        assert "pair_style      soft" in script
        assert "minimize" in script

    def test_conformed_data_valid(self, new_workflow_output):
        data = parse_lammps_data(new_workflow_output / "03_Conformation/system_conformed.data")
        # Coordinates should be set (not all zero after displacement)
        coords = [(a[3], a[4], a[5]) for a in data.atoms.values()]
        nonzero = [c for c in coords if any(abs(v) > 1e-6 for v in c)]
        assert len(nonzero) > 0, "All coordinates are zero after displacement"

    def test_reproducible_with_seed(self, tmp_path):
        """Same seed → same atom count (grafts are deterministic)."""
        from topon.workflows.cg_network import run
        root_a = run(GRAPH_NODES, GRAPH_EDGES, CONFIG, EXPERIMENTAL,
                     tmp_path / "a", seed=99)
        root_b = run(GRAPH_NODES, GRAPH_EDGES, CONFIG, EXPERIMENTAL,
                     tmp_path / "b", seed=99)
        da = parse_lammps_data(root_a / "02_Chemistry/system.data")
        db = parse_lammps_data(root_b / "02_Chemistry/system.data")
        assert da.header["atoms"] == db.header["atoms"], \
            "Same seed produced different atom counts"

    def test_ref_script_matches_new_workflow(self, new_workflow_output):
        """LAMMPS script produced by new workflow matches reference."""
        new_script = new_workflow_output / "04_Simulation/minimize_1_serial.in"
        ref_script = REF_DIR / "04_Simulation/minimize_1_serial.in"
        if not ref_script.exists():
            pytest.skip("Reference LAMMPS script not found")
        assert new_script.read_text() == ref_script.read_text()
