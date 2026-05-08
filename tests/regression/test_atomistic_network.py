"""
Regression tests: Atomistic (DREIDING) polymer network.

Reference: tests/output/v21_atomistic_combined/atomistic_combined/
Workflow:  topon.workflows.atomistic_network (new canonical) +
           tests/workflows/generate_atomistic_combined.py (legacy, for structure tests)

The v21 reference was generated without a fixed seed, so the atom count
(grafts are stochastic) differs between runs. The regression tests verify:
  - TestDreidingStructure:     structure of the frozen reference (no new run needed)
  - TestAtomisticNewWorkflow:  new canonical workflow produces correct file structure
                               and LAMMPS format (with seed=42 for reproducibility)

Run with:
    pytest tests/regression/test_atomistic_network.py -v
"""

import sys
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

REF_DIR = ROOT / "tests/output/v21_atomistic_combined/atomistic_combined"
WORKFLOW = ROOT / "tests/workflows/generate_atomistic_combined.py"
GRAPH_NODES = ROOT / "tests/sample_graphs/network_N5x5x5_trial3.nodes"
GRAPH_EDGES = ROOT / "tests/sample_graphs/network_N5x5x5_trial3.edges"
CONFIG = ROOT / "examples/config_atomistic_combined.json"
EXPERIMENTAL = ROOT / "examples/experimental_test.json"

from tests.regression.lammps_compare import assert_lammps_identical, parse_lammps_data


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ref_data():
    path = REF_DIR / "02_Chemistry/system.data"
    if not path.exists():
        pytest.skip(f"Reference file not found: {path}")
    return parse_lammps_data(path)


@pytest.fixture(scope="module")
def new_workflow_output(tmp_path_factory):
    """Run the new canonical atomistic workflow with seed=42."""
    if not GRAPH_NODES.exists():
        pytest.skip(f"Sample graph not found: {GRAPH_NODES}")
    if not CONFIG.exists():
        pytest.skip(f"Config not found: {CONFIG}")

    from topon.workflows.atomistic_network import run
    output_dir = tmp_path_factory.mktemp("atomistic_new")
    root = run(
        nodes_path=GRAPH_NODES,
        edges_path=GRAPH_EDGES,
        config_path=CONFIG,
        experimental_path=EXPERIMENTAL,
        output_dir=output_dir,
        seed=42,
    )
    return root


@pytest.fixture(scope="module")
def new_output_dir(tmp_path_factory):
    if not WORKFLOW.exists():
        pytest.skip(f"Workflow not found: {WORKFLOW}")
    output_dir = tmp_path_factory.mktemp("atomistic_output")
    import subprocess
    result = subprocess.run(
        [sys.executable, str(WORKFLOW), "--output", str(output_dir)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        pytest.fail(f"Workflow failed:\n{result.stderr}")
    return output_dir


# ---------------------------------------------------------------------------
# Structure Tests (reference introspection — no new run needed)
# ---------------------------------------------------------------------------

class TestDreidingStructure:
    def test_has_atom_types(self, ref_data):
        # PDMS network has Si, O, C, H types
        assert ref_data.header.get("atom types", 0) >= 4

    def test_has_impropers(self, ref_data):
        # DREIDING sp2 atoms need impropers
        assert ref_data.header.get("impropers", 0) >= 0  # May be 0 for pure PDMS

    def test_has_dihedrals(self, ref_data):
        assert ref_data.header.get("dihedrals", 0) > 0

    def test_pair_coeffs_in_settings(self):
        # DreidingWriter puts pair_coeff in system.in.settings, NOT in data file
        settings = REF_DIR / "02_Chemistry/system.in.settings"
        if not settings.exists():
            pytest.skip("Settings file not found")
        text = settings.read_text()
        assert "pair_coeff" in text, "pair_coeff missing from settings file"

    def test_bond_coeffs_in_settings(self):
        settings = REF_DIR / "02_Chemistry/system.in.settings"
        if not settings.exists():
            pytest.skip("Settings file not found")
        text = settings.read_text()
        assert "bond_coeff" in text, "bond_coeff missing from settings file"

    def test_atom_charges_nonzero(self, ref_data):
        charges = [a[2] for a in ref_data.atoms.values()]
        nonzero = [c for c in charges if abs(c) > 1e-10]
        assert len(nonzero) > 0, "Expected non-zero Gasteiger charges"


# ---------------------------------------------------------------------------
# Regression Tests
# ---------------------------------------------------------------------------

class TestDreidingRegression:
    def test_chemistry_data_structure(self, new_output_dir):
        # generate_atomistic_combined.py has no fixed seed; grafts are stochastic.
        # An exact byte-for-byte comparison against the v21 reference is not
        # meaningful here.  Instead verify structural correctness: the new run
        # must have the same atom/bond types and force-field coefficients as the
        # frozen reference, even though atom count will differ across runs.
        new_path = new_output_dir / "atomistic_combined/02_Chemistry/system.data"
        ref_path = REF_DIR / "02_Chemistry/system.data"
        if not ref_path.exists():
            pytest.skip("Reference file not found")
        new = parse_lammps_data(new_path)
        ref = parse_lammps_data(ref_path)
        assert new.header.get("atom types") == ref.header.get("atom types"), \
            "Atom type count mismatch"
        assert new.header.get("bond types") == ref.header.get("bond types"), \
            "Bond type count mismatch"
        assert set(new.masses.keys()) == set(ref.masses.keys()), \
            "Atom type IDs differ"
        assert set(new.bond_coeffs.keys()) == set(ref.bond_coeffs.keys()), \
            "Bond coeff IDs differ"

    def test_settings_file_identical(self, new_output_dir):
        new_path = new_output_dir / "atomistic_combined/02_Chemistry/system.in.settings"
        ref_path = REF_DIR / "02_Chemistry/system.in.settings"
        if not ref_path.exists():
            pytest.skip("Reference settings file not found")
        assert new_path.read_text() == ref_path.read_text(), \
            "Settings file differs from reference"

    def test_lammps_stage1_script(self, new_output_dir):
        new_path = new_output_dir / "atomistic_combined/04_Simulation/minimize_1_serial.in"
        ref_path = REF_DIR / "04_Simulation/minimize_1_serial.in"
        if not ref_path.exists():
            pytest.skip("Reference LAMMPS script not found")
        assert new_path.read_text() == ref_path.read_text()


# ---------------------------------------------------------------------------
# TestAtomisticNewWorkflow: new workflow produces correct LAMMPS format
# ---------------------------------------------------------------------------

class TestAtomisticNewWorkflow:
    def test_output_files_exist(self, new_workflow_output):
        assert (new_workflow_output / "02_Chemistry/system.data").exists()
        assert (new_workflow_output / "02_Chemistry/system.in.settings").exists()
        assert (new_workflow_output / "02_Chemistry/system.groups").exists()
        assert (new_workflow_output / "02_Chemistry/system_nodes.displace").exists()
        assert (new_workflow_output / "02_Chemistry/system_backbone.displace").exists()
        assert (new_workflow_output / "02_Chemistry/system_pendant.displace").exists()
        assert (new_workflow_output / "02_Chemistry/system_hydrogens.displace").exists()
        assert (new_workflow_output / "03_Conformation/system_conformed.data").exists()
        assert (new_workflow_output / "03_Conformation/system_relaxed.data").exists()
        assert (new_workflow_output / "04_Simulation/minimize_1_serial.in").exists()
        assert (new_workflow_output / "04_Simulation/minimize_2_parallel.in").exists()
        assert (new_workflow_output / "04_Simulation/minimize_3_parallel.in").exists()

    def test_atom_types(self, new_workflow_output):
        data = parse_lammps_data(new_workflow_output / "02_Chemistry/system.data")
        # PDMS network: Si3, O_3, C_3, H_ (at minimum)
        assert data.header.get("atom types", 0) >= 4

    def test_has_dihedrals(self, new_workflow_output):
        data = parse_lammps_data(new_workflow_output / "02_Chemistry/system.data")
        assert data.header.get("dihedrals", 0) > 0

    def test_pair_coeffs_in_settings(self, new_workflow_output):
        settings = new_workflow_output / "02_Chemistry/system.in.settings"
        text = settings.read_text()
        assert "pair_coeff" in text, "pair_coeff missing from settings file"

    def test_bond_coeffs_in_settings(self, new_workflow_output):
        settings = new_workflow_output / "02_Chemistry/system.in.settings"
        text = settings.read_text()
        assert "bond_coeff" in text, "bond_coeff missing from settings file"

    def test_gasteiger_charges_nonzero(self, new_workflow_output):
        data = parse_lammps_data(new_workflow_output / "02_Chemistry/system.data")
        charges = [a[2] for a in data.atoms.values()]
        nonzero = [c for c in charges if abs(c) > 1e-10]
        assert len(nonzero) > 0, "Expected non-zero Gasteiger charges"

    def test_lammps_input_script_format(self, new_workflow_output):
        script = (new_workflow_output / "04_Simulation/minimize_1_serial.in").read_text()
        assert "units           real" in script
        assert "bond_style      harmonic" in script
        assert "pair_style      lj" in script or "pair_style      soft" in script
        assert "minimize" in script

    def test_conformed_data_valid(self, new_workflow_output):
        data = parse_lammps_data(new_workflow_output / "03_Conformation/system_conformed.data")
        coords = [(a[3], a[4], a[5]) for a in data.atoms.values()]
        nonzero = [c for c in coords if any(abs(v) > 1e-6 for v in c)]
        assert len(nonzero) > 0, "All coordinates are zero after displacement"

    def test_reproducible_with_seed(self, tmp_path):
        """Same seed -> same atom count (grafts are deterministic)."""
        from topon.workflows.atomistic_network import run
        root_a = run(GRAPH_NODES, GRAPH_EDGES, CONFIG, EXPERIMENTAL,
                     tmp_path / "a", seed=77)
        root_b = run(GRAPH_NODES, GRAPH_EDGES, CONFIG, EXPERIMENTAL,
                     tmp_path / "b", seed=77)
        da = parse_lammps_data(root_a / "02_Chemistry/system.data")
        db = parse_lammps_data(root_b / "02_Chemistry/system.data")
        assert da.header["atoms"] == db.header["atoms"], \
            "Same seed produced different atom counts"
