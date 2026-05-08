# Topon CLI Reference

Topon provides a command-line interface installed as `topon` when the package is installed.

```
topon [--version] [--help] <command> [options]
```

---

## Commands

### `topon generate`

Run the full pipeline from a JSON configuration file.

```bash
topon generate CONFIG_PATH [--output DIR] [--dry-run]
```

| Argument / Option | Description |
|---|---|
| `CONFIG_PATH` | Path to the JSON configuration file (required) |
| `--output`, `-o` | Override the `study.output_dir` from config |
| `--dry-run` | Validate the config and exit without running |

**Examples:**

```bash
# Run full pipeline
topon generate examples/config_cg.json

# Override output directory
topon generate examples/config_cg.json --output ./my_run

# Validate only
topon generate examples/config_cg.json --dry-run
```

The pipeline runs 6 stages in sequence: Topology → Analysis → Assignment → Chemistry → Conformation → Output. LAMMPS data files and input scripts are written to `output_dir/study_name/`.

---

### `topon validate`

Validate a configuration file without running the pipeline.

```bash
topon validate CONFIG_PATH
```

Prints `Configuration is valid!` or lists all validation errors. Useful for catching schema errors before submitting to HPC.

---

### `topon init`

Create a new configuration file populated with defaults.

```bash
topon init [--output FILE] [--full]
```

| Option | Default | Description |
|---|---|---|
| `--output`, `-o` | `config.json` | Output file path |
| `--full` | off | Include all options with default values |

**Example:**

```bash
topon init --output my_config.json
```

Edit the generated file to set your topology source, DP values, entanglements, etc.

---

### `topon simbox`

Pack a crosslink simulation box and write LAMMPS input files. Builds Epoxy-PDMS, Amino-PDMS, and AM0270-POSS molecules, packs them at target density, and writes DREIDING-parameterised LAMMPS data + input scripts.

```bash
topon simbox [--output DIR] [--n-epoxy N] [--n-amino N] [--n-poss N]
             [--density FLOAT] [--seed INT]
```

| Option | Default | Description |
|---|---|---|
| `--output`, `-o` | `simbox_output` | Output directory |
| `--n-epoxy` | `50` | Number of Epoxy-PDMS molecules |
| `--n-amino` | `25` | Number of Amino-PDMS molecules |
| `--n-poss` | `10` | Number of AM0270-POSS molecules |
| `--density` | `0.85` | Target packing density (g/cm³) |
| `--seed` | `42` | Random seed for reproducible packing |

**Examples:**

```bash
# Default system (~85 molecules)
topon simbox

# Large production run, no POSS
topon simbox --output pdms_box --n-epoxy 600 --n-amino 300 --n-poss 0

# 50% POSS system
topon simbox --output poss50 --n-epoxy 60 --n-amino 15 --n-poss 15

# 100% POSS system
topon simbox --output poss100 --n-epoxy 60 --n-amino 0 --n-poss 30
```

**Output files:**

```
simbox_output/
├── system.data          # LAMMPS data file
├── ff_coeffs.in         # Force field coefficients (PairCoeffs, BondCoeffs, …)
├── 1_push_off.in        # Stage 1: soft-core push-off
├── 2_minimize.in        # Stage 2: energy minimization
├── 3_nvt_equilibrate.in # Stage 3: NVT equilibration
└── 4_npt_compress.in    # Stage 4: NPT compression to target density
```

Then run LAMMPS:

```bash
cd simbox_output && lmp -in 1_push_off.in
```

---

### `topon chain`

Build a single polymer chain in solvent and write DREIDING LAMMPS files.

```bash
topon chain --chain-smiles SMILES --dp N [options]
```

| Option | Default | Description |
|---|---|---|
| `--output`, `-o` | `chain_output` | Output directory |
| `--chain-smiles` | *(required)* | SMILES for the polymer repeat unit |
| `--dp` | *(required)* | Degree of polymerization |
| `--solvent-smiles` | `None` (toluene) | SMILES for single solvent molecule. Ignored if `--solvent-mixture` is set. |
| `--n-solvent` | `None` (auto) | Number of solvent molecules (auto-calculated if omitted) |
| `--solvent-mixture` | `None` | Multi-solvent JSON: `'[{"smiles":"...","weight_fraction":0.5}, ...]'` |
| `--graft-density` | `0.0` | Graft attachment probability per backbone unit (0–1) |
| `--graft-smiles` | `None` | SMILES for graft repeat unit (required if `--graft-density > 0`) |
| `--graft-dp` | `5` | Repeat units per side chain |
| `--density` | `0.85` | Target packing density (g/cm³) |
| `--seed` | `42` | Random seed |

**Examples:**

```bash
# PDMS chain in toluene
topon chain --chain-smiles "[Si](C)(C)O" --dp 20 \
            --solvent-smiles "Cc1ccccc1" --n-solvent 200

# Fluorinated PDMS chain in THF
topon chain --chain-smiles "[Si](C)(CCC(F)(F)F)O" --dp 15 \
            --solvent-smiles "C1CCOC1" --n-solvent 100 --output fpdms_in_thf

# Chain with grafts
topon chain --chain-smiles "[Si](C)(C)O" --dp 30 \
            --graft-density 0.1 --graft-smiles "[Si](C)(C)O" --graft-dp 5 \
            --solvent-smiles "Cc1ccccc1" --n-solvent 150
```

**Output files:**

```
chain_output/
├── system.data          # LAMMPS data file (chain + solvent, DREIDING)
├── ff_coeffs.in         # Force field coefficient commands
├── settings.in          # Pair coefficients for re-application after soft push-off
├── groups.txt           # Reactive group definitions
├── 1_minimize.in        # Stage 1: soft push-off + energy minimization
├── 2_nvt.in             # Stage 2: NVT equilibration
└── 3_npt.in             # Stage 3: NPT equilibration
```

Then run LAMMPS:

```bash
cd chain_output && lmp -in 1_minimize.in
```

---

### `topon analyze`

Analyze a topology graph file and print statistics. *(Currently prints placeholder — not yet wired to analysis module.)*

```bash
topon analyze GRAPH_PATH [--format text|json]
```

---

### `topon gui`

Launch the Streamlit web GUI for interactive configuration. *(Not yet implemented — requires `pip install topon[gui]`.)*

```bash
topon gui [--port PORT]
```

---

## Global Options

| Option | Description |
|---|---|
| `--version` | Show version and exit |
| `--help` | Show help message and exit |

Each command also accepts `--help`:

```bash
topon generate --help
topon simbox --help
```

---

## Python API (alternative to CLI)

All CLI commands are thin wrappers around importable functions:

```python
# generate equivalent
from topon.config import load_config
from topon.pipeline import Pipeline

config = load_config("examples/config_cg.json")
pipe = Pipeline(config)
pipe.run()

# simbox equivalent
from topon.simbox.workflow import run_workflow

run_workflow(
    output_dir="simbox_output",
    n_epoxy=600,
    n_amino=300,
    n_poss=0,
    density=0.85,
    seed=42,
)

# validate equivalent
from topon.config import load_config, validate_config

config = load_config("config.json")
errors = validate_config(config)

# chain equivalent
from topon.singlechain.workflow import run_workflow as chain_workflow

chain_workflow(
    output_dir="chain_output",
    chain_smiles="[Si](C)(C)O",
    dp=20,
    solvent_smiles="Cc1ccccc1",
    n_solvent=200,
    density=0.85,
    seed=42,
)
```

---

### `python -m topon.protein_network`

MARTINI 3 protein-network generator (sequence-string -> LAMMPS). Standalone
CLI -- not part of `topon ...` because the protein-network workflow is a peer
of `simbox` and `singlechain`, not a pipeline run.

```bash
python -m topon.protein_network <subcommand> [options]
```

Subcommands: `generate` (single run), `sweep` (water-content sweep into `wXX/`
subdirs), `topology` (BFM topology JSON only).

```bash
python -m topon.protein_network generate \
    --block-seq GGRPSDSYGAPGGGN --n-repeats 6 --n-chains 4 \
    --equil-steps 5000 --water-density 0 --output runs/resilin_dry --seed 42
```

Full flag reference is in [`docs/protein_network.md`](protein_network.md).
