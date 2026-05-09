# topon — Usage Guide

This is the reference for **running topon end-to-end**: CLI flags, sub-system APIs, recipes for the most common system types, and the JSON config schema (appendix).

For the package layout and design rationale, see [ARCHITECTURE.md](ARCHITECTURE.md) first.

---

## 1. Install

```bash
pip install -e .
```

Runtime dependencies: `numpy`, `networkx`, `rdkit`, `pydantic`, `scipy`. LAMMPS (`lmp`) must be on `PATH` only if you want to *run* the generated systems — generation itself does not need it.

---

## 2. Quick start

The fastest path from a checked-out repo to a runnable LAMMPS system:

```bash
# 1. Generate a starter config
topon init --output my_run.json

# 2. Validate it
topon validate my_run.json

# 3. Run the full six-stage pipeline
topon generate my_run.json --output ./runs
```

The output lands in `./runs/<study.name>/{topology, 02_Chemistry, 03_Conformation, 04_Simulation}/`. Stage 4 directory contains `system.data`; stage 6 contains the LAMMPS `.in` scripts.

To run the simulation (optional):

```bash
cd ./runs/<study.name>/04_Simulation/
lmp -in minimize_1_serial.in
```

For ready-to-use config files, see `examples/config_cg.json`, `examples/config_atomistic.json`, etc.

---

## 3. CLI reference

The package installs a `topon` console script that dispatches to several sub-commands. All sub-commands accept `--help`.

```
topon [--version] [--help] <command> [options]
```

### 3.1 `topon generate` — run the full pipeline

```bash
topon generate CONFIG_PATH [--output DIR] [--dry-run]
```

| Argument / Option | Description |
|---|---|
| `CONFIG_PATH` | Path to the JSON config file (required) |
| `--output`, `-o` | Override `study.output_dir` from config |
| `--dry-run` | Validate the config and exit without running |

Runs the six-stage pipeline (Topology → Analysis → Assignment → Chemistry → Conformation → Output). LAMMPS data files and input scripts are written to `output_dir/study_name/`.

```bash
topon generate examples/config_cg.json
topon generate examples/config_cg.json --output ./my_run
topon generate examples/config_cg.json --dry-run
```

### 3.2 `topon validate`

```bash
topon validate CONFIG_PATH
```

Prints `Configuration is valid!` or lists every validation error. Cheap pre-flight check before submitting a long run to HPC.

### 3.3 `topon init`

```bash
topon init [--output FILE] [--full]
```

| Option | Default | Description |
|---|---|---|
| `--output`, `-o` | `config.json` | Output path |
| `--full` | off | Include all options with default values (otherwise minimal) |

Edit the generated file to set topology source, DP, entanglements, etc.

### 3.4 `topon simbox` — pack a crosslink box

Independent sub-system. Builds Epoxy-PDMS / Amino-PDMS / AM0270-POSS molecules, packs them at target density, and writes DREIDING LAMMPS data + input scripts.

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

```bash
# Default ~85-molecule system
topon simbox

# Production: 600 epoxy + 300 amino, no POSS
topon simbox --output pdms_box --n-epoxy 600 --n-amino 300 --n-poss 0

# 50% POSS
topon simbox --output poss50 --n-epoxy 60 --n-amino 15 --n-poss 15
```

Output (single flat directory):

```
simbox_output/
├── system.data           # LAMMPS data file (atom_style full)
├── ff_coeffs.in          # Force-field coefficients
├── settings.in           # pair_coeff / bond_coeff / angle_coeff / dihedral_coeff
├── groups.txt            # group definitions by reactive-group type
├── 1_minimize.in         # Stage 1: soft push-off + CG minimisation
├── 2_nvt.in              # Stage 2: NVT thermalisation
├── 3_npt.in              # Stage 3: NPT density equilibration
└── 4b_crosslink.in       # Stage 4: crosslink template (fix bond/react or bond/create)
```

Then run LAMMPS:

```bash
cd simbox_output && lmp -in 1_minimize.in
```

See §4.2 for the simbox Python API.

### 3.5 `topon chain` — single chain in solvent

Build a single polymer chain in solvent and emit DREIDING LAMMPS files.

```bash
topon chain --chain-smiles SMILES --dp N [options]
```

| Option | Default | Description |
|---|---|---|
| `--output`, `-o` | `chain_output` | Output directory |
| `--chain-smiles` | *(required)* | SMILES for the polymer repeat unit |
| `--dp` | *(required)* | Degree of polymerization |
| `--solvent-smiles` | `None` (toluene fallback) | Single-solvent SMILES; ignored if `--solvent-mixture` set. If both are unset, `run_workflow` falls back to toluene. |
| `--n-solvent` | auto | Solvent molecules; auto-calculated from density if omitted |
| `--solvent-mixture` | `None` | Multi-solvent JSON: `'[{"smiles":"...","weight_fraction":0.5}, ...]'` |
| `--graft-density` | `0.0` | Graft attachment probability per backbone unit (0–1) |
| `--graft-smiles` | `None` | SMILES for graft repeat unit (required if `--graft-density > 0`) |
| `--graft-dp` | `5` | Repeat units per side chain |
| `--density` | `0.85` | Target packing density (g/cm³) |
| `--seed` | `42` | Random seed |

```bash
# PDMS in toluene
topon chain --chain-smiles "[Si](C)(C)O" --dp 20 \
            --solvent-smiles "Cc1ccccc1" --n-solvent 200

# Fluorinated PDMS in THF
topon chain --chain-smiles "[Si](C)(CCC(F)(F)F)O" --dp 15 \
            --solvent-smiles "C1CCOC1" --n-solvent 100 --output fpdms_in_thf

# With grafts
topon chain --chain-smiles "[Si](C)(C)O" --dp 30 \
            --graft-density 0.1 --graft-smiles "[Si](C)(C)O" --graft-dp 5 \
            --solvent-smiles "Cc1ccccc1" --n-solvent 150
```

Output:

```
chain_output/
├── system.data       # LAMMPS data (chain + solvent, DREIDING)
├── ff_coeffs.in
├── settings.in       # pair coefficients (re-applied after soft push-off)
├── groups.txt
├── 1_minimize.in
├── 2_nvt.in
└── 3_npt.in
```

### 3.6 `topon analyze` — graph statistics

Analyze a topology graph file and print statistics.

```bash
topon analyze GRAPH_PATH [--format text|json] [--nodes NODES_PATH]
```

| Argument / Option | Description |
|---|---|
| `GRAPH_PATH` | Path to `.gpickle`, `.nodes`, or `.edges` file |
| `--format`, `-f` | Output format: `text` (default) or `json` |
| `--nodes` | Companion `.nodes` file (required when `GRAPH_PATH` is a `.edges` file) |

```bash
topon analyze network.gpickle
topon analyze network.nodes
topon analyze network.edges --nodes network.nodes
topon analyze network.gpickle --format json
```

The CLI dispatches to `topon.analysis.report.analyze_graph()` and prints degree distribution, connectivity, and topology statistics. Available as both CLI and Python import.

### 3.7 `topon gui` — Streamlit GUI (not implemented)

```bash
topon gui [--port PORT]
```

> **Status: not yet implemented.** Reserved for a future `pip install topon[gui]` extra.

### 3.8 `python -m topon.protein_network` — topro

Protein-network generator (sequence-string → MARTINI 3 LAMMPS). **Standalone CLI** because the protein-network workflow is a peer of `simbox` and `singlechain`, not a `Pipeline` run. See §4.1 for the full topro reference.

```bash
python -m topon.protein_network <subcommand> [options]
```

Sub-commands: `generate`, `sweep`, `topology`. Common flags listed in §4.1.

### Global options

| Option | Description |
|---|---|
| `--version` | Show version and exit |
| `--help` | Show help message and exit |

---

## 4. Sub-systems

### 4.1 topro — protein networks

`topon.protein_network` (the user-facing name is **topro** — *topological protein network*) turns a residue-sequence string into a coarse-grained MARTINI 3 polymer-network LAMMPS input set. It mirrors the legacy CHARMM-atomistic protein-network workflow but emits ~28 MARTINI beads per 15-residue resilin block (vs ~250 atoms in the all-atom equivalent).

**At a glance:**

```
  block_seq + n_repeats + n_chains
              ↓
    bfm.generate_topology  →  protein_network_topology.json (snapshots)
              ↓
   builder.build_protein_system  + residues.py table
              ↓
        water.pack_water  (optional)
              ↓
   lammps_writer.write_lammps  (.data + .in.settings + .in.groups + relaxation/*.in)
```

**Inputs:**
- A one-letter repeat block such as the resilin consensus `GGRPSDSYGAPGGGN`
- Coverage today is 8 amino acids: `GLY`, `ALA`, `ARG`, `PRO`, `SER`, `ASP`, `TYR`, `ASN` (the resilin reference set). To extend, re-run `tools/extract_residues_from_itp.py` against a polyply ITP that contains the new residues; the table regenerates without code edits.
- TYR positions are interpreted as crosslinker sites (dityrosine); the BFM topology stage stochastically links lattice-adjacent TYR pairs (SC4–SC4, TN6 bead).
- Optional W water beads packed on a voxel grid around the protein.

**CLI:**

```bash
python -m topon.protein_network generate \
    --block-seq GGRPSDSYGAPGGGN \
    --n-repeats 6 --n-chains 4 \
    --equil-steps 5000 \
    --water-density 0 \
    --output runs/resilin_dry/ --seed 42
```

Sub-commands:
- `generate` — single run; writes data + settings + groups + input scripts + topology JSON
- `sweep` — repeat across a comma-separated list of water densities into `wXX/` sub-directories (XX = `density × 10` rounded; mirrors the legacy convention)
- `topology` — run only the BFM lattice stage; emit the snapshot JSON

**Common flags** (all sub-commands):

| Flag | Default | Meaning |
|---|---|---|
| `--block-seq` | `GGRPSDSYGAPGGGN` | One-letter repeat block |
| `--n-repeats` | `6` | Repeats per chain |
| `--n-chains` | `4` | Number of chains |
| `--segs-per-block` | `2` | BFM segments per repeat (2 = end + Y; 3 = +1 NC node) |
| `--equil-steps` | `5000` | Monte-Carlo BFM equilibration steps (`0` = skip) |
| `--target-packing` | `0.45` | Volume fraction used to size the lattice |
| `--min-intrachain-sep` | `2` | Minimum Y-index gap for intra-chain crosslinks |
| `--lattice-scale-ang` | `None` (auto) | Å per BFM lattice unit. Default auto-scales so BB–BB equilibrium length ≈ MARTINI 3.6 Å. |
| `--sc-jitter-ang` | `1.5` | Random offset of sidechain beads from BB (Å) |
| `--snapshot-label` | `gel_point` | Which BFM snapshot to build from (`gel_point`, `post_gel_1`, ...) |
| `--seed` | `42` | RNG seed |
| `--quiet` | off | Suppress per-stage status prints |
| `--hierarchical-stage1` | off | Use core-topon-style progressive freeze/unfreeze in stage 1 (safer for BFM-derived topology) |
| `--output` | *required* | Output directory (or JSON path for `topology` sub-command) |

**Water and ion flags** (`generate` and `sweep` only):

| Flag | Default | Used in | Meaning |
|---|---|---|---|
| `--water-density` | `0.0` | `generate` | Water beads per nm³ (`0` = dry, `~10` = bulk MARTINI water) |
| `--water-densities` | `"0,4,8"` | `sweep` | Comma-separated densities to sweep (subdirs `wXX/`) |
| `--water-exclusion` | `4.0` | `generate`, `sweep` | Min protein–water distance (Å) |
| `--water-bead` | `W` | `generate`, `sweep` | Water bead type: `W` (4 H₂O/bead, bulk), `SW` (3 H₂O, confined), `TW` (2 H₂O, very tight) |
| `--n-na-ions` | `0` | `generate`, `sweep` | NA⁺ ions to pack |
| `--n-cl-ions` | `0` | `generate`, `sweep` | CL⁻ ions to pack |

**Python API:**

```python
from topon.protein_network.workflow import run_protein_network

paths = run_protein_network(
    block_seq="GGRPSDSYGAPGGGN",
    n_repeats=6,
    n_chains=4,
    output_dir="runs/resilin_dry",
    equil_steps=5_000,
    water_density_w_per_nm3=0.0,   # ~10 for bulk MARTINI water
    seed=42,
)
# paths is a dict mapping artifact kind -> path:
#   topology_json                          BFM snapshots
#   data, settings, groups                 LAMMPS data + settings + groups files
#   stage1, stage2, stage3                 the three relaxation input scripts
```

Lower-level entry points (each module is independently usable):

- `topon.protein_network.bfm.generate_topology(...)` — BFM cubic-lattice topology
- `topon.protein_network.topology_io.{save_topology, load_topology, get_snapshot}` — JSON round-trip; format matches the legacy `topo_*.json`
- `topon.protein_network.builder.build_protein_system(snapshot, sequence_3letter, library, ...)` — chain chemistry from a snapshot
- `topon.protein_network.water.pack_water(system, library, ...)` — voxel-grid W packer
- `topon.protein_network.lammps_writer.write_lammps(system, library, output_dir, ...)` — LAMMPS file emitter
- `topon.protein_network.martini_ff.MartiniLibrary.from_package_data()` — vendored pruned MARTINI 3 protein FF + W water bead definition

**Output files:**

```
<output_dir>/
├── protein_network_topology.json   # BFM snapshots (gel_point, post_gel_1, ...)
├── protein_network.data            # LAMMPS Atoms (full), Bonds, Angles, Dihedrals, Impropers, Masses
├── protein_network.in.settings     # explicit pair_coeff / bond_coeff / etc.
├── protein_network.in.groups       # group protein/water molecule definitions
└── relaxation/
    ├── protein_network_stage1.in   # soft-push overlap removal
    ├── protein_network_stage2.in   # LJ epsilon ramp 0.001 → 1.0 (nve/limit)
    └── protein_network_stage3.in   # tight CG min + NVT/NPT @ 310 K
                                    #   → ../system_equilibrated.data
```

The 3-stage relaxation protocol mirrors the legacy CHARMM stages, ported to MARTINI 3 (`lj/cut/coul/cut`, `dielectric=15`, no PPPM, T=310 K). Each stage reads the previous stage's `write_data` output via `system_after_soft.data` / `system_ramped.data` handoff files.

**Running:**

```bash
# Serial
cd <output_dir>/relaxation/
lmp -in protein_network_stage1.in
lmp -in protein_network_stage2.in
lmp -in protein_network_stage3.in

# Parallel (HPC)
mpirun -np <N> lmp -in protein_network_stage1.in
# ... etc.
```

**Known approximations** (port of a GROMACS force field to LAMMPS — none affect topology, all affect energetics):

1. **Reaction-field electrostatics.** GROMACS `coulombtype = reaction-field, ε_r=15` is approximated with `pair_style lj/cut/coul/cut` + `dielectric 15.0`; the RF correction term is dropped. For quantitative electrostatic agreement, recompile LAMMPS with USER-MISC `pair_style coul/diel`.
2. **Restricted-bending angle** (GROMACS `funct=10`). MARTINI 3 IDP backbone uses `½ K (cos θ − cos θ₀)² / sin²θ`. Approximated with `angle_style cosine/squared`, dropping the `1/sin²θ` factor. Acceptable for IDP relaxation; review for folded chains.
3. **Multi-term proper dihedrals** (GROMACS `funct=9`). Each GROMACS quartet is emitted as a separate LAMMPS `dihedral_style charmm` coefficient set at the same atom indices. Energetically exact.
4. **Constraints** (GROMACS `[ constraints ]`). Emitted as ULTRA-stiff harmonic bonds (`K = 1e6 kJ/mol/nm²`, matching the reference's `#ifdef FLEXIBLE`). For true rigidity, switch to `fix shake` — the input script header notes how.
5. **Per-atom exclusions** (`[ exclusions ]`). Beyond-1-2 TYR-ring exclusions are not emitted as data-file rows; the input script uses `special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0` instead. More aggressive than the reference's `nrexcl=1` plus selective ring exclusions, but yields correct intra-ring behaviour with no per-pair plumbing.

**Not yet supported** (tracked in `internal/DEVELOPMENT_INTERNAL.md`):
- Folded-protein parameterisations (virtual_sites_2/3/4, Go contact maps)
- Elastic network restraints
- NaCl and other ions (vendored bead types are kept; no packer yet)
- MARTINI inside `simbox/` (simbox stays DREIDING-only)
- Residues beyond the 8-AA resilin set (extend by re-running `tools/extract_residues_from_itp.py`)

**Citations** for the vendored MARTINI 3 force field:

- Souza et al., *Nature Methods* 2021, 10.1038/s41592-021-01098-3
- Kroon et al., *eLife* 2024, 10.7554/elife.90627.2
- Grunewald et al., *Nature Communications* 2022, 10.1038/s41467-021-27627-4

### 4.2 simbox — molecule packing

`topon.simbox` packs individual molecules into a periodic simulation box and emits LAMMPS input scripts for crosslinking studies. Independent of the polymer-network pipeline.

**Workflow:**

```
MoleculeLibrary       →  Molecule objects (RDKit mol + reactive-site annotations)
       ↓
BoxPacker.pack()      →  PackedBox (placed molecules with 3D coordinates)
       ↓
assemble(packed)      →  AssembledSystem (merged RDKit mol, reactive-site registry)
       ↓
write_lammps(system)  →  system.data, settings.in, groups.txt, ff_coeffs.in
       ↓
write_inputs(system)  →  1_minimize.in, 2_nvt.in, 3_npt.in, 4b_crosslink.in
```

Quickest path: `topon simbox` (§3.4) or `topon.simbox.workflow.run_workflow()`.

**`Molecule`** (`topon/simbox/molecule.py`) — RDKit Mol with explicit H, ETKDGv3 + MMFF-optimised 3D conformer, and reactive-site annotations auto-detected via SMARTS:

| Site | SMARTS |
|---|---|
| `epoxide` | `[C]1[O][C]1` |
| `primary_amine` | `[NX3;H2;!$([NH2]C=O)]` |
| `secondary_amine` | `[NX3;H1]([#6])[#6]` |

```python
from topon.simbox.molecule import Molecule
mol = Molecule.from_smiles("EpoxyPDMS", "C1OC1COCCC[Si](C)(C)O...")
mol = Molecule.from_pdb("MyMol", "path/to/file.pdb")
mol = Molecule.from_mol("MyMol", rdkit_mol_object)
```

**`MoleculeLibrary`** (`topon/simbox/library.py`) — pre-built siloxane molecules:

```python
from topon.simbox.library import MoleculeLibrary
lib = MoleculeLibrary()
epoxy  = lib.epoxy_pdms(n_dms=2)    # Glycidoxypropyl-PDMS, ~500 g/mol
amino  = lib.amino_pdms(n_dms=8)    # Aminopropyl-PDMS, ~850 g/mol
poss   = lib.am0270_poss()           # AminopropylIsooctyl POSS, ~1267 g/mol
custom = lib.custom("C1OC1", name="MyEpoxide")
```

Structures:
- **Epoxy-PDMS**: `Epoxide-CH₂-O-CH₂CH₂CH₂-Si(Me)-[O-Si(Me)₂]ₙ-O-Si(Me)-CH₂CH₂CH₂-O-CH₂-Epoxide`
- **Amino-PDMS**: `H₂N-CH₂CH₂CH₂-Si(Me)-[O-Si(Me)₂]ₙ-O-Si(Me)-CH₂CH₂CH₂-NH₂`
- **AM0270 POSS**: Si₈O₁₂ cube cage, corner 0 with `-CH₂CH₂CH₂-NH₂`, corners 1–7 with 2,4,4-trimethylpentyl (isooctyl, inert)

**`BoxPacker`** — grid-based spatial hashing for O(N) overlap detection:

```python
from topon.simbox.packer import BoxPacker

packer = BoxPacker(
    density=0.85,        # g/cm³
    min_dist=2.0,        # Å
    seed=42,
    max_attempts=1000,   # placement attempts per molecule
    growth_factor=1.05,  # box expansion when packing fails
)
packed = packer.pack([(epoxy, 100), (amino, 50), (poss, 10)])
```

Algorithm: compute initial box from total mass and target density → shuffle insertion order → for each molecule random rotation (Shoemake quaternion) + random translation, with min-image overlap check → if `max_attempts` exceeded, grow box by `growth_factor` and retry (up to 20 rounds).

**`AssembledSystem`** — merged Mol with global bookkeeping:

```python
from topon.simbox.system import assemble
system = assemble(packed)
# system.mol               — merged RDKit Mol
# system.box_lengths       — ndarray([Lx, Ly, Lz]) in Å
# system.molecule_ids      — per-atom LAMMPS molecule ID (1-based)
# system.species_names     — per-molecule species name
# system.reactive_sites    — list of ReactiveSiteEntry (global atom index + group name)
```

**Writers:**

```python
from topon.simbox.writer import write_lammps
from topon.simbox.inputs import write_inputs

write_lammps(system, output_dir="output/simbox")
write_inputs(system, output_dir="output/simbox", temperature=300.0, pressure=1.0)
```

Stage 1 — soft push-off + minimisation:
- Phase A: `pair_style soft` with ramped prefactor (0→60) + brief NVT to resolve overlaps
- Phase B: switch to `lj/cut` DREIDING potentials + conjugate-gradient minimisation

Stage 4b — crosslink template (user must configure):
- **Option A** (`fix bond/react`): template-based reactions with molecule pre/post files
- **Option B** (`fix bond/create`): simple distance-based bond formation

**One-call workflow** (canonical entry point):

```python
from topon.simbox.workflow import run_workflow

files = run_workflow(
    output_dir="output/simbox_run",
    n_epoxy=600, n_amino=300, n_poss=0,
    density=0.85, seed=42,
)
```

`run_workflow` also activates `UniversalTypeMapper`, a context manager that patches `topon.forcefield.dreiding` at write time to enforce stable DREIDING type IDs across all compositions, keeping pre-defined `fix bond/react` templates compatible.

### 4.3 singlechain — solubility utility

`topon.singlechain` builds a single polymer chain in a solvent box for solubility studies. Use the CLI (§3.5) or the Python entry point:

```python
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

## 5. Recipes

Each recipe shows a config + the command. Most knobs live in the JSON config; see Appendix A for the full schema.

### 5.1 CG network with entanglements

```json
{
  "study": { "name": "cg_entangled", "output_dir": "./runs" },
  "topology": {
    "source": "generate",
    "generator": {
      "lattice_size": "6x6x6",
      "lattice_type": "SC",
      "max_functionality": 4,
      "degree_distribution": "0:13,1:25"
    }
  },
  "assignment": {
    "node_types": { "method": "degree", "degree": { "mapping": {"1": "end", "2": "A", "3": "A", "4": "A"} } },
    "edge_types": { "method": "uniform", "uniform": { "type": "A" } },
    "dp_distribution": { "default": { "mean": 25, "pdi": 1.0 } },
    "entanglements": {
      "enabled": true, "target": 5, "target_type": "count",
      "kink_params": { "overshoot": 0.2, "z_amp": 0.5, "sigma": 0.15 }
    }
  },
  "chemistry": {
    "model_type": "coarse_grained",
    "node_type_map": {
      "end": { "molecule": "Si", "is_end_cap": true },
      "A":   { "molecule": "Si", "is_end_cap": false }
    },
    "edge_type_map": { "A": { "monomer": "PDMS" } }
  }
}
```

```bash
topon generate cg_entangled.json --output ./runs
```

### 5.2 Atomistic network with POSS junctions

```json
{
  "study": { "name": "atomistic_poss", "output_dir": "./runs" },
  "topology": { "source": "load",
    "existing_files": { "nodes_file": "output/network.nodes",
                        "edges_file": "output/network.edges" } },
  "assignment": {
    "node_types": { "method": "degree",
      "degree": { "mapping": {"1": "end", "2": "A", "3": "A", "4": "POSS"} } },
    "edge_types": { "method": "uniform", "uniform": { "type": "A" } },
    "dp_distribution": { "default": { "mean": 10, "pdi": 1.0 } }
  },
  "chemistry": {
    "model_type": "atomistic", "target_density": 1.1,
    "node_type_map": {
      "end":  { "molecule": "[Si](C)(C)C", "is_end_cap": true },
      "A":    { "molecule": "Si",          "is_end_cap": false },
      "POSS": { "molecule": "POSS_AM0270", "is_end_cap": false }
    },
    "edge_type_map": { "A": { "monomer": "PDMS" } }
  }
}
```

```bash
topon generate atomistic_poss.json --output ./runs
```

### 5.3 Atomistic network with grafts and entanglements

Combine the entanglement recipe (5.1) with a `grafts` block in `assignment`:

```json
"grafts": {
  "enabled": true,
  "per_edge_type": {
    "A": { "graft_density": 0.05, "side_chain_monomer": "PDMS", "side_chain_dp": 5 }
  }
}
```

A complete combined-features workflow is also available as a script: `tests/workflows/generate_atomistic_combined.py`.

### 5.4 simbox crosslink workflow

```bash
topon simbox --output runs/simbox_crosslink \
             --n-epoxy 600 --n-amino 300 --n-poss 0 --seed 42
```

```bash
cd runs/simbox_crosslink
lmp -in 1_minimize.in
lmp -in 2_nvt.in
lmp -in 3_npt.in
# Edit 4b_crosslink.in to choose Option A (fix bond/react) or B (fix bond/create), then:
lmp -in 4b_crosslink.in
```

### 5.5 topro resilin network (dry, then wet)

```bash
# Dry: 16 chains × 6 repeats, no water
python -m topon.protein_network generate \
    --block-seq GGRPSDSYGAPGGGN --n-repeats 6 --n-chains 16 \
    --water-density 0 --output runs/resilin_dry --seed 42

# Sweep across water densities
python -m topon.protein_network sweep \
    --block-seq GGRPSDSYGAPGGGN --n-repeats 6 --n-chains 16 \
    --water-densities 0,1,4 --output runs/resilin_sweep --seed 42
```

---

## 6. Python API (alternative to CLI)

The CLI is a thin wrapper. Equivalent Python:

```python
# topon generate equivalent — full pipeline
from topon.config import load_config
from topon.pipeline import Pipeline

config = load_config("examples/config_cg.json")
pipe = Pipeline(config)
pipe.run()

# topon validate equivalent
from topon.config import load_config, validate_config
errors = validate_config(load_config("config.json"))

# topon simbox equivalent
from topon.simbox.workflow import run_workflow
run_workflow("simbox_output", n_epoxy=600, n_amino=300, n_poss=0,
             density=0.85, seed=42)

# topon chain equivalent
from topon.singlechain.workflow import run_workflow as chain_workflow
chain_workflow("chain_output", chain_smiles="[Si](C)(C)O", dp=20,
               solvent_smiles="Cc1ccccc1", n_solvent=200,
               density=0.85, seed=42)

# topro equivalent
from topon.protein_network.workflow import run_protein_network
run_protein_network(block_seq="GGRPSDSYGAPGGGN", n_repeats=6, n_chains=4,
                    output_dir="runs/resilin_dry", equil_steps=5_000,
                    water_density_w_per_nm3=0.0, seed=42)
```

For lower-level entry points (each stage individually), see ARCHITECTURE.md §2 — the per-stage `Module` line tells you which import path drives that stage.

---

## 7. Workflow scripts

`tests/workflows/` contains executable Python scripts that wrap the Pipeline + sub-systems for end-to-end studies. These are *not* part of the package API; they are the project's own driver scripts (run them with `python <path>`).

| Script | Purpose |
|---|---|
| `generate_atomistic_combined.py` | Atomistic with grafts + entanglements; uses the Pipeline |
| `generate_cg_combined.py` | CG with combined features |
| `generate_simbox_crosslink.py` | Drives `simbox.workflow.run_workflow` end-to-end |
| `run_v41_matrix.py`, `run_v42_matrix.py`, `run_v43_core_topon.py` | Versioned sweep drivers — parameterise the workflows above and write into `tests/output/v<NN>/<cell>/` |
| `analyze_v41_collapse.py` | Post-processing analysis for the v41 sweep |

Output goes to `tests/output/v<NN>/<cell>/` (gitignored). Don't introduce a parallel `runs/` folder — the convention is `tests/output/`.

---

## 8. Running tests

```bash
# Unit tests (~5 s)
pytest tests/unit/ -v

# Regression tests (~1.5 h — generates full systems and byte-compares against frozen references)
pytest tests/regression/ -v
```

Regression references are frozen with `seed=42`. Before modifying anything in `topon/writers/` or `topon/simbox/writer.py`, run regression first to confirm the baseline; make the change; re-run. For `topon/protein_network/lammps_writer.py` the minimum is `pytest tests/unit/`.

---

## Appendix A — JSON config schema

`topon` is configured by a single JSON file. Use `topon init` to generate a starter; edit as needed.

Top-level sections:

```json
{
  "study":      { ... },
  "topology":   { ... },
  "assignment": { ... },
  "chemistry":  { ... },
  "output":     { ... }
}
```

### `study`

| Key | Type | Default | Description |
|---|---|---|---|
| `name` | string | `"my_network"` | Study name; used as sub-directory under `output_dir` |
| `output_dir` | string | `"./output"` | Root output directory |

### `topology`

| Key | Type | Default | Description |
|---|---|---|---|
| `source` | `"generate"` \| `"load"` | `"load"` | Generate a new topology or load an existing one |
| `generator` | object | — | Settings for the C / Python generator (when `source="generate"`) |
| `existing_files` | object | — | File paths (when `source="load"`) |

#### `topology.generator`

| Key | Type | Default | Description |
|---|---|---|---|
| `exe_path` | string \| null | `null` | Path to `generator.exe`; `null` → use Python generator |
| `lattice_size` | string | `"6x6x6"` | Lattice dimensions, e.g. `"8x8x8"` |
| `lattice_type` | `"SC"` \| `"BCC"` \| `"FCC"` | `"SC"` | Lattice type |
| `periodicity` | string | `"111"` | Periodicity per axis (`1`=periodic, `0`=open) |
| `max_functionality` | int | `6` | Maximum crosslink degree per node |
| `max_trials` | int | `1000000` | Trials before giving up |
| `max_saves` | int | `1` | Number of networks to save |
| `degree_distribution` | string | `"0:0,1:0"` | Target degree distribution |

Degree distribution format: `"d:N"` requires N nodes of degree d; `"e:N"` requires N edges total; omitted degrees are unconstrained. Example: `"0:15,1:30,e:371"`.

#### `topology.existing_files`

Provide either `gpickle_file` OR both `nodes_file` + `edges_file`.

| Key | Type | Default | Description |
|---|---|---|---|
| `nodes_file` | string \| null | `null` | Path to `.nodes` file |
| `edges_file` | string \| null | `null` | Path to `.edges` file |
| `gpickle_file` | string \| null | `null` | Path to NetworkX `.gpickle` file |

### `assignment`

Controls how graph attributes (types, DP, defects, entanglements, grafts, copolymers) are written before chemistry is built.

#### `assignment.node_types`

`method` ∈ `"degree"` / `"positional"` / `"random"` / `"explicit"`.

```json
"node_types": { "method": "degree", "degree": {"mapping": {"1": "end", "2": "A", "3": "A", "4": "A"}} }
"node_types": { "method": "positional", "positional": {"dimension": "z", "num_layers": 2, "layer_types": ["A","B"]} }
"node_types": { "method": "random", "random": {"type_ratios": {"A": 70, "B": 30}} }
"node_types": { "method": "explicit", "explicit": {"0": "POSS", "1": "Si"} }
```

#### `assignment.edge_types`

`method` ∈ `"uniform"` / `"random"` / `"composite"`.

```json
"edge_types": { "method": "uniform", "uniform": {"type": "A"} }
"edge_types": { "method": "random", "random": {"type_ratios": {"A": 60, "B": 40}} }
"edge_types": { "method": "composite", "composite": {"dimension": "z", "num_layers": 3, "layer_types": ["A","B","A"]} }
```

#### `assignment.dp_distribution`

```json
"dp_distribution": {
  "default": { "mean": 25, "pdi": 1.0 },
  "per_edge_type": {
    "A": { "mean": 20, "pdi": 1.2 },
    "B": { "mean": 40, "pdi": 1.5 }
  }
}
```

`pdi` = polydispersity index (Schulz-Zimm distribution). `1.0` = monodisperse.

#### `assignment.defects`

```json
"defects": {
  "primary_loops": { "enabled": true, "target": 10, "target_type": "count" }
}
```

#### `assignment.entanglements`

```json
"entanglements": {
  "enabled": true, "target": 5, "target_type": "count",
  "kink_params": { "overshoot": 0.2, "z_amp": 0.5, "sigma": 0.15 }
}
```

Or distribution mode (average per chain):

```json
"entanglements": {
  "enabled": true,
  "avg_crosslinks_per_chain": 2.0,
  "kink_params": { "overshoot": 0.2, "z_amp": 0.5, "sigma": 0.15 }
}
```

| `kink_params` key | Default | Description |
|---|---|---|
| `overshoot` | `0.2` | How far the kink extends past the midpoint (0–1) |
| `z_amp` | `0.5` | Out-of-plane amplitude of the Gaussian kink |
| `sigma` | `0.15` | Width of the Gaussian kink |

#### `assignment.grafts`

```json
"grafts": {
  "enabled": true,
  "per_edge_type": {
    "A": { "graft_density": 0.05, "side_chain_monomer": "PDMS", "side_chain_dp": 5 }
  }
}
```

#### `assignment.copolymer`

```json
"copolymer": {
  "enabled": true,
  "per_edge_type": {
    "A": {
      "arrangement": "block",
      "composition": [
        { "monomer": "A", "fraction": 0.5 },
        { "monomer": "B", "fraction": 0.5 }
      ]
    }
  }
}
```

`arrangement` ∈ `"block"` / `"alternating"` / `"random"` / `"gradient"`.

### `chemistry`

| Key | Type | Default | Description |
|---|---|---|---|
| `model_type` | `"coarse_grained"` \| `"atomistic"` | `"coarse_grained"` | Force-field resolution |
| `target_density` | float | `0.9` | Target density in g/cm³ |

#### `chemistry.node_type_map`

```json
"node_type_map": {
  "end": { "molecule": "[Si](C)(C)C", "is_end_cap": true },
  "A":   { "molecule": "Si",          "is_end_cap": false },
  "B":   { "molecule": "POSS",        "is_end_cap": false }
}
```

Built-in molecule names: `"Si"`, `"POSS"` (Si₈O₁₂ cage), `"POSS_AM0270"` (AM0270 aminopropyl POSS). Any SMILES string is also accepted.

#### `chemistry.edge_type_map`

```json
"edge_type_map": { "A": { "monomer": "PDMS" }, "B": { "monomer": "FPDMS" } }
```

#### `chemistry.monomers`

Built-in defaults:

| Name | SMILES | Description |
|---|---|---|
| `PDMS` | `[Si](C)(C)O` | Polydimethylsiloxane |
| `FPDMS` | `[Si](C)(CCC(F)(F)F)O` | Fluorinated PDMS |
| `Phenyl` | `[Si](C)(c1ccccc1)O` | Phenyl-PDMS |

Add custom monomers:

```json
"monomers": {
  "MyMonomer": {
    "smiles": "[Si](CC)(CC)O",
    "chain_head": "Si",
    "chain_tail": "O"
  }
}
```

#### `chemistry.connection`

```json
"connection": { "auto_bridge": true, "default_bridge_atom": "O" }
```

`auto_bridge`: when the chain head and node atom are the same element (e.g. both Si), automatically inserts a bridge atom. Set `false` to always use direct bonds.

### `output`

| Key | Default | Description |
|---|---|---|
| `lammps_data` | `true` | Write LAMMPS data file |
| `lammps_inputs` | `true` | Write LAMMPS input scripts |
| `visualization` | `true` | Write HTML visualization |
| `analysis_report` | `true` | Write analysis report |
| `save_attributed_graph` | `true` | Save attributed graph as `.gpickle` |
