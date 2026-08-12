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

### 2.1 Interactive shell (recommended)

Run `topon` (or `python -m topon`) with no arguments on a real terminal. You land in the `topon>` REPL, where every subcommand is one token away:

```text
$ python -m topon
   +================ ... banner ... ================+
   Type `help` for the command list, or `exit` to leave.

topon> init --preset cg_kg --output my_run.json
Wrote my_run.json (preset: cg_kg, copied from config.json)

topon> doctor my_run.json
[ok]    schema_gap_extras: ...
Summary: 0 error / 0 warn / 1 ok

topon> generate my_run.json
... runs the 6-stage pipeline ...

topon> inspect output_my_run
... atom counts, box, next LAMMPS commands ...

topon> exit
bye.
```

Shell built-ins: `help`, `help <cmd>`, `exit | quit | q | Ctrl-D`. Arrow keys for history if `readline` is installed. `help <cmd>` shows the same click `--help` page you'd see in one-shot mode.

The shell auto-launches when stdin is a TTY. To force it from a non-TTY context use `topon shell`; to skip it and just print the banner, use `topon --no-shell`.

### 2.2 One-shot mode

Same commands, but each one re-invokes the CLI:

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

### 3.3 `topon init` — starter config that runs as-is

```bash
topon init                              # fastest: write atomistic_pdms preset
topon init --preset cg_kg               # different preset
topon init --interactive                # prompt-driven walk through 6 knobs
topon init --preset martini_resilin     # print the right MARTINI CLI invocation
```

| Option | Default | Description |
|---|---|---|
| `--output`, `-o` | `config.json` | Path for the new config file |
| `--preset` | `atomistic_pdms` | One of `atomistic_pdms`, `cg_kg`, `poss`, `martini_resilin`, `charmm_resilin`. The first three copy a bundled demo `config.json`; the last two print the right `python -m` invocation (those paths use a separate CLI). |
| `--interactive`, `-i` | off | Prompt for the 5–6 knobs that actually vary (study name, output dir, model type, lattice type+size, max functionality, DP, density) and write the result. |

The non-interactive default copies `examples/demos/polymer/atomistic/basic/config.json`. Every preset-produced file passes `topon validate` immediately.

### 3.3a `topon doctor` — semantic lint

```bash
topon doctor my_run.json           # informational + warns
topon doctor my_run.json --strict  # warns also exit 1
```

Where `validate` is a Pydantic schema check, `doctor` runs a small rule registry sourced from known footguns (`internal/DEVELOPMENT_INTERNAL.md` issues + things new users trip on). Current rules:

| Rule | Level | Catches |
|---|---|---|
| `lattice_size_format` | error | `"lattice_size": 5` instead of `"5x5x5"` |
| `unknown_node_type` | warn | `assignment.node_types.degree.mapping` references a type that's not in `chemistry.node_type_map` (would silently fall through to Si — P0-2) |
| `poss_at_internal_junction` | warn | POSS mapped to degree >= 2 (hits known bug P1-H) |
| `atomistic_graft_non_pdms` | warn | Graft density set on a non-PDMS atomistic monomer (currently silently skipped) |
| `dp_below_kuhn` | warn | DP < 5 — conformation/entanglement edge cases |
| `defects_endcap_safe` | ok | Reminder that primary-loop defects safely skip end-caps post-2026-05-10 |
| `schema_gap_extras` | ok | Config has `conformation`/`simulation`/`execution` (not Pydantic-validated; CLI handles via `load_config_full`) |

Adding a new rule: write `check_<name>(cfg, raw) -> list[Issue]` in `topon/diagnostics/rules.py` and append to `RULE_REGISTRY`.

### 3.3b `topon inspect <run_dir>` — post-run summary

```bash
topon inspect runs/my_study
topon inspect examples/demos/polymer/atomistic/graft/expected_output
```

Replaces hand-grepping `system.data` headers after a long pipeline. Parses each stage directory (Pipeline layout: `02_Chemistry/`, `03_Conformation/`, `04_Simulation/`) or the flat `expected_output/`-style layout, and prints:

- atom count, atom-type count, box dimensions
- per-stage status (which files landed, what they say)
- the next LAMMPS commands to run

### 3.3c `topon recipes` — common use cases

```bash
topon recipes
```

Prints a "I want X -> run Y" cheatsheet covering all sub-systems (polymer networks via Pipeline, MARTINI/CHARMM protein networks, simbox, single-chain, batch workflows). Edit `topon/cli.py:recipes()` to add rows.

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

#### CHARMM36m atomistic builder

`topon.protein_network.charmm` is the **all-atom** alternative to the MARTINI path
above: same BFM topology JSON, same sequence string, but ~250 atoms per 15-residue
block instead of ~28 beads, with the CHARMM36m force field (`lj/charmm/coul/long`,
PPPM, CMAP, TIP3P water + NaCl).

```bash
python -m topon.protein_network.charmm.build_systems \
    --topology topo.json --snapshot gel_point \
    --block_seq GGRPSDSYGAPGGGN --n_repeats 18 \
    --water_contents 0,35,55,65,75 \
    --output runs/resilin_atomistic/
```

Writes one `w<XX>/` per water content, each with `protein_network.data`,
`.in.settings`, `.in.groups`, `charmm36m.cmap`, and `relaxation/protein_network_stage{1,2,3}.in`.

| Flag | Default | Meaning |
|---|---|---|
| `--topology` | *required* | BFM topology JSON (from `bfm.generate_topology`) |
| `--snapshot` | `gel_point` | Snapshot label or integer index |
| `--block_seq` | `GGRPSDSYGAPGGGN` | One-letter repeat block |
| `--n_repeats` | *from topology* | Override repeats per chain |
| `--water_contents` | `0,35,55,65,75` | Comma-separated water **weight percent** values |
| `--salt_conc` | `0.15` | Background NaCl (mol/L) |
| `--target_density` | `0.85` | Initial density (g/cm³) used to size the box |
| `--lattice_scale` | auto | Å per BFM unit (auto-sized from `--target_density`) |
| `--no-image-flags` | off | Emit legacy **7-column** Atoms and keep **all** crosslinks (no winding drops). Exact topology, but **single-rank only**. Default emits 10-column image flags and drops winding-cycle crosslinks (MPI-safe). |
| `--physical-backbone` | off | Build physically correct geometry — see below |
| `--xpro-cis-fraction` | `0.0` | With `--physical-backbone`, fraction of X-Pro peptide bonds seeded **cis** (~`0.05` = the folded-protein value, i.e. 95–96 % trans) |

**Stage 1 must run serial.** It switches to `pair_style soft 1.0`, giving a ~3 Å
communication cutoff — shorter than the ~5.5 Å longest bond — so under MPI domain
decomposition a bond straddling a domain boundary loses its partner
(*"Bond atom missing in image check"*). Run stage 1 on one rank; stages 2/3 are
MPI-safe (stage 2 sets `comm_modify cutoff 14`):

```bash
cd <out>/w55/relaxation/
mpirun -n 1  lmp -in protein_network_stage1.in      # serial by design
mpirun -n 104 lmp -in protein_network_stage2.in
mpirun -n 104 lmp -in protein_network_stage3.in     # -> ../system_equilibrated.data
```

**`--physical-backbone`** (opt-in; the default placement is unchanged). The default
builder drops each residue's atoms at its lattice anchor with a small random jitter
and lets minimisation sort out the geometry. That leaves the *barrier-locked*
degrees of freedom to chance: minimisation falls ~50/50 into the cis and trans
basins, and the ~20 kcal/mol omega barrier then freezes the result — measured
~12 % cis on non-proline peptide bonds (physical: <0.1 %) and ~50/50 D/L CA
chirality. With the flag, every atom is placed from the CHARMM RTF internal
-coordinate tables (real bonds/angles, planar impropers, real rotamers, 100 % L
-chirality), the backbone is coiled to ~3.8 Å CA–CA so minimisation needs no
violent expansion, and the writer adds `fix restrain` blocks that hold omega trans
(minus the `--xpro-cis-fraction` X-Pro subset) and the CA chirality improper at L
through stages 1–3, released before each stage's dynamics.

Verified on a 25 × 18 natpro system through stages 1–3: **non-Pro 0.03 % cis,
X-Pro 8 % cis, 98.9 % trans overall, chirality 100 % L, bond median 1.23 Å.**

**Crosslink methods** (`bfm.generate_topology(crosslink_method=...)`, API-only):

| Value | Meaning |
|---|---|
| `"adjacent"` (default) | Lattice-adjacent Y pairs react; snapshots at the gel point and beyond |
| `"winding_safe"` | Same, but rejects crosslinks that would close a periodic winding cycle → the writer drops **zero** bonds |
| `"distance"` | Distance-based candidate search (honours `pre_gel_conversions`) |
| `"none"` | **No crosslinking.** Emits a single conv=0 snapshot labelled `uncrosslinked` with `reactions=[]`. The starting point for **in-situ** crosslinking (form the dityrosine bonds during MD with LAMMPS `fix bond/react`) rather than stitching them at build time. |

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
| `verify_lattice_box.py` | Audits the recorded periodic cell on all four lattices, then runs a BCC network through the pipeline + LAMMPS stage 1 |
| `verify_mixed_lattice.py` | Builds a mixture end-to-end and reports the edge-length shells and bond-length tail against an SC baseline |
| `compare_generators.py` | Sweeps the C searcher against the Python generator over lattices, sizes, mixtures and distribution modes; `--lammps` also builds and minimises a subset |

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
| `lattice_type` | `"SC"` \| `"BCC"` \| `"FCC"` \| `"Diamond"` \| `"MIX"` | `"SC"` | Lattice type; `MIX` overlays SC/BCC/FCC (see below) |
| `mix_fractions` | object | `{"SC":1,"BCC":0,"FCC":0}` | Sublattice fractions for `MIX`; must sum to 1 |
| `mix_cutoff` | float | `1.0` | Neighbour cutoff for `MIX`, in cell units |
| `periodicity` | string | `"111"` | Periodicity per axis (`1`=periodic, `0`=open); see below |
| `max_functionality` | int | `6` | Maximum crosslink degree per node |
| `max_trials` | int | `1000000` | Trials before giving up |
| `max_saves` | int | `1` | Number of networks to save |
| `degree_distribution` | string | `"0:0,1:0"` | Target degree distribution |

Degree distribution format: `"d:N"` requires N nodes of degree d; `"e:N"` requires N edges total; omitted degrees are unconstrained. Example: `"0:15,1:30,e:371"`.

##### Diamond (`lattice_type: "Diamond"`)

Two interpenetrating FCC sublattices offset by ¼ along the body diagonal:
8 sites per cubic cell, **every site exactly 4-coordinated by
construction**. A `max_functionality: 4` network therefore needs no
pruning at all, which makes it the cleanest backbone for a tetrafunctional
network and much faster to generate than sculpting SC or FCC down to 4.

```json
"generator": { "lattice_type": "Diamond", "lattice_size": "6x6x6",
               "max_functionality": 4, "degree_distribution": "" }
```

An `NxNxN` Diamond has `8N³` sites and `16N³` bonds at a nearest-neighbour
distance of `√3/4 ≈ 0.433` cells. Both generators build it identically.

##### Boundaries (`periodicity`)

One digit per axis, `1` periodic and `0` open. An open axis omits its
wrap-around bonds, so the lattice grows a **free surface** there and the
sites on it lose coordination. The site set is unchanged either way.

```json
"generator": { "lattice_size": "6x6x6", "periodicity": "110" }
```

That builds a slab: periodic in x and y, open in z. On a 4x4x4 SC lattice
the bond count goes 192 → 176 → 160 → 144 as you open one, two and three
axes, and the surface sites drop from degree 6 to 5.

**Open boundaries interact with `degree_distribution`.** Corner and edge
sites on a free surface have very low coordination, and the centred
lattices lose the most:

| lattice (4x4x4) | min degree, `"111"` | `"110"` | `"000"` |
|---|---|---|---|
| SC | 6 | 5 | 3 |
| BCC | 8 | 4 | 1 (2 such sites) |
| FCC | 12 | 8 | 3 |
| Diamond | 4 | 2 | 1 (22 such sites) |

So the usual `"0:0,1:0"` (no isolated nodes, no dangling ends) is
**unsatisfiable** on a fully open BCC or Diamond: the only way to clear a
degree-1 site is to cut its last bond, which makes it degree 0, and that
is forbidden too. Both generators decline rather than claim success. On
SC the minimum stays at 3, so the same request is fine. Either drop the
`1:0` term on open lattices, or leave `degree_distribution` empty and let
`max_functionality` do the work.

`max_functionality` still applies on top, so a partially open lattice
reaches the ceiling with less pruning than a closed one.

**What an open axis does to the data file.** Coordinates are wrapped into
the box only on periodic axes. An open axis keeps its atoms where they
were placed and the box grows to contain them, so a junction on the free
surface stays next to the chains bonded to it instead of being split
across the cell. The `.nodes` file records the boundaries in a
`# PERIODICITY 100` header (written only when an axis is open), and the
conformation stage reads it.

An open axis also gets **12 Å of vacuum** between the outermost atom and
the box face, matching the pair cutoff the generated scripts use. That is
not cosmetic: LAMMPS *deletes* atoms that leave a non-periodic (`f`) face,
and the geometry handed to stage 1 is strained enough that surface atoms
move several Å in the first few dozen steps. With only 1 Å of clearance a
bonded atom was lost at step 49. Override with `open_axis_pad` if a run
needs more, or less when the extra volume matters.

The generated LAMMPS scripts still say `boundary p p p` — they are not
periodicity-aware, and changing them is out of scope here. Set
`boundary p f f` yourself to match; the data file is already correct for
it. Verified on a Diamond `100` network: `p f f` and `p p p` both
complete stage-1 minimization, with no bond crossing an open face.

##### Mixed lattices (`lattice_type: "MIX"`)

All three cubic lattices share the cell corner and each adds sites on top
of it: BCC one body centre, FCC three face centres. `MIX` puts the corner
in every cell, the body centre with probability `mix_fractions.BCC`, and
each face centre with probability `mix_fractions.FCC`. The `SC` entry is
the remainder and places no site of its own, which is what makes the
three a partition summing to 1. Expected site count is
`Nx*Ny*Nz * (1 + f_bcc + 3*f_fcc)`.

```json
"generator": {
  "lattice_size": "6x6x6",
  "lattice_type": "MIX",
  "mix_fractions": {"SC": 0.2, "BCC": 0.4, "FCC": 0.4},
  "max_functionality": 4
}
```

The point of mixing is more neighbour distances. A pure SC lattice offers
a single edge length; the mixture above offers four (0.5, 0.707, 0.866,
1.0 cell units), which smooths the distribution of strand end-to-end
distances. Three things to know before using it:

- **`MIX` at `{"SC": 1}` reproduces `SC` exactly**, down to node ids. That
  is *not* true at the other two corners. `MIX` connects by distance
  cutoff rather than by a fixed neighbour pattern, so at `{"BCC": 1}` the
  1.0 cutoff also admits the corner-corner shell and every node carries 14
  neighbours instead of BCC's 8 (18 instead of 12 at `{"FCC": 1}`). Use
  `lattice_type: "BCC"` or `"FCC"` when you want the canonical
  coordination.
- **Bond lengths spread.** A body centre and a face centre can land 0.5
  cells apart, half the SC spacing. DP is assigned independently of edge
  length, so strands of the same DP get built at bond lengths differing by
  up to 2x. Watch for FENE strain on the long edges.
- **The split is a coarse dial.** Per the strand-realism analysis the
  SC/BCC/FCC percentages are a weak, ill-conditioned knob: many splits fit
  a given target comparably well. Site jitter and a Gaussian-weighted edge
  rule move the strand statistics much more. Treat the fractions as
  SC-heavy for short strands shifting toward BCC/FCC as strand length
  grows, and verify by measurement rather than by tuning percentages.

Lowering `mix_cutoff` below 1.0 drops the corner-corner shell, which
disconnects the always-present corner sublattice from itself. 1.0 is the
default for that reason.

#### `topology.existing_files`

Provide either `gpickle_file` OR both `nodes_file` + `edges_file`.

| Key | Type | Default | Description |
|---|---|---|---|
| `nodes_file` | string \| null | `null` | Path to `.nodes` file |
| `edges_file` | string \| null | `null` | Path to `.edges` file |
| `gpickle_file` | string \| null | `null` | Path to NetworkX `.gpickle` file |

##### `.nodes` file format

Whitespace-separated `NodeID X Y Z Degree`, with `#` starting a comment.
An optional `# BOX Lx Ly Lz` header records the periodic cell in lattice
units:

```
# BOX 6 6 6
# NodeID X Y Z Degree
0 0.000000 0.000000 0.000000 3
1 1.000000 0.000000 0.000000 3
```

The header is optional and files without it load exactly as before, but
**write it for any lattice whose sites are not integer-spaced.** Without
it topon estimates the cell as `max - min + 1` over the coordinates,
which is exact for SC but overshoots BCC, FCC and Diamond because their
basis sites sit at fractional offsets and never reach the cell edge. A
4x4x4 BCC or FCC is estimated at 4.5, and since that value drives every
minimum-image calculation, about a third of BCC edges (a quarter of FCC)
get built at twice their true bond length. Anything topon generates
records the header for you; the caveat applies to hand-written or
externally-produced files.

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

##### Choosing pairs on a conformation instead of on crosslink distance

By default candidates are ranked by the distance between their crosslinks,
which is a property of the network rather than of the chains. Two chains can
be nearest neighbours by crosslink and never come near each other, and a kink
placed there aims one chain at a partner that is not present.

`select_entanglements` accepts a `chain_paths` argument — a mapping of
`frozenset((u, v))` to that chain's bead path — and ranks candidates by how
much of the two chains actually lies alongside. The assignment stage does not
draw the conformation itself; the caller supplies one, in the same units as
the node positions. `tests/workflows/entangle_by_proximity.py` shows the
sequence: draw a provisional conformation with no entanglements, rank on it,
select, then draw the final one with the kinks.

`topon/conformation/paths.py` provides `bridging_walk` for that first pass —
a random walk of fixed bond length that closes exactly on its far junction.
A straight chain will not do, since it lies on its chord and so cannot say
anything about which chains meet.

Measured on a 354-chain network, 281 candidates, 0.20 entanglements per chain:

| ranking | median proximity of chosen pairs | chosen pairs whose chains never touch |
|---|---|---|
| crosslink distance (the default) | 39 | 8 of 33 |
| on a conformation | **156** | **0 of 35** |

The pool median is 50, so the default ranking is slightly worse than choosing
at random.

##### Drawing a path around what is already there

`topon/conformation/paths.py` also carries the pieces used to place a chain
into an occupied box. None of it is wired into `Pipeline`; it is called from
the workflow scripts under `tests/workflows/`.

| name | does |
|---|---|
| `Clearance(points, box, radius)` | the beads already present, as a minimum-image nearest-neighbour query. `near`, `worst`, `ok` |
| `bridging_walk(..., avoid=)` | the same fixed-bond random walk, keeping each step clear of `avoid` |
| `loop_around(target, i, radius, n_pts, phase, avoid, span)` | waypoints encircling a strand. `span` is turns, so 0.5 is a hook and 2.0 is two turns |
| `taut_leg(start, end, n_bonds, bond, avoid, placed)` | a deterministic leg: exact bonds, lands on its end, avoids `avoid`, its own earlier beads, and `placed` |
| `route_through(start, end, waypoints, n_bonds, bond, avoid)` | visits every waypoint in order using `taut_leg` for each leg |

Two things about this matter more than they look.

**Draw around, not through.** A path placed with no regard for the beads
already there lands on top of them. Measured on a relaxed melt at density
0.85, routing one chain took the closest pair in the system from 0.502 σ to
0.195 and put 153 beads inside 0.5 σ where there had been none. At 0.195 σ the
WCA energy is of order 10⁵ kT, so the next minimisation does not relax that
contact, it shoves — hard enough to drag chains through each other, which
rewrites whatever topology was just built. `Clearance` is what avoids making
the overlap in the first place; through a real relaxed melt it takes the
tightest contact a routed path makes from 0.081 σ to 0.822.

**Spend the slack deliberately.** A chain carries far more contour than its
route needs, 77 σ for a route of about 21 in one measured case. A random walk
disposes of the rest by wandering, and the wandering crosses the target again
on its own account, so the entanglement count stops being a property of the
design: the same pair, same site, same winding, drawn on three seeds, measured
4, 7 and 0. `route_through` spends it deterministically instead, which is what
makes a requested count repeatable.

`walk_through` (random legs) and `route_through` (deterministic legs) take the
same arguments and give the same guarantees about bonds and junctions. Use the
first when a melt-like conformation is wanted and the second when a specific
topology is.

##### Choosing which neighbour shell to entangle

`shell_weights` biases the draw toward particular neighbour shells. Shells are
numbered from 1, closest first, and are read off the lattice rather than
assumed: the closest approach between two strands takes a handful of discrete
values (0.20, 0.35, 0.41, 0.50 lattice units on a mixed SC/BCC/FCC network),
and those bands are what "first neighbour" and "second neighbour" mean.

```json
"entanglements": {
  "enabled": true,
  "avg_crosslinks_per_chain": 2.0,
  "shell_weights": { "1": 0.7, "2": 0.3 }
}
```

Naming a shell restricts the draw to the shells named and weights them in
proportion. Omitting the key entirely, which is the default, draws from every
shell equally and is the behaviour of every earlier version. It multiplies
into `placement_bias_kind` rather than replacing it, so spatial and shell
biasing compose.

**Only the first shell reliably produces an entanglement.** Measured with a
primitive-path analysis, each pair checked on its own after the full
three-stage protocol:

| shell | gap | realised as asked |
|---|---|---|
| 1 | 12.3 σ | 5 of 7 |
| 2 | 21.4 σ | 2 of 16 |
| 3 | 24.7 σ | 0 of 16 |

The reason is the pair's gap divided by the chain's chord: 0.29 in the first
shell, 0.50 in the second. A chain has to spend contour reaching its partner,
and past about a third of a chord it runs out. That ratio is scale-free and
so a property of the lattice: it is 0.29 / 0.50 / 0.58 for every SC/BCC/FCC
mixture whatever the fractions, unchanged by box size, and worse for the pure
lattices (FCC 0.71, BCC 0.82). Weighting the outer shells up is allowed and
will not give you more entanglements there.

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
> **`gradient` is broken:** it ignores the requested composition and emits a
> hard 50:50 split (ask for A=0.1 and you still get A=0.50). For two monomers
> at equal fractions and even DP it is byte-identical to `block`. See
> [JOURNAL](JOURNAL.md).

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
