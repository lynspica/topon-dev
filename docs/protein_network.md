# MARTINI 3 protein-network generator

`topon.protein_network` turns a residue sequence into a coarse-grained MARTINI 3
polymer-network LAMMPS input set. It is a peer of `topon.simbox` and
`topon.singlechain` -- it does **not** plug into the central `topon.pipeline`,
because the central pipeline starts from a lattice topology that does not
naturally describe a long protein chain with sparse crosslinks.

The generator mirrors the topro CHARMM-atomistic protein-network workflow but
emits MARTINI 3 coarse-grained beads (~28 beads per 15-residue resilin block,
versus ~250 atoms in the all-atom CHARMM equivalent).

## At a glance

```
  block_seq + n_repeats + n_chains
              v
    bfm.generate_topology   ----->  protein_network_topology.json  (snapshots)
              v
   builder.build_protein_system  + residues.py table
              v
        water.pack_water         (optional)
              v
   lammps_writer.write_lammps    (.data + .in.settings + .in.groups + .in)
```

Inputs the generator understands:

* a one-letter repeat block such as the resilin consensus `GGRPSDSYGAPGGGN`
* the 8 residues currently covered by the auto-extracted residue table:
  `GLY`, `ALA`, `ARG`, `PRO`, `SER`, `ASP`, `TYR`, `ASN`
* TYR positions are interpreted as crosslinker sites (dityrosine);
  the BFM topology stage stochastically links lattice-adjacent TYR pairs
* optional W water beads packed on a voxel grid around the protein

## CLI

```
python -m topon.protein_network generate \
    --block-seq GGRPSDSYGAPGGGN \
    --n-repeats 6 \
    --n-chains 4 \
    --equil-steps 5000 \
    --water-density 0 \
    --output runs/resilin_dry/ \
    --seed 42
```

Subcommands:

* `generate`  -- single run; writes data + settings + groups + input script
  + the topology JSON.
* `sweep`     -- repeat the run across a comma-separated list of water
  densities into `wXX/` subdirectories (mirroring the topro convention where
  XX is `density * 10` rounded).
* `topology`  -- run only the BFM lattice stage; emit the snapshot JSON.

Common flags (all subcommands):

| Flag | Default | Meaning |
|---|---|---|
| `--block-seq` | `GGRPSDSYGAPGGGN` | One-letter repeat block (resilin consensus). |
| `--n-repeats` | 6 | Repeats per chain. |
| `--n-chains` | 4 | Number of chains. |
| `--segs-per-block` | 2 | BFM segments per repeat (2 = end + Y; 3 = +1 NC node). |
| `--equil-steps` | 5000 | Monte-Carlo equilibration steps (0 = skip). |
| `--target-packing` | 0.45 | Volume fraction used to size the lattice. |
| `--min-intrachain-sep` | 2 | Minimum Y-index gap for intra-chain crosslinks. |
| `--lattice-scale-ang` | 4.7 | Angstroms per BFM lattice unit (~MARTINI sigma). |
| `--sc-jitter-ang` | 1.5 | Random offset of sidechain beads from BB (A). |
| `--snapshot-label` | `gel_point` | Which BFM snapshot to build from (`gel_point`, `post_gel_1`...). |
| `--seed` | 42 | RNG seed. |
| `--quiet` | off | Suppress per-stage status prints. |

## Python API

```python
from topon.protein_network.workflow import run_protein_network

paths = run_protein_network(
    block_seq="GGRPSDSYGAPGGGN",
    n_repeats=6,
    n_chains=4,
    output_dir="runs/resilin_dry",
    equil_steps=5_000,
    water_density_w_per_nm3=0.0,   # set ~10 for bulk MARTINI water
    seed=42,
)
# paths is a dict: {data, settings, groups, in, topology_json}
```

Lower-level entry points (each module is independently usable):

* `topon.protein_network.bfm.generate_topology(...)` -- BFM cubic-lattice topology
* `topon.protein_network.topology_io.{save_topology, load_topology, get_snapshot}`
  -- JSON round-trip; format is byte-compatible with topro's `topo_*.json`.
* `topon.protein_network.builder.build_protein_system(snapshot, sequence_3letter, library, ...)`
  -- chain chemistry from a snapshot.
* `topon.protein_network.water.pack_water(system, library, ...)` -- voxel-grid W packer.
* `topon.protein_network.lammps_writer.write_lammps(system, library, output_dir, ...)`
  -- LAMMPS file emitter.
* `topon.protein_network.martini_ff.MartiniLibrary.from_package_data()` -- loads
  the vendored pruned MARTINI 3 protein FF + W water bead definition.

## Output files

```
<output_dir>/
  protein_network_topology.json   # BFM snapshots (gel_point, post_gel_1, ...)
  protein_network.data            # LAMMPS data: Atoms (full), Bonds, Angles,
                                   #             Dihedrals, Impropers, Masses
  protein_network.in.settings     # explicit pair_coeff / bond_coeff / etc.
  protein_network.in.groups       # group protein molecule N,
                                   # group water molecule N definitions
  relaxation/
    protein_network_stage1.in     # soft-push overlap removal
    protein_network_stage2.in     # LJ epsilon ramp 0.001 -> 1.0 (nve/limit)
    protein_network_stage3.in     # tight CG min + brief NVT/NPT @ 310 K
                                   #   -> ../system_equilibrated.data
```

The 3-stage relaxation protocol mirrors topro's CHARMM stages, swapped to
MARTINI 3 (lj/cut/coul/cut, dielectric=15, no PPPM, T=310 K). Each stage reads
the previous stage's `write_data` output via the `system_after_soft.data` and
`system_ramped.data` handoff files. After stage 3, `system_equilibrated.data`
sits in the parent directory ready to feed into your supercomputer annealing
protocol.

### Running

Both serial and parallel execution work. Serial is the default; the only
parallel-flavored directive in the stages (`comm_modify mode single cutoff 14.0`
in stage 2) is a no-op on one processor.

Serial (one-machine workstation runs):
```
cd <output_dir>/relaxation/
lmp -in protein_network_stage1.in
lmp -in protein_network_stage2.in
lmp -in protein_network_stage3.in
```

Parallel (HPC):
```
mpirun -np <N> lmp -in protein_network_stage1.in
mpirun -np <N> lmp -in protein_network_stage2.in
mpirun -np <N> lmp -in protein_network_stage3.in
```

### Pre-generated demo

`tests/output/v33_protein_network_demo_16chain_w4/` contains a freshly
generated 16-chain x 6-repeat resilin system at 4 W/nm^3 hydration -- handy
for inspecting the data file format and the 3 relaxation scripts without
re-running the generator.

## How the residue table is built

`tools/extract_residues_from_itp.py` parses
`tests/_martini_extracted/Martini_Ahmet/itp_files/nat_pro.itp` (a
polyply-generated MARTINI 3 ITP for the 270-residue resilin chain) and writes
`topon/protein_network/residues.py`. Re-running the extractor after dropping a
new polyply ITP into the source directory regenerates the table without code
edits. The extractor also prunes the 16 MB master MARTINI 3 ITP down to just
the bead types referenced by the protein (~30 KB) and copies the W water
solvent ITP into the package data folder.

## Known approximations

This is a first-cut LAMMPS port of a force field designed for GROMACS. Five
things are approximated rather than reproduced exactly. None affect topology;
all affect energy:

1. **Reaction-field electrostatics**. GROMACS uses `coulombtype = reaction-field
   epsilon_r = 15`. We approximate with `pair_style lj/cut/coul/cut` plus
   `dielectric 15.0`. The RF correction term is omitted. For quantitative
   electrostatic agreement on charged systems, recompile LAMMPS with the
   USER-MISC `pair_style coul/diel` and edit the `.in` script.
2. **Restricted-bending angle (GROMACS funct=10)**. MARTINI 3's IDP backbone
   uses `U = 0.5 K (cos t - cos t0)^2 / sin^2 t`. We approximate with
   `angle_style cosine/squared` (`U = K (cos t - cos t0)^2`), omitting the
   `1/sin^2 t` factor.
3. **Proper dihedrals (GROMACS funct=9, multi-term)**. Each term in the
   GROMACS quartet is emitted as a separate LAMMPS `dihedral_style charmm`
   coefficient set with the same atom indices. This matches GROMACS energetics
   exactly.
4. **Constraints (`[ constraints ]`)**. GROMACS treats them as rigid
   distances (SHAKE/LINCS). We emit them as ULTRA-stiff harmonic bonds
   (`K = 1e6 kJ/mol/nm^2`, matching the reference's `#ifdef FLEXIBLE` block).
   For true rigidity, switch to `fix shake` -- the input script header notes
   how.
5. **Per-atom exclusions (`[ exclusions ]`)**. The TYR ring's
   beyond-1-2 exclusions are not emitted as LAMMPS-data-file rows; instead the
   input script uses `special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0` (1-2,
   1-3, and 1-4 are all excluded). This is more aggressive than the reference's
   `nrexcl=1` plus selective ring exclusions, but yields correct intra-ring
   behavior with no per-pair plumbing.

## Citations

Vendored MARTINI 3 force field tables:

* Souza, P C T et al., "Martini 3: a general purpose force field for
  coarse-grained molecular dynamics", *Nature Methods* 2021,
  10.1038/s41592-021-01098-3.
* Kroon, P C et al., "Martinize2 and Vermouth: Unified framework for
  topology generation", *eLife* 2024, 10.7554/elife.90627.2.
* Grunewald, F et al., "Polyply; a python suite for facilitating simulations
  of macromolecules and nanomaterials", *Nature Communications* 2022,
  10.1038/s41467-021-27627-4.

## What is NOT included (yet)

* **Folded-protein parameterizations**: this writer assumes IDP backbone
  chemistry. Folded MARTINI 3 proteins use `[ virtual_sites_2/3/4 ]` and
  optional Go contact maps -- neither is supported.
* **Elastic network restraints**: not emitted. Add manually if you need
  one for a folded chain.
* **Ions other than W water**: NaCl and friends are deferred. The extractor
  retains their bead types in the pruned FF for future use.
* **Simbox MARTINI support**: `topon/simbox/` remains DREIDING-only.
* **Sequence-residue coverage beyond the 8 amino acids in the resilin
  reference**: extending requires re-running the extractor against a polyply
  ITP that contains the new residues, then committing the updated
  `residues.py` and pruned FF data files.
