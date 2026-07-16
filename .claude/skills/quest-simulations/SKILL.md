---
name: quest-simulations
description: >
  How to run topon / Python / LAMMPS jobs on Northwestern's Quest HPC
  cluster — login, the p32566 allocation, file transfer, the conda env,
  SLURM short-partition templates, and the gotchas (mamba activate fails
  in batch; rsync absent in local git-bash). Use whenever a task involves
  generating datasets or running simulations on Quest rather than locally.
---

# Running simulations on Quest (Northwestern HPC)

> Owner-maintained. Paths under `pym5012` / project `p32566` are specific
> to this account — update them if the account or allocation changes.

## 1. Login & access

```bash
ssh pym5012@quest.it.northwestern.edu
```

Key-based auth is configured, so non-interactive ssh works (`ssh -o
BatchMode=yes ...` succeeds without a password). Agents can run remote
commands directly:

```bash
ssh -o BatchMode=yes pym5012@quest.it.northwestern.edu "squeue -u \$USER"
```

Requires Northwestern VPN / on-campus network. If ssh hangs, the user is
off-VPN — stop and tell them, don't retry in a loop.

## 2. Allocation & where we work

- **SLURM account / allocation:** `p32566`  (use `#SBATCH -A p32566`)
- **Project storage (large, shared):** `/gpfs/projects/p32566/`
- **Main working tree for topon dataset work:**
  `/gpfs/projects/p32566/Studies/npz/`
  - `topon/`            — the topon Python package (keep in sync with local)
  - `internal/`         — workflow scripts (dataset builders, etc.)
  - `generator_debug13` — the compiled C topology generator
  - `funcsweep_20k/`    — example study dir (SLURM scripts + output)
- **Home:** `/home/pym5012` (also reachable as `/gpfs/home/pym5012`)

## 3. Transferring files local → Quest

`rsync` is **not** available in the local Windows git-bash shell, so use
`scp`:

```bash
# single file
scp -o BatchMode=yes path/to/file.py \
    pym5012@quest.it.northwestern.edu:/gpfs/projects/p32566/Studies/npz/...

# directory
scp -r -o BatchMode=yes localdir/ \
    pym5012@quest.it.northwestern.edu:/gpfs/projects/p32566/Studies/npz/dest/
```

After editing any topon source locally, re-`scp` just the changed files;
Quest's copy does **not** auto-sync. When changing one module, also check
its imports are present on Quest (e.g. a builder that imports a sibling
`dataset_20k_builder.py` needs that file shipped too).

## 4. Python environment (IMPORTANT gotcha)

There is a dedicated conda env `topon` with numpy / networkx / pydantic:

```
/home/pym5012/.conda/envs/topon/bin/python
```

**Do not `mamba activate topon` (or `conda activate`) inside a SLURM /
non-interactive shell** — it fails with *"Run 'mamba init' ..."* because
the shell hasn't been initialized. Instead **call the env's python
binary directly**:

```bash
PY=/home/pym5012/.conda/envs/topon/bin/python
$PY -c "import topon, numpy, networkx, pydantic; print('ok')"
```

The base login shell loads no modules by default (`module list` is
empty); the cluster's bare `python` has **no** numpy. Always use `$PY`.

## 5. Running topon dataset builders

The builders need the C generator exe and `topon` on the path:

```bash
export TOPON_GENERATOR_EXE=/gpfs/projects/p32566/Studies/npz/generator_debug13
export PYTHONPATH=/gpfs/projects/p32566/Studies/npz
PY=/home/pym5012/.conda/envs/topon/bin/python

cd /gpfs/projects/p32566/Studies/npz/internal/workflows/topology_dataset
$PY dataset_funcsweep_20k_builder.py --out <dir> --workers 52 \
    --shard <i> --n_shards <n> --only_N <N>
```

- Builders are **resume-aware**: NPZs already on disk are detected and
  only missing IDs are retried. Safe to resubmit a killed job verbatim.
- **Shard** a big dataset across jobs with `--shard i --n_shards n`; each
  shard owns a disjoint ID slice and writes its own
  `metadata.shardNN.csv` + `build_log.shardNN.txt` into a shared `--out`.
- Concatenate per-shard CSVs when done:
  ```bash
  head -1 metadata.shard00.csv > metadata.csv
  for f in metadata.shard*.csv; do tail -n+2 "$f" >> metadata.csv; done
  ```

## 6. SLURM short-partition template (4 h, 1 node, 52 cores)

```bash
#!/bin/bash
#SBATCH -A p32566
#SBATCH -p short                 # 4 h max
#SBATCH -N 1
#SBATCH --ntasks-per-node=52
#SBATCH --mem-per-cpu=512M
#SBATCH -t 3:59:59
#SBATCH --job-name=myjob
#SBATCH --output=slurm-%x-%j.out

module purge all
PY=/home/pym5012/.conda/envs/topon/bin/python
export TOPON_GENERATOR_EXE=/gpfs/projects/p32566/Studies/npz/generator_debug13
export PYTHONPATH=/gpfs/projects/p32566/Studies/npz

cd /gpfs/projects/p32566/Studies/npz/internal/workflows/topology_dataset
$PY <builder>.py --workers 52 ...
```

Submit / monitor:

```bash
sbatch run.sh
squeue -u $USER
scancel <jobid>
```

Sizing notes (topology generation, observed): a 52-core node sustains
~200 networks/min for small/medium systems; large low-functionality
N=20-style lattices drop to ~25/min. Budget shards so each finishes
inside 4 h; if one hits the wall, just resubmit (resume handles it).

## 7. Running LAMMPS (MD: minimization / equilibration / UTS)

- **Binary:** `/home/pym5012/lammps-22Jul2025/src/lmp_mpi`
- **MPI module:** `mpi/openmpi-4.1.1-gcc.10.2.0`

```bash
#SBATCH -A p32566
#SBATCH -p short
#SBATCH -N 4
#SBATCH --ntasks-per-node=52
#SBATCH -t 3:59:59

module purge all
module load mpi/openmpi-4.1.1-gcc.10.2.0
mpirun -n 208 /home/pym5012/lammps-22Jul2025/src/lmp_mpi -in input.in
```

(`-n 208` = 4 nodes × 52 ranks.) For chained stages with dependencies
(e.g. anneal → extended-equil → UTS) use `sbatch --parsable` +
`--dependency=afterok:$JOBID` to serialize. topon's generated
`06_Equilibration/submit_all.sh` is the canonical example.

## 8. Quick checklist before submitting

1. `ssh ... "squeue -u \$USER"` works → on VPN, key OK.
2. `$PY -c "import topon"` succeeds on Quest → package + env in place.
3. `$TOPON_GENERATOR_EXE` exists and is executable (for builders).
4. Output dir is under `/gpfs/projects/p32566/...` (project quota, not home).
5. Job time ≤ 3:59:59 for `short`; shard or use a longer partition if not.
