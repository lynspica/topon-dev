# MARTINI resilin — 2-step workflow

End-to-end MARTINI 3 protein-network build for the **natpro** resilin sequence at user-chosen system size + crosslink conversion + water content. Designed so you inspect the topology sweep before committing to a multi-hour LAMMPS run.

System parameters baked in:
- **Sequence**: `GGRPSDSYGAPGGGN` × 18 = 270 residues per chain (the vendored MARTINI 3 nat_pro template covers exactly this; 504 beads per chain)
- **segs_per_block**: 2 (TYR-eligible "Y nodes" every other backbone bead — matches topro's V42 convention)
- **System sizes**: 50 chains and 100 chains
- **Conversion targets**: pre-gel snapshots at 10 %, 15 %, 20 %, **25 %** (the user's target), 30 %, 40 %, 50 %, plus gel point + post-gel if reached
- **seed**: 42 (deterministic)

## Files

| Path | Purpose |
|---|---|
| `01_topology_sweep.py` | Step 1: generate BFM topologies + snapshots + CSV + plot |
| `02_build_martini.py`  | Step 2: pick a checkpoint, build MARTINI LAMMPS system |
| `topologies/`          | Output of step 1 — one JSON per snapshot, plus `summary.csv` and `summary.png` |
| `runs/`                | Output of step 2 — one folder per build |

## Step 1 — topology sweep

```bash
python examples/workflows/martini_resilin/01_topology_sweep.py
```

Wall: ~5 seconds. Produces per-snapshot topology JSONs in `topologies/` and a `summary.csv` listing every snapshot's conversion, crosslink count, and absolute path. The matching `summary.png` shows crosslink growth vs. conversion for both system sizes with a vertical line at the 25 % target, plus a bar chart of gel-point conversions.

Look at the CSV / plot, decide which `(system, snapshot)` you want.

### What the columns mean

| Column | Meaning |
|---|---|
| `conversion` | fraction of TYR sites reacted = `n_crosslinks * 2 / total_TYR` |
| `active_frac` | **reacted TYR in the largest connected component / total TYR** — the load-carrying / "active" crosslinker fraction |
| `n_lcc_chains` | how many chains have at least one TYR in the LCC |
| `gel_conv` | conversion at which all chains first joined one cluster, or `no_gel` |

`active_frac` is the metric to target if you want a fixed *load-bearing* crosslinker density. Below the gel point it's smaller than `conversion` (some crosslinks are stuck in tiny isolated sub-clusters and don't transmit load); near the gel point the two converge.

### Example summary output

```
system           snapshot                conv    xlinks   active   chains_in_LCC
natpro_50chain   pre_gel_conv0250        0.251   113/450  0.247    40/50      <- ~25% active
natpro_50chain   pre_gel_conv0500        0.500   225/450  0.500    48/50
natpro_50chain   no_gel                  0.793   357/450  0.784    49/50      (gel not reached)
natpro_100chain  pre_gel_conv0250        0.250   225/900  0.243    87/100     <- ~25% active
natpro_100chain  gel_point               0.711   640/900  0.711    100/100    (gel)
natpro_100chain  post_gel_10             0.788   709/900  0.788    100/100
```

At these system sizes the standard topological **gel point sits around 71 %** conversion (100 chains) or **isn't reached** (50 chains). The `pre_gel_conv0250` snapshot gives a sub-percolated load-carrying network with ~80–87 % of chains already in the LCC — a meaningful working point even though it's below gel.

## Step 2 — build a MARTINI LAMMPS system

```bash
python examples/workflows/martini_resilin/02_build_martini.py \
    --system   natpro_100chain \
    --snapshot pre_gel_conv0250 \
    --water-density 4.0 \
    --water-bead    W
```

| Flag | Default | Notes |
|---|---|---|
| `--system` | (required) | `natpro_50chain` or `natpro_100chain` |
| `--snapshot` | (required) | One of the snapshot labels from step 1's CSV (`pre_gel_conv0250`, `gel_point`, etc.) |
| `--water-density` | `4.0` | W beads per nm³. **4.0 ≈ 50 % of bulk water** when using the regular W bead (4 H₂O/bead, bulk ≈ 8.25 W/nm³). Matches the v42 "W04nm3" medium-water reference cell. |
| `--water-bead` | `W` | `W` (4 H₂O/bead), `SW` (3), or `TW` (2). "Regular size" = W. |

The script re-runs BFM with the same seed + parameters as step 1, then drives chemistry → water packing → LAMMPS writing through `topon.protein_network.workflow.run_protein_network`. Output lands at `runs/<system>__<snapshot>__W<density>/`.

## Step 3 — run LAMMPS (optional)

The build script prints the exact commands. For a 100-chain × 0.25-conversion × W=4 system that's roughly:

```bash
cd examples/workflows/martini_resilin/runs/natpro_100chain__pre_gel_conv0250__W4/relaxation
lmp -in protein_network_stage1.in   # ~10–30 s soft-push overlap removal
lmp -in protein_network_stage2.in   # ~10–30 min LJ-epsilon ramp (system-size dependent)
lmp -in protein_network_stage3.in   # ~1–5 min tight CG min + brief NVT/NPT @ 310 K
```

After stage 3, `system_equilibrated.data` sits in the parent folder ready for production.

## Sanity defaults — keep the two scripts in sync

Both scripts depend on the same BFM knobs (`target_packing`, `equil_steps`, `crosslink_method`, `max_crosslink_distance_ang`, `lattice_scale_ang`, `pre_gel_conversions`). The values are duplicated literally in both files — if you tune one for a different chain length or sequence, update the other to match. The `seed=42` is what guarantees the topology re-generated in step 2 matches the one indexed in step 1.

## Reference

- `topon.protein_network.bfm.generate_topology` — the BFM Monte Carlo + crosslinking generator
- `topon.protein_network.workflow.run_protein_network` — end-to-end MARTINI workflow used by step 2
- Historical runs the natpro template was validated against: `tests/output/v42/natpro__*` (8 chains × 270 residues)
