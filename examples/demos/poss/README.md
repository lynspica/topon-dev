# POSS-junction networks

Atomistic network with POSS (Polyhedral Oligomeric Silsesquioxane, Si₈O₁₂ cage) at the four-functional junctions instead of bare Si. The provided config maps degree-4 nodes to `POSS_AM0270` (aminopropyl POSS, AM0270) and degree ≤ 3 nodes to a bare Si atom; chains are PDMS.

## Run

```bash
topon generate examples/demos/poss/config.json
```

Or via Python:

```bash
python examples/run_via_api.py demos/poss/config.json
```

## What the config encodes

```json
"node_types.degree.mapping": {
  "1": "end", "2": "A", "3": "A", "4": "POSS"
}
"node_type_map": {
  "end":  { "molecule": "[Si](C)(C)C",  "is_end_cap": true },
  "A":    { "molecule": "Si",            "is_end_cap": false },
  "POSS": { "molecule": "POSS_AM0270",   "is_end_cap": false }
}
```

`POSS_AM0270` is a built-in molecule string handled by `ChemistryBuilder._place_poss_am0270()` — places the Si₈O₁₂ cage with one aminopropyl arm at corner 0 and seven 2,4,4-trimethylpentyl (isooctyl, inert) arms at corners 1–7. Explicit hydrogens are placed (V30–V31).

## Output

LAMMPS data + input scripts written to `output_atomistic_poss/atomistic_poss/`. Atom counts scale steeply with POSS fraction — a 5×5×5 lattice with degree-4 nodes mapped to POSS easily reaches ~10 k atoms.

## Variants worth trying

- **POSS fraction sweep** — change `assignment.node_types.degree.mapping` so degrees 3 and 4 both map to POSS, or use `method: random` with a POSS ratio. See V32 changelog (POSS sweeps with seven compositions).
- **POSS-only network** — set `chemistry.target_density = 1.2` and replace all junctions with POSS.
- **simbox crosslink + POSS** — for a different POSS workflow (packed AM0270-POSS molecules with epoxy / amine crosslinking templates), see `topon simbox` in [`docs/USAGE.md`](../../../docs/USAGE.md) §3.4.
