# Topon Configuration Reference

All pipeline behaviour is controlled by a single JSON file. Use `topon init` to generate a starter file, then edit as needed.

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

---

## `study`

| Key | Type | Default | Description |
|---|---|---|---|
| `name` | string | `"my_network"` | Study name; used as sub-directory under `output_dir` |
| `output_dir` | string | `"./output"` | Root output directory |

```json
"study": {
  "name": "pdms_6x6x6",
  "output_dir": "./runs"
}
```

---

## `topology`

| Key | Type | Default | Description |
|---|---|---|---|
| `source` | `"generate"` \| `"load"` | `"load"` | Generate a new topology or load an existing one |
| `generator` | object | — | Settings for the C / Python generator (used when `source="generate"`) |
| `existing_files` | object | — | File paths (used when `source="load"`) |

### `topology.generator`

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

**Degree distribution format:**

```
"0:15,1:30,e:371"
```

- `d:N` → require exactly N nodes of degree d
- `e:N` → require exactly N edges total
- Omitted degrees are unconstrained

**Examples:**

```json
"generator": {
  "exe_path": null,
  "lattice_size": "6x6x6",
  "lattice_type": "SC",
  "max_functionality": 4,
  "degree_distribution": "0:13,1:25"
}
```

```json
"generator": {
  "lattice_type": "BCC",
  "lattice_size": "4x4x4",
  "max_functionality": 6,
  "degree_distribution": "e:248"
}
```

### `topology.existing_files`

| Key | Type | Default | Description |
|---|---|---|---|
| `nodes_file` | string \| null | `null` | Path to `.nodes` file |
| `edges_file` | string \| null | `null` | Path to `.edges` file |
| `gpickle_file` | string \| null | `null` | Path to `.gpickle` (NetworkX) file |

Provide either `gpickle_file` OR both `nodes_file` + `edges_file`.

```json
"topology": {
  "source": "load",
  "existing_files": {
    "nodes_file": "output/network.nodes",
    "edges_file": "output/network.edges"
  }
}
```

---

## `assignment`

Controls how graph attributes (types, DP, defects, entanglements, grafts, copolymers) are assigned to nodes and edges before chemistry is built.

### `assignment.node_types`

| Key | Type | Default | Description |
|---|---|---|---|
| `method` | `"degree"` \| `"positional"` \| `"random"` \| `"explicit"` | `"degree"` | Assignment method |

**method = `"degree"`** — assign type by node degree:
```json
"node_types": {
  "method": "degree",
  "degree": { "mapping": {"1": "end", "2": "A", "3": "A", "4": "A"} }
}
```

**method = `"positional"`** — assign type by layer along an axis:
```json
"node_types": {
  "method": "positional",
  "positional": { "dimension": "z", "num_layers": 2, "layer_types": ["A", "B"] }
}
```

**method = `"random"`** — assign type by weighted random:
```json
"node_types": {
  "method": "random",
  "random": { "type_ratios": {"A": 70, "B": 30} }
}
```

**method = `"explicit"`** — assign per-node-ID:
```json
"node_types": {
  "method": "explicit",
  "explicit": {"0": "POSS", "1": "Si"}
}
```

---

### `assignment.edge_types`

| Key | Type | Default | Description |
|---|---|---|---|
| `method` | `"uniform"` \| `"random"` \| `"composite"` | `"uniform"` | Assignment method |

**method = `"uniform"`** — all edges get the same type:
```json
"edge_types": {
  "method": "uniform",
  "uniform": { "type": "A" }
}
```

**method = `"random"`** — edges assigned by weighted random:
```json
"edge_types": {
  "method": "random",
  "random": { "type_ratios": {"A": 60, "B": 40} }
}
```

**method = `"composite"`** — assign by spatial layer:
```json
"edge_types": {
  "method": "composite",
  "composite": { "dimension": "z", "num_layers": 3, "layer_types": ["A", "B", "A"] }
}
```

---

### `assignment.dp_distribution`

Controls the degree of polymerization (chain length) per edge.

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

---

### `assignment.defects`

Inject primary loops (parallel edges between the same two nodes).

```json
"defects": {
  "primary_loops": {
    "enabled": true,
    "target": 10,
    "target_type": "count"
  }
}
```

| Key | Values | Description |
|---|---|---|
| `enabled` | bool | Enable defect injection |
| `target` | int | Count or percentage value |
| `target_type` | `"count"` \| `"percentage"` | How `target` is interpreted |

---

### `assignment.entanglements`

Select pairs of topologically entangled chains (physical knots).

```json
"entanglements": {
  "enabled": true,
  "target": 5,
  "target_type": "count",
  "kink_params": {
    "overshoot": 0.2,
    "z_amp": 0.5,
    "sigma": 0.15
  }
}
```

Alternatively, use distribution mode to specify average crosslinks per chain:
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

---

### `assignment.grafts`

Attach side chains to backbone beads at a given density.

```json
"grafts": {
  "enabled": true,
  "per_edge_type": {
    "A": {
      "graft_density": 0.05,
      "side_chain_monomer": "PDMS",
      "side_chain_dp": 5
    }
  }
}
```

| Key | Default | Description |
|---|---|---|
| `graft_density` | `0.5` | Probability of attachment at each backbone bead (0–1) |
| `side_chain_monomer` | `"PDMS"` | Monomer name for side-chain beads |
| `side_chain_dp` | `5` | Number of beads per side chain |

---

### `assignment.copolymer`

Generate per-bead monomer sequences along each chain.

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

| `arrangement` | Description |
|---|---|
| `"block"` | All of one monomer, then all of the next |
| `"alternating"` | Strict A-B-A-B… pattern |
| `"random"` | Randomly sampled at each position by fraction |
| `"gradient"` | Linear composition gradient along the chain |

---

## `chemistry`

Maps abstract graph types → concrete molecules and controls the force field model.

| Key | Type | Default | Description |
|---|---|---|---|
| `model_type` | `"coarse_grained"` \| `"atomistic"` | `"coarse_grained"` | Force field resolution |
| `target_density` | float | `0.9` | Target density in g/cm³ |

### `chemistry.node_type_map`

Maps node type labels → molecule definitions.

```json
"node_type_map": {
  "end": { "molecule": "[Si](C)(C)C", "is_end_cap": true },
  "A":   { "molecule": "Si",          "is_end_cap": false },
  "B":   { "molecule": "POSS",        "is_end_cap": false }
}
```

Built-in molecule names: `"Si"`, `"POSS"` (Si8O12 cage), `"POSS_AM0270"` (AM0270 aminopropyl POSS). Any SMILES string is also accepted.

### `chemistry.edge_type_map`

Maps edge type labels → monomer names.

```json
"edge_type_map": {
  "A": { "monomer": "PDMS" },
  "B": { "monomer": "FPDMS" }
}
```

### `chemistry.monomers`

Monomer library. Built-in defaults:

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

### `chemistry.connection`

```json
"connection": {
  "auto_bridge": true,
  "default_bridge_atom": "O"
}
```

`auto_bridge`: when the chain head and node atom are the same element (e.g., both Si), automatically inserts a bridge atom between them. Set `false` to always use direct bonds.

---

## `output`

| Key | Default | Description |
|---|---|---|
| `lammps_data` | `true` | Write LAMMPS data file |
| `lammps_inputs` | `true` | Write LAMMPS input scripts |
| `visualization` | `true` | Write HTML visualization |
| `analysis_report` | `true` | Write analysis report |
| `save_attributed_graph` | `true` | Save attributed graph as `.gpickle` |

---

## Complete Example: CG Network with Entanglements

```json
{
  "study": {
    "name": "cg_entangled",
    "output_dir": "./runs"
  },
  "topology": {
    "source": "generate",
    "generator": {
      "exe_path": null,
      "lattice_size": "6x6x6",
      "lattice_type": "SC",
      "max_functionality": 4,
      "degree_distribution": "0:13,1:25"
    }
  },
  "assignment": {
    "node_types": {
      "method": "degree",
      "degree": { "mapping": {"1": "end", "2": "A", "3": "A", "4": "A"} }
    },
    "edge_types": {
      "method": "uniform",
      "uniform": { "type": "A" }
    },
    "dp_distribution": {
      "default": { "mean": 25, "pdi": 1.0 }
    },
    "entanglements": {
      "enabled": true,
      "target": 5,
      "target_type": "count",
      "kink_params": { "overshoot": 0.2, "z_amp": 0.5, "sigma": 0.15 }
    }
  },
  "chemistry": {
    "model_type": "coarse_grained",
    "target_density": 0.9,
    "node_type_map": {
      "end": { "molecule": "Si", "is_end_cap": true },
      "A":   { "molecule": "Si", "is_end_cap": false }
    },
    "edge_type_map": { "A": { "monomer": "PDMS" } },
    "monomers": {
      "PDMS": { "smiles": "[Si](C)(C)O", "chain_head": "Si", "chain_tail": "O" }
    }
  }
}
```

## Complete Example: Atomistic POSS Network

```json
{
  "study": { "name": "atomistic_poss", "output_dir": "./runs" },
  "topology": {
    "source": "load",
    "existing_files": {
      "nodes_file": "output/network.nodes",
      "edges_file": "output/network.edges"
    }
  },
  "assignment": {
    "node_types": {
      "method": "degree",
      "degree": { "mapping": {"1": "end", "2": "A", "3": "A", "4": "POSS"} }
    },
    "edge_types": { "method": "uniform", "uniform": { "type": "A" } },
    "dp_distribution": { "default": { "mean": 10, "pdi": 1.0 } }
  },
  "chemistry": {
    "model_type": "atomistic",
    "target_density": 1.1,
    "node_type_map": {
      "end":  { "molecule": "[Si](C)(C)C", "is_end_cap": true },
      "A":    { "molecule": "Si",           "is_end_cap": false },
      "POSS": { "molecule": "POSS_AM0270",  "is_end_cap": false }
    },
    "edge_type_map": { "A": { "monomer": "PDMS" } },
    "monomers": {
      "PDMS": { "smiles": "[Si](C)(C)O", "chain_head": "Si", "chain_tail": "O" }
    }
  }
}
```
