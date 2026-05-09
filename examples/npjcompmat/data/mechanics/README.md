# Mechanics network files

One subfolder per case in `data/csv/mechanics.csv`. Each subfolder contains the
three network files that define the case's topology:

| File | Format | Notes |
|---|---|---|
| `network_N6x6x6_trial*.edges` | plain text | Edge list (source, target) per line. |
| `network_N6x6x6_trial*.nodes` | plain text | Node list with attributes. |
| `network_N6x6x6_trial*.gpickle` | Python pickle | `networkx` graph — redundant with the two above, but convenient to load. |

The trial index in the filename varies per case (it's the seed used during
generation). The `graph_file` column in `data/csv/mechanics.csv` holds the exact
relative path for each row (from the demo root), so the mapping is unambiguous:

```python
import pandas as pd, pickle
df = pd.read_csv('data/csv/mechanics.csv')
row = df.iloc[0]
with open(row['graph_file'], 'rb') as f:    # e.g. data/mechanics/10_0_26_75_43_53_9/network_N6x6x6_trial6.gpickle
    G = pickle.load(f)                      # networkx Graph
```

The `.edges` + `.nodes` pair is the portable, language-agnostic representation;
use `.gpickle` only if you trust this source (standard pickle security caveat).

327 cases total.
