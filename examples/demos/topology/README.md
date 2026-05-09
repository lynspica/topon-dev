# Topology demos

Topology generation in isolation — build the network graph without driving it through chemistry / conformation / writers. Useful for inspecting graph properties (degree distribution, connectivity, cycles) before committing to a full pipeline run, and for benchmarking.

Two sub-folders:

- [`end_linking/`](end_linking/) — SC-lattice end-linking topology generation. Two paths to the same output: the C generator (`generator.exe`) and the pure-Python port (`topon.topology.generator_python`). Speed comparison in [`end_linking/speed_logs/benchmark.md`](end_linking/speed_logs/benchmark.md).
- [`bfm/`](bfm/) — BFM (Bond-Fluctuation Model) topology used by topro for protein-network workflows.

The topologies produced here can be fed back into the polymer or protein chemistry demos by switching `topology.source` to `"load"` in any of those configs.
