# Benchmark workflows

## CausFate simulation

After installing CausFate, run from the repository root:

```sh
Rscript benchmark/simulation_causfate.R
# Root-free reconstruction
CAUSFATE_ROOT=none Rscript benchmark/simulation_causfate.R
```

The example covers 3–9 nodes and disconnected, sparse, tree and dense
reference networks.

Default settings: 20 resamples, resampling fraction 0.2, skeleton
thresholds 0.05–0.15 in steps of 0.01, no edge-count filtering, pruning
threshold 0.8, and at most two parents per node.

## Other methods

Source `benchmark/benchmark_functions.R` to access the comparator functions.
The Python scripts are `notears_root.py`, `golem_root.py` and
`dag_gnn_root.py`. Each supports `--help` and `--root none`.

### Requirements

- NOTEARS: `numpy`, `pandas`, `scipy`, `networkx`.
- GOLEM and DAG-GNN: the above packages plus `torch` and `gcastle`.
- PC-stable: `bnlearn`.
- FGES: `causalDisco`, `caugi`, `rJava`, and Tetrad.

### PC-stable output

`pcExtensionOrRaw()` uses `bnlearn::cextend()` to obtain a DAG. When a
consistent extension cannot be found, it returns the raw partially directed
graph and records the reason in `pc_extension_error`. The
`pc_output_type` attribute identifies the selected output.

For directional edge counting, use `nrow(bnlearn::arcs(net))`: an unresolved
adjacency is represented by two records, A→B and B→A.
