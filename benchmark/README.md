# Simulation example and comparator helpers

From the package root, after installing the updated local package:

```sh
Rscript benchmark/simulation_causfate.R
CAUSFATE_ROOT=none Rscript benchmark/simulation_causfate.R
```

Alternatively, edit `root` to `"A"` or `NULL`. Results, cached candidate DAGs,
and plots are separated under `results/simulation_causfate/root_A/` and
`results/simulation_causfate/root_free/`. Existing results are reused by
default. Use a new `CAUSFATE_OUTPUT_DIR` when testing a different code or
parameter version; cache reuse does not validate provenance.

The example explicitly passes seed 42 to `BNLearning()`. Earlier wrappers
called `set.seed(42)` but omitted the function's seed argument, whose default
uses the current time. Their independent resampling runs are therefore not
guaranteed to reproduce exactly; use saved candidate ensembles for exact
aggregation comparisons.

The CausFate example uses 20 resamples, fraction 0.2, skeleton thresholds
0.05–0.15 in steps of 0.01, no candidate-DAG edge-count filtering, trim threshold
0.8, and at most two parents. `discrete` is a legacy filename key for
disconnected (edge-free) reference graphs. The example runs CausFate only.

## Additional methods

`benchmark_functions.R` contains classical, PC-stable/FGES and continuous-
optimization helpers. The Python wrappers `notears_root.py`, `golem_root.py`
and `dag_gnn_root.py` must be available in the supplied `script_dir`. All
accept `--help`; `--root none` disables the root constraint. NOTEARS requires
NumPy, pandas, SciPy and NetworkX. GOLEM and DAG-GNN additionally require
PyTorch and gCastle. PC-stable uses bnlearn; FGES uses causalDisco, caugi,
rJava and a configured Tetrad installation (see `configureTetrad()`).
These optional dependencies are not needed to run the CausFate-only example.
Keep versions, parameters, input matrices and seeds with results when using
these helpers. This folder is not a complete runner for every manuscript
analysis or a guarantee of numerical equality across software versions.

## PC-stable extension-first fallback

`pcExtensionOrRaw()` first attempts a consistent DAG extension with
`bnlearn::cextend()`. If no extension is possible, it returns the original
partially directed graph, records `pc_output_type = "raw_fallback"` and the
extension error, and emits a warning. Learning failures remain failures;
they are not replaced with empty graphs.

An extension may fail when the estimated directions and v-structures cannot
all be preserved in an acyclic orientation of the same skeleton. Undirected
edges alone do not imply failure. Saved outputs reproduced extension failure
for IntestD with and without root information and HHR without root information.

Plot and evaluate the selected output, retaining raw output in the `cpdag`
attribute for inspection. In the benchmark's directional-record convention,
each unresolved adjacency contributes both A→B and B→A records. Use
`bnlearn::arcs(net)` and its row count when implementing this convention, not
`bnlearn::narcs(net)`, which counts an undirected adjacency once. The raw graph
may display unresolved adjacencies without arrowheads; this does not imply
two established causal directions. Existing manuscript result files are not
changed by this wrapper update.
