# CausFate
CausFate is a computational framework that leverages causal inference to dissect cellular state transitions and to identify cell fate-determining features. 
- It reconstructs model-based cell-state relationships without requiring densely sampled intermediate cells.
- It is applicable to single-cell omics, bulk omics and microarray data.
- It prioritizes candidate fate determinants complementary to differential analysis.

## Requirements
```
R-4.3.3, Matrix-1.6-5, Seurat-5.0.3, bnlearn-4.9.4, doParallel-1.0.17, dplyr-1.1.4, foreach-1.5.2, MASS-7.3-60.0.1, corpcor-1.6.10, igraph-2.0.3, graph-1.80.0, magrittr-2.0.3, matrixStats-1.3.0, infotheo-1.2.0.1, parallel-4.3.3, Rgraphviz-2.46.0, rlang-1.1.3, tibble-3.2.1, tidyr-1.3.1, tidyselect-1.2.1, space-0.1-1.1, statmod-1.5.0
```

## Installation
First, install `tidyverse` and the required Bioconductor dependencies:
```
install.packages("tidyverse")
install.packages("BiocManager")
BiocManager::install(c("graph", "Rgraphviz"))
```
Then install the `space` package:
```
install.packages("space")
```
If installation of `space` fails, install it from the CRAN Archive instead:
```
install.packages("https://cran.r-project.org/src/contrib/Archive/space/space_0.1-1.1.tar.gz",
repos = NULL, type = "source")
```
Finally, install `CausFate`:
```
install.packages("devtools")
devtools::install_github("Lan-lab/CausFate")
```


## Tutorial
Here we provide demos for reconstructing causal cell-state networks and ranking potential fate-determining features using both bulk and single-cell datasets: https://github.com/Lan-lab/CausFate/tree/main/tutorials.

We also provide a [CausFate simulation example](benchmark/simulation_causfate.R) using predefined reference networks with 3–9 nodes (disconnected, sparse, tree and dense). The historical filename key `discrete` denotes the disconnected, edge-free setting. Additional comparator helpers and Python wrappers are included. See [benchmark instructions](benchmark/README.md).

## Perturbation and effect metrics

The updated `PerturbResult()` defaults to zeroing (`deletion = FALSE`,
`perturb_ratio = 0`). Ratios strictly between 0 and 1 scale values down
(knockdown), ratios above 1 scale values up (knockup), and 1 leaves values
unchanged. These operations scale the input values; they do not themselves
establish biological loss or gain of function. Use `deletion = TRUE` to remove
the feature before SEM refitting.

This software update is distinct from the reported analyses: the main HMR,
HHATAC and zebrafish single-feature analyses used permutation, whereas the
HMR hidden-node analysis and bulk/microarray analyses used feature deletion.
Graded scaling was not used in those analyses. The notebooks explain their
chosen operation; the new zeroing example is not the original HMR ranking.

```r
effMat <- EffectMatrix(perturbRes, dist_metric = "diff_mean")
score <- diffScore(effMat, edgeSet, abs = TRUE)
# Other metrics: "W1", "W2", "energy", "mmd"
```

`diff_mean` preserves the legacy signed sum of reference-minus-perturbed
coefficient differences across runs. Absolute values are applied by
`diffScore(abs = TRUE)`, not by `EffectMatrix()`. For equal repeat counts,
the signed sum and signed mean shift differ by a common scale factor.
The existing quadratic `mmd` implementation returns MMD squared with a fixed
edge-specific Gaussian-kernel bandwidth. The legacy `mode = "mean"` and other
`mode` calls remain supported; conflicting `mode` and `dist_metric` values
produce an error. Wasserstein metrics require the R package `transport`.

`combineDAGsmpl(..., model_averaging = "joint")` is the default. The
`model_averaging = "two-tier"` option is retained for method comparisons.
GRN-based perturbation additionally requires a configured Python environment
with CellOracle, accessed via `reticulate`.

## Overview of CausFate
![CausFate workflow](CausFate%20workflow.png)
