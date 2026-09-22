# CausFate
CausFate is a computational framework that leverages causal inference to dissect cellular state transitions and to identify cell fate-determining features. 
- It reconstructs model-based cell-state relationships without requiring densely sampled intermediate cells.
- It is applicable to single-cell omics, bulk omics and microarray data.
- It prioritizes candidate fate determinants complementary to differential analysis.

## Requirements
```
R-4.3.3, Matrix-1.6-5, Seurat-5.0.3, bnlearn-4.9.4, doParallel-1.0.17, dplyr-1.1.4, foreach-1.5.2, MASS-7.3-60.0.1, corpcor-1.6.10, igraph-2.0.3, graph-1.80.0, magrittr-2.0.3, matrixStats-1.3.0, infotheo-1.2.0.1, parallel-4.3.3, Rgraphviz-2.46.0, rlang-1.1.3, tibble-3.2.1, tidyr-1.3.1, tidyselect-1.2.1, space-0.1-1.1, statmod-1.5.0, SeuratObject-5.3.0, Rcpp-1.1.1, transport-0.15-4, reticulate-1.45.0
```

Optional Python dependency for GRN-based perturbation: Python 3 with
[CellOracle](https://github.com/morris-lab/CellOracle).

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

We also provide a [CausFate simulation example](benchmark/simulation_causfate.R) using predefined reference networks with 3–9 nodes (disconnected, sparse, tree and dense). See [benchmark instructions](benchmark/README.md).

## Perturbation and effect metrics

CausFate supports both single-feature (GRN-free) and GRN-based perturbation.
For GRN-free perturbation of single-cell data:

```r
perturbRes <- PerturbResult(
  net_struc = dag_struc,
  data = HMR,
  meta = meta,
  index = index,
  n_sample = 10,
  mode = "single_cell",             # Input data mode
  perturbation = "single_feature",  # GRN-free perturbation
  perturb_ratio = 0,                # Zeroing (default)
  ncores = 1
)
```

Other GRN-free perturbation types can be selected as follows:

```r
PerturbResult(..., perturb_ratio = 0.5)  # Knockdown: any value between 0 and 1
PerturbResult(..., perturb_ratio = 2)    # Knockup: any value greater than 1
PerturbResult(..., perturb_ratio = "deletion") # Feature deletion
```

For CellOracle-based GRN perturbation:

```r
perturbRes <- PerturbResult(
  net_struc = dag_struc,
  data = seurat_object,
  mode = "single_cell",               # Input data mode
  perturbation = "GRN",               # GRN-based perturbation
  save_dir = "celloracle_results",
  base_GRN = base_GRN,
  group.by = "celltype",
  genes = NULL,                       # NULL perturbs all available regulators
  n_jobs = 1
)
```

The effect matrix can then be calculated using:

```r
effMat <- EffectMatrix(perturbRes, dist_metric = "diff_mean")
# Other metrics: "W1", "W2", "energy", "mmd"
```

## Overview of CausFate
![CausFate workflow](CausFate%20workflow.png)
