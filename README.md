# CausFate
CausFate is a computational framework that leverages causal inference to dissect cellular state transitions and to identify cell fate-determining features. 
- It can infer causality in the absence of continuity in cell states.
- It is applicable to single-cell omics, bulk omics and microarray data.
- It outperforms differential analyses in enriching crucial genes.

## Requirements
```
R-4.3.3, Matrix-1.6-5, Seurat-5.0.3, SeuratObject-5.3.0, bnlearn-4.9.4, doParallel-1.0.17, dplyr-1.1.4, foreach-1.5.2, MASS-7.3-60.0.1, corpcor-1.6.10, igraph-2.0.3, graph-1.80.0, magrittr-2.0.3, matrixStats-1.3.0, infotheo-1.2.0.1, parallel-4.3.3, reticulate-1.45.0, Rgraphviz-2.46.0, rlang-1.1.3, str2str-1.0.0, tibble-3.2.1, tidyr-1.3.1, tidyselect-1.2.1, space-0.1-1.1, statmod-1.5.0, transport-0.15-4
Python CellOracle (required only for GRN mode)
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

We also provide benchmarking scripts for simulated datasets generated from pre-defined reference causal networks with varying numbers of nodes (3–9) and structural densities (discrete, sparse, tree and dense): https://github.com/Lan-lab/CausFate/tree/main/benchmark.

## Overview of CausFate
<img width="3090" height="4084" alt="CausFate Overview" src="https://github.com/user-attachments/assets/c8e4dc64-0481-479f-9569-5af100a97a3c" />
