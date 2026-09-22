#' Prepare CellOracle objects
#'
#' @param object A Seurat object.
#' @param base_GRN A base GRN data frame or file path.
#' @param save_dir Directory in which to save the Oracle and Links files.
#' @param group.by Cell-type column in the Seurat metadata.
#' @param assay Assay used by CellOracle.
#' @param counts_layer Raw-count layer.
#' @param embedding Seurat dimensional reduction.
#' @param features Genes used by CellOracle; `NULL` uses the variable features.
#' @param n_pca_dims Number of principal components used for KNN imputation.
#' @param k Number of nearest neighbours.
#' @param alpha Ridge coefficient used for GRN inference.
#' @param bagging_number Number of bagging estimators.
#' @param n_jobs Number of CPU cores used by CellOracle.
#' @param verbal Logical; whether to print progress messages.
#'
#' @return A list containing the Oracle and Links objects and their file paths.
#' @export
PrepareCellOracle <- function(
    object,
    base_GRN,
    save_dir,
    group.by = "celltype",
    assay = "RNA",
    counts_layer = "counts",
    embedding = "umap",
    features = NULL,
    n_pca_dims = 50,
    k = NULL,
    alpha = 10,
    bagging_number = 20,
    n_jobs = -1,
    verbal = TRUE
) {
  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object.")
  }
  
  if (is.null(base_GRN)) {
    stop("'base_GRN' is required when Oracle and Links are not provided.")
  }
  
  if (!group.by %in% colnames(object@meta.data)) {
    stop("Metadata column '", group.by, "' was not found.")
  }
  
  if (!embedding %in% names(object@reductions)) {
    stop("Reduction '", embedding, "' was not found.")
  }
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  counts <- SeuratObject::LayerData(
    object,
    assay = assay,
    layer = counts_layer
  )
  
  if (!inherits(counts, "matrix")) {
    counts <- as.matrix(counts)
  }
  
  if (is.null(features)) {
    features <- SeuratObject::VariableFeatures(object[[assay]])
    
    if (length(features) == 0) {
      features <- rownames(counts)
    }
  }
  
  features <- intersect(features, rownames(counts))
  counts <- counts[features, , drop = FALSE]
  
  meta <- as.character(object@meta.data[colnames(counts), group.by])
  embedding_mat <- SeuratObject::Embeddings(
    object[[embedding]]
  )[colnames(counts), , drop = FALSE]
  
  celloracle <- reticulate::import(
    "celloracle",
    convert = FALSE
  )
  
  reticulate::py_run_string("
import anndata as ad
import pandas as pd
import numpy as np

def _ciber_make_anndata(x, cells, genes, groups, embedding, groupby):
    adata = ad.AnnData(X=x)
    adata.obs_names = cells
    adata.var_names = genes
    adata.obs[groupby] = pd.Categorical(groups)
    adata.obsm['X_ciber'] = np.asarray(embedding)
    return adata
")
  
  py <- reticulate::import_main(convert = FALSE)
  
  adata <- py$`_ciber_make_anndata`(
    reticulate::r_to_py(t(counts)),
    reticulate::r_to_py(as.list(colnames(counts))),
    reticulate::r_to_py(as.list(rownames(counts))),
    reticulate::r_to_py(as.list(meta)),
    reticulate::r_to_py(embedding_mat),
    group.by
  )
  
  oracle <- celloracle$Oracle()
  
  oracle$import_anndata_as_raw_count(
    adata = adata,
    cluster_column_name = group.by,
    embedding_name = "X_ciber"
  )
  
  if (is.character(base_GRN)) {
    oracle$import_TF_data(
      TF_info_matrix_path = normalizePath(
        base_GRN,
        winslash = "/",
        mustWork = TRUE
      )
    )
  } else {
    oracle$import_TF_data(
      TF_info_matrix = reticulate::r_to_py(base_GRN)
    )
  }
  
  oracle$perform_PCA()
  
  n_pca_dims <- min(
    as.integer(n_pca_dims),
    nrow(counts),
    ncol(counts) - 1L
  )
  
  if (is.null(k)) {
    k <- max(2L, floor(0.025 * ncol(counts)))
  }
  
  k <- min(as.integer(k), ncol(counts) - 1L)
  
  oracle$knn_imputation(
    n_pca_dims = n_pca_dims,
    k = k,
    balanced = TRUE,
    b_sight = min(k * 8L, ncol(counts) - 1L),
    b_maxl = min(k * 4L, ncol(counts) - 1L),
    n_jobs = as.integer(n_jobs)
  )
  
  links <- oracle$get_links(
    cluster_name_for_GRN_unit = group.by,
    alpha = alpha,
    bagging_number = as.integer(bagging_number),
    verbose_level = if (verbal) 1L else 0L,
    n_jobs = as.integer(n_jobs)
  )
  
  oracle_path <- file.path(
    save_dir,
    "CIBER.celloracle.oracle"
  )
  
  links_path <- file.path(
    save_dir,
    "CIBER.celloracle.links"
  )
  
  oracle$to_hdf5(oracle_path)
  links$to_hdf5(links_path)
  
  list(
    oracle = oracle,
    links = links,
    oracle_path = oracle_path,
    links_path = links_path
  )
}


#' CellOracle-based GRN perturbation
#'
#' @param object A Seurat object.
#' @param net_struc A cell-type DAG.
#' @param oracle A CellOracle Oracle object or file path.
#' @param links A CellOracle Links object or file path.
#' @param save_dir Directory for CellOracle and intermediate SCM results.
#' @param base_GRN Base GRN used when `oracle` or `links` is `NULL`.
#' @param group.by Cell-type column in the Seurat metadata.
#' @param assay Assay used for SCM calculation.
#' @param layer Expression layer used for SCM calculation.
#' @param counts_layer Raw-count layer used to build the CellOracle object.
#' @param embedding Seurat dimensional reduction used by CellOracle.
#' @param features Genes used when building the CellOracle object.
#' @param genes Genes to perturb; `NULL` uses all available regulatory genes.
#' @param perturb_value Target expression value passed to CellOracle.
#' @param n_propagation Number of CellOracle propagation steps.
#' @param n_pca_dims Number of principal components used for KNN imputation.
#' @param k Number of nearest neighbours.
#' @param alpha Ridge coefficient.
#' @param bagging_number Number of bagging estimators.
#' @param n_jobs Number of CPU cores used by CellOracle.
#' @param filter_links Logical; whether to filter the Links object.
#' @param overwrite Logical; whether to overwrite cached SCM results.
#' @param verbal Logical; whether to print progress messages.
#'
#' @return A list compatible with the output of `PerturbResult()`.
#' @export
GRNPerturbResult <- function(
    object,
    net_struc,
    oracle = NULL,
    links = NULL,
    save_dir,
    base_GRN = NULL,
    group.by = "celltype",
    assay = "RNA",
    layer = "data",
    counts_layer = "counts",
    embedding = "umap",
    features = NULL,
    genes = NULL,
    perturb_value = 0,
    n_propagation = 3,
    n_pca_dims = 50,
    k = NULL,
    alpha = 10,
    bagging_number = 20,
    n_jobs = -1,
    filter_links = TRUE,
    overwrite = FALSE,
    verbal = TRUE
) {
  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object.")
  }
  
  if (!group.by %in% colnames(object@meta.data)) {
    stop("Metadata column '", group.by, "' was not found.")
  }
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  scm_dir <- file.path(save_dir, "scm")
  dir.create(scm_dir, showWarnings = FALSE)
  
  celloracle <- reticulate::import(
    "celloracle",
    convert = FALSE
  )
  
  if (is.null(oracle) || is.null(links)) {
    prepared <- PrepareCellOracle(
      object = object,
      base_GRN = base_GRN,
      save_dir = file.path(save_dir, "celloracle"),
      group.by = group.by,
      assay = assay,
      counts_layer = counts_layer,
      embedding = embedding,
      features = features,
      n_pca_dims = n_pca_dims,
      k = k,
      alpha = alpha,
      bagging_number = bagging_number,
      n_jobs = n_jobs,
      verbal = verbal
    )
    
    oracle <- prepared$oracle
    links <- prepared$links
  } else {
    if (is.character(oracle)) {
      oracle <- celloracle$load_hdf5(
        normalizePath(oracle, winslash = "/", mustWork = TRUE)
      )
    }
    
    if (is.character(links)) {
      links <- celloracle$load_hdf5(
        normalizePath(links, winslash = "/", mustWork = TRUE)
      )
    }
  }
  
  data <- SeuratObject::LayerData(
    object,
    assay = assay,
    layer = layer
  )
  
  data <- as.matrix(data)
  
  meta <- droplevels(
    factor(object@meta.data[colnames(data), group.by])
  )
  ctypes <- bnlearn::nodes(net_struc)
  
  mem.ref <- gem2mem(data, meta, "mean")
  
  if (!all(ctypes %in% colnames(mem.ref))) {
    stop("Some DAG nodes are absent from '", group.by, "'.")
  }
  
  mem.ref <- mem.ref[, ctypes, drop = FALSE]
  scm.ref <- get_scm(mem.ref, net_struc, "ref")
  
  if (filter_links) {
    links$filter_links()
  }
  
  oracle$get_cluster_specific_TFdict_from_Links(
    links_object = links
  )
  
  oracle_genes <- reticulate::py_to_r(
    oracle$adata$var_names$tolist()
  )
  
  regulatory_genes <- reticulate::py_to_r(
    oracle$all_regulatory_genes_in_TFdict
  )
  
  available_genes <- Reduce(
    intersect,
    list(regulatory_genes, oracle_genes, rownames(mem.ref))
  )
  
  if (is.null(genes)) {
    genes <- available_genes
  } else {
    missing_genes <- setdiff(genes, available_genes)
    
    if (length(missing_genes) > 0) {
      warning(
        "Skipped unavailable genes: ",
        paste(missing_genes, collapse = ", ")
      )
    }
    
    genes <- intersect(genes, available_genes)
  }
  
  if (length(genes) == 0) {
    stop("No valid genes are available for perturbation.")
  }
  
  cache_files <- file.path(
    scm_dir,
    paste0(make.names(genes, unique = TRUE), ".rds")
  )
  
  scm.list <- vector("list", length(genes))
  cached <- file.exists(cache_files) & !overwrite
  scm.list[cached] <- lapply(cache_files[cached], readRDS)
  
  if (verbal && any(cached)) {
    message(sum(cached), " cached perturbations loaded.")
  }
  
  if (any(!cached)) {
    oracle$fit_GRN_for_simulation(
      GRN_unit = "cluster",
      alpha = alpha,
      use_cluster_specific_TFdict = TRUE
    )
    
    reticulate::py_run_string("
import numpy as np

def _ciber_group_mean(oracle, groupby, levels):
    matrix = oracle.adata.layers['simulated_count']
    labels = oracle.adata.obs[groupby].astype(str).to_numpy()
    result = []

    for level in levels:
        mean = matrix[labels == str(level)].mean(axis=0)
        mean = mean.A1 if hasattr(mean, 'A1') else np.asarray(mean).ravel()
        result.append(mean)

    return np.stack(result, axis=1)
")
    
    py <- reticulate::import_main(convert = FALSE)
    
    oracle_group.by <- reticulate::py_to_r(
      oracle$cluster_column_name
    )
    
    common_genes <- intersect(
      rownames(mem.ref),
      oracle_genes
    )
    
    for (i in which(!cached)) {
      gene <- genes[i]
      
      if (verbal) {
        message("[", i, "/", length(genes), "] Perturbing ", gene)
      }
      
      condition <- reticulate::dict()
      condition[[gene]] <- as.numeric(perturb_value)
      
      oracle$simulate_shift(
        perturb_condition = condition,
        GRN_unit = "cluster",
        n_propagation = as.integer(n_propagation)
      )
      
      simulated.mem <- reticulate::py_to_r(
        py$`_ciber_group_mean`(
          oracle,
          oracle_group.by,
          reticulate::r_to_py(as.list(ctypes))
        )
      )
      
      rownames(simulated.mem) <- oracle_genes
      colnames(simulated.mem) <- ctypes
      
      perturb.mem <- mem.ref
      perturb.mem[common_genes, ] <- simulated.mem[
        common_genes,
        ,
        drop = FALSE
      ]
      
      scm.list[[i]] <- get_scm(
        perturb.mem,
        net_struc,
        gene
      )
      
      saveRDS(
        scm.list[[i]],
        cache_files[i]
      )
    }
  }
  
  scm.inter <- do.call(cbind, scm.list)
  colnames(scm.inter) <- genes
  
  list(
    n_sample = 1L,
    n_permutation = 1L,
    ref = list(ref = scm.ref),
    perturb = list(perturb = scm.inter)
  )
}

# ============================================================
# Internal helpers
# ============================================================

#' Calculate the one-dimensional energy distance
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @return A non-negative numeric scalar.
#' @keywords internal
.effect_energy1d <- function(x, y) {
  
  x <- sort(as.numeric(x))
  y <- sort(as.numeric(y))
  
  nx <- length(x)
  ny <- length(y)
  
  if (nx == 0 || ny == 0) {
    return(NA_real_)
  }
  
  within_x <- 2 * sum(
    (2 * seq_len(nx) - nx - 1) * x
  ) / nx^2
  
  within_y <- 2 * sum(
    (2 * seq_len(ny) - ny - 1) * y
  ) / ny^2
  
  cy <- c(0, cumsum(y))
  k <- findInterval(x, y)
  
  cross <- sum(
    k * x - cy[k + 1] +
      cy[ny + 1] - cy[k + 1] -
      (ny - k) * x
  ) / (nx * ny)
  
  sqrt(max(
    2 * cross - within_x - within_y,
    0
  ))
}


#' Calculate a linear-time MMD estimate
#'
#' @param x A numeric vector containing the reference values.
#' @param y A numeric vector containing the perturbed values.
#'
#' @return A non-negative numeric scalar.
#' @keywords internal
.effect_mmd_linear <- function(x, y) {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  if (length(x) == 0 || length(y) < 2) {
    return(0)
  }
  
  sigma <- stats::sd(c(x, y))
  
  if (!is.finite(sigma) || sigma <= 0) {
    return(0)
  }
  
  # Repeat the reference empirical distribution to match y
  x <- rep(x, length.out = length(y))
  
  n <- floor(length(y) / 2) * 2
  
  if (n < 2) {
    return(0)
  }
  
  x <- x[seq_len(n)]
  y <- y[seq_len(n)]
  
  id1 <- seq(1, n, 2)
  id2 <- id1 + 1
  
  kernel <- function(a, b) {
    exp(-(a - b)^2 / (2 * sigma^2))
  }
  
  mmd2 <- mean(
    kernel(x[id1], x[id2]) +
      kernel(y[id1], y[id2]) -
      kernel(x[id1], y[id2]) -
      kernel(x[id2], y[id1])
  )
  
  sqrt(max(mmd2, 0))
}


#' Calculate squared MMD with a Gaussian kernel
#'
#' @param x A numeric vector containing the reference values.
#' @param y A numeric vector containing the perturbed values.
#' @param sigma Gaussian-kernel bandwidth.
#' @param kxx Optional precomputed within-reference kernel mean.
#'
#' @return A non-negative numeric scalar.
#' @keywords internal
.effect_mmd2_rbf <- function(x, y, sigma, kxx = NULL) {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  if (!is.finite(sigma) || sigma <= 0) {
    return(0)
  }
  
  if (is.null(kxx)) {
    kxx <- mean(
      exp(
        -outer(x, x, "-")^2 /
          (2 * sigma^2)
      )
    )
  }
  
  kyy <- mean(
    exp(
      -outer(y, y, "-")^2 /
        (2 * sigma^2)
    )
  )
  
  kxy <- mean(
    exp(
      -outer(x, y, "-")^2 /
        (2 * sigma^2)
    )
  )
  
  max(kxx + kyy - 2 * kxy, 0)
}


#' Estimate edge-specific MMD bandwidths
#'
#' @param diffBN_ref A matrix of reference coefficients, with edges in rows.
#' @param diffBN_raw A three-dimensional array of perturbation coefficients.
#' @param n_sample Maximum number of perturbation values sampled per edge.
#' @param seed Random seed used for sampling.
#'
#' @return A numeric vector containing one bandwidth per edge.
#' @keywords internal
.effect_get_mmd_bandwidth <- function(
    diffBN_ref,
    diffBN_raw,
    n_sample = 2000,
    seed = 1) {
  
  set.seed(seed)
  
  n_edge <- dim(diffBN_raw)[1]
  n_gene <- dim(diffBN_raw)[2]
  n_rep <- dim(diffBN_raw)[3]
  
  sigma <- vapply(
    seq_len(n_edge),
    FUN.VALUE = numeric(1),
    FUN = function(j) {
      
      n <- min(n_sample, n_gene * n_rep)
      
      gene_id <- sample.int(
        n_gene,
        n,
        replace = TRUE
      )
      
      rep_id <- sample.int(
        n_rep,
        n,
        replace = TRUE
      )
      
      raw_sample <- diffBN_raw[
        cbind(
          rep(j, n),
          gene_id,
          rep_id
        )
      ]
      
      z <- c(
        as.numeric(diffBN_ref[j, ]),
        raw_sample
      )
      
      s <- stats::median(stats::dist(z))
      
      if (!is.finite(s) || s <= 0) {
        s <- stats::sd(z)
      }
      
      if (!is.finite(s) || s <= 0) {
        s <- 0
      }
      
      s
    }
  )
  
  sigma
}


#' Calculate effect matrix
#'
#' Calculate an effect matrix from perturbation results using coefficient
#' differences or distribution-based distance metrics.
#'
#' @param diffBN A list containing perturbation results. The object must contain
#'   \code{ref} and either \code{perturb} or \code{raw}. Depending on the
#'   selected metric, \code{n_sample}, \code{n_permutation}, or a precomputed
#'   \code{diff} component may also be used.
#'
#' @param dist_metric Effect metric (case-insensitive): \code{"diff_mean"}
#'   (default), \code{"W1"}, \code{"W2"}, \code{"energy"}, or \code{"mmd"}.
#'
#' @param mmd_n_sample Integer specifying the maximum number of perturbed
#'   coefficient values sampled per edge when estimating the Gaussian kernel
#'   bandwidth for \code{dist_metric = "mmd"}. Default is 2000.
#'
#' @param mmd_seed Integer random seed used for bandwidth estimation when
#'   \code{dist_metric = "mmd"}. Default is 1.
#'
#' @details
#' The returned effect matrix has perturbed genes or features in rows and
#' network edges in columns.
#'
#' \code{dist_metric = "diff_mean"} calculates the signed sum of paired
#' reference-minus-perturbed coefficient differences. No absolute value
#' or division by the number of runs is applied here. With equal repeat counts,
#' this differs from the signed mean shift by a common scale factor.
#' \code{diffScore(..., abs = TRUE)} takes absolute edge effects before summing
#' over the selected edges. Use that default for magnitude-based rankings.
#'
#' \code{dist_metric = "diff_mean"} calculates effects using
#' \code{get_diffCoeff(..., mode = "mean")} for each paired
#' reference and perturbation result, followed by summation across runs.
#'
#' \code{dist_metric = "W1"} and \code{dist_metric = "W2"} calculate the first- and
#' second-order Wasserstein distances, respectively, between the reference and
#' perturbed coefficient distributions for each gene-edge pair.
#'
#' \code{dist_metric = "energy"} calculates the one-dimensional energy distance
#' between the reference and perturbed coefficient distributions.
#'
#' \code{dist_metric = "mmd"} calculates the biased quadratic Gaussian-kernel MMD
#' squared statistic. A single edge-specific bandwidth is estimated using the
#' median heuristic and is shared across all perturbed genes for that edge.
#'
#' @return A data frame containing the effect matrix, with perturbed genes or
#'   features as rows and network edges as columns.
#'
#' @examples
#' \dontrun{
#' effect_mean <- EffectMatrix(diffBN, dist_metric = "diff_mean")
#' effect_w2 <- EffectMatrix(diffBN, dist_metric = "w2")
#' effect_energy <- EffectMatrix(diffBN, dist_metric = "energy")
#' effect_mmd <- EffectMatrix(diffBN, dist_metric = "MMD")
#' }
#'
#' @export
EffectMatrix <- function(
    diffBN,
    dist_metric = "diff_mean",
    mmd_n_sample = 2000,
    mmd_seed = 1) {
  normalize_metric <- function(value) {
    metrics <- c(diff_mean = "diff_mean", w1 = "W1", w2 = "W2",
                 energy = "energy", mmd = "mmd")
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !tolower(value) %in% names(metrics)) {
      stop("Invalid effect metric. Use diff_mean, W1, W2, energy, or mmd.")
    }
    unname(metrics[[tolower(value)]])
  }
  metric <- normalize_metric(dist_metric)
  
  if (is.null(diffBN$ref) || length(diffBN$ref) == 0) {
    stop("diffBN$ref is missing or empty.")
  }
  
  # ==========================================================
  # Signed coefficient-difference metric
  # ==========================================================
  
  if (metric == "diff_mean") {
    
    if (is.null(diffBN$perturb) || length(diffBN$perturb) == 0) {
      stop(
        "dist_metric = 'diff_mean' requires diffBN$perturb."
      )
    }
    
    n_sample <- diffBN$n_sample
    
    if (is.null(n_sample)) {
      n_sample <- length(diffBN$ref)
    }
    
    n_permutation <- diffBN$n_permutation
    
    if (is.null(n_permutation)) {
      
      if (length(diffBN$perturb) %% n_sample != 0) {
        stop(
          "Cannot infer n_permutation from diffBN$perturb."
        )
      }
      
      n_permutation <- length(diffBN$perturb) / n_sample
    }
    
    if (
      n_sample * n_permutation !=
      length(diffBN$perturb)
    ) {
      stop(
        "n_sample * n_permutation does not match ",
        "the number of perturbation results."
      )
    }
    
    diffBN_result <- vector(
      "list",
      length(diffBN$perturb)
    )
    
    k <- 1L
    
    for (i in seq_len(n_sample)) {
      
      for (j in seq_len(n_permutation)) {
        
        diffBN_result[[k]] <- get_diffCoeff(
          scm.inter = diffBN$perturb[[k]],
          scm.ref = diffBN$ref[[i]],
          mode = "mean"
        )
        
        k <- k + 1L
      }
    }
    
    result <- as.matrix(diffBN_result[[1]])
    
    if (length(diffBN_result) > 1) {
      
      for (k in 2:length(diffBN_result)) {
        result <- result +
          as.matrix(diffBN_result[[k]])
      }
    }
    
    return(as.data.frame(result))
  }
  
  # ==========================================================
  # New benchmark metrics
  # ==========================================================
  
  raw_list <- if (!is.null(diffBN$raw)) {
    diffBN$raw
  } else {
    diffBN$perturb
  }
  
  if (is.null(raw_list) || length(raw_list) == 0) {
    stop(
      "diffBN must contain either $raw or $perturb."
    )
  }
  
  n_edge <- nrow(diffBN$ref[[1]])
  n_features <- ncol(raw_list[[1]])
  
  if (is.null(n_edge) || is.null(n_features)) {
    stop(
      "diffBN$ref and perturbation results must be ",
      "matrix-like objects."
    )
  }
  
  if (
    any(vapply(
      diffBN$ref,
      nrow,
      integer(1)
    ) != n_edge)
  ) {
    stop(
      "All reference matrices must contain the same number of edges."
    )
  }
  
  raw_dim <- vapply(
    raw_list,
    function(x) {
      paste(dim(as.matrix(x)), collapse = "x")
    },
    character(1)
  )
  
  if (length(unique(raw_dim)) != 1) {
    stop(
      "All perturbation matrices must have identical dimensions."
    )
  }
  
  # ==========================================================
  # Assemble empirical coefficient distributions
  #
  # ref: edge x reference replicate
  # raw: edge x feature x perturbation replicate
  # ==========================================================
  
  diffBN_ref <- do.call(
    cbind,
    lapply(
      diffBN$ref,
      as.matrix
    )
  )
  
  diffBN_raw <- simplify2array(
    lapply(
      raw_list,
      as.matrix
    )
  )
  
  # simplify2array drops the third dimension for a single run
  if (length(dim(diffBN_raw)) == 2) {
    
    diffBN_raw <- array(
      diffBN_raw,
      dim = c(
        nrow(diffBN_raw),
        ncol(diffBN_raw),
        1L
      )
    )
  }
  
  if (
    length(dim(diffBN_raw)) != 3 ||
    dim(diffBN_raw)[1] != n_edge ||
    dim(diffBN_raw)[2] != n_features
  ) {
    stop(
      "Unable to construct the perturbation coefficient array."
    )
  }
  
  # ==========================================================
  # Fixed bandwidth per edge for quadratic MMD
  # ==========================================================
  
  if (metric == "mmd") {
    
    mmd_sigma <- .effect_get_mmd_bandwidth(
      diffBN_ref = diffBN_ref,
      diffBN_raw = diffBN_raw,
      n_sample = mmd_n_sample,
      seed = mmd_seed
    )
    
    mmd_kxx <- vapply(
      seq_len(n_edge),
      FUN.VALUE = numeric(1),
      FUN = function(j) {
        
        sigma <- mmd_sigma[j]
        
        if (!is.finite(sigma) || sigma <= 0) {
          return(0)
        }
        
        x <- as.numeric(
          diffBN_ref[j, ]
        )
        
        mean(
          exp(
            -outer(x, x, "-")^2 /
              (2 * sigma^2)
          )
        )
      }
    )
  }
  
  # ==========================================================
  # Calculate feature x edge effect matrix
  # ==========================================================
  
  result <- matrix(
    NA_real_,
    nrow = n_features,
    ncol = n_edge
  )
  
  for (i in seq_len(n_features)) {
    
    for (j in seq_len(n_edge)) {
      
      x <- as.numeric(
        diffBN_ref[j, ]
      )
      
      y <- as.numeric(
        diffBN_raw[j, i, ]
      )
      
      result[i, j] <- switch(
        metric,
        
        W1 = transport::wasserstein1d(
          x,
          y,
          p = 1
        ),
        
        W2 = transport::wasserstein1d(
          x,
          y,
          p = 2
        ),
        
        energy = .effect_energy1d(
          x,
          y
        ),
        
        mmd = .effect_mmd2_rbf(
          x,
          y,
          sigma = mmd_sigma[j],
          kxx = mmd_kxx[j]
        )
      )
    }
  }
  
  rownames(result) <-
    colnames(raw_list[[1]])
  
  colnames(result) <-
    rownames(diffBN$ref[[1]])
  
  as.data.frame(result)
}

#' Calculate graph coefficients after gene perturbation
#'
#' `PerturbResult()` calculates graph coefficients from reference and perturbed
#' expression data. By default, each gene is perturbed in turn by multiplying
#' its expression by `perturb_ratio`. Set `perturb_ratio = "deletion"` to use
#' feature deletion instead.
#'
#' @param net_struc A `bn.fit` network structure.
#' @param data A data frame or Seurat object containing expression data.
#' @param meta A character or factor vector of cell-type labels corresponding to
#'   the columns of `data` in single-cell mode.
#' @param index A list of column indices, one element per cell sample.
#' @param n_sample Number of cell samples.
#' @param n_permutation Retained for backward compatibility; ratio perturbation
#'   and deletion each produce one perturbation result per cell sample.
#' @param perturb_ratio A non-negative numeric scalar by which expression of
#'   the perturbed feature is multiplied, or `"deletion"` to remove each
#'   feature in turn. A value of `0` performs zeroing, values between `0` and
#'   `1` perform knockdown, and values greater than `1` perform knockup.
#' @param mode Character string specifying the input data mode:
#'   `"single_cell"` or `"bulk"`.
#' @param perturbation Character string specifying the perturbation mode:
#'   `"single_feature"` for GRN-free perturbation or `"GRN"` for
#'   CellOracle-based GRN perturbation.
#' @param ncores Number of cores used for parallel computation.
#' @param verbal Logical; whether to print progress messages.
#' @param replace Retained for backward compatibility with permutation-based
#'   perturbation; it is not used by ratio perturbation or deletion.
#' @param oracle A CellOracle Oracle object or file path used in GRN mode.
#' @param links A CellOracle Links object or file path used in GRN mode.
#' @param save_dir Directory for CellOracle and intermediate SCM results.
#' @param group.by Cell-type column in the Seurat metadata.
#' @param assay Assay used for SCM calculation in GRN mode.
#' @param layer Expression layer used for SCM calculation in GRN mode.
#' @param genes Genes to perturb in GRN mode; `NULL` uses all available
#'   regulatory genes.
#' @param perturb_value Target expression value passed to CellOracle in GRN
#'   mode.
#' @param n_propagation Number of CellOracle propagation steps.
#' @param alpha Ridge coefficient used by CellOracle.
#' @param filter_links Logical; whether to filter the CellOracle Links object.
#' @param overwrite Logical; whether to overwrite cached SCM results.
#' @param n_pca_dims Number of principal components used for KNN imputation.
#' @param k Number of nearest neighbours.
#' @param bagging_number Number of bagging estimators.
#' @param n_jobs Number of CPU cores used by CellOracle.
#' @param base_GRN Base GRN used when `oracle` or `links` is `NULL`.
#' @param counts_layer Raw-count layer used to build the CellOracle object.
#' @param embedding Seurat dimensional reduction used by CellOracle.
#' @param features Genes used when building the CellOracle object.
#' @param ... Additional arguments passed to `GRNPerturbResult()` in GRN mode.
#'
#' @return A list containing reference and perturbed graph coefficients.
#' @export
PerturbResult <- function(
    net_struc,
    data,
    meta = NULL,
    index = NULL,
    n_sample = 1,
    n_permutation = 1,
    perturb_ratio = 0,
    mode = "single_cell",
    perturbation = "single_feature",
    ncores = 1,
    verbal = FALSE,
    replace = FALSE,
    oracle = NULL,
    links = NULL,
    save_dir = NULL,
    group.by = "celltype",
    assay = "RNA",
    layer = "data",
    genes = NULL,
    perturb_value = 0,
    n_propagation = 3,
    alpha = 10,
    filter_links = TRUE,
    overwrite = FALSE,
    n_pca_dims = 50,
    k = NULL,
    bagging_number = 20,
    n_jobs = -1,
    base_GRN = NULL,
    counts_layer = "counts",
    embedding = "umap",
    features = NULL,
    ...
) {
  mode <- match.arg(mode, c("single_cell", "bulk"))
  perturbation <- match.arg(perturbation, c("single_feature", "GRN"))
  dots <- list(...)
  if ("deletion" %in% names(dots)) {
    stop("'deletion' is not a separate argument; use perturb_ratio = \"deletion\".")
  }
  if (perturbation == "GRN" && mode != "single_cell") {
    stop("GRN perturbation is only supported with mode = 'single_cell'.")
  }
  is_deletion <- identical(perturb_ratio, "deletion")
  is_numeric_ratio <-
    length(perturb_ratio) == 1 &&
    is.numeric(perturb_ratio) &&
    is.finite(perturb_ratio) &&
    perturb_ratio >= 0
  if (perturbation == "single_feature" && !is_deletion && !is_numeric_ratio) {
    stop(
      "'perturb_ratio' must be a single non-negative finite number or ",
      "'deletion'."
    )
  }

  doParallel::registerDoParallel(ncores)
  if (perturbation == "single_feature" && mode == "single_cell") {
    if (is.null(n_sample)) {
      stop("Sample number not indicated!")
    }
    if (verbal) {
      print(paste0(
        "Permutation method: ",
        ifelse(
          is_deletion,
          "deletion.",
          paste0("expression ratio ", perturb_ratio, ".")
        )
      ))
    }
    if (class(data)[1] == "Seurat") {
      data <- as.data.frame(as.matrix(data@assays$RNA@data))
    }
    refs <- vector(mode = "list", length = n_sample)
    raws <- vector(mode = "list", length = n_sample)
    for (sample_id in seq_len(n_sample)) {
      data_tmp <- gem2mem(
        data[, index[[sample_id]], drop = FALSE],
        meta[index[[sample_id]]],
        "mean"
      )
      scm.ref <- get_scm(mem = data_tmp, graph = net_struc, id = "ref")
      Features <- rownames(data_tmp)
      scm.inter <- foreach::foreach(
        gene_id = seq_len(nrow(data_tmp)),
        .combine = cbind,
        .export = "get_scm"
      ) %dopar% {
        perturb_tmp <- data_tmp
        if (is_deletion) {
          perturb_tmp <- perturb_tmp[-gene_id, , drop = FALSE]
        } else {
          perturb_tmp[gene_id, ] <-
            perturb_tmp[gene_id, ] * perturb_ratio
        }
        get_scm(
          mem = perturb_tmp,
          graph = net_struc,
          id = Features[gene_id]
        )
      }
      refs[[sample_id]] <- scm.ref
      raws[[sample_id]] <- scm.inter
      if (verbal) {
        print(paste0("Sample ", sample_id, " calculation done."))
      }
    }
    names(raws) <- paste0("sample_", seq_len(n_sample), "/", n_sample)
    names(refs) <- paste0("sample_", seq_len(n_sample), "/", n_sample)
    diffCoeff_final_result <- list(
      perturb = raws,
      ref = refs,
      n_sample = n_sample,
      n_permutation = 1
    )
  } else if (perturbation == "single_feature" && mode == "bulk") {
    scm.ref <- get_scm(mem = data, graph = net_struc, id = "ref")
    Features <- rownames(data)
    scm.inter <- foreach::foreach(
      gene_id = seq_len(nrow(data)),
      .combine = cbind,
      .export = "get_scm"
    ) %dopar% {
      perturb_tmp <- data
      if (is_deletion) {
        perturb_tmp <- perturb_tmp[-gene_id, , drop = FALSE]
      } else {
        perturb_tmp[gene_id, ] <-
          perturb_tmp[gene_id, ] * perturb_ratio
      }
      get_scm(
        mem = perturb_tmp,
        graph = net_struc,
        id = Features[gene_id]
      )
    }
    diffCoeff_final_result <- list(
      perturb = list(scm.inter),
      ref = list(scm.ref),
      n_sample = 1,
      n_permutation = 1
    )
  } else if (perturbation == "GRN") {
    grn_args <- list(
      object = data,
      net_struc = net_struc,
      oracle = oracle,
      links = links,
      save_dir = save_dir,
      base_GRN = base_GRN,
      group.by = group.by,
      assay = assay,
      layer = layer,
      genes = genes,
      perturb_value = perturb_value,
      n_propagation = n_propagation,
      alpha = alpha,
      filter_links = filter_links,
      overwrite = overwrite,
      n_pca_dims = n_pca_dims,
      k = k,
      bagging_number = bagging_number,
      n_jobs = n_jobs,
      counts_layer = counts_layer,
      embedding = embedding,
      features = features,
      verbal = verbal
    )
    
    grn_args <- modifyList(
      grn_args,
      dots
    )
    
    diffCoeff_final_result <- do.call(
      GRNPerturbResult,
      grn_args
    )
    
  }
  return(diffCoeff_final_result)
}

#' Calculate coefficient differences
#'
#' Calculate coefficient differences from a set of `bn.fit` results and an
#' unperturbed reference.
#'
#' @param scm.inter A data frame of perturbed `bn.fit` coefficients.
#' @param scm.ref A data frame of reference `bn.fit` coefficients.
#' @param mode Character string specifying how coefficient differences are
#'   calculated. Currently, only `"mean"` is supported.
#'
#' @return A data frame of coefficient differences.
get_diffCoeff <- function(scm.inter = NULL, scm.ref = NULL, mode = "mean") {
  if (mode == "mean") {
    nGene <- ncol(scm.inter)
    diffCoeff_trans <- rep(scm.ref, nGene) - scm.inter
    diffCoeff <- diffCoeff_trans %>%
      t() %>%
      as.data.frame()
  } else {
    diffCoeff <- NULL
    warning("Invalid mode!")
  }
  return(diffCoeff)
}

#' Calculate coefficients after permutation-based perturbation
#'
#' Calculate graph coefficients after gene-wise permutation of single-cell
#' expression data.
#'
#' @param net_struc A `bn.fit` network structure.
#' @param data A data frame of single-cell expression data.
#' @param meta A character or factor vector of cell-type labels.
#' @param index A list of column indices, one element per cell sample.
#' @param n_sample Number of cell samples.
#' @param n_permutation Number of gene permutations per cell sample.
#' @param diffCoeff_mode Character string specifying the coefficient-difference
#'   calculation mode.
#' @param verbal Logical; whether to print progress messages.
#' @param replace Logical; whether permutation sampling is performed with
#'   replacement.
#'
#' @return A list containing reference and perturbed graph coefficients.
run_diffCoeff <- function(
    net_struc,
    data,
    meta,
    index,
    n_sample = 20,
    n_permutation = 20,
    diffCoeff_mode = "mean",
    verbal = F,
    replace = F
) {
  diffCoeff_ref <- vector(mode = "list", length = n_sample)
  meta <- as.matrix(meta)
  Features <- rownames(data)
  ds <- vector(mode = "list", length = n_sample)
  ds_per_all <- vector(mode = "list", length = n_sample)
  data_tmp <- vector(mode = "list", length = n_sample)
  meta_tmp <- vector(mode = "list", length = n_sample)
  if (verbal) {
    print("Begin diffCoeff reference calculation...")
  }
  for (i in 1:n_sample) {
    data_tmp[[i]] <- data[, index[[i]]]
    if (class(data_tmp[[i]])[1] == "Seurat") {
      data_tmp[[i]] <- as.data.frame(as.matrix(data_tmp[[i]]@assays$RNA@data))
    }
    meta_tmp[[i]] <- meta[index[[i]]]
    ds[[i]] <- gem2mem(data_tmp[[i]], meta_tmp[[i]], "mean")
    scm.ref <- get_scm(mem = ds[[i]], graph = net_struc, id = "ref")
    diffCoeff_ref[[i]] <- scm.ref
    names(diffCoeff_ref)[i] <- paste0("sample:", i, "_ref")
    
    ds_per_all[[i]] <- vector(mode = "list", length = n_permutation)
    
    for (j in 1:n_permutation) {
      data_smpl_per <- data_tmp[[i]]
      data_smpl_per <- foreach::foreach(k = 1:nrow(data_smpl_per)) %dopar%
        {
          sample(data_smpl_per[k, ], replace = replace)
        } # as a list
      data_smpl_per <- mapply(c, data_smpl_per) %>%
        t() %>%
        as.data.frame()
      rownames(data_smpl_per) <- rownames(data_tmp[[i]])
      colnames(data_smpl_per) <- colnames(data_tmp[[i]])
      data_smpl_per[] <- lapply(data_smpl_per, as.numeric)
      ds_per_all[[i]][[j]] <- gem2mem(
        gem = data_smpl_per,
        meta = meta_tmp[[i]],
        "mean"
      )
    }
    if (verbal) print(paste0("Reference #", i, " calculation complete."))
  }
  if (verbal) {
    ("Begin gene permutation...")
  }
  paral_index <- data.frame(
    i = rep(1:n_sample, each = n_permutation),
    h = rep(1:n_permutation, times = n_sample)
  )
  diffCoeff_final_result <- foreach::foreach(
    k = 1:(n_sample * n_permutation),
    .combine = rbind
  ) %dopar%
    {
      i <- paral_index[k, "i"]
      h <- paral_index[k, "h"]
      ds_per <- ds[[i]]
      scm.inter <- foreach::foreach(j = 1:nrow(ds[[i]]), .combine = cbind) %do%
        {
          if (j > 1) {
            ds_per[j - 1, ] <- ds[[i]][j - 1, ]
          }
          ds_per[j, ] <- ds_per_all[[i]][[h]][j, ]
          get_scm(ds_per, graph = net_struc, id = Features[j])
        }
      if (verbal) {
        print(paste0("n_sample:", i, "/n_permutation:", h, " done."))
      }
      list(k = k, perturb = scm.inter)
    }
  if (verbal) {
    print("All calculation done successfully!")
  }
  if (n_sample * n_permutation == 1) {
    l_name <- paste0(
      "sample_",
      paral_index[1, "i"],
      "/",
      n_sample,
      "_permute_",
      paral_index[1, "h"],
      "/",
      n_permutation
    )
    diffCoeff_final_result <- list(
      perturb = list(diffCoeff_final_result$perturb),
      ref = diffCoeff_ref,
      n_permutation = 1,
      n_sample = 1
    )
    names(diffCoeff_final_result$perturb) <- l_name
    return(diffCoeff_final_result)
  }
  row_index <- as.vector(unlist(diffCoeff_final_result[, "k"]))
  rownames(diffCoeff_final_result) <- paste0(
    "sample_",
    paral_index[row_index, "i"],
    "/",
    n_sample,
    "_permute_",
    paral_index[row_index, "h"],
    "/",
    n_permutation
  )
  diffCoeff_final_result <- list(
    perturb = diffCoeff_final_result[, "perturb"],
    ref = diffCoeff_ref,
    n_permutation = n_permutation,
    n_sample = n_sample
  )
  return(diffCoeff_final_result)
}

#' Permute expression values for a specific gene
#'
#' @param data_per A data frame containing expression data.
#' @param gene_pos Integer index of the gene to permute.
#'
#' @return The data frame with the selected gene permuted.
g_per <- function(data_per, gene_pos) {
  data_per[gene_pos, ] <- as.numeric(sample(
    data_per[gene_pos, ],
    length(data_per[gene_pos, ])
  ))
  return(data_per)
}

#' Calculate structural causal model coefficients
#'
#' Calculate `bn.fit` coefficients for a network and a cell-type expression
#' data frame.
#'
#' @param mem A data frame of cell-type expression values.
#' @param graph A `bn.fit` network structure.
#' @param id Character identifier for the sampling and perturbation result.
#'
#' @return A data frame of fitted edge coefficients.
get_scm <- function(mem = NULL, graph = NULL, id = NULL) {
  if (is.null(mem)) {
    stop("MEM data is missing")
  }
  if (is.null(graph)) {
    stop("DAG used for linear regression is missing")
  }
  if (is.null(id)) {
    stop(
      "id should be given to label the information of sampling and permutation, best in the form of s1p1"
    )
  }
  
  fit_result <- bnlearn::bn.fit(graph, mem)
  scm <- as.data.frame(bnlearn::arcs(graph))
  graph_coef <- stats::coef(fit_result)
  scm[id] <- sapply(
    rownames(scm),
    function(x) graph_coef[[scm[x, "to"]]][scm[x, "from"]]
  )
  scm_result <- scm[, id, drop = FALSE]
  rownames(scm_result) <- paste(scm[["from"]], scm[["to"]], sep = "~")
  return(scm_result)
}
