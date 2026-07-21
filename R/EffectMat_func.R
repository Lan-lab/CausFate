#' CellOracle-based GRN perturbation
#'
#' @param object Seurat object
#' @param net_struc Cell-type DAG
#' @param oracle Preprocessed CellOracle Oracle object
#' @param links Path to a CellOracle Links file
#' @param save_dir Directory for intermediate SCM results
#' @param group.by Cell-type column in Seurat metadata
#' @param assay Assay used for SCM calculation
#' @param layer Expression layer used for SCM calculation
#' @param genes Genes to perturb; NULL uses all available regulatory genes
#' @param perturb_value Target expression value
#' @param n_propagation Number of CellOracle propagation steps
#' @param alpha Ridge coefficient used for GRN fitting
#' @param filter_links Whether to filter Links before GRN fitting
#' @param overwrite Whether to overwrite cached results
#' @param verbal Whether to print progress
#'
#' @return A PerturbResult-compatible list
#' @export
GRNPerturbResult <- function(
    object,
    net_struc,
    oracle,
    links,
    save_dir,
    group.by = "celltype",
    assay = "RNA",
    layer = "data",
    genes = NULL,
    perturb_value = 0,
    n_propagation = 3,
    alpha = 10,
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

  data <- SeuratObject::LayerData(
    object = object,
    assay = assay,
    layer = layer
  )

  meta <- object@meta.data[[group.by]]
  ctypes <- bnlearn::nodes(net_struc)

  mem.ref <- gem2mem(data, meta, "mean")

  if (!all(ctypes %in% colnames(mem.ref))) {
    stop("Some DAG nodes are absent from '", group.by, "'.")
  }

  mem.ref <- mem.ref[, ctypes, drop = FALSE]
  scm.ref <- get_scm(mem.ref, net_struc, "ref")

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
    save_dir,
    paste0(make.names(genes, unique = TRUE), ".rds")
  )

  scm.list <- vector("list", length(genes))
  cached <- file.exists(cache_files) & !overwrite

  scm.list[cached] <- lapply(cache_files[cached], readRDS)

  if (verbal && any(cached)) {
    message(sum(cached), " cached perturbations loaded.")
  }

  if (any(!cached)) {
    celloracle <- reticulate::import(
      "celloracle",
      convert = FALSE
    )

    links.object <- celloracle$load_hdf5(
      normalizePath(links, winslash = "/", mustWork = TRUE)
    )

    if (filter_links) {
      links.object$filter_links()
    }

    oracle$get_cluster_specific_TFdict_from_Links(
      links_object = links.object
    )

    oracle$fit_GRN_for_simulation(
      GRN_unit = "cluster",
      alpha = alpha,
      use_cluster_specific_TFdict = TRUE
    )

    reticulate::py_run_string(
      "
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
"
    )

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

      saveRDS(scm.list[[i]], cache_files[i])
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

#' Calculate effect matrix
#'
#' Calculate effect matrix from perturbation results.
#'
#' @param diffBN list of diffBN raw result
#' @param mode character, determine the diffBN function
#'
#' @return Effect matrix
#' @export
EffectMatrix <- function(diffBN, mode = "mean") {
  diffBN_result <- list(NULL)
  n_sample <- diffBN$n_sample
  n_permutation <- diffBN$n_permutation
  if (mode == "mean") {
    k <- 1
    for (i in 1:(n_sample)) {
      for (j in 1:(n_permutation)) {
        diffBN_result[[k]] <- get_diffCoeff(
          scm.inter = diffBN$perturb[[k]],
          scm.ref = diffBN$ref[[i]],
          mode = mode
        )
        k <- k + 1
      }
    }
    names(diffBN_result) <- names(diffBN$perturb)
    tmp <- as.matrix(diffBN_result[[1]])
    if (n_sample * n_permutation > 1) {
      for (k in 2:(n_sample * n_permutation)) {
        tmp <- tmp + as.matrix(diffBN_result[[k]])
      }
    }
    diffBN_result <- tmp
  } else if (mode == "OT") {
    n_edge <- nrow(diffBN$ref[[1]])
    n_features <- ncol(diffBN$perturb[[1]])
    diffBN_ref <- str2str::ld2a(diffBN$ref) %>% as.data.frame()
    diffBN_raw <- str2str::ld2a(diffBN$perturb)
    # row=gene col=edge
    diffBN_result <- foreach::foreach(i = 1:n_features, .combine = rbind) %do%
      {
        foreach::foreach(j = 1:n_edge, .combine = c) %do%
          {
            transport::wasserstein1d(
              as.numeric(diffBN_ref[j, ]),
              diffBN_raw[j, i, ]
            )
          }
      }
    colnames(diffBN_result) <- rownames(diffBN$ref[[1]])
    rownames(diffBN_result) <- colnames(diffBN$perturb[[1]])
  } else {
    stop("Invalid diffBN calculation mode!")
  }

  return(as.data.frame(diffBN_result))
}

#' PerturbResult
#'
#' PerturbResult is a function that calculate coefficients of the graph structure based on a given expression data
#'
#' @param net_struc bn.fit structure
#' @param data data.frame of the expression data
#' @param meta character of metadata corresponding to expression data
#' @param index list of sampling results
#' @param n_permutation integer, times of gene permutation, usually 20-50
#' @param n_sample integer, times of cell sampling, usually 20-50
#' @param deletion, logical, whether the permutation uses deletion, if this is set to true, n_permutation and replace will not be used
#' @param mode character, determine the diffCoeff function, usually "mean"
#' @param ncores integer, sets the number of cores used in the parallel computing
#' @param mode character, can be 'single_cell' or 'bulk'
#' @param verbal logical
#' @param replace logical, whether the sampling will be performed with replacing
#'
#' @return diffCoeff results
#' @export
PerturbResult <- function(
    net_struc,
    data,
    meta = NULL,
    index = NULL,
    n_sample = 1,
    n_permutation = 1,
    deletion = T,
    mode = "single_cell",
    ncores = 1,
    verbal = F,
    replace = F,
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
    overwrite = FALSE
) {
  doParallel::registerDoParallel(ncores)
  if (mode == "single_cell") {
    if (is.null(n_sample)) {
      stop("Sample number not indicated!")
    }
    if (verbal) {
      print(paste0(
        "Permutation method: ",
        ifelse(
          deletion,
          "deletion.",
          paste0(
            "sampling ",
            ifelse(replace, "with replacing.", "without replacing.")
          )
        )
      ))
    }
    if (deletion) {
      if (class(data)[1] == "Seurat") {
        data <- as.data.frame(as.matrix(data@assays$RNA@data))
      }
      refs <- vector(mode = "list", length = n_sample)
      raws <- vector(mode = "list", length = n_sample)
      for (k in 1:n_sample) {
        data_tmp <- gem2mem(data[, index[[k]]], meta[index[[k]]], "mean")
        scm.ref <- get_scm(mem = data_tmp, graph = net_struc, id = "ref") # get reference
        Features <- rownames(data_tmp)
        scm.inter <- foreach::foreach(
          i = 1:nrow(data_tmp),
          .combine = cbind
        ) %dopar%
          {
            get_scm(mem = data_tmp[-i, ], graph = net_struc, id = Features[i])
          } # get diffCoeff for every gene
        refs[[k]] <- scm.ref
        raws[[k]] <- scm.inter
        if (verbal) print(paste0("Sample ", k, " calculation done."))
      }
      names(raws) <- paste0("sample_", 1:n_sample, "/", n_sample)
      names(refs) <- paste0("sample_", 1:n_sample, "/", n_sample)
      diffCoeff_final_result <- list(
        perturb = raws,
        ref = refs,
        n_sample = n_sample,
        n_permutation = 1
      )
    } else {
      diffCoeff_final_result <- run_diffCoeff(
        net_struc,
        data,
        meta,
        index,
        n_sample,
        n_permutation,
        verbal,
        replace
      )
    }
  } else if (mode == "bulk") {
    scm.ref <- get_scm(mem = data, graph = net_struc, id = "ref") # get reference
    Features <- rownames(data)
    scm.inter <- foreach::foreach(i = 1:nrow(data), .combine = cbind) %dopar%
      {
        get_scm(mem = data[-i, ], graph = net_struc, id = Features[i])
      } # get diffCoeff for every gene
    diffCoeff_final_result <- list(
      perturb = list(scm.inter),
      ref = scm.ref,
      n_sample = 1,
      n_permutation = 1
    )
  } else if (mode == "GRN") {
    diffCoeff_final_result <- GRNPerturbResult(
      object = data,
      net_struc = net_struc,
      oracle = oracle,
      links = links,
      save_dir = save_dir,
      group.by = group.by,
      assay = assay,
      layer = layer,
      genes = genes,
      perturb_value = perturb_value,
      n_propagation = n_propagation,
      alpha = alpha,
      filter_links = filter_links,
      overwrite = overwrite,
      verbal = verbal
    )
  } else {
    diffCoeff_final_result <- NULL
    print("Wrong mode!")
  }
  return(diffCoeff_final_result)
}

#' get_diffCoeff
#'
#' get_diffCoeff: based on a set of bn.fit results and a bn.fir reference, calculate the diffCoeff results
#'
#' @param scm.inter data.frame of bn.fit results
#' @param scm.ref data.frame of bn.fit reference, which don't have any gene deletion or permutation
#' @param mode character, determine how the diffCoeff is calculate, usually use "mean"
#'
#' @return data.frame of diffCoeff
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

#' run_diffCoeff
#'
#' run_diffCoeff is a function that calculate diffCoeff results based on a single cell expression data
#'
#' @param net_struc bn.fit structure
#' @param data data.frame of single cell expression data
#' @param meta character of metadata corresponding to expression data
#' @param index list of sampling results
#' @param n_permutation integer, times of gene permutation, usually 20-50
#' @param n_sample integer, times of cell sampling, usually 20-50
#' @param diffCoeff_mode character, determine the diffCoeff function, usually "mean"
#'
#' @return diffCoeff results
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

#' Permutate specific gene data
#'
#' @param data_per data.frame
#' @param gene_pos interger
#'
#' @return data.frame
g_per <- function(data_per, gene_pos) {
  data_per[gene_pos, ] <- as.numeric(sample(
    data_per[gene_pos, ],
    length(data_per[gene_pos, ])
  ))
  return(data_per)
}

#' get_scm
#'
#' get_scm:calculate bn.fit results according to a certain neetwork structure and a cell type expression data.frame
#'
#' @param mem data.frame of cell type expression
#' @param graph bn.fit structure
#' @param id character, id should be given to label the information of sampling and permutation, best in the form of s1p1
#'
#' @return bn.fit results
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
  scm <- bnlearn::arcs(graph) %>% as.data.frame()
  graph_coef <- stats::coef(fit_result)
  scm[id] <- sapply(
    rownames(scm),
    function(x) graph_coef[[scm[x, "to"]]][scm[x, "from"]]
  )
  scm_result <- scm[, id, drop = FALSE]
  rownames(scm_result) <- paste(scm[["from"]], scm[["to"]], sep = "~")
  return(scm_result)
}
