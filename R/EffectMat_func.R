#' Calculate effect matrix
#'
#' When function `PerturbResult` is done, the ouput of `PerturbResult` can be the input of this function to produce effect matrix
#'
#' @param diffCoeff list of diffCoeff raw result
#' @param mode character, determine the diffCoeff function, usually "mean"
#'
#' @return Effect matrix
#' @export
EffectMatrix <- function(diffCoeff, mode = "mean") {
  diffCoeff_result <- list(NULL)
  n_sample <- diffCoeff$n_sample
  n_permutation <- diffCoeff$n_permutation
  if (mode == "mean") {
    k <- 1
    for (i in 1:(n_sample)) {
      for (j in 1:(n_permutation)) {
        diffCoeff_result[[k]] <- get_diffCoeff(
          scm.inter = diffCoeff$perturb[[k]],
          scm.ref = diffCoeff$ref[[i]],
          mode = mode
        )
        k <- k + 1
      }
    }
    names(diffCoeff_result) <- names(diffCoeff$perturb)
    tmp <- as.matrix(diffCoeff_result[[1]])
    if (n_sample * n_permutation > 1) {
      for (k in 2:(n_sample * n_permutation)) {
        tmp <- tmp + as.matrix(diffCoeff_result[[k]])
      }
    }
    diffCoeff_result <- tmp
  } else {
    stop("Invalid diffCoeff calculation mode!")
  }

  return(as.data.frame(diffCoeff_result))
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
  replace = F
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
