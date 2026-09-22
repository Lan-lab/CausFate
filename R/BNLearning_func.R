#' Learn an undirected graph
#'
#' Use an undirected-graph algorithm to identify connections between nodes.
#'
#' @param mem A data matrix with genes in rows and cell types in columns.
#' @param param Numeric tuning parameter for the undirected-graph method,
#'   generally ranging from 0 to 1.
#' @param plot Logical; whether to plot the learned graph.
#' @param weight Logical; whether to display edge weights in the plot.
#' @param ugMethod Character string specifying the undirected-graph method.
#'
#' @return The result returned by the selected undirected-graph method.
UG_methods <- function(
  mem,
  param = 0.2,
  plot = FALSE,
  weight = FALSE,
  ugMethod = "cmi2ni"
) {
  if (ugMethod == "GeneNet") {
    ug <- network.GeneNet(expr.data = (mem), fdr = param)
  } else if (ugMethod == "ns") {
    ug <- network.ns(expr.data = (mem), alpha = param)
  } else if (ugMethod == "glasso") {
    ug <- network.glasso(expr.data = (mem), lambda = param)
  } else if (ugMethod == "glassosf") {
    ug <- network.glassosf(expr.data = (mem), alpha = param)
  } else if (ugMethod == "pcacmi") {
    ug <- network.pcacmi(expr.data = (mem), lambda = param)
  } else if (ugMethod == "cmi2ni") {
    ug <- network.cmi2ni(expr.data = (mem), lambda = param)
  } else if (ugMethod == "space") {
    ug <- network.space(expr.data = (mem), alpha = param)
  } else if (ugMethod == "bayesianglasso") {
    ug <- network.bayesianglasso(expr.data = (mem), prob = param)
  }

  if (plot) {
    ugPlot(ug, main = paste0(ugMethod, ", para = ", param), weight = weight)
  }
  return(ug)
}


#' Learn directed graphs over a parameter grid
#'
#' Learn a collection of network structures for all values in `params` using
#' the selected undirected- and directed-graph methods.
#'
#' @param mem A data matrix with genes in rows and cell types in columns.
#' @param params Numeric vector of tuning parameters for the undirected-graph
#'   method, generally ranging from 0 to 1.
#' @param whiteList A data frame or character vector specifying arcs that must
#'   be present in the network.
#' @param blackList A data frame or character vector specifying arcs that must
#'   not be present in the network.
#' @param root Character vector of root nodes, which cannot have incoming arcs.
#' @param ugMethod Character string specifying the undirected-graph method.
#' @param dagMethod Character string specifying the method used to determine arc
#'   directions.
#' @param ncores Number of cores available for computation.
#'
#' @return A list of network structures, one for each parameter value.
DG_grid <- function(
  mem,
  params = c(1:100) / 100,
  whiteList = NULL,
  blackList = NULL,
  root = NULL,
  ugMethod = "cmi2ni",
  dagMethod = "hc",
  ncores = 1
) {
  # doParallel::registerDoParallel(cores = ncores)
  dags <- foreach::foreach(param = params, .combine = c) %do%
    {
      tmp <- list()
      tmp[[1]] <- DG_methods(
        mem = (mem),
        param = param,
        root = root,
        whiteList = whiteList,
        blackList = blackList,
        ugMethod = ugMethod,
        dagMethod = dagMethod
      )
      tmp
    }
  names(dags) <- paste0("param_", params)
  dags
}

#' Learn a directed graph
#'
#' Learn a directed network using one tuning parameter and the selected graph
#' learning methods.
#'
#' @param mem A data matrix with genes in rows and cell types in columns.
#' @param param Numeric tuning parameter for the undirected-graph method.
#' @param root Character vector of root nodes, which cannot have incoming arcs.
#' @param whiteList A data frame or character vector specifying arcs that must
#'   be present in the network.
#' @param blackList A data frame or character vector specifying arcs that must
#'   not be present in the network.
#' @param plot Logical; whether to plot the learned network.
#' @param ugMethod Character string specifying the undirected-graph method.
#' @param dagMethod Character string specifying the method used to determine arc
#'   directions.
#'
#' @return The network structure learned with the specified parameter and
#'   methods.
DG_methods <- function(
  mem,
  param = 0.2,
  root = NULL,
  whiteList = NULL,
  blackList = NULL,
  plot = FALSE,
  ugMethod = "cmi2ni",
  dagMethod = "hc"
) {
  "Edges in whiteList and blackList should be in the form of node1~node2"
  ug <- UG_methods((mem), param = param, ugMethod = ugMethod)
  nodes <- colnames((mem))
  blackList <- genBList(ug, nodes = nodes, root = root, blackList = blackList)
  whiteList <- parse_edge_list(whiteList)
  dag <- get(dagMethod)(
    as.data.frame(mem),
    blacklist = blackList,
    whitelist = whiteList
  )
  main <- paste0(
    "method = ",
    ugMethod,
    ", ng = ",
    nrow((mem)),
    ", param = ",
    param
  )
  if (plot) {
    bnlearn::graphviz.plot(dag, main = main, shape = "circle")
  }
  dag
}


#' Learn directed graphs from single-cell samples
#'
#' Learn network structures across parameter values and cell samples from a
#' Seurat object.
#'
#' @param dat A Seurat object.
#' @param frac Fraction of cells included in each sample.
#' @param N_smpl Number of cell samples.
#' @param params Numeric vector of tuning parameters for the undirected-graph
#'   method.
#' @param whiteList A data frame or character vector specifying arcs that must
#'   be present in the network.
#' @param blackList A data frame or character vector specifying arcs that must
#'   not be present in the network.
#' @param root Character vector of root nodes, which cannot have incoming arcs.
#' @param ugMethod Character string specifying the undirected-graph method.
#' @param dagMethod Character string specifying the method used to determine arc
#'   directions.
#' @param ncores Number of cores used for computation.
#' @param seed Numeric random seed.
#'
#' @return A list of network structures for each parameter value and cell
#'   sample.
DG_smpl <- function(
  dat,
  frac = 0.3,
  N_smpl = 10,
  params = c(1:100) / 100,
  whiteList = NULL,
  blackList = NULL,
  root = NULL,
  ugMethod = "cmi2ni",
  dagMethod = "hc",
  ncores = 64,
  seed = Sys.time()
) {
  set.seed(as.numeric(seed))
  meta <- dat@meta.data %>% tibble::rownames_to_column("CellID")
  ctypes <- names(table(meta$celltype))
  random_number <- stats::runif(N_smpl, 1, as.numeric(seed))

  doParallel::registerDoParallel(cores = ncores)
  dags.smpl <- foreach::foreach(i_smpl = 1:N_smpl, .combine = c) %dopar%
    {
      set.seed(random_number[i_smpl])
      smpl.ttl <- 50
      meta.smpl <- tibble::tibble()
      while (TRUE) {
        meta.smpl <- meta %>%
          dplyr::group_by(celltype) %>%
          dplyr::sample_frac(frac)
        if (length(unique(meta.smpl$celltype)) == length(ctypes)) {
          break
        }
        smpl.ttl <- smpl.ttl - 1
        if (smpl.ttl == 0) {
          stop(
            "Failed sampling attempts reach limit: consider if the sampling fraction is too low or wrong dataset!"
          )
        }
      }
      dat.smpl <- subset(dat, cells = meta.smpl$CellID)
      meta_tmp <- dat.smpl$celltype
      tmp <- data.frame(dat.smpl@assays$RNA@data)
      mem.smpl <- data.frame(matrix(0, nrow = nrow(tmp), ncol = length(ctypes)))
      colnames(mem.smpl) <- ctypes
      rownames(mem.smpl) <- rownames(tmp)
      for (c_t in ctypes) {
        mem.smpl[c_t] <- tmp[, meta_tmp == c_t, drop = FALSE] %>%
          rowMeans()
      }
      dags <- list()
      dags[[1]] <- DG_grid(
        mem.smpl,
        params = params,
        whiteList = whiteList,
        blackList = blackList,
        root = root,
        ugMethod = ugMethod,
        dagMethod = dagMethod,
        ncores = ncores
      )
      dags
    }
  dags.smpl
}

#' Learn Bayesian network structures
#'
#' Learn network structures across parameter values, with optional cell
#' sampling for single-cell data.
#'
#' @param dat A Seurat object in single-cell mode or an expression matrix in
#'   bulk mode.
#' @param frac Fraction of cells included in each sample.
#' @param N_smpl Number of cell samples.
#' @param params Numeric vector of tuning parameters for the undirected-graph
#'   method.
#' @param whiteList A data frame or character vector specifying arcs that must
#'   be present in the network.
#' @param blackList A data frame or character vector specifying arcs that must
#'   not be present in the network.
#' @param root Character vector of root nodes, which cannot have incoming arcs.
#' @param ugMethod Character string specifying the undirected-graph method.
#' @param dagMethod Character string specifying the method used to determine arc
#'   directions.
#' @param ncores Number of cores used for computation.
#' @param mode Character string specifying `"single_cell"` or `"bulk"` mode.
#' @param seed Numeric random seed.
#'
#' @return A list of learned network structures.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @import bnlearn
#'
BNLearning <- function(
  dat,
  frac = 0.3,
  N_smpl = 10,
  params = c(1:100) / 100,
  whiteList = NULL,
  blackList = NULL,
  root = NULL,
  ugMethod = "cmi2ni",
  dagMethod = "hc",
  ncores = 1,
  mode = "single_cell",
  seed = Sys.time()
) {
  if (mode == "single_cell") {
    BNLearn_result <- DG_smpl(
      dat,
      frac,
      N_smpl,
      params,
      whiteList,
      blackList,
      root,
      ugMethod,
      dagMethod,
      ncores,
      seed
    )
  } else if (mode == "bulk") {
    BNLearn_result <- DG_grid(
      dat,
      params,
      whiteList,
      blackList,
      root,
      ugMethod,
      dagMethod,
      ncores
    )
  } else {
    BNLearn_result <- NULL
    print("Wrong mode!")
  }
  BNLearn_result
}

#' Learn the principal DAG from data
#'
#' `learnDAG` uses `BNLearning` to generate a set of possible DAGs representing
#' the differentiation network and analyses these graphs to give the final DAG.
#'
#' @param dat A Seurat object in single-cell mode or an expression matrix in
#'   bulk mode.
#' @param frac Fraction of cells included in each sample.
#' @param N_smpl Number of cell samples.
#' @param params Numeric vector of tuning parameters for the undirected-graph
#'   method.
#' @param whiteList A data frame or character vector specifying arcs that must
#'   be present in the network.
#' @param blackList A data frame or character vector specifying arcs that must
#'   not be present in the network.
#' @param root Character vector of root nodes, which cannot have incoming arcs.
#' @param ugMethod Character string specifying the undirected-graph method.
#' @param dagMethod Character string specifying the method used to determine arc
#'   directions.
#' @param ncores Number of cores used for computation.
#' @param mode Character string specifying `"single_cell"` or `"bulk"` mode.
#' @param Emin Minimum number of arcs required for a DAG to be retained.
#' @param Emax Maximum number of arcs allowed for a DAG to be retained.
#'
#' @return A `bn` object representing the principal directed acyclic graph.
#'
#' @export
#'
learnDAG <- function(
  dat,
  frac = 0.3,
  N_smpl = 10,
  params = c(1:100) / 100,
  whiteList = NULL,
  blackList = NULL,
  root = NULL,
  ugMethod = "cmi2ni",
  dagMethod = "hc",
  ncores = 1,
  mode = "single_cell",
  Emin = 1,
  Emax = 50
) {
  DAGs <- BNLearning(
    dat,
    frac = frac,
    N_smpl = N_smpl,
    params = params,
    whiteList = whiteList,
    blackList = blackList,
    root = root,
    ugMethod = ugMethod,
    dagMethod = dagMethod,
    ncores = ncores,
    mode = mode
  )
  if (mode == "single_cell") {
    ctypes <- names(table(dat$celltype))
    # Old method:
    # tmp <- lapply(DAGs,function(x){df2mat(combineDAGs(x,Emin = 1,Emax = 50))})
    # b = norm_mat(tmp[[1]][ctypes,ctypes]);for (i in 2:length(tmp)) {b = b + norm_mat(tmp[[i]][ctypes,ctypes])}
    # print(b)
    # return(b)

    # crudeUG <-foreach::foreach(DAG=DAGs,.combine=rbind)%dopar%{rmCyc_sorted(combineDAGs(DAG, Emin = Emin, Emax = Emax))}
    # crudeUG <- crudeUG%>%
    #   dplyr::group_by(from, to) %>%
    #   dplyr::summarise(freq = sum(freq)) %>%
    #   as.data.frame() %>%
    #   dplyr::mutate(edge = paste0(from, "~", to)) %>%
    #   tibble::remove_rownames() %>%
    #   tibble::column_to_rownames("edge")

    crudeUG <- combineDAGsmpl(DAGs, Emin = Emin, Emax = Emax)
    crudeDAG <- rmCyc_sorted(crudeUG)
  } else if (mode == "bulk") {
    ctypes <- colnames(dat)
    crudeUG <- combineDAGs(DAGs, Emin = Emin, Emax = Emax)
    crudeDAG <- rmCyc_sorted(crudeUG)
  } else {
    error("Invalid mode!")
  }
  outputDAG <- bnlearn::empty.graph(ctypes)
  strc1 <- crudeDAG[, 1:2]
  strc1[, 1] <- as.character(strc1[, 1])
  strc1[, 2] <- as.character(strc1[, 2])
  bnlearn::arcs(outputDAG) <- strc1
  graphviz.plot(outputDAG, shape = "ellipse")
  return(outputDAG)
}

#' Combine DAGs and remove cycles
#'
#' Combine learned DAGs and remove cycles from the resulting network.
#'
#' @param DAGs A list of learned DAGs.
#' @param mode Character string specifying `"single_cell"` or `"bulk"` mode.
#' @param Emin Minimum number of arcs required for a DAG to be retained.
#' @param Emax Maximum number of arcs allowed for a DAG to be retained.
#'
#' @return A `bn` object representing the combined acyclic graph.
#' @export
combineAndRmCyc <- function(
  DAGs,
  mode = "single_cell",
  Emin = NULL,
  Emax = NULL
) {
  if (mode == "single_cell") {
    ctypes <- names(table(dat$celltype))
    # Old method:
    # tmp <- lapply(DAGs,function(x){df2mat(combineDAGs(x,Emin = 1,Emax = 50))})
    # b = norm_mat(tmp[[1]][ctypes,ctypes]);for (i in 2:length(tmp)) {b = b + norm_mat(tmp[[i]][ctypes,ctypes])}
    # print(b)
    # return(b)

    crudeUG <- foreach::foreach(DAG = DAGs, .combine = rbind) %dopar%
      {
        rmCyc_sorted(combineDAGs(DAG, Emin = Emin, Emax = Emax))
      }
    crudeUG <- crudeUG %>%
      dplyr::group_by(from, to) %>%
      dplyr::summarise(freq = sum(freq)) %>%
      as.data.frame() %>%
      dplyr::mutate(edge = paste0(from, "~", to)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("edge")

    # crudeUG <- combineDAGsmpl(DAGs, Emin = Emin,Emax=Emax)
    crudeDAG <- rmCyc_sorted(crudeUG)
  } else if (mode == "bulk") {
    ctypes <- colnames(dat)
    crudeUG <- combineDAGs(DAGs, Emin = Emin, Emax = Emax)
    crudeDAG <- rmCyc_sorted(crudeUG)
  } else {
    error("Invalid mode!")
  }
  outputDAG <- bnlearn::empty.graph(ctypes)
  strc1 <- crudeDAG[, 1:2]
  strc1[, 1] <- as.character(strc1[, 1])
  strc1[, 2] <- as.character(strc1[, 2])
  bnlearn::arcs(outputDAG) <- strc1
  graphviz.plot(outputDAG, shape = "ellipse")
  return(outputDAG)
}

#' Remove cycles in descending edge-weight order
#'
#' @param dS A summarized network structure.
#'
#' @return A summarized acyclic network structure.
#' @keywords internal
rmCyc_sorted <- function(dS) {
  dS <- rmCyc2(dS)
  ctypes <- unique(union(dS[, "from"], dS[, "to"]))
  dS <- dS[order(dS[, "freq"], decreasing = TRUE), ]
  edge_tune <- logical(nrow(dS))
  tmp <- data.frame(from = character(), to = character(), freq = numeric())
  for (e in 1:nrow(dS)) {
    if ((hasCyc(rbind(tmp, dS[e, ])))) {
      edge_tune[e] <- FALSE
    } else {
      tmp <- rbind(tmp, dS[e, ])
      edge_tune[e] <- TRUE
    }
  }
  dS[edge_tune, ]
}

#' Test a network for cycles
#'
#' @param dS A summarized network structure.
#'
#' @return `TRUE` if the network contains a cycle; otherwise, `FALSE`.
#' @keywords internal
hasCyc <- function(dS) {
  Es <- dS[c("from", "to")]
  Es[, 1] <- as.character(Es[, 1])
  Es[, 2] <- as.character(Es[, 2])
  Ys <- Es %>%
    as.matrix() %>%
    (graph::ftM2graphNEL) %>%
    (Rgraphviz::layoutGraph) %>%
    (graph::nodeRenderInfo) %>%
    .$nodeY %>%
    sort() %>%
    rev()
  coords <- Es %>% dplyr::mutate(from = Ys[from], to = Ys[to])
  return(!all(coords$from > coords$to))
}

#' Generate bootstrap indices
#'
#' Create sampling indices for `run_diffCoeff()`.
#'
#' @param meta A character or factor vector of cell-type labels.
#' @param bootstrap_times Number of bootstrap samples.
#' @param ratio Fraction of cells included in each sample.
#'
#' @return A list of sampling-index vectors.
#' @export
#'
bootstrap_index <- function(meta, bootstrap_times, ratio) {
  celltype <- names(table(meta))
  index <- list(NULL)
  for (j in 1:bootstrap_times) {
    index[[j]] <- NA
    for (i in 1:length(table(meta))) {
      tmp <- which(meta == celltype[i])
      index[[j]] <- c(index[[j]], sample(tmp, length(tmp) * ratio, replace = F))
    }
    index[[j]] <- index[[j]][-1]
  }

  return(index)
}

#' Combine DAGs across cell samples
#'
#' Combine the network structures in `dags.smpl` after filtering them by their
#' numbers of arcs.
#'
#' @param dags.smpl A list of network structures returned by `DG_smpl()`.
#' @param Emin Minimum number of arcs required for a network to be retained.
#' @param Emax Maximum number of arcs allowed for a network to be retained.
#' @param ncores Number of cores used for computation.
#' @param model_averaging Character string specifying the model-averaging
#'   strategy: `"joint"` or `"two-tier"`.
#'
#' @return A data frame containing the summarized network structure.
#' @export
combineDAGsmpl <- function(
    dags.smpl,
    Emin = NULL,
    Emax = NULL,
    ncores = 64,
    model_averaging = c("joint", "two-tier")
) {
  
  model_averaging <- match.arg(model_averaging)
  
  doParallel::registerDoParallel(cores = ncores)
  
  if (model_averaging == "joint") {
    
    Efreqs.tmp <- foreach::foreach(
      i_smpl = 1:length(dags.smpl),
      .combine = rbind
    ) %dopar% {
      combineDAGs(
        dags.smpl[[i_smpl]],
        Emin = Emin,
        Emax = Emax
      )
    }
    
  } else if (model_averaging == "two-tier") {
    
    Efreqs.tmp <- foreach::foreach(
      i_smpl = 1:length(dags.smpl),
      .combine = rbind
    ) %dopar% {
      
      x <- combineDAGs(
        dags.smpl[[i_smpl]],
        Emin = Emin,
        Emax = Emax
      )
      
      if (is.null(x) ||
          nrow(x) == 0 ||
          !all(c("from", "to", "freq") %in% colnames(x))) {
        return(NULL)
      }
      
      x <- x[
        x$freq != 0,
        ,
        drop = FALSE
      ]
      
      if (nrow(x) == 0) {
        return(NULL)
      }
      
      rmCyc(x)
    }
  }
  
  Efreqs.tmp %>%
    dplyr::group_by(from, to) %>%
    dplyr::summarise(freq = sum(freq), .groups = "drop") %>%
    as.data.frame() %>%
    dplyr::mutate(edge = paste0(from, "~", to)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("edge")
}


#' Combine DAGs
#'
#' Combine network structures after filtering them by their numbers of arcs.
#'
#' @param dags A list of network structures returned by `DG_grid()`.
#' @param Emin Minimum number of arcs required for a network to be retained.
#' @param Emax Maximum number of arcs allowed for a network to be retained.
#'
#' @return A data frame containing the summarized network structure.
#' @export
combineDAGs <- function(dags, Emin = NULL, Emax = NULL) {
  Efreqs <- foreach::foreach(i = 1:length(dags), .combine = rbind) %do%
    {
      bool1 <- ifelse(is.null(Emin), TRUE, bnlearn::narcs(dags[[i]]) >= Emin)
      bool2 <- ifelse(is.null(Emax), TRUE, bnlearn::narcs(dags[[i]]) <= Emax)
      bool <- bool1 && bool2
      # bool <- (bnlearn::narcs(dags[[i]]) >= Emin) && (bnlearn::narcs(dags[[i]]) <= Emax)
      if (isTRUE(bool)) {
        Efreq <- bnlearn::arcs(dags[[i]]) %>%
          as.data.frame() %>%
          dplyr::mutate(freq = 1)
        Efreq
      } else {
        NULL
      }
    }
  if (is.null(Efreqs)) {
    return(data.frame(from = character(), to = character(), freq = numeric()))
  }
  Efreqs %>%
    dplyr::group_by(from, to) %>%
    dplyr::summarise(freq = sum(freq)) %>%
    as.data.frame() %>%
    dplyr::mutate(edge = paste0(from, "~", to)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("edge")
}

#' Trim a DAG using fitted coefficients
#'
#' Modify a network structure according to fitted `bn.fit` coefficients.
#'
#' @param dat_tmp A matrix of gene expression values across cell types.
#' @param e A `bn.fit` network structure.
#' @param min_arc Minimum number of arcs pointing to a node.
#' @param max_arc Maximum number of arcs pointing to a node.
#' @param threshold_value Numeric threshold used to filter arcs according to
#'   fitted coefficients.
#' @param plot Logical; whether to plot the trimmed DAG.
#'
#' @return The filtered network structure.
#' @export
trimDAG <- function(
  dat_tmp,
  e,
  min_arc = 2,
  max_arc = 4,
  threshold_value = 0.9,
  plot = TRUE
) {
  ref <- bnlearn::bn.fit(e, dat_tmp)
  celltype <- names(ref)
  arc_modified <- bnlearn::arcs(ref)

  # cutoff & min
  for (i in celltype) {
    cut_off <- sum(abs(unlist(ref[[i]][4])[-1])) * (1 - threshold_value)
    if (length(unlist(ref[[i]][4])[-1]) == 0) {
      i <- i
    } else if (length(unlist(ref[[i]][4])) < (min_arc + 2)) {
      i <- i
    } else {
      for (j in 1:length(unlist(ref[[i]][4])[-1])) {
        if (abs(unlist(ref[[i]][4])[-1])[j] < cut_off) {
          tmp <- strsplit(
            names(unlist(ref[[i]][4])[-1])[j],
            split = "coefficients."
          )[[1]][2]
          target <- intersect(
            which(arc_modified[, "to"] == i),
            which(arc_modified[, "from"] == tmp)
          )
          arc_modified <- arc_modified[-target, ]
        }
      }
    }
  }
  bnlearn::arcs(e) <- arc_modified
  ref <- bnlearn::bn.fit(e, dat_tmp)

  # max
  for (i in celltype) {
    if (length(unlist(ref[[i]][4])[-1]) > max_arc) {
      cut_off <- as.numeric(sort(abs(unlist(ref[[i]][4])[-1]), decreasing = T)[
        max_arc
      ])
      for (j in 1:length(unlist(ref[[i]][4])[-1])) {
        if (abs(unlist(ref[[i]][4])[-1])[j] < cut_off) {
          tmp <- strsplit(
            names(unlist(ref[[i]][4])[-1])[j],
            split = "coefficients."
          )[[1]][2]
          target <- intersect(
            which(arc_modified[, "to"] == i),
            which(arc_modified[, "from"] == tmp)
          )
          arc_modified <- arc_modified[-target, ]
        }
      }
    }
  }
  bnlearn::arcs(e) <- arc_modified
  if (plot) {
    bnlearn::graphviz.plot(e, shape = "ellipse")
  }
  return(e)
}

#' Remove cycles using graph layout
#'
#' Remove cycles from a network according to its hierarchical layout.
#'
#' @param dS A summarized network structure.
#'
#' @return A summarized network structure with cycles removed.
#' @export
rmCycL <- function(dS) {
  Es <- dS[c("from", "to")]
  Es[, 1] <- as.character(Es[, 1])
  Es[, 2] <- as.character(Es[, 2])
  # gp <- Es %>% as.matrix %>% graph::ftM2graphNEL
  # lgp <- Rgraphviz::layoutGraph(gp)
  # Ys <- graph::nodeRenderInfo(lgp)$nodeY
  Ys <- Es %>%
    as.matrix() %>%
    (graph::ftM2graphNEL) %>%
    (Rgraphviz::layoutGraph) %>%
    (graph::nodeRenderInfo) %>%
    .$nodeY %>%
    sort() %>%
    rev()
  coords <- Es %>% dplyr::mutate(from = Ys[from], to = Ys[to])
  dS.fltr <- dS[coords$from > coords$to, ]
  dS.fltr
}

#' Reduce undirected edges to directed edges
#'
#' @param dS A summarized network structure.
#'
#' @return A summarized directed network structure.
rmCyc2 <- function(dS) {
  dS <- dS %>% df2mat()
  dS[dS - t(dS) < 0] <- 0
  dS %>% mat2df()
}

#' Reduce undirected edges and normalize edge weight
#'
#' @param N A numeric adjacency matrix.
#'
#' @return A normalized adjacency matrix.
norm_mat <- function(N) {
  M <- N
  for (c in 1:nrow(M)) {
    for (d in 1:ncol(M)) {
      if (c < d) {
        if (M[c, d] + M[d, c] > 0) {
          if (M[c, d] > M[d, c]) {
            M[c, d] <- 1
            M[d, c] <- 0
          } else {
            M[c, d] <- 0
            M[d, c] <- 1
          }
        }
      }
    }
  }
  return(M)
}

#' Remove cycles
#'
#' Remove cycles according to the hierarchical structure and network matrix.
#'
#' @param dS A summarized network structure.
#' @param method Character string specifying the cycle-removal method.
#'
#' @return A summarized network structure with cycles removed.
#' @export
rmCyc <- function(dS, method = "sorted") {
  if (method == "sorted") {
    return(rmCyc_sorted(dS))
  }
  tmp <- dS %>% rmCyc2()
  dS.fltr <- tmp %>% rmCycL()
  edge_rev <- dS.fltr %>%
    dplyr::mutate(edge_rev = paste0(to, "~", from)) %>%
    .$edge_rev
  edge_tune <- setdiff(rownames(dS), rownames(tmp)) %>%
    setdiff(., edge_rev) %>%
    union(., rownames(dS.fltr))
  dS.tune <- dS[edge_tune, ]
  dS.tune
}


#' Generate edge set
#'
#' @param fromSet Character vector of starting nodes.
#' @param toSet Character vector of ending nodes.
#' @param sep Character string used to separate node names.
#'
#' @return A character vector of edge names.
#' @export
setEdges <- function(fromSet, toSet, sep = "~") {
  edges <- paste(rep(fromSet, each = length(toSet)), toSet, sep = sep)
  selfs <- paste(fromSet, fromSet, sep = sep)
  return(setdiff(edges, selfs))
}


#' Summarize a gene-expression matrix
#'
#' Summarize a single-cell expression matrix into a cell-type expression matrix.
#'
#' @param gem A data frame or matrix of single-cell expression values.
#' @param meta A character or factor vector of cell-type labels.
#' @param FUN Character string specifying `"mean"` or `"median"`.
#'
#' @return A data frame of summarized expression values.
#' @export
gem2mem <- function(gem = NULL, meta = NULL, FUN = c("mean", "median")) {
  if (is.null(meta)) {
    warning(
      "Celltypes information is missing, MEM is identical to input GEM ..."
    )
    mem <- gem
  } else {
    celltype <- names(table(meta))
    mem <- data.frame(matrix(0, nrow = nrow(gem), ncol = length(celltype)))
    colnames(mem) <- celltype
    rownames(mem) <- rownames(gem)
    for (c in celltype) {
      if (FUN == "mean") {
        mem[, c] <- rowMeans(gem[, meta == c])
      } else if (FUN == "median") {
        mem[, c] <- matrixStats::rowMedians(gem[, meta == c])
      }
    }
  }
  return(mem)
}


#' Generate blacklist
#'
#' @param ug An undirected graph represented as a data frame.
#' @param nodes Character vector of node names.
#' @param root Character vector of root nodes.
#' @param blackList An optional data frame or character vector of excluded arcs.
#'
#' @return A data frame containing excluded arcs.
genBList <- function(ug, nodes, root = NULL, blackList = NULL) {
  edgeAll <- setEdges(nodes, nodes, sep = "~")
  edgeCor <- paste0(c(ug$node1, ug$node2), "~", c(ug$node2, ug$node1))
  blacklist <- setdiff(edgeAll, edgeCor) %>%
    union(., setEdges(nodes, root)) %>%
    as.data.frame() %>%
    rlang::set_names("edge") %>%
    tidyr::separate(edge, c("from", "to"), sep = "~")
  blacklist <- rbind(blacklist, parse_edge_list(blackList))
  blacklist
}

#' Combine values into a unique vector
#'
#' @param ... Vectors to combine.
#'
#' @return A vector of unique values.
#' @keywords internal
alter <- function(...) {
  vec <- unique(c(rbind(...)))
  return(vec)
}

#' Get N edges
#'
#' @param dags A list of network structures.
#'
#' @return A named numeric vector containing the number of arcs in each
#'   network.
getNEdges <- function(dags) {
  nEdges <- sapply(
    1:length(dags),
    function(x) bnlearn::narcs(dags[[x]])
  ) %>%
    rlang::set_names(names(dags))
  nEdges
}

#' Convert a network data frame to a matrix
#'
#' Transform a summarized network structure into an adjacency matrix.
#'
#' @param DAG A summarized network structure.
#' @param ctypes Optional character vector of cell types.
#' @export
#'
#' @return An adjacency matrix corresponding to the summarized network.
df2mat <- function(DAG, ctypes = NULL) {
  if (class(DAG)[[1]] != "data.frame") {
    DAG <- DAG %>%
      as.data.frame() %>%
      mutate(freq = 1)
  }
  DAG.dm <- tidyr::spread(DAG, key = to, value = freq) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("from") %>%
    dplyr::mutate(dplyr::across(
      tidyselect::everything(),
      .fns = ~ replace_na(., 0)
    ))

  if (is.null(ctypes)) {
    ctypes <- union(rownames(DAG.dm), colnames(DAG.dm))
  }
  DAG.dm[setdiff(ctypes, rownames(DAG.dm)), ] <- 0
  DAG.dm[, setdiff(ctypes, colnames(DAG.dm))] <- 0
  DAG.dm <- DAG.dm[ctypes, ctypes] %>% as.matrix()
  DAG.dm
}

#' Convert a network matrix to a data frame
#'
#' Transform an adjacency matrix into a summarized network structure.
#'
#' @param DAG.dm A matrix representing a network structure.
#' @export
#'
#' @return A data frame containing the summarized network structure.
mat2df <- function(DAG.dm) {
  DAG <- DAG.dm %>%
    t() %>%
    as.table() %>%
    as.data.frame() %>%
    rlang::set_names(c("to", "from", "freq")) %>%
    subset(freq != 0) %>%
    .[, c("from", "to", "freq")] %>%
    dplyr::mutate(
      from = as.character(from),
      to = as.character(to),
      edge = paste0(from, "~", to)
    ) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("edge")
  DAG
}

#' Convert a graph data frame to a Bayesian network
#'
#' @param dS A data frame with columns named `from` and `to`.
#' @param node_names Character vector of node names.
#' @param plot Logical; whether to plot the network.
#'
#' @return A `bn` object.
#' @export
df2bn <- function(dS, node_names, plot = TRUE) {
  e <- bnlearn::empty.graph(node_names)
  arcs(e) <- dS[, c("from", "to")]
  if (plot) {
    graphviz.plot(e, shape = "ellipse")
  }
  return(e)
}

#' Extract bn coefficients
#'
#' @param dag A `bn` object.
#' @param data A data frame or matrix used to fit the network.
#'
#' @return A data frame containing arcs and their absolute coefficients.
#' @export
getCoef <- function(dag, data) {
  fit_res <- bnlearn::bn.fit(dag, data)
  Efreq <- bnlearn::arcs(dag) %>%
    as.data.frame() %>%
    dplyr::mutate(freq = 0)
  for (i in 1:nrow(Efreq)) {
    Efreq[i, "freq"] <- abs(fit_res[[Efreq[i, "to"]]][["coefficients"]][[Efreq[
      i,
      "from"
    ]]])
  }
  return(Efreq)
}

#' Filter edges by bn coefficients
#'
#' @param dag A `bn` object.
#' @param data A data frame or matrix used to fit the network.
#' @param ctypes Character vector of node names.
#' @param threshold Minimum absolute coefficient required to retain an edge.
#'
#' @return A filtered `bn` object.
#' @export
edgeFilter <- function(dag, data, ctypes, threshold = .1) {
  return(
    dag %>%
      getCoef(data) %>%
      .[.[, "freq"] > threshold, 1:2] %>%
      df2bn(ctypes, plot = F)
  )
}

#' Parse a whitelist or blacklist
#'
#' @param edges A data frame or character vector of edges.
#'
#' @return A data frame with `from` and `to` columns, or `NULL`.
#' @keywords internal
parse_edge_list <- function(edges) {
  if (is.null(edges)) {
    return(NULL)
  }

  if (is.data.frame(edges)) {
    colnames_needed <- c("from", "to")

    if (!all(colnames_needed %in% names(edges))) {
      # Assume first two columns are from and to
      if (ncol(edges) < 2) {
        stop("Data frame must have at least two columns.")
      }
      edges <- edges[, 1:2]
      names(edges) <- colnames_needed
    }

    return(edges)
  } else if (is.character(edges)) {
    # Expecting strings like "A~B"
    split_edges <- strsplit(edges, "~")
    if (any(sapply(split_edges, length) != 2)) {
      stop(
        "Each string must contain exactly one '~' character separating 'from' and 'to'."
      )
    }
    from_to <- do.call(rbind, split_edges)
    return(data.frame(
      from = from_to[, 1],
      to = from_to[, 2],
      stringsAsFactors = FALSE
    ))
  } else {
    stop("Input must be either a data.frame or a character vector.")
  }
}
