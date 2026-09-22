str2arc <- function(arc_list, sep = "~") {
  if (is.null(arc_list)) return(data.frame(from = character(), to = character()))
  do.call(rbind, lapply(strsplit(arc_list, sep, fixed = TRUE), function(x) {
    if (length(x) != 2L) stop("Invalid arc specification: ", paste(x, collapse = sep))
    data.frame(
      from = x[1],
      to = strsplit(x[2], "", fixed = TRUE)[[1]],
      stringsAsFactors = FALSE
    )
  }))
}

generateBayesianParams <- function(ref_net, intercept_coef = 5, parent_coef = 2, sd = 1) {
  node_names <- bnlearn::nodes(ref_net)
  out <- lapply(node_names, function(node) {
    parents <- bnlearn::parents(ref_net, node)
    coef <- c(
      as.integer(stats::runif(1) * 10 * intercept_coef) / 10,
      vapply(
        parents,
        function(x) as.integer(stats::runif(1) * 10 * parent_coef) / 10,
        numeric(1)
      )
    )
    names(coef) <- c("(Intercept)", parents)
    list(coef = coef, sd = as.integer(stats::runif(1) * 10 * sd) / 10)
  })
  names(out) <- node_names
  out
}

generateExpressionSeurat <- function(
    ref_net,
    n_feature = 500,
    n_node_smpl = 50,
    n_edge_smpl = 100,
    edge_span = 1.1
) {
  node_names <- bnlearn::nodes(ref_net)
  
  expressions <- lapply(seq_len(n_feature), function(i) {
    fitted <- bnlearn::custom.fit(
      ref_net,
      dist = generateBayesianParams(ref_net)
    )
    as.matrix(bnlearn::rbn(fitted, n = n_node_smpl))
  })
  
  data <- do.call(cbind, lapply(seq_len(n_node_smpl), function(i) {
    t(vapply(
      expressions,
      function(x) x[i, node_names],
      numeric(length(node_names))
    ))
  }))
  
  meta_nodes <- rep(node_names, times = n_node_smpl)
  meta <- meta_nodes
  ref_arcs <- bnlearn::arcs(ref_net)
  
  if (nrow(ref_arcs) > 0) {
    for (i in seq_len(nrow(ref_arcs))) {
      df_a <- data[
        ,
        sample(
          which(meta_nodes == ref_arcs[i, 1]),
          n_edge_smpl,
          replace = TRUE
        ),
        drop = FALSE
      ]
      
      df_b <- data[
        ,
        sample(
          which(meta_nodes == ref_arcs[i, 2]),
          n_edge_smpl,
          replace = TRUE
        ),
        drop = FALSE
      ]
      
      frac <- stats::runif(
        n_edge_smpl,
        min = -(edge_span - 1) / 2,
        max = 1 + (edge_span - 1) / 2
      )
      
      edge_data <- vapply(seq_len(n_edge_smpl), function(j) {
        df_a[, j] * frac[j] + df_b[, j] * (1 - frac[j])
      }, numeric(n_feature))
      
      data <- cbind(data, edge_data)
      meta <- c(
        meta,
        ifelse(frac < 0.5, ref_arcs[i, 1], ref_arcs[i, 2])
      )
    }
  }
  
  rownames(data) <- paste0("Gene_", seq_len(n_feature))
  colnames(data) <- paste0("Cell_", seq_len(ncol(data)))
  
  seu <- Seurat::CreateSeuratObject(
    counts = data,
    data = data,
    meta.data = data.frame(
      celltype = meta,
      row.names = colnames(data)
    )
  )
  
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
  seu[["RNA"]] <- methods::as(seu[["RNA"]], "Assay")
  seu
}

getSimulationArcs <- function(n_nodes, data_type) {
  if (data_type == "discrete") return(NULL)
  
  arcs <- list(
    sparse = list(
      `3` = c("A~B"),
      `4` = c("A~B"),
      `5` = c("A~B", "A~C"),
      `6` = c("A~B", "A~C", "B~D"),
      `7` = c("A~B", "A~C", "B~D", "C~E"),
      `8` = c("A~B", "A~C", "B~D", "C~E", "B~F"),
      `9` = c("A~B", "A~C", "B~D", "C~E", "B~F", "C~G")
    ),
    tree = list(
      `3` = c("A~B", "A~C"),
      `4` = c("A~B", "A~C", "B~D"),
      `5` = c("A~B", "A~C", "B~D", "C~E"),
      `6` = c("A~B", "A~C", "B~D", "C~E", "B~F"),
      `7` = c("A~B", "A~C", "B~D", "C~E", "B~F", "C~G"),
      `8` = c("A~B", "A~C", "B~D", "C~E", "B~F", "C~G", "E~H"),
      `9` = c("A~B", "A~C", "B~D", "C~E", "B~F", "C~G", "E~H", "G~I")
    ),
    dense = list(
      `3` = c(
        "A~B",
        "B~C",
        "A~C"
      ),
      `4` = c(
        "A~B",
        "B~C",
        "B~D",
        "A~C",
        "A~D"
      ),
      `5` = c(
        "A~B",
        "A~C",
        "B~D",
        "C~E",
        "A~D",
        "A~E"
      ),
      `6` = c(
        "A~BDEF",
        "A~C",
        "B~D",
        "C~E",
        "B~F"
      ),
      `7` = c(
        "A~B",
        "A~C",
        "B~DF",
        "C~EG",
        "A~DFEG"
      ),
      `8` = c(
        "A~B",
        "A~C",
        "B~D",
        "C~E",
        "B~F",
        "C~G",
        "E~H",
        "A~DEFG"
      ),
      `9` = c(
        "A~B",
        "A~C",
        "B~D",
        "C~E",
        "B~F",
        "C~G",
        "E~H",
        "G~I",
        "A~DEFG"
      )
    )
  )
  
  if (!data_type %in% names(arcs)) {
    stop("Unknown data_type: ", data_type)
  }
  
  out <- arcs[[data_type]][[as.character(n_nodes)]]
  
  if (is.null(out)) {
    stop(
      "No simulation structure for n_nodes = ",
      n_nodes,
      ", data_type = ",
      data_type
    )
  }
  
  out
}

generateSimulationData <- function(n_nodes, data_type, seed = 42) {
  node_names <- LETTERS[seq_len(n_nodes)]
  ref <- bnlearn::empty.graph(node_names)
  arc_short <- getSimulationArcs(n_nodes, data_type)
  
  if (!is.null(arc_short)) {
    bnlearn::arcs(ref) <- str2arc(arc_short)
  }
  
  set.seed(seed)
  edge_span <- if (data_type == "dense") 1.1 else 1.2
  
  list(
    Seurat_obj = generateExpressionSeurat(
      ref,
      edge_span = edge_span
    ),
    ref_net = ref
  )
}

getSimulationConfig <- function(n_nodes, data_type) {
  cfg <- list(
    params = seq(0.05, 0.15, 0.01),
    N_smpl = 20,
    freq_cutoff = NULL,
    trim_threshold = 0.8,
    max_arc = 2
  )
  cfg
}

getCachedDagSamples <- function(file, overwrite, run_fun) {
  if (file.exists(file) && !overwrite) {
    return(readRDS(file))
  }
  
  dag_smpl <- run_fun()
  saveRDS(dag_smpl, file)
  dag_smpl
}

buildCausFateNet <- function(
    dag_smpl,
    node_names,
    freq_cutoff = NULL
) {
  net <- causfate::combineDAGsmpl(
    dag_smpl,
    Emin = NULL,
    Emax = NULL,
    ncores = 1
  )
  
  if (is.null(net) || nrow(net) == 0) {
    return(bnlearn::empty.graph(node_names))
  }
  
  if (!is.null(freq_cutoff)) {
    mat <- causfate::df2mat(net)
    net <- causfate::mat2df(
      mat * (mat > freq_cutoff)
    )
  }
  
  if (is.null(net) || nrow(net) == 0) {
    return(bnlearn::empty.graph(node_names))
  }
  
  net <- causfate::rmCyc(net)
  out <- bnlearn::empty.graph(node_names)
  
  if (nrow(net) > 0) {
    edge <- net[, 1:2, drop = FALSE]
    edge[] <- lapply(edge, as.character)
    bnlearn::arcs(out) <- edge
  }
  
  out
}

benchmarkNets <- function(
    mem,
    root = NULL,
    methods = c("mmhc", "rsmax2", "h2pc"),
    seed = 42
) {
  blacklist <- if (is.null(root)) {
    NULL
  } else {
    data.frame(
      from = setdiff(colnames(mem), root),
      to = root,
      stringsAsFactors = FALSE
    )
  }
  
  out <- lapply(methods, function(method) {
    set.seed(seed)
    getExportedValue("bnlearn", method)(
      mem,
      blacklist = blacklist
    )
  })
  
  names(out) <- methods
  out
}

configureTetrad <- function(
    java_home = Sys.getenv("JAVA_HOME", unset = ""),
    heap = "2g"
) {
  if (nzchar(java_home)) {
    java_home <- normalizePath(
      path.expand(java_home),
      mustWork = TRUE
    )
    
    Sys.setenv(
      JAVA_HOME = java_home,
      PATH = paste(
        file.path(java_home, "bin"),
        Sys.getenv("PATH"),
        sep = .Platform$path.sep
      )
    )
  }
  
  for (pkg in c("causalDisco", "caugi", "rJava")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required.")
    }
  }
  
  status <- causalDisco::verify_tetrad()
  
  if (!isTRUE(status$installed)) {
    stop(
      "Tetrad is not installed. ",
      "Run causalDisco::install_tetrad() first."
    )
  }
  
  if (!isTRUE(status$java_ok)) {
    stop(
      "Java JDK 21 or later is required. Detected: ",
      status$java_version
    )
  }
  
  tetrad_dir <- Sys.getenv(
    "TETRAD_DIR",
    unset = ""
  )
  
  if (!nzchar(tetrad_dir)) {
    tetrad_dir <- getOption(
      "causalDisco.tetrad_cache",
      ""
    )
  }
  
  jars <- if (
    nzchar(tetrad_dir) &&
    dir.exists(tetrad_dir)
  ) {
    list.files(
      tetrad_dir,
      pattern = "^tetrad.*\\.jar$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  } else {
    character()
  }
  
  if (!length(jars)) {
    stop(
      "Cannot locate the Tetrad JAR. ",
      "Run causalDisco::install_tetrad() first."
    )
  }
  
  jars <- normalizePath(
    jars,
    mustWork = TRUE
  )
  
  if (!rJava::.jniInitialized) {
    rJava::.jinit(
      parameters = paste0("-Xmx", heap),
      classpath = jars
    )
  } else {
    current <- normalizePath(
      rJava::.jclassPath(),
      mustWork = FALSE
    )
    
    for (jar in setdiff(jars, current)) {
      rJava::.jaddClassPath(jar)
    }
  }
  
  rJava::.jfindClass(
    "edu/cmu/tetrad/data/Knowledge"
  )
  
  invisible(status)
}

pcalgAmatToBnlearn <- function(amat) {
  amat <- as.matrix(amat)
  labels <- colnames(amat)
  
  if (is.null(labels)) {
    stop("amat must contain column names.")
  }
  
  net <- bnlearn::empty.graph(labels)
  
  if (length(labels) < 2) {
    return(net)
  }
  
  for (i in seq_len(length(labels) - 1L)) {
    for (j in seq.int(i + 1L, length(labels))) {
      if (
        amat[i, j] == 0 &&
        amat[j, i] == 1
      ) {
        net <- bnlearn::set.arc(
          net,
          labels[i],
          labels[j],
          check.cycles = FALSE
        )
      }
      
      if (
        amat[i, j] == 1 &&
        amat[j, i] == 0
      ) {
        net <- bnlearn::set.arc(
          net,
          labels[j],
          labels[i],
          check.cycles = FALSE
        )
      }
      
      if (
        amat[i, j] == 1 &&
        amat[j, i] == 1
      ) {
        net <- bnlearn::set.edge(
          net,
          labels[i],
          labels[j],
          check.cycles = FALSE
        )
      }
    }
  }
  
  net
}

.validateBenchmarkMatrix <- function(mem, root = NULL) {
  X <- as.data.frame(mem)
  X[] <- lapply(X, as.numeric)
  
  if (nrow(X) < 2 || ncol(X) < 2) {
    stop(
      "mem must contain at least two observations ",
      "and two nodes."
    )
  }
  
  if (
    is.null(colnames(X)) ||
    anyDuplicated(colnames(X))
  ) {
    stop("mem must have unique node names.")
  }
  
  if (
    anyNA(X) ||
    any(!is.finite(as.matrix(X)))
  ) {
    stop("mem contains NA, NaN, or infinite values.")
  }
  
  zero_var <- vapply(
    X,
    stats::sd,
    numeric(1)
  ) == 0
  
  if (any(zero_var)) {
    stop(
      "Zero-variance nodes: ",
      paste(
        colnames(X)[zero_var],
        collapse = ", "
      )
    )
  }
  
  if (
    !is.null(root) &&
    !root %in% colnames(X)
  ) {
    stop(
      "root is not present in mem: ",
      root
    )
  }
  
  X
}

# A raw estimated PDAG may not admit a consistent DAG extension. This means
# its remaining edges cannot be oriented while preserving the skeleton,
# existing directions and v-structures without cycles; it is not a plotting
# failure and is not caused merely by the presence of undirected edges.
# Do not invent directions or remove edges to force an extension.
pcExtensionOrRaw <- function(pc_cpdag) {
  extension_error <- NULL
  result <- tryCatch(bnlearn::cextend(pc_cpdag), error = function(e) {
    extension_error <<- conditionMessage(e)
    NULL
  })
  if (is.null(result)) {
    result <- pc_cpdag
    attr(result, "pc_output_type") <- "raw_fallback"
    attr(result, "pc_extension_error") <- extension_error
    warning("PC-stable: using raw output because DAG extension failed: ",
            extension_error, call. = FALSE)
  } else {
    attr(result, "pc_output_type") <- "dag_extension"
  }
  result
}

benchmarkPCFGES <- function(
    mem,
    root = "A",
    alpha = 0.01,
    ges_lambda = NULL,
    num_cores = 1,
    verbose = FALSE,
    continue_on_error = TRUE
) {
  X <- .validateBenchmarkMatrix(
    mem,
    root
  )
  
  labels <- colnames(X)
  
  if (is.null(ges_lambda)) {
    ges_lambda <- 0.5 * log(nrow(X))
  }
  
  penalty_discount <- 2 * ges_lambda / log(nrow(X))
  
  blacklist <- if (is.null(root)) {
    NULL
  } else {
    data.frame(
      from = setdiff(labels, root),
      to = root,
      stringsAsFactors = FALSE
    )
  }
  
  knowledge <- causalDisco::knowledge() |>
    causalDisco::add_vars(labels)
  
  if (!is.null(root)) {
    knowledge <- knowledge |>
      causalDisco::add_exogenous(
        tidyselect::all_of(root)
      )
  }
  
  cluster <- if (num_cores > 1) {
    parallel::makeCluster(num_cores)
  } else {
    NULL
  }
  
  if (!is.null(cluster)) {
    on.exit(
      parallel::stopCluster(cluster),
      add = TRUE
    )
  }
  
  pc_error <- NULL
  pc_cpdag <- NULL
  
  pc <- tryCatch(
    {
      pc_cpdag <- bnlearn::pc.stable(
        X,
        cluster = cluster,
        blacklist = blacklist,
        test = "zf",
        alpha = alpha,
        max.sx = NULL,
        debug = verbose,
        undirected = FALSE
      )
      
      # Prefer a consistent DAG extension when one exists. cextend() is not
      # just a plotting conversion: it must preserve the skeleton, existing
      # arc directions and v-structures while orienting the remaining edges
      # without introducing directed cycles. An estimated partially directed
      # graph need not admit such an extension; undirected edges alone are
      # not a reason for failure. Do not force directions or delete edges to
      # make this step succeed.
      # Saved real-data outputs reproduce "no consistent extension" for
      # IntestD (root-informed and root-free) and HHR (root-free).
      # Historical fallback figures were generated separately from the raw
      # PC output. The wrapper now applies the same extension-first fallback
      # automatically and records pc_output_type / pc_extension_error.
      # For a raw fallback, draw the raw graph and score its arcs as stored:
      # an undirected adjacency has two directional records, counted
      # separately in the benchmark, not two established causal directions.
      dag <- pcExtensionOrRaw(pc_cpdag)
      
      if (
        !is.null(root) &&
        any(bnlearn::arcs(dag)[, "to"] == root)
      ) {
        stop(
          "PC-stable violated the root constraint."
        )
      }
      
      dag
    },
    error = function(e) {
      pc_error <<- e
      
      if (!continue_on_error) {
        stop(e)
      }
      
      NULL
    }
  )
  
  ges_error <- NULL
  ges_cpdag <- NULL
  
  ges <- tryCatch(
    {
      method <- causalDisco::ges(
        engine = "tetrad",
        score = "sem_bic",
        penalty_discount = penalty_discount,
        symmetric_first_step = FALSE,
        max_degree = -1,
        parallelized = num_cores > 1,
        faithfulness_assumed = FALSE
      )
      
      fit <- causalDisco::disco(
        data = X,
        method = method,
        knowledge = knowledge
      )
      
      amat <- t(
        caugi::as_adjacency(
          fit$caugi
        )
      )
      
      storage.mode(amat) <- "integer"
      dimnames(amat) <- list(
        labels,
        labels
      )
      
      ges_cpdag <- pcalgAmatToBnlearn(
        amat
      )
      
      dag <- bnlearn::cextend(
        ges_cpdag
      )
      
      if (
        !is.null(root) &&
        any(bnlearn::arcs(dag)[, "to"] == root)
      ) {
        stop(
          "FGES violated the root constraint."
        )
      }
      
      dag
    },
    error = function(e) {
      ges_error <<- e
      
      if (!continue_on_error) {
        stop(e)
      }
      
      NULL
    }
  )
  
  out <- list(
    pc_stable = pc,
    ges = ges
  )
  
  attr(out, "status") <- c(
    pc_stable = if (is.null(pc)) "failed" else attr(pc, "pc_output_type"),
    ges = if (is.null(ges)) "failed" else "completed"
  )
  
  attr(out, "errors") <- list(
    pc_stable = if (is.null(pc_error)) {
      NULL
    } else {
      conditionMessage(pc_error)
    },
    ges = if (is.null(ges_error)) {
      NULL
    } else {
      conditionMessage(ges_error)
    }
  )
  
  attr(out, "cpdag") <- list(
    pc_stable = pc_cpdag,
    ges = ges_cpdag
  )
  attr(out, "pc_extension_error") <- attr(pc, "pc_extension_error")
  
  out
}

configureGradientBenchmark <- function(
    python,
    script_dir
) {
  if (!nzchar(python)) {
    stop(
      "Set CAUSFATE_GRADIENT_PYTHON to the ",
      "gradient environment's Python executable."
    )
  }
  
  python <- normalizePath(
    path.expand(python),
    mustWork = TRUE
  )
  
  script_dir <- normalizePath(
    path.expand(script_dir),
    mustWork = TRUE
  )
  
  Sys.setenv(
    CAUSFATE_GRADIENT_PYTHON = python,
    CAUSFATE_GRADIENT_SCRIPT_DIR = script_dir,
    PYTHONNOUSERSITE = "1",
    PYTHONPATH = "",
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    CUDA_VISIBLE_DEVICES = ""
  )
  
  list(
    python = python,
    script_dir = script_dir
  )
}

defaultGradientConfig <- function() {
  list(
    threshold = 0.3,
    notears = list(
      lambda1 = 0.1,
      loss_type = "l2",
      max_iter = 100,
      h_tol = 1e-8,
      rho_max = 1e16,
      standardize = TRUE
    ),
    golem = list(
      mode = "two-stage",
      lambda1_ev = 0.02,
      lambda1_nv = 0.002,
      lambda2 = 5,
      learning_rate = 1e-3,
      num_iter_ev = 20000,
      num_iter_nv = 20000,
      checkpoint_iter = 5000,
      device = "cpu",
      standardize = TRUE
    ),
    dag_gnn = list(
      encoder_type = "mlp",
      decoder_type = "mlp",
      encoder_hidden = 64,
      decoder_hidden = 64,
      epochs = 300,
      k_max_iter = 100,
      batch_size = 100,
      lr = 3e-3,
      lr_decay = 200,
      gamma = 1,
      tau_a = 0,
      h_tolerance = 1e-8,
      device = "cpu",
      device_id = "0",
      standardize = TRUE
    )
  )
}

.readGradientMatrix <- function(
    file,
    labels,
    binary = FALSE
) {
  if (!file.exists(file)) {
    stop(
      "Expected Python output was not created: ",
      file
    )
  }
  
  x <- as.matrix(
    utils::read.csv(
      file,
      row.names = 1,
      check.names = FALSE
    )
  )
  
  storage.mode(x) <- "double"
  
  if (
    !setequal(rownames(x), labels) ||
    !setequal(colnames(x), labels)
  ) {
    stop(
      "Python output nodes do not match mem."
    )
  }
  
  x <- x[
    labels,
    labels,
    drop = FALSE
  ]
  
  if (binary) {
    x <- 1L * (x != 0)
  }
  
  x
}

.adjacencyToBnlearnDag <- function(
    adjacency,
    root = NULL
) {
  adjacency <- 1L * (
    as.matrix(adjacency) != 0
  )
  
  diag(adjacency) <- 0L
  labels <- colnames(adjacency)
  
  if (
    is.null(labels) ||
    !identical(
      rownames(adjacency),
      labels
    )
  ) {
    stop(
      "adjacency must have matching row ",
      "and column names."
    )
  }
  
  if (
    !is.null(root) &&
    any(adjacency[, root] != 0)
  ) {
    stop(
      "Gradient method violated ",
      "the root constraint."
    )
  }
  
  dag <- bnlearn::empty.graph(
    labels
  )
  
  idx <- which(
    adjacency != 0,
    arr.ind = TRUE
  )
  
  if (nrow(idx) > 0) {
    for (i in seq_len(nrow(idx))) {
      dag <- bnlearn::set.arc(
        dag,
        rownames(adjacency)[idx[i, 1]],
        colnames(adjacency)[idx[i, 2]],
        check.cycles = TRUE
      )
    }
  }
  
  dag
}

.gradientArgs <- function(
    method,
    cfg,
    threshold,
    seed
) {
  switch(
    method,
    notears = c(
      "--lambda1", cfg$lambda1,
      "--loss-type", cfg$loss_type,
      "--max-iter", cfg$max_iter,
      "--h-tol", cfg$h_tol,
      "--rho-max", cfg$rho_max,
      "--threshold", threshold
    ),
    golem = c(
      "--mode", cfg$mode,
      "--lambda1-ev", cfg$lambda1_ev,
      "--lambda1-nv", cfg$lambda1_nv,
      "--lambda2", cfg$lambda2,
      "--learning-rate", cfg$learning_rate,
      "--num-iter-ev", cfg$num_iter_ev,
      "--num-iter-nv", cfg$num_iter_nv,
      "--checkpoint-iter", cfg$checkpoint_iter,
      "--threshold", threshold,
      "--seed", seed,
      "--device", cfg$device
    ),
    dag_gnn = c(
      "--encoder-type", cfg$encoder_type,
      "--decoder-type", cfg$decoder_type,
      "--encoder-hidden", cfg$encoder_hidden,
      "--decoder-hidden", cfg$decoder_hidden,
      "--epochs", cfg$epochs,
      "--k-max-iter", cfg$k_max_iter,
      "--batch-size", cfg$batch_size,
      "--lr", cfg$lr,
      "--lr-decay", cfg$lr_decay,
      "--gamma", cfg$gamma,
      "--tau-a", cfg$tau_a,
      "--h-tolerance", cfg$h_tolerance,
      "--threshold", threshold,
      "--seed", seed,
      "--device", cfg$device,
      "--device-id", cfg$device_id
    )
  ) |>
    as.character()
}

.runGradientMethod <- function(
    mem,
    method,
    root,
    threshold,
    seed,
    cfg,
    python,
    script_dir,
    output_dir,
    verbose = TRUE
) {
  X <- .validateBenchmarkMatrix(
    mem,
    root
  )
  
  script <- c(
    notears = "notears_root.py",
    golem = "golem_root.py",
    dag_gnn = "dag_gnn_root.py"
  )[[method]]
  
  script_file <- normalizePath(
    file.path(
      script_dir,
      script
    ),
    mustWork = TRUE
  )
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  input_file <- file.path(
    output_dir,
    "mem.csv"
  )
  
  utils::write.csv(
    X,
    input_file,
    row.names = FALSE,
    quote = FALSE
  )
  
  args <- c(
    shQuote(script_file),
    "--input",
    shQuote(input_file),
    "--output-dir",
    shQuote(output_dir),
    "--root",
    if (is.null(root)) "none" else root,
    .gradientArgs(
      method,
      cfg,
      threshold,
      seed
    )
  )
  
  if (!isTRUE(cfg$standardize)) {
    args <- c(
      args,
      "--no-standardize"
    )
  }
  
  log <- system2(
    python,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )
  
  status <- attr(
    log,
    "status"
  )
  
  if (
    !is.null(status) &&
    status != 0L
  ) {
    stop(
      toupper(method),
      " failed:\n",
      paste(
        tail(log, 100),
        collapse = "\n"
      )
    )
  }
  
  binary <- .readGradientMatrix(
    file.path(
      output_dir,
      paste0(
        method,
        "_binary.csv"
      )
    ),
    colnames(X),
    binary = TRUE
  )
  
  if (verbose) {
    message(
      "  ",
      toupper(
        gsub(
          "_",
          "-",
          method
        )
      ),
      ": ",
      sum(binary),
      " arcs"
    )
  }
  
  .adjacencyToBnlearnDag(
    binary,
    root
  )
}

benchmarkGradientNets <- function(
    mem,
    root = "A",
    seed = 42,
    config = defaultGradientConfig(),
    python,
    script_dir,
    output_dir,
    verbose = TRUE
) {
  methods <- c(
    "notears",
    "golem",
    "dag_gnn"
  )
  
  out <- lapply(methods, function(method) {
    .runGradientMethod(
      mem = mem,
      method = method,
      root = root,
      threshold = config$threshold,
      seed = seed,
      cfg = config[[method]],
      python = python,
      script_dir = script_dir,
      output_dir = file.path(
        output_dir,
        method
      ),
      verbose = verbose
    )
  })
  
  names(out) <- methods
  out
}

colorEdges <- function(
    net,
    ref,
    node_col,
    main = ""
) {
  node_col <- node_col[
    bnlearn::nodes(net)
  ]
  
  node_col[is.na(node_col)] <- "#D9D9D9"
  
  if (bnlearn::narcs(net) == 0) {
    g <- igraph::make_empty_graph(
      n = length(
        bnlearn::nodes(net)
      ),
      directed = TRUE
    )
    
    igraph::V(g)$name <- bnlearn::nodes(net)
    
    plot(
      g,
      layout = igraph::layout_in_circle(g),
      vertex.color = node_col[igraph::V(g)$name],
      vertex.label = NA,
      vertex.size = 35,
      vertex.frame.color = "black",
      edge.arrow.size = 0.5,
      main = main
    )
    
    return(
      invisible(NULL)
    )
  }
  
  ref_arcs <- as.data.frame(
    bnlearn::arcs(ref)
  )
  
  ref_names <- paste0(
    ref_arcs$from,
    "~",
    ref_arcs$to
  )
  
  dist <- igraph::distances(
    bnlearn::as.igraph(ref),
    mode = "out"
  )
  
  learned <- as.data.frame(
    bnlearn::arcs(net)
  )
  
  learned_names <- paste0(
    learned$from,
    "~",
    learned$to
  )
  
  learned_dist <- dist[
    cbind(
      match(
        learned$from,
        rownames(dist)
      ),
      match(
        learned$to,
        colnames(dist)
      )
    )
  ]
  
  edge_col <- ifelse(
    learned_names %in% ref_names,
    "black",
    ifelse(
      is.finite(learned_dist),
      "#C9C9C9",
      "red"
    )
  )
  
  names(edge_col) <- learned_names
  
  graph <- bnlearn::graphviz.plot(
    net,
    shape = "circle",
    render = FALSE,
    main = main
  )
  
  graph::edgeRenderInfo(graph) <- list(
    col = edge_col
  )
  
  graph::nodeRenderInfo(graph) <- list(
    fill = node_col,
    label = setNames(
      rep(
        "",
        length(node_col)
      ),
      names(node_col)
    )
  )
  
  Rgraphviz::renderGraph(
    graph
  )
}

benchmarkSimulationDataset <- function(
    n_nodes,
    data_type,
    cache_dir,
    gradient_output_dir,
    python,
    script_dir,
    root = "A",
    seed = 42,
    overwrite = FALSE,
    gradient_config = defaultGradientConfig()
) {
  sim <- generateSimulationData(
    n_nodes,
    data_type,
    seed
  )
  
  seu <- sim$Seurat_obj
  ref <- sim$ref_net
  
  mem <- causfate::gem2mem(
    as.data.frame(
      seu@assays$RNA@data
    ),
    seu$celltype,
    "mean"
  )
  
  cfg <- getSimulationConfig(
    n_nodes,
    data_type
  )
  
  dataset <- paste0(
    n_nodes,
    "_",
    data_type
  )
  
  dag_smpl <- getCachedDagSamples(
    file.path(
      cache_dir,
      paste0(
        dataset,
        ".rds"
      )
    ),
    overwrite,
    function() {
      set.seed(seed)
      
      causfate::BNLearning(
        seu,
        frac = 0.2,
        N_smpl = cfg$N_smpl,
        params = cfg$params,
        root = root,
        mode = "single_cell",
        ncores = 1,
        dagMethod = "hc",
        ugMethod = "cmi2ni",
        seed = seed
      )
    }
  )
  
  cf <- buildCausFateNet(
    dag_smpl,
    bnlearn::nodes(ref),
    cfg$freq_cutoff
  )
  
  cf <- causfate::trimDAG(
    mem,
    cf,
    threshold_value = cfg$trim_threshold,
    max_arc = cfg$max_arc,
    plot = FALSE
  )
  
  classical <- benchmarkNets(
    mem,
    root,
    seed = seed
  )
  
  pc_fges <- benchmarkPCFGES(
    mem,
    root,
    num_cores = 1,
    verbose = FALSE,
    continue_on_error = TRUE
  )
  
  gradient <- benchmarkGradientNets(
    mem = mem,
    root = root,
    seed = seed,
    config = gradient_config,
    python = python,
    script_dir = script_dir,
    output_dir = file.path(
      gradient_output_dir,
      dataset
    ),
    verbose = TRUE
  )
  
  result <- list(
    ref = ref,
    CIBER = cf,
    mmhc = classical$mmhc,
    rsmax2 = classical$rsmax2,
    h2pc = classical$h2pc,
    pc_stable = pc_fges$pc_stable,
    ges = pc_fges$ges,
    notears = gradient$notears,
    golem = gradient$golem,
    dag_gnn = gradient$dag_gnn
  )
  
  attr(result, "tetrad_status") <- attr(
    pc_fges,
    "status"
  )
  
  attr(result, "tetrad_errors") <- attr(
    pc_fges,
    "errors"
  )
  
  attr(result, "tetrad_cpdag") <- attr(
    pc_fges,
    "cpdag"
  )
  
  result
}

plotBenchmarkStructures <- function(
    result,
    file,
    dataset,
    method_labels,
    node_col
) {
  grDevices::pdf(
    file,
    width = 7,
    height = 6,
    onefile = TRUE
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  for (id in names(method_labels)) {
    net <- result[[id]]
    
    if (!inherits(net, "bn")) {
      next
    }
    
    colorEdges(
      net = net,
      ref = result$ref,
      node_col = node_col,
      main = paste(
        dataset,
        method_labels[[id]],
        sep = " - "
      )
    )
  }
  
  invisible(NULL)
}
