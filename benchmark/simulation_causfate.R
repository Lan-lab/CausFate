library(causfate)
library(tidyverse)
# Run from the package root: Rscript benchmark/simulation_causfate.R
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
benchmark_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else if (file.exists("benchmark_functions.R")) "." else "benchmark"

source(file.path(benchmark_dir, "benchmark_functions.R"))

# Set root <- NULL for root-free reconstruction. The environment override
# also supports: CAUSFATE_ROOT=none Rscript benchmark/simulation_causfate.R
root <- Sys.getenv("CAUSFATE_ROOT", "A")
if (tolower(root) %in% c("none", "null")) root <- NULL
seed <- 42
overwrite <- FALSE

root_tag <- if (is.null(root)) "root_free" else paste0("root_", root)
out_dir <- file.path(Sys.getenv("CAUSFATE_OUTPUT_DIR",
                                file.path(benchmark_dir, "results")),
                     "simulation_causfate", root_tag)
cache_dir <- file.path(out_dir, "dag_smpl")
plot_dir <- file.path(out_dir, "structures")

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

node_col <- setNames(
  c("#999898", "#C69C6D", "#FCEE21", "#29ABE2", "#D88176",
    "#6896E9", "#8CC63F", "#838BC5", "#FCCF2E"),
  LETTERS[1:9]
)

for (n_nodes in 3:9) {
  # Internal key for edge-free (disconnected) reference networks.
  for (data_type in c("discrete", "sparse", "tree", "dense")) {
    dataset <- paste0(n_nodes, "_", data_type)
    result_file <- file.path(out_dir, paste0(dataset, ".rds"))
    
    message("Running ", dataset)
    
    if (file.exists(result_file) && !overwrite) {
      result <- readRDS(result_file)
    } else {
      sim <- generateSimulationData(n_nodes, data_type, seed)
      seu <- sim$Seurat_obj
      ref <- sim$ref_net
      
      mem <- causfate::gem2mem(
        as.data.frame(seu@assays$RNA@data),
        seu$celltype,
        "mean"
      )
      
      cfg <- getSimulationConfig(n_nodes, data_type)
      
      dag_smpl <- getCachedDagSamples(
        file.path(cache_dir, paste0(dataset, ".rds")),
        overwrite,
        function() {
          set.seed(seed)
          causfate::BNLearning(
            seu,
            resampling_fraction = 0.2,
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
        freq_cutoff = cfg$freq_cutoff
      )
      
      cf <- causfate::trimDAG(
        mem,
        cf,
        threshold_value = cfg$trim_threshold,
        max_arc = cfg$max_arc,
        plot = FALSE
      )
      
      result <- list(
        ref = ref,
        CIBER = cf
      )
      attr(result, "example_config") <- list(root = root, seed = seed, settings = cfg)
      
      saveRDS(result, result_file)
    }
    
    grDevices::pdf(
      file.path(plot_dir, paste0(dataset, ".pdf")),
      width = 7,
      height = 6
    )
    
    colorEdges(
      net = result$ref,
      ref = result$ref,
      node_col = node_col,
      main = paste(dataset, "Reference", sep = " - ")
    )
    
    colorEdges(
      net = result$CIBER,
      ref = result$ref,
      node_col = node_col,
      main = paste(dataset, "CausFate", sep = " - ")
    )
    
    grDevices::dev.off()
    
    message(
      dataset,
      " | Reference arcs = ", bnlearn::narcs(result$ref),
      " | CausFate arcs = ", bnlearn::narcs(result$CIBER)
    )
    
    rm(result)
    gc()
  }
}
