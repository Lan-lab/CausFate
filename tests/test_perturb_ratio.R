suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
})
source("R/EffectMat_func.R")
source("R/BNLearning_func.R")

graph <- bnlearn::set.arc(bnlearn::empty.graph(c("A", "B")), "A", "B")
expression <- data.frame(
  A = c(1, 2, 3, 5), B = c(1, 4, 2, 9),
  row.names = paste0("g", 1:4)
)

run_perturbation <- function(ratio) {
  PerturbResult(
    net_struc = graph,
    data = expression,
    mode = "bulk",
    perturbation = "single_feature",
    perturb_ratio = ratio,
    ncores = 1
  )
}

stopifnot(!"deletion" %in% names(formals(PerturbResult)))
deleted <- run_perturbation("deletion")
stopifnot(isTRUE(all.equal(
  as.numeric(deleted$perturb[[1]][1, ]),
  c(1.92857142857143, 2, 1.92307692307692, 0.5)
)))

for (case in list(
  list(ratio = 0, expected = c(1.65384615384615, 1.69491525423729,
                               1.85714285714286, 0.9)),
  list(ratio = 0.5, expected = c(1.73684210526316, 1.72727272727273,
                                 2.05161290322581, 0.828571428571429)),
  list(ratio = 2, expected = c(2, 2.17142857142857,
                               0.882352941176471, 1.9))
)) {
  actual <- run_perturbation(case$ratio)
  stopifnot(isTRUE(all.equal(
    as.numeric(actual$perturb[[1]][1, ]), case$expected,
    tolerance = 1e-6
  )))
}

cells <- as.data.frame(cbind(
  matrix(expression$A, nrow = 4, ncol = 4),
  matrix(expression$B, nrow = 4, ncol = 4)
))
rownames(cells) <- rownames(expression)
cell_labels <- c(rep("A", 4), rep("B", 4))
for (ratio in list("deletion", 0, 0.5, 2)) {
  bulk <- run_perturbation(ratio)
  single_cell <- PerturbResult(
    net_struc = graph,
    data = cells,
    meta = cell_labels,
    index = list(seq_len(8)),
    n_sample = 1,
    mode = "single_cell",
    perturbation = "single_feature",
    perturb_ratio = ratio,
    ncores = 1
  )
  stopifnot(isTRUE(all.equal(
    as.matrix(single_cell$perturb[[1]]),
    as.matrix(bulk$perturb[[1]])
  )))
}

for (invalid in list("DELETE", -1, Inf, NA_real_)) {
  error <- tryCatch(run_perturbation(invalid), error = identity)
  stopifnot(inherits(error, "error"))
}

legacy_error <- tryCatch(
  PerturbResult(
    net_struc = graph,
    data = expression,
    mode = "bulk",
    perturbation = "single_feature",
    deletion = TRUE,
    ncores = 1
  ),
  error = identity
)
stopifnot(inherits(legacy_error, "error"))

cat("Perturbation-ratio interface checks passed.\n")
