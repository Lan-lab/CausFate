library(magrittr)
source(file.path("R", "EffectMat_func.R"))

effect_formals <- names(formals(EffectMatrix))
stopifnot("dist_metric" %in% effect_formals)
stopifnot(!"mode" %in% effect_formals)

ref_fit <- c(`A->B` = 0.8, `A->C` = 0.4)
perturb_fit <- matrix(c(0.3, 0.1), ncol = 1,
                      dimnames = list(c("A->B", "A->C"), "feature_1"))
diff_bn <- list(
  ref = list(ref_fit),
  perturb = list(perturb_fit),
  n_sample = 1
)

effect <- EffectMatrix(diff_bn, dist_metric = "diff_mean")
stopifnot(is.data.frame(effect))

distribution_bn <- list(
  ref = list(
    matrix(c(0.8, 0.4), ncol = 1,
           dimnames = list(c("A->B", "A->C"), "reference")),
    matrix(c(0.7, 0.5), ncol = 1,
           dimnames = list(c("A->B", "A->C"), "reference"))
  ),
  raw = list(
    matrix(c(0.3, 0.1, 0.6, 0.2), nrow = 2,
           dimnames = list(c("A->B", "A->C"), c("feature_1", "feature_2"))),
    matrix(c(0.4, 0.2, 0.5, 0.3), nrow = 2,
           dimnames = list(c("A->B", "A->C"), c("feature_1", "feature_2")))
  )
)

available_metrics <- c("energy", "mmd", "MMD")
if (requireNamespace("transport", quietly = TRUE)) {
  available_metrics <- c("W1", "w1", "W2", "w2", available_metrics)
}

for (metric in available_metrics) {
  value <- EffectMatrix(distribution_bn, dist_metric = metric)
  stopifnot(identical(dim(value), c(2L, 2L)))
}

for (legacy_metric in c("mean", "OT", "MMD_linear")) {
  error <- try(EffectMatrix(distribution_bn, dist_metric = legacy_metric), silent = TRUE)
  stopifnot(inherits(error, "try-error"))
}

mode_error <- try(EffectMatrix(diff_bn, mode = "mean"), silent = TRUE)
stopifnot(inherits(mode_error, "try-error"))
