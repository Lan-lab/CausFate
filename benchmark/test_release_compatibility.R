# Run from the package root with Rscript benchmark/test_release_compatibility.R.
# These fixtures catch changed signs/scales, broken aliases, and loss of raw
# PC output on extension failure. Expected numbers are calculated by hand.
suppressPackageStartupMessages({library(dplyr); library(foreach)})
source("R/EffectMat_func.R")
source("R/findFeatures_func.R")
source("benchmark/benchmark_functions.R")
failures <- character()
check <- function(name, expr) {
  tryCatch({force(expr); cat("PASS:", name, "\n")}, error = function(e) {
    failures <<- c(failures, paste(name, conditionMessage(e), sep = ": "))
    cat("FAIL:", tail(failures, 1), "\n")
  })
}
ref <- function(a, b) matrix(c(a,b), 2, 1, dimnames=list(c("A~B","B~C"),"ref"))
pert <- function(a,b,c,d) matrix(c(a,b,c,d),2,2,dimnames=list(c("A~B","B~C"),c("g1","g2")))
x <- list(ref=list(ref(4,2),ref(6,4)),
          perturb=list(pert(1,5,7,0),pert(3,7,9,2)),
          n_sample=2,n_permutation=1)
want <- matrix(c(6,-6,-6,4),2,2,byrow=TRUE,
               dimnames=list(c("g1","g2"),c("A~B","B~C")))
for (metric in c("diff_mean","DIFF_MEAN")) check(paste("signed",metric), {
  stopifnot(isTRUE(all.equal(as.matrix(EffectMatrix(x,dist_metric=metric)),want)))
})
check("legacy mean and default preserve signed sums", {
  stopifnot(isTRUE(all.equal(as.matrix(EffectMatrix(x,mode="mean")),want)),
            isTRUE(all.equal(as.matrix(EffectMatrix(x)),want)))
})
check("abs is applied only by diffScore", {
  z<-EffectMatrix(x,dist_metric="diff_mean")
  stopifnot(identical(unname(diffScore(z,colnames(z))),c(12,10)),
            identical(unname(diffScore(z,colnames(z),abs=FALSE)),c(0,-2)))
})
check("MMD alias preserves existing estimator", {
  stopifnot(isTRUE(all.equal(EffectMatrix(x,dist_metric="mmd"),EffectMatrix(x,mode="MMD"))),
            isTRUE(all.equal(EffectMatrix(x,dist_metric="MMD"),EffectMatrix(x,mode="MMD"))))
})
check("equivalent simultaneous aliases accepted", {
  stopifnot(isTRUE(all.equal(as.matrix(EffectMatrix(x,mode="mean",dist_metric="DIFF_MEAN")),want)))
})
check("conflicting metrics rejected", {
  e<-tryCatch(EffectMatrix(x,mode="mean",dist_metric="energy"),error=identity)
  stopifnot(inherits(e,"error"),grepl("conflict",conditionMessage(e),ignore.case=TRUE))
})
check("unknown metric rejected", {
  e<-tryCatch(EffectMatrix(x,dist_metric="unknown"),error=identity)
  stopifnot(inherits(e,"error"))
})
check("energy alias", {
  stopifnot(isTRUE(all.equal(EffectMatrix(x,dist_metric="ENERGY"),EffectMatrix(x,mode="energy"))))
})
for(metric in c("w1","w2")) check(paste("Wasserstein",metric), {
  z<-EffectMatrix(x,dist_metric=metric)
  stopifnot(isTRUE(all.equal(as.matrix(z),matrix(c(3,3,3,2),2,2,byrow=TRUE,
                  dimnames=list(c("g1","g2"),c("A~B","B~C"))))),
            isTRUE(all.equal(z,EffectMatrix(x,mode=toupper(metric)))))
})
check("PC successful extension", {
  p<-bnlearn::set.edge(bnlearn::empty.graph(c("A","B")),"A","B")
  z<-pcExtensionOrRaw(p)
  stopifnot(bnlearn::narcs(z)==1,identical(attr(z,"pc_output_type"),"dag_extension"))
})
check("PC failed extension retains all raw arcs", {
  p<-bnlearn::empty.graph(LETTERS[1:4])
  for(e in list(c("A","B"),c("B","C"),c("C","D"),c("D","A"))) p<-bnlearn::set.edge(p,e[1],e[2])
  z<-suppressWarnings(pcExtensionOrRaw(p))
  stopifnot(identical(bnlearn::arcs(z),bnlearn::arcs(p)),nrow(bnlearn::arcs(z))==8,
            identical(attr(z,"pc_output_type"),"raw_fallback"),
            grepl("consistent extension",attr(z,"pc_extension_error")))
})
if(length(failures)) stop(paste(failures,collapse="\n"))
cat("All compatibility checks passed.\n")
