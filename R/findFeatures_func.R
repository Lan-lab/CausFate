#' Generate edge set
#'
#' @param fromSet Character vector of starting nodes.
#' @param toSet Character vector of ending nodes.
#' @param sep Character string used to separate node names.
#'
#' @return A character vector of edge names.
#' @export
#'
#' @examples
#' ctypes <- c("MPPa", "MPPb")
#' s.diffScore <- setEdges(ctypes, ctypes, sep = "~")
setEdges <- function(fromSet, toSet, sep = "~") {
  edges <- paste(rep(fromSet, each = length(toSet)), toSet, sep = sep)
  selfs <- paste(fromSet, fromSet, sep = sep)
  return(setdiff(edges, selfs))
}

#' Get diffScore
#'
#' @param diffCoeff A data frame of previously calculated coefficient
#'   differences.
#' @param edgeSet A character vector of edge names.
#' @param abs Logical; whether to sum the absolute coefficient differences.
#'
#' @return A numeric vector of differential scores.
#' @export
#'
diffScore <- function(diffCoeff, edgeSet, abs = TRUE) {
  if (abs) {
    diffCoeff <- abs(diffCoeff)
  }
  dm <- diffCoeff[intersect(colnames(diffCoeff), edgeSet)]
  return(rowSums(dm))
}

#' Rank genes according to diffScores
#'
#' @param diffScores A named numeric vector of previously calculated
#'   differential scores.
#'
#' @return A character vector of ranked gene names.
#' @export
#'
dsRank <- function(diffScores) {
  diffScores %>%
    sort(decreasing = TRUE) %>%
    names()
}
