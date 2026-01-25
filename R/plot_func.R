#' plot_brweight
#'
#' @param brweight numeric
#' @param title character
#' @param font_size interger
#'
#' @return plot
#' @export
plot_brweight <- function(brweight, title = NULL, font_size = 40) {
  gp <- graph::ftM2graphNEL(as.matrix(brweight$branch))
  eAttrs <- list()
  nAttrs <- list()
  attrs <- list()
  eAttrs$label <- brweight$weight
  labels <- graph::nodes(gp)
  fsize <- rep(font_size, length(labels))
  names(fsize) <- labels
  nAttrs$fontsize <- fsize
  nAttrs$fillcolor <- c("HSC" = "grey")
  attrs$edge$fontsize <- font_size
  # attrs$edge$fontcolor <- 'red'
  plot(gp, edgeAttrs = eAttrs, nodeAttrs = nAttrs, attrs = attrs, main = title)
}

#' dagPlot
#'
#' @param DAG list
#' @param weight logic
#' @param ... other
#'
#' @return plot
#' @export
dagPlot <- function(DAG, weight = FALSE, ...) {
  gp <- DAG[, c("from", "to")] %>%
    as.matrix() %>%
    graph::ftM2graphNEL(., edgemode = "directed")
  eAttrs <- list()
  if (isTRUE(weight)) {
    eAttrs$lwd <- DAG$freq
    names(eAttrs$lwd) <- rownames(DAG)
  }
  plot(gp, edgeAttrs = eAttrs, ...)
}

#' ugPlot
#'
#' @param ug data.frame
#' @param weight logic
#' @param ... other
#'
#' @return plot
#' @export
ugPlot <- function(ug, weight = FALSE, ...) {
  gp <- ug[, c("node1", "node2")] %>%
    as.matrix() %>%
    graph::ftM2graphNEL(., edgemode = "undirected")
  eAttrs <- list()
  if (isTRUE(weight)) {
    ug <- ug %>% dplyr::mutate(
      edge1 = paste0(ug$node1, "~", ug$node2),
      edge2 = paste0(ug$node2, "~", ug$node1)
    )
    edge <- alter(paste0(ug$node1, "~", ug$node2), paste0(ug$node2, "~", ug$node1))
    rownames(ug) <- edge[edge %in% graph::edgeNames(gp)]
    eAttrs$lwd <- ug[graph::edgeNames(gp), "strength"]
    names(eAttrs$lwd) <- graph::edgeNames(gp)
  }
  plot(gp, edgeAttrs = eAttrs)
}
