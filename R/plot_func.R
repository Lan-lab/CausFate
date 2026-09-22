#' Plot branch weights
#'
#' @param brweight A list containing branch and weight information.
#' @param title Optional character string used as the plot title.
#' @param font_size Numeric font size.
#'
#' @return The plotted graph, invisibly.
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

#' Plot a directed acyclic graph
#'
#' @param DAG A data frame containing `from`, `to`, and optionally `freq`
#'   columns.
#' @param weight Logical; whether to display edge weights.
#' @param ... Additional arguments passed to `plot()`.
#'
#' @return The plotted graph, invisibly.
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

#' Plot an undirected graph
#'
#' @param ug A data frame containing `node1`, `node2`, and optionally `strength`
#'   columns.
#' @param weight Logical; whether to display edge weights.
#' @param ... Additional arguments passed to `plot()`.
#'
#' @return The plotted graph, invisibly.
#' @export
ugPlot <- function(ug, weight = FALSE, ...) {
  gp <- ug[, c("node1", "node2")] %>%
    as.matrix() %>%
    graph::ftM2graphNEL(., edgemode = "undirected")
  eAttrs <- list()
  if (isTRUE(weight)) {
    ug <- ug %>%
      dplyr::mutate(
        edge1 = paste0(ug$node1, "~", ug$node2),
        edge2 = paste0(ug$node2, "~", ug$node1)
      )
    edge <- alter(
      paste0(ug$node1, "~", ug$node2),
      paste0(ug$node2, "~", ug$node1)
    )
    rownames(ug) <- edge[edge %in% graph::edgeNames(gp)]
    eAttrs$lwd <- ug[graph::edgeNames(gp), "strength"]
    names(eAttrs$lwd) <- graph::edgeNames(gp)
  }
  plot(gp, edgeAttrs = eAttrs)
}

#' Plot a network with reference-based edge colours
#'
#' @param net A `bn` network to plot.
#' @param ref A reference `bn` network.
#' @param main Character string used as the plot title.
#' @param node_col Named vector of node colours.
#' @param fontsize Font size for node labels.
#' @param node_size Node diameter in inches.
#' @param sep Spacing between nodes (`ranksep` and `nodesep`).
#'
#' @return The rendered graph, invisibly.
#' @export
colorEdges <- function(
  net,
  ref,
  main = "",
  node_col = NULL,
  fontsize = 12,
  node_size = 0.7,
  sep = 1.2
) {
  # 1. Edge Color Logic (using explicit bnlearn namespace)
  ref_arcs <- bnlearn::arcs(ref) %>% as.data.frame()
  ref_arcs <- paste0(ref_arcs$from, '~', ref_arcs$to)

  df <- igraph::distances(bnlearn::as.igraph(ref), mode = "out") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("from") %>%
    tidyr::pivot_longer(-from, names_to = "to", values_to = "dist") %>%
    mutate(arc = paste0(from, '~', to)) %>%
    mutate(
      col = ifelse(
        arc %in% ref_arcs,
        'black',
        ifelse(is.infinite(dist), 'red', '#C9C9C9')
      )
    )

  edge_colors <- df$col
  names(edge_colors) <- df$arc

  # 2. Convert to graphNEL
  gNEL <- bnlearn::as.graphNEL(net)
  # CRITICAL FIX: Use graph::nodes for graphNEL objects to avoid bnlearn conflict
  node_names <- graph::nodes(gNEL)

  # 3. Setup Layout Attributes
  attrs <- list(
    graph = list(
      ranksep = as.character(sep),
      nodesep = as.character(sep),
      overlap = "false"
    )
  )

  # 4. Define Node Attributes for Layout
  node_attrs <- list(
    shape = setNames(rep("circle", length(node_names)), node_names),
    fixedsize = setNames(rep(TRUE, length(node_names)), node_names),
    width = setNames(
      rep(as.character(node_size), length(node_names)),
      node_names
    ),
    height = setNames(
      rep(as.character(node_size), length(node_names)),
      node_names
    )
  )

  # 5. Perform Layout Calculation
  graph <- Rgraphviz::layoutGraph(gNEL, attrs = attrs, nodeAttrs = node_attrs)

  # 6. Final Drawing Overrides
  nodeRenderInfo(graph) <- list(
    fill = if (!is.null(node_col)) node_col[node_names] else "transparent",
    fontsize = fontsize,
    # Standardizing points for rendering
    width = node_size * 72,
    height = node_size * 72
  )

  if (nrow(bnlearn::arcs(net)) > 0) {
    # Ensure edge colors are applied to the calculated graph
    edgeRenderInfo(graph) <- list(col = edge_colors)
  }

  # 7. Final Render
  Rgraphviz::renderGraph(graph)
  if (main != "") title(main)
}
