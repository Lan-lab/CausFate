str2arc <- function(arc_list, sep = '~'){
  lapply(strsplit(arc_list, sep), function(fs){
    assert_that(length(fs)==2)
    father <- fs[[1]]
    sons <- strsplit(fs[[2]],'') %>% unlist
    data.frame(from=father,to=sons)
  }) %>% do.call(rbind,.)
}

generateBayesianParams <- function(ref_net, intercept_coef = 5, parent_coef = 2, sd = 1){
  nodes <- nodes(ref_net)
  node_params <- lapply(nodes, function(node){
    coef <- c(as.integer(runif(1)*10*intercept_coef)/10,
              (sapply(parents(ref_net, node),function(parent){
                as.integer(runif(1)*10*parent_coef)/10
              }))) %>% as.numeric %>% set_names(c("(Intercept)",parents(ref_net, node)))
    list(coef = coef, sd = as.integer(runif(1)*10*sd)/10)
  }) %>% set_names(nodes)
  node_params
}

generateExpressionSeurat <- function(ref_net, n_feature=500, n_node_smpl=500, n_edge_smpl =1000, edge_span = 1.1){
  expressions <- vector(mode="list",length=n_feature)
  
  for(i in 1:n_feature){
    node_params <- generateBayesianParams(ref_net)
    fitted <- custom.fit(ref_net, dist = node_params)
    expressions[[i]]=rbn(fitted,n=n_node_smpl)
  }
  nodes <- nodes(ref_net)
  data <- str2str::ld2a(expressions)
  data <- str2str::a2ld(data, along = 1)
  data <- do.call(rbind, data) %>% t
  meta <- rep(nodes, times=n_node_smpl)
  meta_nodes <- meta
  
  ref_net_arcs <- arcs(ref_net)
  if(nrow(ref_net_arcs)>0){for(i in 1:nrow(ref_net_arcs)){
    df_a <- data[,sample(which(meta_nodes==ref_net_arcs[i,1]), n_edge_smpl, replace = T)]
    df_b <- data[,sample(which(meta_nodes==ref_net_arcs[i,2]), n_edge_smpl, replace = T)]
    frac <- runif(n_edge_smpl, min=-(edge_span-1)/2, max=(1+(edge_span-1)/2))
    df <- do.call(cbind,lapply(1:n_edge_smpl, function(j){
      df_a[,j]*frac[j] + df_b[,j]*(1-frac[j])
    }))
    data <- cbind(data,df)
    meta <- c(meta, ifelse(frac<0.5,ref_net_arcs[i,1],ref_net_arcs[i,2]))
  }}
  
  rownames(data) <- paste0("Gene_",1:n_feature)
  colnames(data) <- paste0("Cell_",1:ncol(data))
  data <- as.matrix(data)
  
  Seurat_obj <- CreateSeuratObject(counts = data, data = data, meta.data = list(celltype=meta))
  
  Seurat_obj <- ScaleData(Seurat_obj)
  Seurat_obj <- FindVariableFeatures(Seurat_obj)
  
  # Need this if Seurat version > 5
  Seurat_obj[["RNA"]] <- as(object = Seurat_obj[["RNA"]], Class = "Assay")
  
  Seurat_obj
}

colorEdges <- function(net, ref, main=""){
  ref_arcs <- arcs(ref) %>% as.data.frame
  ref_arcs <- paste0(ref_arcs$from,'~',ref_arcs$to)
  df <- distances(as.igraph(ref), mode="out") %>% as.data.frame %>% 
    rownames_to_column("from") %>% pivot_longer(-from, names_to = "to", values_to = "dist") %>%
    mutate(arc = paste0(from,'~',to)) %>% 
    mutate(col = ifelse(arc %in% ref_arcs, 'black', ifelse(is.infinite(dist), 'red', '#C9C9C9')))
  col <- df$col
  names(col) <- df$arc
  draw_para_arc <- list(col=col)
  graph <- graphviz.plot(net, shape='circle', render=F, main=main)
  if(nrow(arcs(net)>0)){edgeRenderInfo(graph) <- draw_para_arc}
  node_col <- c('#999898', '#C69C6D', '#FCEE21', '#29ABE2', '#D88176', '#6896E9', '#8CC63F', '#838BC5', '#FCCF2E')
  names(node_col) <-  c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I')
  node_col <- node_col[nodes(net)]
  nodeRenderInfo(graph) <- list(fill = node_col)
  tryCatch(renderGraph(graph),  error = function(e) {
    message("Caught error at iteration ",": ", e$message)
  })
  
}

benchmarkNets <- function(mem,root="A", methods = c("hc","mmhc","h2pc","rsmax2")){
    res <- lapply(methods, function(m){
    get(m)(mem, blacklist = data.frame(from = nodes(ref_net),to=root))})
    names(res) <- methods
    res
}

str2bn <- function(arcs_short, plot = T){
  arcs <- str2arc(arcs_short)
  nodes <- union(arcs$from,arcs$to)
  df2bn(arcs, nodes, plot = plot)
}
# 
# colorEdgesGrob <- function(net, ref, main=""){
#   ref_arcs <- arcs(ref) %>% as.data.frame
#   ref_arcs <- paste0(ref_arcs$from,'~',ref_arcs$to)
#   df <- distances(as.igraph(ref_net), mode="out") %>% as.data.frame %>% 
#     rownames_to_column("from") %>% pivot_longer(-from, names_to = "to", values_to = "dist") %>%
#     mutate(arc = paste0(from,'~',to)) %>% 
#     mutate(col = ifelse(arc %in% ref_arcs, 'black', ifelse(is.infinite(dist), 'red', '#C9C9C9')))
#   col <- df$col
#   names(col) <- df$arc
#   draw_para_arc <- list(col=col)
#   graph <- graphviz.plot(net, shape='circle', render=F, main=main)
#   edgeRenderInfo(graph) <- draw_para_arc
#   
#   grid::grid.newpage()
#   renderGraph(graph)
#   gridGraphics::grid.echo() 
#   grid::grid.grab()  
# }

netDiff <- function(net, ref, main=""){
  if(nrow(arcs(ref))!=0){ref_arcs <- arcs(ref) %>% as.data.frame
  ref_arcs <- paste0(ref_arcs$from,'~',ref_arcs$to)}else{ref_arcs = vector(mode="character")}
  df <- distances(as.igraph(ref), mode="out") %>% as.data.frame %>% 
    rownames_to_column("from") %>% pivot_longer(-from, names_to = "to", values_to = "dist") %>%
    mutate(arc = paste0(from,'~',to)) %>% 
    mutate(type = ifelse(arc %in% ref_arcs, 'direct', ifelse(is.infinite(dist), 'wrong', 'indirect')))
  
  if(nrow(arcs(net))!=0){
    net_arcs <- arcs(net) %>% as.data.frame
  net_arcs <- paste0(net_arcs$from,'~',net_arcs$to)
  }else{
    net_arcs <- vector(mode="character")
  }
  
  
  arc_types <- df %>% filter(arc %in% net_arcs) %>% .$type
  counts <- c(sum(arc_types=='direct'),sum(arc_types=='wrong'),sum(arc_types=='indirect'))
  names(counts) <- c("direct","wrong","indirect")
  counts["Hamming"] = (length(ref_arcs) - counts["direct"]) + counts["wrong"]
  counts["recall"] = counts["direct"] / length(ref_arcs)
  counts["precision"] = counts["direct"] / length(net_arcs)
  counts["relaxed_precision"] = (counts["direct"] + counts["indirect"]) / length(net_arcs)
  counts["Jaccard"] = counts["direct"] / (length(ref_arcs) + counts["wrong"] + counts["indirect"])
  counts
}

plotStat <- function(df, by, title, y_sign_pos = c(1.1,1.2,1.3), y_scale_breaks=seq(0, 1, 0.2)){
  df %>% ggplot(., aes(x=factor(Algorithm, levels=c("CausFate","HC","MMHC","RSMAX2","H2PC")), y=!!sym(by), fill=Algorithm)) +
    stat_summary(fun=mean, geom="bar", col="black", width=0.8) +
    stat_summary(fun.data=mean_sdl, 
                 fun.args = list(mult=1), 
                 geom="errorbar", width=.2) +
    geom_point(position = position_jitter(width = 0.3, height = 0), color = "black", size = 0.6) + 
    geom_signif(comparisons=list(c("CausFate", "MMHC")),
                y_position = y_sign_pos[1], tip_length = 0.02, vjust=0.5, map_signif_level = T) +
    geom_signif(comparisons=list(c("CausFate", "RSMAX2")),
                y_position = y_sign_pos[2], tip_length = 0.02, vjust=0.5, map_signif_level = T) +
    geom_signif(comparisons=list(c("CausFate", "H2PC")),
                y_position = y_sign_pos[3], tip_length = 0.02, vjust=0.5, map_signif_level = T) +
    ggtitle(title) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position="none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text.x = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 2.5) +
    scale_y_continuous(breaks=y_scale_breaks) 
  
}

generate_data_only <- function(n_nodes, structure_type,
                               ref_net_arcs_short = NULL, edge_span = 1.2,
                               return_net = TRUE) {
  dataset_name <- glue("{n_nodes}_{structure_type}")
  message("Generating data for ", dataset_name)
  
  node_names <- LETTERS[1:n_nodes]
  
  # Build the reference network
  if (structure_type == "discrete") {
    ref_net <- empty.graph(nodes = node_names)
  } else if (!is.null(ref_net_arcs_short)) {
    ref_net <- str2bn(ref_net_arcs_short, plot = FALSE)
    tmp <- empty.graph(nodes = node_names)
    arcs(tmp) <- arcs(ref_net)
    ref_net <- tmp
  } else {
    stop("You must provide ref_net_arcs_short for non-discrete structures.")
  }
  
  # Generate synthetic expression data
  set.seed(42)
  Seurat_obj <- generateExpressionSeurat(ref_net, edge_span = edge_span)
  
  # Return as needed
  if (return_net) {
    return(list(Seurat_obj = Seurat_obj, ref_net = ref_net))
  } else {
    return(Seurat_obj)
  }
}

plotAllStatsByDataType <- function(data_type){
  all_stats <- lapply(3:9,function(n_nodes){
    nets <- readRDS(glue("{dataset_path}/{n_nodes}_{data_type}.rds"))
    stat <- lapply(c("CausFate", "mmhc", "h2pc", "rsmax2"), function(i){
      netDiff(nets[[i]],nets[["ref"]])
    }) %>% do.call(rbind, .) %>% as.data.frame
    stat[,"n_node"] <- n_nodes
    stat[,"Algorithm"] <- c("CausFate", "MMHC", "H2PC", "RSMAX2")
    stat
  }) %>% do.call(rbind, .) 
  
  df <- all_stats %>% pivot_longer(cols = c("Hamming","recall","precision","relaxed_precision","Jaccard"),names_to = "index") %>%
    group_by(index) %>%
    mutate(
      vmin = suppressWarnings(min(value, na.rm = TRUE)),
      vmax = suppressWarnings(max(value, na.rm = TRUE)),
      all_na = all(is.na(value)),
      value_scaled = case_when(
        is.na(value) ~ 0,                            # all NA → 0
        vmin == vmax & vmax != 0 ~ 1,                # all equal and nonzero → 1
        vmin == vmax & vmax == 0 ~ 0,               # all equal and zero → 0.1
        TRUE ~ (value) / (vmax)        # normal scaling
      )
      
    ) %>%
    ungroup() %>%
    mutate(Algorithm = factor(Algorithm, levels = c("RSMAX2","MMHC","H2PC","CausFate")),
           index = factor(index, levels = c("precision","relaxed_precision","recall","Jaccard","Hamming")),
           n_node = factor(n_node, levels = 9:3)) %>%
    mutate(
      tile_x = as.numeric(as.factor(index)),
      tile_y = as.numeric(as.factor(Algorithm)),
      # Subdivide tile vertically by n_node
      bar_height = 0.8 / length(unique(n_node)),
      bar_ymin = tile_y - 0.4 + (as.numeric((n_node))-1) * bar_height,
      bar_ymax = bar_ymin + bar_height * 0.9,  # small gap
      bar_xmin = tile_x - 0.4,
      bar_xmax = tile_x - 0.4 + value_scaled * 0.8,  # scaled per index
      print_text = case_when(
        is.na(value) ~ "NaN",
        TRUE ~ as.character(round(value,2))
      )
    ) 
  df %>%
    ggplot(aes(x=(index),y=(Algorithm))) +
    geom_tile(fill="white") +
    geom_rect(
      aes(xmin = bar_xmin, xmax = bar_xmax,
          ymin = bar_ymin, ymax = bar_ymax,),
      fill = "skyblue"
    ) +
    geom_text(aes(x = tile_x + 0.4,
                  y = (bar_ymin + bar_ymax) / 2,label = print_text), hjust=1) +
    geom_text(
      data = df %>% filter(index == "precision"),  # only last column
      aes(
        x = 0,                       # or bar_xmin + offset
        y = (bar_ymin + bar_ymax)/2,
        label = n_node
      ),
      hjust = -0.1,  # slightly to the right of bar end
      size = 2.5,
      color = "black"
    )+
    theme_minimal() + xlab("") + ylab("Algorithm") +
    ggtitle(label=glue("Stats for {data_type} data")) +
    scale_x_discrete(labels = c(
      "Hamming" = "Hamming Dist.",
      "recall" = "Recall Rate",
      "precision" = "Precision Rate",
      "relaxed_precision" = "Relaxed Prec.",
      "Jaccard" = "Jaccard Index"
    )) +
    theme(axis.ticks = element_blank(), aspect.ratio = 1.2, panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(color="black"))
  
}

bn_to_svg <- function(net, ref_net, main = NULL,
                      width = 6, height = 6,
                      file = tempfile(fileext = ".svg")) {
  # Render the BN to an SVG file
  svglite::svglite(file, width = width, height = height, standalone = TRUE)
  colorEdges(net, ref_net, main = main) # your BN plotting function
  grDevices::dev.off()
  
  file
}

getAllStatsByDataType <- function(data_type){
  all_stats <- lapply(3:9,function(n_nodes){
    nets <- readRDS(glue("{dataset_path}/{n_nodes}_{data_type}.rds"))
    stat <- lapply(c("CausFate", "mmhc", "h2pc", "rsmax2"), function(i){
      netDiff(nets[[i]],nets[["ref"]])
    }) %>% do.call(rbind, .) %>% as.data.frame
    stat[,"n_node"] <- n_nodes
    stat[,"Algorithm"] <- c("CausFate", "MMHC", "H2PC", "RSMAX2")
    stat
  }) %>% do.call(rbind, .) 
  all_stats
}
