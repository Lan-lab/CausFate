library(assertthat)
library(tidyverse)
library(causfate)
library(bnlearn)
library(Seurat)
library(igraph)
library(glue)
library(ggplot2)
library(ggsignif)
library(Rgraphviz)
library(ggpubr)
library(zeallot)
library(glue)

source("~/CausFate/benchmark/simulation_func.R")
dataset_path <- "~/CausFate/benchmark/result"

# 3-discrete --------------------------------------------------------------

dataset_name = "3_discrete"
set.seed(42)
c(Seurat_obj, ref_net) %<-% generate_data_only(n_nodes=3, structure_type = "discrete")

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat

outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 3-sparse --------------------------------------------------------------

dataset_name = "3_sparse"
set.seed(42)
c(Seurat_obj, ref_net) %<-% generate_data_only(n_nodes=3, structure_type = "sparse",
                                               ref_net_arcs_short = c("A~B"))

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 3-tree --------------------------------------------------------------

dataset_name = "3_tree"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=3, structure_type = "tree",
                                                          ref_net_arcs_short = c("A~B","A~C"))) 

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 3-dense --------------------------------------------------------------

dataset_name = "3_dense"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=3, structure_type = "dense",
                                                          ref_net_arcs_short = c("A~B","B~C","A~C"))) 


set.seed(42)
Seurat_obj<-generateExpressionSeurat(ref_net)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = .6, max_arc = 2, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 4-discrete --------------------------------------------------------------

dataset_name = "4_discrete"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=4, structure_type = "discrete"))

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat

outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 4-sparse --------------------------------------------------------------

dataset_name = "4_sparse"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=4, structure_type = "sparse",
                                                          ref_net_arcs_short = c("A~B")))


set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 4-tree --------------------------------------------------------------

dataset_name = "4_tree"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=4, structure_type = "tree",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D"))) 

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 4-dense --------------------------------------------------------------

dataset_name = "4_dense"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=4, structure_type = "dense",
                                                          ref_net_arcs_short = c("A~B","B~C","B~D","A~C","A~D"))) 


set.seed(42)
Seurat_obj<-generateExpressionSeurat(ref_net)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.1, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 15, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = .6, max_arc = 2, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 5-discrete --------------------------------------------------------------

dataset_name = "5_discrete"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=5, structure_type = "discrete"))

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat

outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 5-sparse --------------------------------------------------------------

dataset_name = "5_sparse"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=5, structure_type = "sparse",
                                                          ref_net_arcs_short = c("A~B","A~C")))


set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 5-tree --------------------------------------------------------------

dataset_name = "5_tree"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=5, structure_type = "tree",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D","C~E"))) 

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 5-dense --------------------------------------------------------------

dataset_name = "5_dense"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=5, structure_type = "dense",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D","C~E","A~D","A~E"))) 


set.seed(42)
Seurat_obj<-generateExpressionSeurat(ref_net)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.1, 0.2, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 10, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
 titles <- paste0("Param: ", params)
 for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
min_Hamming <- Inf
thres <- max_arc <-0
for(i in seq(0.1,0.9,0.05)){
  for (j in 1:3){
    tmp <- trimDAG(mem, outputDAG, threshold_value = i, max_arc = j, plot = F)
    tmp2 <- netDiff(tmp,ref_net)['Hamming']
    if(tmp2 <= min_Hamming){
      min_Hamming = tmp2
      thres = i
      max_arc =j
    }
  }
}
tmp <- trimDAG(mem, outputDAG, threshold_value = .9, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 6-discrete --------------------------------------------------------------

dataset_name = "6_discrete"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=6, structure_type = "discrete"))

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat

outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 6-sparse --------------------------------------------------------------

dataset_name = "6_sparse"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=6, structure_type = "sparse",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D")))


set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 6-tree --------------------------------------------------------------

dataset_name = "6_tree"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=6, structure_type = "tree",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D","C~E","B~F"))) 

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 6-dense --------------------------------------------------------------

dataset_name = "6_dense"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=6, structure_type = "dense",
                                                          ref_net_arcs_short = c("A~BDEF","A~C","B~D","C~E","B~F"))) 


set.seed(42)
Seurat_obj<-generateExpressionSeurat(ref_net)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 10, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
min_Hamming <- Inf
thres <- max_arc <-0
for(i in seq(0.1,0.9,0.05)){
  for (j in 1:3){
    tmp <- trimDAG(mem, outputDAG, threshold_value = i, max_arc = j, plot = F)
    tmp2 <- netDiff(tmp,ref_net)['Hamming']
    if(tmp2 <= min_Hamming){
      min_Hamming = tmp2
      thres = i
      max_arc =j
    }
  }
}
tmp <- trimDAG(mem, outputDAG, threshold_value = .75, max_arc = 3, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 7-discrete --------------------------------------------------------------

dataset_name = "7_discrete"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=7, structure_type = "discrete"))

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat

outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 7-sparse --------------------------------------------------------------

dataset_name = "7_sparse"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=7, structure_type = "sparse",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D","C~E")))


set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 7-tree --------------------------------------------------------------

dataset_name = "7_tree"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=7, structure_type = "tree",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D","C~E","B~F","C~G"))) 

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 7-dense --------------------------------------------------------------

dataset_name = "7_dense"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=7, structure_type = "dense",
                                                          ref_net_arcs_short = c("A~B","A~C","B~DF","C~EG","A~DFEG"))) 


set.seed(42)
Seurat_obj<-generateExpressionSeurat(ref_net)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 10, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
min_Hamming <- Inf
thres <- max_arc <-0
for(i in seq(0.1,0.9,0.05)){
  for (j in 1:3){
    tmp <- trimDAG(mem, outputDAG, threshold_value = i, max_arc = j, plot = F)
    tmp2 <- netDiff(tmp,ref_net)['Hamming']
    if(tmp2 <= min_Hamming){
      min_Hamming = tmp2
      thres = i
      max_arc =j
    }
  }
}
tmp <- trimDAG(mem, outputDAG, threshold_value = .9, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 8-discrete --------------------------------------------------------------

dataset_name = "8_discrete"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=8, structure_type = "discrete"))

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat

outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 8-sparse --------------------------------------------------------------

dataset_name = "8_sparse"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=8, structure_type = "sparse",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D","C~E","B~F")))


set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 8-tree --------------------------------------------------------------

dataset_name = "8_tree"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=8, structure_type = "tree",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D","C~E","B~F","C~G","E~H"))) 

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 0.9, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 8-dense --------------------------------------------------------------

dataset_name = "8_dense"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=8, structure_type = "dense",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D","C~E","B~F","C~G","E~H","A~DEFG"))) 


set.seed(42)
Seurat_obj<-generateExpressionSeurat(ref_net)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 10, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 5)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
min_Hamming <- Inf
thres <- max_arc <-0
for(i in seq(0.1,0.9,0.05)){
  for (j in 1:3){
    tmp <- trimDAG(mem, outputDAG, threshold_value = i, max_arc = j, plot = F)
    tmp2 <- netDiff(tmp,ref_net)['Hamming']
    if(tmp2 <= min_Hamming){
      min_Hamming = tmp2
      thres = i
      max_arc =j
    }
  }
}
tmp <- trimDAG(mem, outputDAG, threshold_value = .65, max_arc = 3, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 9-discrete --------------------------------------------------------------

dataset_name = "9_discrete"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=9, structure_type = "discrete"))

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat

outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 9-sparse --------------------------------------------------------------

dataset_name = "9_sparse"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=9, structure_type = "sparse",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D","C~E","B~F","C~G")))


set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 1, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 9-tree --------------------------------------------------------------

dataset_name = "9_tree"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=9, structure_type = "tree",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D","C~E","B~F","C~G","E~H","G~I"))) 

set.seed(42)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 20, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 0)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
tmp <- trimDAG(mem, outputDAG, threshold_value = 0.7, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))

# 9-dense --------------------------------------------------------------

dataset_name = "9_dense"
set.seed(42)
zeallot::`%<-%`(c(Seurat_obj, ref_net),generate_data_only(n_nodes=9, structure_type = "dense",
                                                          ref_net_arcs_short = c("A~B","A~C","B~D","C~E","B~F","C~G","E~H","G~I","A~DEFG"))) 


set.seed(42)
Seurat_obj<-generateExpressionSeurat(ref_net)
data<-Seurat_obj@assays$RNA@data %>% as.data.frame
meta<-Seurat_obj$celltype

set.seed(42)
params <- seq(0.05, 0.15, 0.01)
dag_smpl <- BNLearning(Seurat_obj,
                       frac = 0.2, N_smpl = 10, params = params,
                       root = "A", mode = "single_cell", ncores = 1,
                       dagMethod = "hc", ugMethod = "cmi2ni")
# titles <- paste0("Param: ", params)
# for(i in 1:length(params))graphviz.plot(dag_smpl[[1]][[i]], main = titles[[i]])
net <- combineDAGsmpl(dag_smpl, Emin = NULL, Emax = NULL, ncores = 1)
net_mat <- net %>% df2mat(); net_mat
net <- (net_mat * (net_mat > 5)) %>%
  mat2df() %>%
  rmCyc()
outputDAG <- df2bn(net, nodes(ref_net), plot = TRUE)
mem <- gem2mem(data, meta, "mean")
min_Hamming <- Inf
thres <- max_arc <-0
for(i in seq(0.1,0.9,0.05)){
  for (j in 1:3){
    tmp <- trimDAG(mem, outputDAG, threshold_value = i, max_arc = j, plot = F)
    tmp2 <- netDiff(tmp,ref_net)['Hamming']
    if(tmp2 <= min_Hamming){
      min_Hamming = tmp2
      thres = i
      max_arc =j
    }
  }
}
tmp <- trimDAG(mem, outputDAG, threshold_value = 0.9, max_arc = 1, plot = T)

other_nets <- benchmarkNets(mem)
for(i in 1:length(other_nets))colorEdges(other_nets[[i]], ref_net, main = names(other_nets)[[i]])
all_nets <- c(list(ref=ref_net,CausFate =tmp), other_nets)
saveRDS(all_nets, file=glue("{dataset_path}/{dataset_name}.rds"))


# Output ------------------------------------------------------------------



for(n_nodes in 3:9){
  for(data_type in c("discrete","sparse","tree","dense")){
    all_nets <- readRDS(glue("{dataset_path}/{n_nodes}_{data_type}.rds"))
    for(net_type in c("ref","CausFate","mmhc","h2pc","rsmax2")){
      bn_to_svg(all_nets[[net_type]],all_nets[['ref']],file=glue("{dataset_path}/svg/{n_nodes}_{data_type}_{net_type}.svg"))
    }
  }
}

for(data_type in c("discrete","sparse","tree","dense")){
  plotAllStatsByDataType(data_type)
ggsave(glue("{dataset_path}/fig/{data_type}.pdf"),scale = 2)
}

for(data_type in c("discrete","sparse","tree","dense")){
  full_stats <- getAllStatsByDataType(data_type)
  write.csv(full_stats, glue("{dataset_path}/stats/{data_type}.csv"))
}
