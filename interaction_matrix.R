library(Seurat)
library(psych)
library(dplyr)
## Interaction Matrix ##
setwd('~/Desktop/doublet/bm/Doublets/')
interactions_csv <- read.csv('~/Desktop/doublet/wiring/Leuk_interaction_affinities.csv',header = TRUE)

## Define a function to calculate the average expression of each cell type
table(Idents(bm))
bm.averages <- AverageExpression(bm, slot = 'data', return.seurat = TRUE)
# Three modes: RNA, SCT, integrated
dim(bm.averages@assays$RNA@data)
head(bm.averages@assays$RNA@data)
bm_avg_expr <- bm.averages@assays$RNA@data
write.csv(bm_avg_expr, file = "./bm_avg_expr.csv")

## Define a function to calculate the interaction between cell type m and cell type n
## The basic concept is the sum of all pairs of expression score divided by Kd
## But some are without Kd
## And 
cal_gm_mean <- function(gene_list,cell){
  expr <- c()
  for (i in 1:length(gene_list)){
    expr <- c(expr,bm_avg_expr[gene_list[i],cell])
  }
  gm_mean <- geometric.mean(expr)
  return(gm_mean)
}

cal_expr_score <- function(ligand,receptor,transmitter,receiver,method){
  ligand <- unlist(strsplit(ligand,','))
  receptor <- unlist(strsplit(receptor,','))
  l_expr <- cal_gm_mean(ligand,transmitter)
  r_expr <- cal_gm_mean(receptor,receiver)
  if (method == 'counts'){
    expr_score <- l_expr * r_expr
    return(expr_score)
  }
  else if (method == 'Kd'){
    Kd_nM <- filter(interactions_csv,gene_L == ligand,gene_R == receptor)$Kd_nM
    expr_score <- (l_expr * r_expr) / Kd_nM
    return(expr_score)
  }
}


cal_intxn_score <- function(transmitter_expr,receiver_expr,ligand,receptor,method){

}