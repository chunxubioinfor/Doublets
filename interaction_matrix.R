library(Seurat)
library(psych)
## Interaction Matrix ##
setwd('~/Desktop/doublet/bm/Doublets/')
interactions_csv <- read.csv('~/Desktop/doublet/wiring/Leuk_interaction_affinities.csv',header = TRUE)

## Define a function to calculate the average expression of each cell type
table(Idents(bm))
bm.averages <- AverageExpression(bm, slot = 'data', return.seurat = TRUE)
# Three modes: RNA, SCT, integrated
dim(bm.averages@assays$RNA@data)
head(bm.averages@assays$RNA@data)

## Define a function to calculate the interaction between cell type m and cell type n
## The basic concept is the sum of all pairs of expression score divided by Kd
## But some are without Kd
## And 
lignad_receptor <- c(interactions_csv$Gene_L,interactions_csv$Gene_R)

pair_labeler <- function(){
  
}


cal_intxn_score <- function(transmitter_expr,receiver_expr,ligand,receptor){

}