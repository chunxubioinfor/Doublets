library(Seurat)
library(psych)
library(dplyr)
library(reshape2)
library(ggplot2)
library(Ipaper)
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
    if (gene_list[i] %in% rownames(bm_avg_expr)){
      expr <- c(expr,bm_avg_expr[gene_list[i],cell])
    } else {
      expr <- c(expr,0)
    }
  }
  gm_mean <- geometric.mean(expr,na.rm = TRUE)
  return(gm_mean)
}

cal_expr_score <- function(ligand,receptor,transmitter,receiver,method){
  ligand_list <- unlist(strsplit(ligand,','))
  receptor_list <- unlist(strsplit(receptor,','))
  l_expr <- cal_gm_mean(ligand_list,transmitter)
  r_expr <- cal_gm_mean(receptor_list,receiver)
  if (method == 'counts'){
    expr_score <- l_expr * r_expr
  } else if (method == 'Kd'){
    Kd_nM <- filter(interactions_csv,Gene_L == ligand,Gene_R == receptor)$Kd_nM
    if (is.na(Kd_nM)){
      expr_score <- 0
    } else {
      expr_score <- (l_expr * r_expr) / Kd_nM
    }
  }
  return(expr_score)
}


cal_intxn_score <- function(lr_pairs,cell_m,cell_n,method){
  expr_score <- c()
  for (i in 1:nrow(lr_pairs)){
    print(paste('Hi',i))
    ligand <- lr_pairs[i,]$Gene_L
    receptor <- lr_pairs[i,]$Gene_R
    # if cell_m is the transmitter while cell_n is the receiver
    expr_score_1 <- cal_expr_score(ligand,receptor,cell_m,cell_n,method)
    # else if cell_m is the receiver while cell_n is the transmitter
    expr_score_2 <- cal_expr_score(ligand,receptor,cell_n,cell_m,method)
    expr_score <-c(expr_score,expr_score_1 + expr_score_2)
  }
  intxn_score <- sum(expr_score)
  return(intxn_score)
}

intxn_mtx <- matrix(0,nrow = 16,ncol = 16,dimnames = list(sort(levels(bm)),sort(levels(bm))))
for (i in 1:nrow(intxn_mtx)){
  for(j in 1:ncol(intxn_mtx)){
    intxn_score <- cal_intxn_score(interactions_csv,cell_m = rownames(intxn_mtx)[i],cell_n = colnames(intxn_mtx)[j],method = 'Kd')
    doublets_name <- paste(sort(c(rownames(intxn_mtx)[i],colnames(intxn_mtx)[j])),collapse = '_')
    intxn_mtx[i,j] <- intxn_score
  }
}
write.csv(intxn_mtx,'./intxn_mtx_counts.csv')
write.csv(intxn_mtx,'./intxn_mtx_Kd.csv')

melted_intxn_df <- melt(intxn_mtx)
p1 <-  ggplot(data = melted_intxn_df,aes(x=Var1,y=Var2,fill=value)) +
    geom_tile(height = -1,width = 1) +
    scale_fill_gradient2(low = '#000000',high = '#FFEA45',mid = '#21325A',midpoint = 15) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1,size = 13,face = 'bold'),axis.text.y = element_text(size = 10,face = 'bold')) +
  labs(x=NULL,y=NULL)
write_fig(p1,'./hm_counts.png',width = 11,height = 10,res = 300,show = FALSE)
write_fig(p1,'./hm_Kd.png',width = 11,height = 10,res = 300,show = FALSE)

