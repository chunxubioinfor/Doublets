# Identify possible surface elements 
# Suppose a doublet pair
# NK cells and Monocytes

library(Seurat)
library(CIMseq)
library(plyr)
library(dplyr)
library(BiocManager)
library(scde)
library(stringr)
BiocManager::install("scde")
setwd('~/Desktop/doublet/bm/Doublets/')
options(stringsAsFactors = F)
# Devide the NK cells into NK-M-doublets-specific group and singlets group
bm.origin <- readRDS('./input_old/bm_origin.rds')
output_path <- './output_v3/output/'
file <- list.files(path = output_path)
doublets_df <- data.frame()
for (i in 1:length(file)){
  print(i)
  result <- readRDS(paste(output_path,file[i],sep = ''))
  result_df <- as.data.frame(result[[4]])
  result_df$barcode <- row.names(result_df)
  if (i == 1){
    doublets_df <- result_df
  } else {
    doublets_df <- rbind.fill(doublets_df,result_df)
  }
}
NK_cell_singlets <- WhichCells(bm,idents = 'NK cell')
NK_cell_singlets <- unlist(lapply(NK_cell_singlets,function(x) paste(x,'-1',sep = '')))
NK_mono_doublets <- filter(doublets_df,`NK cell` == 1,`Monocyte` == 1)$barcode

cell_name <- vector(mode = 'numeric',length = nrow(bm.origin@meta.data))
names(cell_name) <- row.names(bm.origin@meta.data)
bm.origin <- AddMetaData(bm.origin,cell_name,col.name = 'cell_name')
cell_name[NK_cell_singlets] <- 'NK_cell_singlets'
cell_name[NK_mono_doublets] <- 'NK_mono_doublets'
bm_matrix <- as.data.frame(bm.origin@assays$RNA@counts)

data(es.mef.small)

