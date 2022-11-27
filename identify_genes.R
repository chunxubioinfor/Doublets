# Identify possible surface elements 
# Suppose a doublet pair
# NK cells and Monocytes

library(Seurat)
library(CIMseq)
library(plyr)
library(dplyr)
library(BiocManager)
library(scde)
BiocManager::install("scde")
setwd('~/Desktop/doublet/bm/Doublets/')
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
NK_mono_doublets <- filter(doublets_df,`NK cell` == 1,`Monocyte` == 1)$barcode



