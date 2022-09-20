## This is a R scripts run on server ##
## Library as less as you can and save the output carefully##

library(Seurat)
library(stringr)
library(dplyr)
library(tidyverse)
library(CIMseq)
library(Cairo)

setwd('data/hancx/project/doublets/doublets_bm/')
samples <- readRDS('./input/samples.rds')
bm_origin <- readRDS('./input/bm_origin.rds')
singlets_list <- readRDS('./input/singlets_list.rds')
doublets_list <- readRDS('./input/doublets_list.rds')
all.genes <- readRDS('./input/all.genes.rds')
bm <- readRDS('./input/bm.rds')


sObj.si_list <- list()
connection_list <- list()
de_result_list <- list()
mult.pred_list <- list()

for (i in 1:length(samples)){
  sample_id <- samples[i]
  bm_sub <- subset(bm_origin,orig.ident == sample_id)
  singlets <- singlets_list[str_detect(singlets_list,sample_id)]
  doublets <- doublets_list[str_detect(doublets_list,sample_id)]
  raw_counts <- as.data.frame(bm_sub@assays$RNA@counts,stringsAsFactors = F)
  counts.sng <- select(raw_counts,singlets)[all.genes,]
  counts.sng[is.na(counts.sng)] <- 0
  counts.sng <- as.matrix(counts.sng)
  counts.sng[1:2, 1:2]
  counts.mul <- select(raw_counts,doublets)[all.genes,]
  counts.mul[is.na(counts.mul)] <- 0
  counts.mul <- as.matrix(counts.mul)
  counts.mul[1:2, 1:2]
  
  bm_singlets <- subset(bm,orig.ident == sample_id)
  classes <- as.character(as.data.frame(Idents(bm_singlets))$`Idents(bm_singlets)`)
  dim <- Embeddings(object = bm_singlets, reduction = "pca")
  features <- match(VariableFeatures(bm_singlets),rownames(bm_singlets))
  
  print(paste('Now the deconvolution is',sample_id,'is ready!'))
  cObjSng.si <-CIMseqSinglets(counts=counts.sng, classification=classes, norm.to=10000)
  cObjMul.si <- CIMseqMultiplets(counts=counts.mul, features=features,norm.to=10000)
  
  sObj.si <- CIMseqSwarm(cObjSng.si, cObjMul.si, maxiter=100, swarmsize=110, nSyntheticMultiplets=200, seed=123)
  sObj.si_list[[i]] <- sObj.si
  deconvolution_result <- sObj.si@fractions
  de_result_list[[i]] <- deconvolution_result
  
  si.edges <- calculateEdgeStats(sObj.si, cObjSng.si, cObjMul.si, multiplet.factor=2, maxCellsPerMultiplet=3)
  print(paste('The top ten connection within',sample_id,'is as below:'))
  connection_list[[i]] <- si.edges %>% filter(pval < 5e-2 & weight > 3) %>% arrange(desc(score)) %>% head(n=10)
  print('The circos plot is ready!')
  Cairo::CairoPNG(filename = paste(sample_id,'.png',sep = ''),
                  width = 7,
                  height = 7,
                  units = "in",
                  dpi = 300
  )
  plotSwarmCircos(sObj.si, cObjSng.si, cObjMul.si, weightCut=3,
                  maxCellsPerMultiplet=3, alpha=0.05, h.ratio=0.3,
                  depleted=F, multiplet.factor=2,legend=F
  )
  dev.off()
  print(paste('The circos plot of',sample_id,'is stored in',paste(sample_id,'.png',sep = ''),sep = ' '))
  print(paste('The deconvolution of',sample_id,'is done!',sep = ' '))
  mult.pred <- adjustFractions(singlets=cObjSng.si, multiplets=cObjMul.si, swarm=sObj.si, binary=T, maxCellsPerMultiplet=3,multiplet.factor = 2)
  mult.pred_list[[i]] <- mult.pred
}

saveRDS(sObj.si_list,'./output/sObj.si_list.rds')
saveRDS(connection_list,'./output/connection_list.rds')
saveRDS(de_result_list,'./output/de_result_list.rds')
