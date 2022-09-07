library(Seurat)
library(dplyr)
library(Matrix)
library(reticulate)
library(DropletUtils)
library(ggplot2)
library(celldex)
library(SingleR)
library(devtools)
library(CIMseq)
library(future)
library(tidyverse)
use_python("/usr/bin/python3")

setwd('~/Desktop/doublet/bm/')
# Merge multiple samples
fs=list.files('~/Desktop/doublet/bm/GSE120221_RAW','^GSM')
fs
samples=str_split(fs,'_',simplify = T)[,1]
lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  folder=paste0("GSE120221_RAW/", str_split(y[1],'_',simplify = T)[,1])
  dir.create(folder,recursive = T)
  file.rename(paste0("GSE120221_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  file.rename(paste0("GSE120221_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("GSE120221_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})


dir <- file.path('./GSE120221_RAW',samples)
names(dir) <- samples
#合并方法1
bm.counts <- Read10X(data.dir = dir)
pbmc.live <- CreateSeuratObject(counts = pbmc.live.data,
                                project = "pbmc_live_10X",
                                min.cells = 3,
                                min.features = 200)

# QC
pbmc.live[["percent.mt"]] <- PercentageFeatureSet(pbmc.live, pattern = "^MT-")
VlnPlot(pbmc.live, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc.live <- subset(pbmc.live, subset = nCount_RNA > 1000 & percent.mt < 10 & nFeature_RNA >500)
pbmc.live <- subset(pbmc.live, subset = nCount_RNA < 25000 & nFeature_RNA < 5000)
VlnPlot(pbmc.live, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dir.create("./pbmc_filtered")
write10xCounts("./pbmc_filtered/outs/",pbmc.live[["RNA"]]@counts,version = "3")

scRNA1 = CreateSeuratObject(counts, min.cells=1)
dim(scRNA1)   #查看基因数和细胞总数
table(scRNA1@meta.data$orig.ident)  #查看每个样本的细胞数