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

samples=list.files("~/Desktop/doublet/bm/GSE120221_RAW/")
dir <- file.path('~/Desktop/doublet/bm/GSE120221_RAW',samples)
names(dir) <- samples
bm.counts <- Read10X(data.dir = dir)
bm = CreateSeuratObject(bm.counts,
                        project = 'bone marrow',
                        min.cells=3,
                        )
dim(scRNA1)   #查看基因数和细胞总数
table(scRNA1@meta.data$orig.ident)  #查看每个样本的细胞数


# Standard workflow
## Normalize the data
pbmc.live <- NormalizeData(object = pbmc.live,
                           normalization.method = "LogNormalize", 
                           scale.factor = 1e4
                           # round(median(pbmc@meta.data$nCount_RNA))
)

## Feature selection
## Identification of highly variable features and plot variable features
all.genes <- rownames(pbmc.live)
pbmc.live <- FindVariableFeatures(pbmc.live, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc.live), 10)
plot1 <- VariableFeaturePlot(pbmc.live)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

## Scale the data 
pbmc.live <- ScaleData(pbmc.live, features = all.genes, vars.to.regress = "percent.mt")

## Run linear dimensional reduction
pbmc.live <- RunPCA(pbmc.live, features = VariableFeatures(object = pbmc.live))

## Determine the ‘dimensionality’ of the dataset
pbmc.live <- JackStraw(pbmc.live, num.replicate = 100)
pbmc.live <- ScoreJackStraw(pbmc.live, dims = 1:20)
JackStrawPlot(pbmc.live, dims = 1:15)
ElbowPlot(pbmc.live)

# We choose 15 as the final
# Run non-linear dimensional reduction
pbmc.live <- FindNeighbors(pbmc.live, dims = 1:15)
pbmc.live <- FindClusters(pbmc.live, resolution = 0.5)
pbmc.live <- RunUMAP(pbmc.live, dims = 1:15)
DimPlot(pbmc.live, reduction = "umap",label = TRUE)
saveRDS(pbmc.live, file = "./pbmc_test2.rds")

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