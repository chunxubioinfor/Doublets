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
                        min.features = 200
                        )
dim(bm)   #check the dimension of the object
# 23480 * 87349
table(bm@meta.data$orig.ident)  #check the cell numbers of each sample


# Standard workflow
# QC
bm[["percent.mt"]] <- PercentageFeatureSet(bm, pattern = "^MT-")
VlnPlot(bm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
bm <- subset(bm, subset = nCount_RNA > 1000 & percent.mt < 10 & nFeature_RNA >500)
bm <- subset(bm, subset = nCount_RNA < 60000 & nFeature_RNA < 6000)
VlnPlot(bm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dir.create("./pbmc_filtered")
write10xCounts("./pbmc_filtered/outs/",bm[["RNA"]]@counts,version = "3")



# Filter the doublets
doublets <- read.csv('./pbmc_filtered/outs/doublets.txt',header = TRUE)
bm@meta.data$DoubletScore <- doublets$doublet_scores
bm@meta.data$DoubletPred <- doublets$predicted_doublets
doublets_list <- rownames(bm@meta.data[bm@meta.data$DoubletPred == 'True',])
singlets_list <- rownames(bm@meta.data[bm@meta.data$DoubletPred == 'False',])
bm <- subset(x = bm, subset = DoubletPred == 'False')


## Normalize the data
bm <- NormalizeData(object = bm,
                    normalization.method = "LogNormalize", 
                    scale.factor = 1e4
                    # round(median(pbmc@meta.data$nCount_RNA))
)
## Feature selection
## Identification of highly variable features and plot variable features
all.genes <- rownames(bm)
bm <- FindVariableFeatures(bm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(bm), 10)
plot1 <- VariableFeaturePlot(bm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


saveRDS(bm,'/Users/kimhan/Desktop/doublet/bm/Doublets/bm_before_scale.rds')
## Scale the data 
bm <- ScaleData(bm, features = all.genes, vars.to.regress = "percent.mt")

## Run linear dimensional reduction
bm <- RunPCA(bm, features = VariableFeatures(object = bm))

## Determine the ‘dimensionality’ of the dataset
bm <- JackStraw(bm, num.replicate = 100)
bm <- ScoreJackStraw(bm, dims = 1:20)
JackStrawPlot(bm, dims = 1:15)
ElbowPlot(bm)

# We choose 15 as the final
# Run non-linear dimensional reduction
bm <- FindNeighbors(bm, dims = 1:15)
bm <- FindClusters(bm, resolution = 0.5)
bm <- RunUMAP(bm, dims = 1:15)
DimPlot(bm, reduction = "umap",label = TRUE)
saveRDS(bm, file = "./pbmc_test2.rds")



