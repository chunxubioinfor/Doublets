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
library(SCOPfunctions)
library(stringr)
library(customLayout)
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
dir.create("./bm_filtered")
write10xCounts("./bm_filtered/outs/",bm[["RNA"]]@counts,version = "3")



# Filter the doublets
doublets <- read.csv('./bm_filtered/outs/doublets.txt',header = TRUE)
bm@meta.data$DoubletScore <- doublets$doublet_scores
bm@meta.data$DoubletPred <- doublets$predicted_doublets
doublets_list <- rownames(bm@meta.data[bm@meta.data$DoubletPred == 'True',])
singlets_list <- rownames(bm@meta.data[bm@meta.data$DoubletPred == 'False',])
bm <- subset(x = bm, subset = DoubletPred == 'False')


## Normalize the data
bm <- NormalizeData(object = bm,
                    normalization.method = "LogNormalize", 
                    scale.factor = 1e4
                    # round(median(bm@meta.data$nCount_RNA))
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
saveRDS(bm, file = "./bm_before_annotation.rds")

# Identify the marker genes
# Automated annotation using SingleR
ref_MI <- celldex::MonacoImmuneData()
ref_HPCA <- celldex::HumanPrimaryCellAtlasData()
bm <- readRDS('./bm_before_annotation.rds')
bm_for_SingleR <- GetAssayData(bm, slot="data")
bm.hesc.main <- SingleR(test = bm_for_SingleR, ref = list(ref_HPCA,ref_MI), labels = list(ref_HPCA$label.main,ref_MI$label.main))
bm.hesc.fine <- SingleR(test = bm_for_SingleR, ref = list(ref_HPCA,ref_MI), labels = list(ref_HPCA$label.fine,ref_MI$label.fine))
bm.hesc.main
bm.hesc.fine
bm_sglr_main <- table(bm.hesc.main$labels,bm@meta.data$seurat_clusters)
write.csv(bm_sglr_main,'./bm_sglr_main.csv')
bm_sglr_fine <- table(bm.hesc.fine$labels,bm@meta.data$seurat_clusters)
write.csv(bm_sglr_fine,'./bm_sglr_fine.scv')

# Manul annotation
top_100_markers <- data.frame(matrix(0,nrow = 100,ncol = length(levels(bm@meta.data$seurat_clusters))))
for (i in 1:length(levels(bm@meta.data$seurat_clusters))){
  cluster_markers <- FindMarkers(bm,ident.1 = i-1,min.pct = 0.25)
  top_100_markers[,i] <- head(rownames(cluster_markers),100)
}
top_100_markers

table(bm@active.ident)
new.cluster.ids <- c("CD4+ naive T cells", "CD4+ memory T cells", 'CD8+ T cells', 'CD14+ monocytes',"B cells",'Late erythroid progenitors',
                     'Early erythrocytes','Hematopoietic stem cell','Early erythroid progenitors','Late erythrocytes','B cell progenitors',
                     'CD16+ monocytes','Granulocyte/monocyte progenitor cells','Myeloid dendritic cells','Late erythrocytes','Plasma cells',
                     'Plasmacytoid dendritic cells','Natural killer cells','B cells')
names(new.cluster.ids) <- levels(bm)
bm <- RenameIdents(bm, new.cluster.ids)
DimPlot(bm, reduction = "umap", label = TRUE, pt.size = 0.6,label.size = 4)
#DimPlot(bm, reduction = "tsne", label = TRUE, pt.size = 0.6,label.size = 4)

# remove the unknown cluster
bm <- subset(bm,idents = 'Unknown',invert = TRUE)
DimPlot(bm, reduction = "umap",label = TRUE,pt.size = 0.5)

sObj.si_list <- list()
circos_plot_list <- list()
lay1 <- lay_new(mat = matrix(1:25, ncol = 5))
lay_set(lay1)
for (i in 1:length(samples)){
  sample_id <- samples[i]
  bm_sub <- subset(bm,subset = orig.ident = i)
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
  
  classes <- as.character(as.data.frame(Idents(bm_sub))$`Idents(bm_sub)`)
  dim <- Embeddings(object = bm_sub, reduction = "pca")
  features <- match(VariableFeatures(bm_sub),rownames(bm_sub))
  
  print('Now the deconvolution of ',sample_id,'is ready!')
  cObjSng.si <-CIMseqSinglets(counts=counts.sng, classification=classes, norm.to=10000)
  cObjMul.si <- CIMseqMultiplets(counts=counts.mul, features=features,norm.to=10000)
  
  sObj.si <- CIMseqSwarm(cObjSng.si, cObjMul.si, maxiter=100, swarmsize=110, nSyntheticMultiplets=200, seed=123)
  sObj.si_list[i] <- sObj.si
  deconvolution_result <- sObj.si@fractions
  
  si.edges <- calculateEdgeStats(sObj.si, cObjSng.si, cObjMul.si, multiplet.factor=2, maxCellsPerMultiplet=3)
  si.edges %>% filter(pval < 5e-2 & weight > 3) %>% arrange(desc(score)) %>% head(n=10)
  plotSwarmCircos(sObj.si, cObjSng.si, cObjMul.si, weightCut=3,
                  maxCellsPerMultiplet=3, alpha=0.05, h.ratio=0.5,
                  depleted=F, multiplet.factor=2)
  print('The circos plot of ',sample_id,'is done!')
  print('The deconvolution of ',sample_id,'is done!')
  mult.pred <- adjustFractions(singlets=cObjSng.si, multiplets=cObjMul.si, swarm=sObj.si.k, binary=T, maxCellsPerMultiplet=3,multiplet.factor = 2)
  
}

# CIMseq deconvolution of doublets
## prepare for the four arguments
raw_counts <- as.data.frame(bm.counts,stringsAsFactors = F)    #get the raw counts matrix
counts.sng <- select(raw_counts,singlets_list)[all.genes,]
counts.sng[is.na(counts.sng)] <- 0
counts.sng <- as.matrix(counts.sng)
counts.sng[1:2, 1:2]
counts.mul <- select(raw_counts,doublets_list)[all.genes,]
counts.mul[is.na(counts.mul)] <- 0
counts.mul <- as.matrix(counts.mul)
counts.mul[1:2, 1:2]

raw_counts <- utils_big_as.matrix(bm.counts, n_slices_init = 1, verbose = T)

classes <- as.character(as.data.frame(Idents(bm))$`Idents(bm)`)
dim <- Embeddings(object = bm, reduction = "pca")
features <- match(VariableFeatures(bm),rownames(bm))

cObjSng.si <-CIMseqSinglets(counts=counts.sng, classification=classes, norm.to=10000)
cObjMul.si <- CIMseqMultiplets(counts=counts.mul, features=features,norm.to=10000)

ercc.s <- grepl("^ERCC\\-[0-9]*$", rownames(counts.sng))
singlets <- counts.sng[!ercc.s, ]
singletERCC <- counts.sng[ercc.s, ]
ercc.m <- grepl("^ERCC\\-[0-9]*$", rownames(counts.mul))
multiplets <- counts.mul[!ercc.m, ]
multipletERCC <- counts.mul[ercc.m, ]

plan(multicore)
options(future.globals.maxSize= 1000000000)
sObj.si <- CIMseqSwarm(cObjSng.si, cObjMul.si, maxiter=100, swarmsize=110, nSyntheticMultiplets=200, seed=123)
deconvolution_result <- sObj.si@fractions
deconvolution_result

si.edges <- calculateEdgeStats(sObj.si, cObjSng.si, cObjMul.si, multiplet.factor=2, maxCellsPerMultiplet=3)
si.edges %>% filter(pval < 5e-2 & weight > 3) %>% arrange(desc(score)) %>% head(n=10)
plotSwarmCircos(sObj.si, cObjSng.si, cObjMul.si, weightCut=3,
                maxCellsPerMultiplet=3, alpha=0.05, h.ratio=0.5,
                depleted=F, multiplet.factor=2)
mult.pred <- adjustFractions(singlets=cObjSng.si, multiplets=cObjMul.si, swarm=sObj.si.k, binary=T, maxCellsPerMultiplet=3,multiplet.factor = 2)


