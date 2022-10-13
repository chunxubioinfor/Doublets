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
library(clustree)
library(ggplot2)
library(ggpubr)
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
bm_origin = CreateSeuratObject(bm.counts,
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
bm <- FindClusters(bm, resolution = 0.9)
bm <- RunUMAP(bm, dims = 1:15)
DimPlot(bm, reduction = "umap",label = TRUE)
saveRDS(bm, file = "./bm_before_annotation_rs.rds")

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
write.csv(bm_sglr_main,'./bm_sglr_main_rs.csv')
bm_sglr_fine <- table(bm.hesc.fine$labels,bm@meta.data$seurat_clusters)
write.csv(bm_sglr_fine,'./bm_sglr_fine_rs.scv')

# Manul annotation
top_150_markers <- data.frame(matrix(0,nrow = 150,ncol = length(levels(bm@meta.data$seurat_clusters))))
for (i in 1:length(levels(bm@meta.data$seurat_clusters))){
  cluster_markers <- FindMarkers(bm,ident.1 = i-1,min.pct = 0.25,only.pos = TRUE)
  if (nrow(cluster_markers) >= 150){
    top_150_markers[,i] <- head(rownames(cluster_markers),150)
  } else {
    top_150_markers[1:nrow(cluster_markers),i] <- row.names(cluster_markers)
  }
  
}
top_150_markers

table(bm@active.ident)
new.cluster.ids <- c("CD4+ memory T cells", "CD4+ memory T cells", 'CD8+ T cells', 'CD4+ naive T cells',"CD14+ monocytes",
                     'B cells','CD8+ T cells','Late erythroid progenitors','NK cells','Early erythrocytes','Neutraphil','Early erythroid progenitors',
                     'Late erythrocytes','CD4+ memory T cells','Megakaryocytes','CD16+ Monocytes','Myeloid dendritic cells','B progenitor cells',
                     'Late erythrocytes','Plasma cells','Plasmacytoid dendritic cells','B progenitor cells','Macrophage',
                     'Early erythroid progenitors','Late erythrocytes','NK T cells','B progenitor cells','Stromal cells')
names(new.cluster.ids) <- levels(bm)
bm <- RenameIdents(bm, new.cluster.ids)
DimPlot(bm, reduction = "umap", label = TRUE, pt.size = 0.6,label.size = 4)

# remove the unknown cluster
bm <- subset(bm,idents = 'Unknown',invert = TRUE)
DimPlot(bm, reduction = "umap",label = TRUE,pt.size = 0.5)

sink('./output.txt')
sObj.si_list <- list()
circos_plot_list <- list()
connection_list <- list()
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
}
sink()

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


library(Seurat)
setwd('/data/hancx/project/doublets/doublets_bm/')
bm <- readRDS('./input/bm.rds')
samples <- readRDS('./input/samples.rds')

for (i in 1:25){
  sample_id <- samples[i]
  bm_singlets <- subset(bm,orig.ident == sample_id)
  file_name <- paste('./singlets_rds/bm',i,'.rds',sep = '')
  saveRDS(bm_singlets,file = file_name)
}

bm_dc_df <- data.frame()
setwd('./output')
file <- list.files(path = './')
for (i in 1:length(file)){
  result <- readRDS(file[i])
  result_df <- result[[3]]
  result_df$merge <- apply(result_df[,c(1:2)],1,function(x) paste(sort(x),collapse ='_'))
  result_df<-result_df[!duplicated(result_df$merge),]
  bm_dc_df <- rbind(bm_dc_df,result_df)
}
count_df <- bm_dc_df %>% select(c('merge','weight'))
count_table <- as.data.frame(table(count_df$merge))
count_table$ymin <- rep(0,22)
count_df <- aggregate(count_df$weight,by = list(count_df$merge),sum)

theme <- theme(panel.background = element_blank(), # 去掉背景格子
               # 显示x平行网格线
               panel.grid.major.x = element_line(colour = "black"), 
               # 显示x轴坐标
               axis.line.x = element_line(colour = "black"),
               axis.title.y = element_blank(),
               axis.text =  element_text(size = 10,face = 'bold'))

ggplot(data = count_table) + geom_segment(aes(x = reorder(Var1,Freq), y = ymin,
                                              xend =reorder(Var1,Freq), yend = Freq)) + ylab('Freqency') +
  geom_point(aes(x = reorder(Var1,Freq), y = Freq), 
             size = 6,colour = '#db5461') +
  scale_y_continuous(breaks=seq(0, 10, 2),expand = c(0,0),limits = c(0,10))+
  coord_flip() + theme 


ggplot(data = count_df) + geom_bar(aes(x = reorder(Group.1,x), y = x),stat = "identity",fill = '#3d5467')  +
  scale_y_continuous(expand = c(0,0),limits = c(0,100))+
  coord_flip() + theme

table(bm$orig.ident)
table(bm$seurat_clusters)
table(Idents(bm))
prop.table(table(Idents(bm)))
cell.prop <- as.data.frame(prop.table(table(Idents(bm))))
ggdonutchart(cell.prop,'Freq',
             label = 'Var1',
             fill = 'Var1')
