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
library(ggplot2)
library(ggrepel)
library(Ipaper)
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
cell_name[NK_cell_singlets] <- 'NK_cell_singlets'
cell_name[NK_mono_doublets] <- 'NK_mono_doublets'
bm.origin <- AddMetaData(bm.origin,cell_name,col.name = 'cell_name')
NK_mono_db_genes <- subset(bm.origin,subset = (cell_name == 'NK_mono_doublets' | cell_name == 'NK_cell_singlets'))
NK_mono_mtx <- as.data.frame(NK_mono_db_genes@assays$RNA@counts)
dim(NK_mono_mtx)
# clean up the dataset
NK_mono_mtx <- clean.counts(NK_mono_mtx, min.lib.size=1000, min.reads = 1, min.detected = 1)
# prepare the classification
NK_mono_df <- data.frame(cell = c(NK_cell_singlets,NK_mono_doublets),cell_type = vector(mode = 'numeric',length = 3072))
NK_mono_df[NK_mono_df$cell %in% NK_cell_singlets,]$cell_type <- 'NK_cell_singlets'
NK_mono_df[NK_mono_df$cell %in% NK_mono_doublets,]$cell_type <- 'NK_mono_doublets'
sg <- as.factor(NK_mono_df$cell_type)
names(sg) <- NK_mono_df$cell
# calculate models
# but it took too much time, so run this step on the. server
o.ifm <- scde.error.models(counts = NK_mono_mtx, groups = sg, n.cores = 1, threshold.segmentation = TRUE, 
                           save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
o.ifm <- readRDS('./o.ifm.rds')
head(o.ifm)
data(es.mef.small)
NK_mono_mtx<-apply(NK_mono_mtx,2,function(x) {storage.mode(x) <- 'integer'; x})

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = NK_mono_mtx, length.out = 400, show.plot = FALSE)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, NK_mono_mtx, o.prior, groups  =  sg, 
                                    n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# top up-regulated genes
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])

# top down-regulated genes
tail(ediff[order(ediff$Z, decreasing  =  TRUE), ])

# # write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "de_results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# convert the z score to p value
for (i in 1:nrow(ediff)){
  Z <- ediff$Z[i]
  if (Z > 0){
    ediff$p[i] <- pnorm(Z,lower.tail=F)
  }
  else {
    ediff$p[i] <- pnorm(Z,lower.tail=T)
    }
}

# volcano plot
ediff$threshold = as.factor(ifelse(ediff$p < 0.001 & abs(ediff$mle) >= 2, ifelse(ediff$mle > 2,'Up','Down'),'NoSignifi'))
p <- ggplot(data = ediff, aes(x = mle, y = -log10(p), colour=threshold)) +
  geom_point(alpha=0.4, size=1) +
  scale_color_manual(values=c("blue", "grey","red")) +
  labs(x="log2(fold change)",y="-log10 (p-value)")
  
write_fig(p,'./v.png')

# highlight the surface ligands and receptors
ligand_receptor <- readRDS('~/Desktop/CCI_DK/ccc/ligand_receptor.rds')
ediff_lr <- ediff[ligand_receptor,]
ediff_lr$gene <- rownames(ediff_lr)
ediff_lr <- na.omit(ediff_lr)
p1 <- ggplot(data = ediff_lr, aes(x = mle, y = -log10(p), colour=threshold)) +
  geom_point(alpha=0.4, size=1,na.rm = TRUE) +
  scale_color_manual(values=c("blue", "grey","red")) +
  labs(x="log2(fold change)",y="-log10 (p-value)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
  geom_text_repel(
    data = subset(ediff_lr, ediff_lr$p < 0.0001 & abs(ediff_lr$mle) >= 4),
    aes(label = gene),
    size = 2,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE,
    max.overlaps = Inf)
write_fig(p1,'./v2.png')
