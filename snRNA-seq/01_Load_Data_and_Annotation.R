library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(scCustomize)

# load snRNA-seq data
# construct a data pathway table first(including the pathway and data type)
data.pathway <- read.csv("../data.detail.csv")
data.list <- list()
data.pathway.multiomes <- subset(data.pathway,data.type=="multiomes")
data.pathway.RNA <- subset(data.pathway,data.type=="snRNA")
for (i in 1:length(data.pathway.multiomes$objectID)) {
  data.list[[i]] <- Read10X(paste(data.pathway.multiomes$pathway[i], "/filtered_feature_bc_matrix",sep = ""),
                            gene.column = 2)
  data.list[[i]] <- data.list[[i]][[1]]
}
for (i in 1:length(data.pathway.RNA$objectID)) {
  data.list[[i+length(data.pathway.multiomes$objectID)]] <- Read10X(paste(data.pathway.RNA$pathway[i], "filtered_feature_bc_matrix",sep = ""),
                                                                    gene.column = 2)
}

# read gene ID and construct a assay with geneID for biomart analysis
data.ID.list <- list()
for (i in 1:length(data.pathway.multiomes$objectID)) {
  data.ID.list[[i]] <- Read10X(paste(data.pathway.multiomes$pathway[i], "/filtered_feature_bc_matrix",sep = ""),
                               gene.column = 1)
  data.ID.list[[i]] <- data.ID.list[[i]][[1]]
}
for (i in 1:length(data.pathway.RNA$objectID)) {
  data.ID.list[[i+length(data.pathway.multiomes$objectID)]] <- Read10X(paste(data.pathway.RNA$pathway[i], "filtered_feature_bc_matrix",sep = ""),
                                                                       gene.column = 1)
}

# find homology gene between human&macaque
# ues emsemble biomart
library(biomaRt)
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl")
macaque <- useMart('ensembl',dataset = "mmulatta_gene_ensembl")
# sample intersect genes
mon.ID <- intersect(rownames(data.ID.list[[5]]),rownames(data.ID.list[[10]]))
hum.ID <- intersect(rownames(data.ID.list[[1]]),rownames(data.ID.list[[8]]))
# mon.hum <- intersect(mon.ID,hum.ID)
# find macaque homology genes in human genelist
mon.hm <- getLDS(attributes = "ensembl_gene_id", values = mon.ID, filters = "ensembl_gene_id",
                 mart = macaque, attributesL = c("external_gene_name","ensembl_gene_id"), martL = human, 
                 uniqueRows = T)
# intersect gene in human&macaque gene
mon.hum.ID <- intersect(mon.hm$Gene.stable.ID.1,hum.ID)
mon.hum.hm <- mon.hm[mon.hm$Gene.stable.ID.1%in%mon.hum.ID, ]
mon.hum.hm <- filter(mon.hum.hm, Gene.name!="")
for (i in c(5,6,7,10,11)) {
  data.ID.list[[i]] <- data.ID.list[[i]][mon.hum.hm$Gene.stable.ID,]
  row.names(data.ID.list[[i]]) <- mon.hum.hm$Gene.name
}
for (i in c(1:4,8,9)) {
  data.ID.list[[i]] <- data.ID.list[[i]][mon.hum.hm$Gene.stable.ID.1,]
  row.names(data.ID.list[[i]]) <- mon.hum.hm$Gene.name
}
rm(human,macaque,mon.hm,mon.hum.hm,hum.ID,mon.hum.ID,mon.ID)

# create object
object.list <- list()
for (i in 1:length(data.pathway$objectID)) {
  object.list[[i]] <- CreateSeuratObject(counts = data.list[[i]],project = data.pathway$objectID[i])
  homo.assay <- CreateAssayObject(counts = data.ID.list[[i]])
  object.list[[i]][["homo"]] <- homo.assay
}
names(object.list) <- data.pathway$objectID
rm(data.list,data.ID.list)

# use scrublets to detect doublets
# output matrix
library(Matrix)
for (i in 1:length(object.list)) {
  sparse.mtx <- Matrix(object.list[[i]]@assays$RNA@counts,sparse = T)
  writeMM(sparse.mtx, file = paste(data.pathway$objectID[i],".matrix.mtx",sep = ""))
  system(paste("gzip ",data.pathway$objectID[i],".matrix.mtx",sep = ""))
}
# run scrublets in python
# read scrublets result to metadata
pathway <- "/sde/wuhaixv/wuhaixv/HM-RM_ACC_scRNA_YJM/2023_newlibrary/scrublets/"
for (i in 1:length(object.list)) {
  ds <- read.csv(paste(pathway,data.pathway$objectID[i],".DoubletScores.csv",sep = ""),header = F)
  pd <- read.csv(paste(pathway,data.pathway$objectID[i],".PredictedDoublets.csv",sep = ""),header = F)
  object.list[[i]]@meta.data[["DoubletScores"]] <- ds
  object.list[[i]]@meta.data[["PredictedDoublets"]] <- pd
  object.list[[i]]@meta.data[["DoubletScores"]] <- unlist(object.list[[i]]@meta.data[["DoubletScores"]])
  object.list[[i]]@meta.data[["PredictedDoublets"]] <- unlist(object.list[[i]]@meta.data[["PredictedDoublets"]])
}

# calculate mt gene percent
# human
for (i in c(1:4,8,9)) {
  object.list[[i]][["percent_mito"]] <- PercentageFeatureSet(object.list[[i]], pattern = "^MT-")
}
macaque.mt <- read.table("/sdf/wuhaixv/2023_Language/macac10_mt.gene.txt",header = T)
# macaque
for (i in c(5,6,7,10,11)) {
  object.list[[i]][["percent_mito"]] <- PercentageFeatureSet(object.list[[i]], features = macaque.mt$Gene)
}

# cluster data per library
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x,verbose = F)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000,verbose = F)
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, npcs = 30, verbose = FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:30,verbose = F)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:30,verbose = F)
  x <- FindClusters(x, resolution = 0.5,verbose = F)
})
plot.genes <- c("SLC17A7","SATB2","SYT1","NRGN","CUX2","RORB","THEMIS","FEZF2","NXPH2","BCL11B",
                "ABO", "SYT6","NPFFR2","VAT1L","ADRA1A", #Ex
                "GAD1","GAD2", "ADARB2","LAMP5","VIP","SST","LHX6","PVALB","SNCG", #In
                "AQP4","SLC1A2","GFAP","ADGRV1","MOG","MOBP","OPALIN","PLP1","MBP","PDGFRA","MYT1","PCDH15","VCAN", #Ast OPC Oli
                "APBB1IP","TYROBP","CX3CR1","CSF1R","CLDN5","FLT1","PDGFRB") #End Per Mic
for (i in 1:length(data.pathway$objectID)) {
  p1 <- DimPlot(object.list[[i]],label = T)+ NoLegend()+
    labs(title = paste("snRNA",data.pathway$objectID[i]))
  p2 <- DotPlot(object.list[[i]], features = unique(plot.genes), dot.scale = 8) + RotatedAxis()+
    labs(title = paste("snRNA",data.pathway$objectID[i]))
  p3 <- FeaturePlot(object.list[[i]],features = "PredictedDoublets")
  p4 <- FeaturePlot(object.list[[i]],features = "DoubletScores")
  p5 <- FeaturePlot(object.list[[i]],features = "nFeature_RNA")
  p6 <- FeaturePlot(object.list[[i]],features = "nCount_RNA")
  print(wrap_plots(p1, p2,p3,p4,p5,p6, ncol = 3))
}

# celltypes annotation 
# check celltype by markers expression first
cell.annotation <- read.csv("../snRNA-annotation.sample.csv",header = T)
for (i in 1:length(data.pathway$objectID)) {
  new.cluster.ids <- cell.annotation[,i+1][1:length(table(object.list[[i]]$seurat_clusters))]
  Idents(object.list[[i]]) <- "seurat_clusters"
  names(new.cluster.ids) <- levels(object.list[[i]])
  object.list[[i]] <- RenameIdents(object.list[[i]], new.cluster.ids)
  object.list[[i]]$sample.celltype <- Idents(object.list[[i]])
  object.list[[i]]$sample.seurat_clusters <- object.list[[i]]$seurat_clusters
  object.list[[i]]$sample <- object.list[[i]]$orig.ident
  object.list[[i]]$sample.celltype.detail <- paste(object.list[[i]]$sample.celltype,
                                                   object.list[[i]]$sample,
                                                   object.list[[i]]$sample.seurat_clusters,sep = "-")
  object.list[[i]]$species <- data.pathway.RNA$species[i]
  object.list[[i]]$donor <- data.pathway.RNA$donor[i]
}

# base on the data quality for each sample, set cutoff
saveRDS(object.list,"202309_analysis/object.list.rds")
object.filter.list <- list()
for (i in 1:length(object.filter.list)) {
  object.filter.list[[i]] <- subset(object.list[[i]],subset = 
                                      nFeature_RNA > 300 & 
                                      nFeature_RNA < 5000 &
                                      percent_mito < 10 & 
                                      PredictedDoublets < 1 &
                                      sample.celltype != "doublets")
}
rm(object.list)

# merge data with homology genes set
for (i in length(object.filter.list)) {
  DefaultAssay(object.filter.list[[i]]) <- "homo"
}
object.filter.list <- lapply(X = object.filter.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = object.filter.list, nfeatures = 2000)
object.filter.list <- PrepSCTIntegration(object.list = object.filter.list, anchor.features = features)
object.filter.list <- lapply(X = object.filter.list, FUN = RunPCA, features = features)
# merge data by cca
anchors <- FindIntegrationAnchors(object.list = object.filter.list, normalization.method = "SCT",anchor.features = features,
                                  reduction = "cca",dims = 1:30)
ACC.new <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(ACC.new) <- "integrated"
rm(anchors,object.filter.list)
# Run the standard workflow for visualization and clustering
ACC.new <- RunPCA(ACC.new, npcs = 30, verbose = FALSE)
ACC.new <- RunUMAP(ACC.new, reduction = "pca", dims = 1:30)
ACC.new <- RunTSNE(ACC.new, reduction = "pca", dims = 1:30)
ACC.new <- FindNeighbors(ACC.new, reduction = "pca", dims = 1:30)
ACC.new <- FindClusters(ACC.new, resolution = 3)

# annotation
# show marker genes expression
DefaultAssay(ACC.new) <- "RNA"
p <- DotPlot(ACC.new, features = unique(plot.genes),group.by="seurat_clusters",
             dot.scale = 8) + RotatedAxis()
print(p)
# show doublets score
FeaturePlot(ACC.new,features = "DoubletScores")
# make the annatation table base on "seurat_clusters"
cell.annotation <- read.csv("202309_analysis/snRNA-annotation.merge.csv",header = T)
# rename subtype
new.cluster.ids <- cell.annotation$subtype
names(new.cluster.ids) <- levels(ACC.new)
ACC.new <- RenameIdents(ACC.new, new.cluster.ids)
ACC.new$subtype <- Idents(ACC.new)
# rename majortype
Idents(ACC.new) <- "seurat_clusters"
new.cluster.ids <- cell.annotation$majortype
names(new.cluster.ids) <- levels(ACC.new)
ACC.new <- RenameIdents(ACC.new, new.cluster.ids)
ACC.new$majortype <- Idents(ACC.new)
# rename cluster
Idents(ACC.new) <- "seurat_clusters"
new.cluster.ids <- cell.annotation$subtype_seurat_cluster
names(new.cluster.ids) <- levels(ACC.new)
ACC.new <- RenameIdents(ACC.new, new.cluster.ids)
ACC.new$subtype_seurat_cluster <- Idents(ACC.new)

saveRDS(ACC.new,file = ".../ACC.merge.filter.rds")

