library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(scCustomize)

# subset EX neurons from ACC data
ACC.new <- readRDS(".../ACC.merge.filter.rds")
Idents(ACC.new) <- "subtype" 
Ex <- subset(ACC.new, idents=levels(ACC.new$subtype)[c(5,8,11,13,14,16,20)])

# VEN related genes (pineline show in Figure S3F)
# find DEG between hum & mon (using Wilicoxen test)
DefaultAssay(Ex) <- "homo"
Ex$subtype_species <- paste(Ex$subtype,Ex$species,sep = "_")
Idents(Ex) <- "subtype_species"
ven.DEG <- FindMarkers(Ex,ident.1 = "EX-ET_Human",ident.2 = "EX-ET_Macaque")
ven.DEG$change = ifelse(ven.DEG$p_val_adj < 0.001 & abs(ven.DEG$avg_log2FC) >= 1, 
                        ifelse(ven.DEG$avg_log2FC> 1 ,'Human','Macaque'), 'Stable')
ven.DEG$label <- ifelse(ven.DEG$p_val_adj < 0.001 & abs(ven.DEG$avg_log2FC) >= 2,
                        as.character(rownames(ven.DEG)), "")
ven.DEG <- ven.DEG[-grep("MT-", rownames(ven.DEG)),]
# find marker genes from species
Idents(Ex) <- "species"
Ex.human <- subset(Ex, idents="Human")
Ex.macaque <- subset(Ex, idents="Macaque")
Idents(Ex.human) <- "subtype"
Idents(Ex.macaque) <- "subtype"
ven.hum.markers <-  FindMarkers(Ex.human,ident.1 = "EX-ET",verbose = F)
ven.hum.markers$gene <- rownames(ven.hum.markers)
ven.mon.markers <-  FindMarkers(Ex.macaque,ident.1 = "EX-ET",verbose = F)
ven.mon.markers$gene <- rownames(ven.mon.markers)

# human and macaque shared marker genes
hm.SMG.hum <- subset(ven.hum.markers,subset = avg_log2FC>1 & p_val_adj <0.05 & pct.1>0.25)
hm.SMG.mon <- subset(ven.mon.markers,subset = avg_log2FC>1 & p_val_adj <0.05 & pct.1>0.25)
hm.SMG <- intersect(rownames(hm.SMG.hum),rownames(hm.SMG.mon))
colnames(hm.SMG.hum) <- paste(colnames(hm.SMG.hum),"_Human",sep = "")
colnames(hm.SMG.mon) <- paste(colnames(hm.SMG.mon),"_Macaque",sep = "")
ven.SMG <- cbind(hm.SMG.hum[hm.SMG,],hm.SMG.mon[hm.SMG,])
write.csv(ven.SMG,"../ven.SMG.csv")

# compare VEN marker gene to each subtype one by one
ven.SMG.eachcelltype <- ven.SMG
for (i in 1:6) {
  temp.human <- FindMarkers(Ex.human,ident.1 = "EX-ET",ident.2 = levels(Ex.human$subtype)[i],verbose = F)
  temp.human$gene <- rownames(temp.human)
  colnames(temp.human) <- paste(colnames(temp.human),"_Human_",levels(Ex.human$subtype)[i],sep = "")
  temp.human <- temp.human[hm.SMG,]
  temp.macaque <- FindMarkers(Ex.macaque,ident.1 = "EX-ET",ident.2 = levels(Ex.macaque$subtype)[i],verbose = F)
  temp.macaque$gene <- rownames(temp.macaque)
  colnames(temp.macaque) <- paste(colnames(temp.macaque),"_Macaque_",levels(Ex.macaque$subtype)[i],sep = "")
  temp.macaque <- temp.macaque[hm.SMG,]
  ven.SMG.eachcelltype <- cbind(ven.SMG.eachcelltype,temp.human)
  ven.SMG.eachcelltype <- cbind(ven.SMG.eachcelltype,temp.macaque)
}
ven.SMG.eachcelltype.filter <- ven.SMG.eachcelltype[complete.cases(ven.SMG.eachcelltype),]
write.csv(ven.SMG.eachcelltype.filter,"../ven.SMG.eachcelltype.filter.csv")
ven.SMG.eachcelltype.filter <- read.csv("2023_newlibrary/202309_analysis/ven.SMG.eachcelltype.filter.csv",row.names = 1)

# VEN DEG between human and macaque
ven.DEG.filter <- subset(ven.DEG,subset = abs(avg_log2FC) >1, p_val_adj <0.05)
hm.DEG <- union(intersect(rownames(hm.SMG.hum),rownames(ven.DEG.filter)),
                intersect(rownames(hm.SMG.mon),rownames(ven.DEG.filter)))
ven.DEG.filter <- ven.DEG[hm.DEG,]
ven.DEG.filter.order <- ven.DEG.filter[order(ven.DEG.filter$avg_log2FC,decreasing = T),]

saveRDS(Ex,"2023_newlibrary/202309_analysis/ACC.EX.rds")

# human and macaque shared marker genes
# set VEN marker genes cutoff 
options(Seurat.object.assay.version = 'v3')
library(Libra)
Ex <- readRDS("2023_newlibrary/202309_analysis/ACC.EX.rds")
Idents(Ex) <- "species"
Ex.human <- subset(Ex, idents="Human")
Ex.macaque <- subset(Ex, idents="Macaque")
Idents(Ex.human) <- "subtype"
Idents(Ex.macaque) <- "subtype"
ven.hum.markers <-  FindMarkers(Ex.human,ident.1 = "EX-ET",verbose = F,
                                min.pct = 0,logfc.threshold = 0)
ven.hum.markers$gene <- rownames(ven.hum.markers)
ven.mon.markers <-  FindMarkers(Ex.macaque,ident.1 = "EX-ET",verbose = F,
                                min.pct = 0,logfc.threshold = 0)
ven.mon.markers$gene <- rownames(ven.mon.markers)
ven.hum.markers <- subset(ven.hum.markers,avg_log2FC>1 & p_val_adj<0.05)
ven.mon.markers <- subset(ven.mon.markers,avg_log2FC>1 & p_val_adj<0.05)

# analysis ET marker genes using Libra (pesudobulk analysis by edgeR)
# human VEN marker genes 
Ex.human$subtype_ET <- droplevels(Ex.human$subtype)
Ex.human$sample <- droplevels(Ex.human$sample)
levels(Ex.human$subtype_ET) <- c(rep("otherEX",6),"EX-ET")
Ex.human$sample.subtype <- paste(Ex.human$sample,Ex.human$subtype_ET)
DE.edgR <- run_de(Ex.human, cell_type_col = "species", label_col = "subtype_ET",
                  replicate_col = "sample.subtype")
write.csv(DE.edgR,"2023_newlibrary/202309_analysis/20240224_ETmarker/human.edgeR.csv")
# macaque VEN marker genes
Ex.macaque$subtype_ET <- droplevels(Ex.macaque$subtype)
Ex.macaque$sample <- droplevels(Ex.macaque$sample)
levels(Ex.macaque$subtype_ET) <- c(rep("otherEX",6),"EX-ET")
Ex.macaque$sample.subtype <- paste(Ex.macaque$sample,Ex.macaque$subtype_ET)
DE.edgR <- run_de(Ex.macaque, cell_type_col = "species", label_col = "subtype_ET",
                  replicate_col = "sample.subtype")
write.csv(DE.edgR,"2023_newlibrary/202309_analysis/20240224_ETmarker/macaque.edgeR.csv")

# set VEN marker genes cutoff 
ven.hum.markers.edgeR <- read.csv(("2023_newlibrary/202309_analysis/20240224_ETmarker/human.edgeR.csv"))
ven.mon.markers.edgeR <- read.csv(("2023_newlibrary/202309_analysis/20240224_ETmarker/macaque.edgeR.csv"))
ven.hum.markers.edgeR <- subset(ven.hum.markers.edgeR,subset = avg_logFC<-1&p_val_adj<0.05)
ven.mon.markers.edgeR <- subset(ven.mon.markers.edgeR,subset = avg_logFC<-1&p_val_adj<0.05)

# intersect wilicoxen test and edgeR result —— human and macaque shared marker genes
hum.markers <- intersect(ven.hum.markers$gene,ven.hum.markers.edgeR$gene)
mon.markers <- intersect(ven.mon.markers$gene,ven.mon.markers.edgeR$gene)

SMG <- intersect(hum.markers,mon.markers)
specific <- setdiff(union(hum.markers,mon.markers),SMG)

# human and macaque DEGs
# wilicoxen test
Idents(Ex) <- "subtype_species"
ven.DEG <- FindMarkers(Ex,ident.1 = "EX-ET_Human",ident.2 = "EX-ET_Macaque",min.pct = 0.25)
ven.DEG$gene <- rownames(ven.DEG)
ven.DEG <- subset(ven.DEG, gene%in%specific & abs(avg_log2FC)>1&p_val_adj<0.05)

# edgeR pesudobulk test
DE.edgR = run_de(Ex, cell_type_col = "subtype", label_col = "species",
                 replicate_col = "sample")
DE.DEseq = run_de(EX, cell_type_col = "subtype", label_col = "species",
                  replicate_col = "sample",de_family = 'pseudobulk', de_method = 'DESeq2')

write.csv(DE.edgR,"../DEG.edgeR.csv")
write.csv(DE.DEseq,"../DEG.DEseq.csv")
ven.DEG.edgeR <- read.csv("../DEG.edgeR.csv")
ven.DEG.edgeR <- subset(ven.DEG.edgeR, abs(avg_logFC)>1&p_val_adj<0.05&cell_type=="EX-ET"&gene%in%specific)

# add mouse data 
mon.acc<- readRDS("../mouse.aca.rds")
mon.acc.meta <- filter(mon.acc@meta.data,subclass_label %in% c("L5 PT CTX","L4/5 IT CTX","L5 NP CTX","L2/3 IT CTX",
                                                               "L6 IT CTX","L5/6 IT CTX","L5 IT CTX","L5 NP CT CTX ",
                                                               "L6 CT CTX","L6b CTX"))
ex.mon.acc <- subset(mon.acc,cells = rownames(mon.acc.meta))
rm(mon.acc)
DefaultAssay(ex.mon.acc) <- "homo"
ex.mon.acc <- NormalizeData(ex.mon.acc)
Idents(ex.mon.acc) <- "subclass_label"
ven.mouse.markers <-  FindMarkers(ex.mon.acc,ident.1 = "L5 PT CTX",verbose = F,
                                  min.pct = 0,logfc.threshold =0)
ven.mouse.markers$gene <- rownames(ven.mouse.markers)
ven.mouse.markers <- subset(ven.mouse.markers,avg_log2FC>1 & p_val_adj<0.05)

# get primate and rodent share marker genes(prSMGs)
RP.SMG <- intersect(ven.mouse.markers$gene,SMG)

# get human specific marker genes(hmSMGs)
Ex.inter <- readRDS("2023_newlibrary/202309_analysis/ACCmergeMouse.rds")
DefaultAssay(Ex.inter) <- "homo"
human_mouse <- FindMarkers(Ex.inter,ident.1 = "Human_EX-ET",ident.2 = "Mouse_L5 PT CTX",
                           logfc.threshold = 0,min.pct = 0)
human_mouse$gene <- rownames(human_mouse)
macaque_mouse <- FindMarkers(Ex.inter,ident.1 = "Macaque_EX-ET",ident.2 = "Mouse_L5 PT CTX",
                             logfc.threshold = 0,min.pct = 0)
macaque_mouse$gene <- rownames(macaque_mouse)
human_mouse_SMG <- subset(human_mouse, gene %in% SMG&avg_log2FC>2 & p_val_adj<0.01)
macaque_mouse_SMG <- subset(macaque_mouse, gene %in% SMG&avg_log2FC>2 & p_val_adj<0.01)
Primate <- intersect(human_mouse_SMG$gene,macaque_mouse_SMG$gene)
human_specific <- subset(human_mouse, gene%in%human.smg$gene & avg_log2FC>2 & p_val_adj<0.01)

# merge human&macaque data with mouse
# test mouse data homology
mon.acc<- readRDS("../mouse.aca.rds")
mon.acc.meta <- filter(mon.acc@meta.data,label %in% c("L5 PT CTX","L4/5 IT CTX","L5 NP CTX","L2/3 IT CTX",
                                                      "L6 IT CTX","L5/6 IT CTX","L5 IT CTX","L5 NP CT CTX ",
                                                      "L6 CT CTX","L6b CTX"))
ex.mon.acc <- subset(mon.acc,cells = rownames(mon.acc.meta))
rm(mon.acc)
DefaultAssay(ex.mon.acc) <- "homo"
DefaultAssay(Ex) <- "homo"
ex.mon.acc$species <- "Mouse"
ex.mon.acc$subtype <- ex.mon.acc$label
ex.mon.acc$sample <- ex.mon.acc$sample_name
Ex.list <- SplitObject(Ex,split.by = "sample")
Ex.list <- lapply(X = Ex.list, FUN = function(x) {
  x <- NormalizeData(x,verbose = F)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000,verbose = F)
})
ex.mon.acc <- NormalizeData(ex.mon.acc,verbose = F)
ex.mon.acc <- FindVariableFeatures(ex.mon.acc, selection.method = "vst", nfeatures = 2000,verbose = F)
features <- SelectIntegrationFeatures(object.list = c(Ex.list,ex.mon.acc))
anchors <- FindIntegrationAnchors(object.list = c(Ex.list,ex.mon.acc), dims = 1:30,
                                  anchor.features = features,reduction = "cca")
# integrated
Ex.inter <- IntegrateData(anchorset = anchors, dims = 1:30,k.weight = 20)
DefaultAssay(Ex.inter) <- "integrated"
# Run the standard workflow for visualization and clustering
Ex.inter <- ScaleData(Ex.inter, verbose = FALSE)
Ex.inter <- RunPCA(Ex.inter, npcs = 30, verbose = FALSE)
Ex.inter <- FindNeighbors(object = Ex.inter, dims = 1:30)
Ex.inter <- FindClusters(object = Ex.inter, resolution = 0.3)
head(x = Idents(object = Ex.inter), 5)
Ex.inter <- RunUMAP(Ex.inter, reduction = "pca", dims = 1:30)
DimPlot(Ex.inter, reduction = "umap", label = T, split.by = "species", group.by = "subtype", label.size = 3)
DimPlot(Ex.inter, reduction = "umap", label = T, group.by = "species", label.size = 3)
DimPlot(Ex.inter, reduction = "umap", label = T, group.by = "seurat_clusters", label.size = 3)
# renaeme subtype name
new.cluster.ids <- c("L2/3 IT","L2/3 IT","L6 CT","L2/3 IT","L6 IT","L5 IT","L6 CT","L5 IT","NP",
                     "L2/3 IT","L5 IT","L6b","L5 IT","L5 ET","L5 IT","L2/3 IT")
names(new.cluster.ids) <- levels(Ex.inter)
Ex.inter <- RenameIdents(Ex.inter, new.cluster.ids)
Ex.inter$cluster.inter <- Idents(Ex.inter)
DimPlot(Ex.inter, reduction = "umap", label = T, label.size = 5,pt.size = 0.5)  + theme(legend.position = 'none') 
Ex.inter$species.celltype <- paste(Ex.inter$species,Ex.inter$subtype,sep = "_")
Idents(Ex.inter) <- "species.celltype"
Cluster_Highlight_Plot(seurat_object = Ex.inter, cluster_name = "Human_EX-ET", highlight_color = "#546de5",
                       background_color = "lightgray") + theme(legend.position = 'none')
Cluster_Highlight_Plot(seurat_object = Ex.inter, cluster_name = "Macaque_EX-ET", highlight_color = "#ff4757",
                       background_color = "lightgray") + theme(legend.position = 'none')
Cluster_Highlight_Plot(seurat_object = Ex.inter, cluster_name = "Mouse_L5 PT CTX", highlight_color = "#87BE4B",
                       background_color = "lightgray") + theme(legend.position = 'none')
saveRDS(Ex.inter,"../ACCmergeMouse.rds")

# compare cell percent in mouse data
library(gplots)
library(RColorBrewer)
Ex.inter$orig.celltype <- paste(Ex.inter$orig.ident,Ex.inter$subtype,sep = "_")
compare = table(Ex.inter@meta.data$species.celltype,Ex.inter@meta.data$cluster.inter)
compare = round(100*compare/rowSums(compare))
roesidecol = c(rep("#546DE5",7),rep("#FF4757",7),rep("#87BE4B",10))
names(roesidecol) = rownames(compare)
heatmap.2(compare, dendrogram=c("none"), trace="none", margins = c(5, 15), colsep=1:100, sepwidth=c(0.05,0.1),
          rowsep=c(4,7,10,14,17,20,24), RowSideColors=roesidecol,col=colorRampPalette(brewer.pal(9, "Reds"))(100))

### FigS
### fig VEN
# VEN marker confirm https://www.nature.com/articles/s41467-020-14952-3
VEN.marker <- read.csv("20211125/ven/paper.marker.F.csv")
Idents(Ex) <- 'subtype'
levels(Ex)[1:7] <- c("EX L6b","EX IT L6","EX IT L5","EX IT L2/3","EX NP","EX CT","EX-ET")
DefaultAssay(Ex) <- "RNA"
VlnPlot(Ex,VEN.marker$gene,stack = T,split.by = "species",cols = c("#546de5","#ff4757"))
VlnPlot(Ex,VEN.marker.YLX,stack = T,split.by = "species",cols = c("#546de5","#ff4757"))
Stacked_VlnPlot(seurat_object = Ex, features = VEN.marker$gene, x_lab_rotate = TRUE,
                split.by = "species",colors_use = c("#546de5","#ff4757"))

# merge with FI data
library(feather)
library(gplots)
library(RColorBrewer)
# devtools::install_github("AllenInstitute/VENcelltypes")
# compare with ven-celltype from Allen
## Read in the data
inputFolder = "/sde/wuhaixv/wuhaixv/HM-RM_ACC_scRNA_YJM/20210707 VENcelltype_Allen/FI/"
Expr.dat <- feather(paste(inputFolder,"data.feather",sep=""))
annoFI   <- read_feather(paste(inputFolder,"anno.feather",sep="")) 
exprData <- as.matrix(Expr.dat[,colnames(Expr.dat)[colnames(Expr.dat)!="sample_id"]])
rownames(exprData) = Expr.dat$sample_id
datFI    <- t(exprData)
datFI    <- log2(datFI+1)
load(paste(inputFolder,"clusterInfo.rda",sep="")) 
infoFI   <- clusterInfo
## Only include excitatory nuclei
clusterType <- annoFI$cluster_type_label 
datFI       <- datFI[,clusterType=="exc"]
annoFI      <- annoFI[clusterType=="exc",]
# use Seurat to merge
brain.data     <- cbind(datFI)
brain.metadata <- data.frame(species=c(rep("humanFI",dim(datFI)[2])),
                             subtype = c(annoFI$cluster_label),
                             sample=rep("humanFI",dim(datFI)[2]))
rownames(brain.metadata) <- colnames(brain.data)
brain <- CreateSeuratObject(brain.data, meta.data = brain.metadata)
brain <- merge(brain,Ex)
brain.list <- SplitObject(object = brain, split.by = "sample")
for (i in 1:length(x = brain.list)) {
  brain.list[[i]] <- NormalizeData(object = brain.list[[i]], verbose = FALSE)
  brain.list[[i]] <- FindVariableFeatures(object = brain.list[[i]], 
                                          selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
brain.anchors    <- FindIntegrationAnchors(object.list = brain.list, dims = 1:30)
brain.integrated <- IntegrateData(anchorset = brain.anchors, dims = 1:30)
## Run the main integration analysis
DefaultAssay(object = brain.integrated) <- "integrated"
brain.integrated <- ScaleData(object = brain.integrated, verbose = FALSE)
brain.integrated <- RunPCA(object = brain.integrated, npcs = 30, verbose = FALSE)
brain.integrated <- RunUMAP(object = brain.integrated, reduction = "pca", dims = 1:30)
resolution       <- 0.2
FindNeighbors.Seurat <- VENcelltypes::FindNeighbors.Seurat  ## REQUIRED FOR CONSISTENT OUTPUT IN WINDOWS AND UNIX
brain.integrated <- FindNeighbors(object = brain.integrated, dims = 1:30)
brain.integrated <- FindClusters(object = brain.integrated, resolution = resolution)
DimPlot(brain.integrated,group.by = "subtype",split.by = "species")

Idents(brain.integrated) <- "species"
p1 <- DimPlot(object = subset(brain.integrated,idents = "Human"), group.by = "subtype",
              label.size = 3, label = T,pt.size = 0.5) +NoLegend()
p2<- DimPlot(object = subset(brain.integrated,idents = "Macaque"), group.by = "subtype",
             label.size = 3, label = T,pt.size = 0.5) +NoLegend()
p3 <- DimPlot(object = subset(brain.integrated,idents = "humanFI"), group.by = "subtype",
              label.size = 3, label = T,pt.size = 0.5)  +NoLegend()
wrap_plots(p1, p2,p3, ncol = 3)

saveRDS(brain.integrated,"2023_newlibrary/202309_analysis/ACCmergeFI.rds")


# test VEN marker to split ET neurons
# ET cluster by different feature
ven.yang <- read.csv("../VENmarkerFromLCM.csv")
Idents(Ex) <- "subtype"
EX.ET <- subset(Ex, idents = "EX-ET")
DefaultAssay(Ex) <- "homo"
EX.ET <- NormalizeData(EX.ET)
all.genes <- rownames(EX.ET)
EX.ET <- ScaleData(EX.ET, features = all.genes)
EX.ET <- RunPCA(EX.ET, features = ven.yang$gene_names)
EX.ET <- FindNeighbors(EX.ET, dims = 1:30)
EX.ET <- FindClusters(EX.ET, resolution = 0.5)
EX.ET <- RunUMAP(EX.ET, dims = 1:15)
DimPlot(EX.ET, reduction = "umap")
DimPlot(EX.ET, reduction = "umap",group.by = "species",cols = c("#546DE5","#FF4757"))
FeaturePlot(EX.ET, features = c("VAT1L","ADRA1A","CHST8","SULF2","BCL11B"),min.cutoff = 0,max.cutoff = 2,ncol = 3)