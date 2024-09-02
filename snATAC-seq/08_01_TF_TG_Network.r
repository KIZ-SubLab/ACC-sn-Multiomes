library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(JASPAR2020)
library(readr)
library(stats)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(Signac)
library(Seurat)

# Running Human chromVAR (Get Peak Motif Colocation Information)
input_dir = "../Human_ATAC/Peak_df/Union/"
temp_df <- read_delim(paste0(input_dir, 'Count_df.txt'), "\t", escape_double = FALSE, trim_ws = TRUE)
peaknames <- read_delim(paste0(input_dir, 'Count_Peaks_Signac.txt'), "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
cellnames <- read_delim(paste0(input_dir, 'Count_Cells_Signac.txt'), "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)

datafr = sparseMatrix(i=temp_df$Peaks, j=temp_df$Cells, x=1)
colnames(datafr) = cellnames$X1
rownames(datafr) = peaknames$X1

brain_assay <- CreateChromatinAssay(
  counts = datafr,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = '../Data/Human/Human_fragments.tsv.gz',
  min.cells = 1
)

adata <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'peaks',
  project = 'ATAC')#,
#  meta.data = metadata
#)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
# add motif information
adata <- AddMotifs(
  object = adata,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
adata <- RunChromVAR(
  object = adata,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Get Peak Motif colocation Information
Motif_Info = Motifs(adata)
Motif_Name = data.frame(Motif_Info@motif.names)
Motif_Name = t(Motif_Name)
Motif_Name = data.frame(Motif_Name)

write.table(Motif_Name, file='Motif_Trans.txt', sep='\t', col.names =FALSE, quote=FALSE)
saveRDS(Motif_Info@data, 'Peak_Motif_colocation.rds')


# Save to .txt
Peak_hg38_file = '../Human_ATAC/Peak_Motif_colocation.rds'
Peak_rheMac10_file = '../Monkey_ATAC/Peak_Motif_colocation.rds'

write.table(as.matrix(Peak_hg38), 'Motif_Colocation/Peak_Motif_hg38.txt', sep='\t', quote=FALSE)
write.table(as.matrix(Peak_rheMac10), 'Motif_Colocation/Peak_Motif_rheMac10.txt', sep='\t', quote=FALSE)