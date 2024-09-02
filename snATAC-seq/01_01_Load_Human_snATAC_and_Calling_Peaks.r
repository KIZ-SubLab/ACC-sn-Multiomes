library(ArchR)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

# Create ArchR Arrow Files
pathFragments = '../../Data/Human'
inputFiles =  c('../../Data/Human/Human_fragments.tsv.gz')
names(inputFiles) = "Human_ATAC"
ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = names(inputFiles),
    minTSS = 3, #Dont set this too high because you can always increase later
    minFrags = 5000, 
    maxFrags=1e7,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
)


# Preprocessing (Filtering Low-quality Cells)
ArrowFiles = 'Human_ATAC.arrow'
adata <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HumanBrain",
  copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
meta = read.delim('../Annotation/ATAC_metadata_ArchR.txt', sep='\t', row.names=1)
Overlap_cell = intersect(rownames(meta), rownames(adata@cellColData))
adata = adata[Overlap_cell,]

# Plotting TSS Enrichment score and #Fragments
df <- getCellColData(adata, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enricdwqhment",
    xlim = c(3, quantile(df[,1], probs = 0.999)),
    ylim = c(2, quantile(df[,2], probs = 0.999))
) + geom_hline(yintercept = 3, lty = "dashed") + geom_vline(xintercept = log10(5000), lty = "dashed")
p


# Add Annotation (Obtained from integration of scRNA-seq data)
adata@cellColData$Annotation = meta[rownames(adata@cellColData), 'Annotation']


# Calling Peaks
pathToMacs2 <- findMacs2()
adata <- addGroupCoverages(ArchRProj = adata, maxCells =1000, minCells=200, 
                               groupBy ="Annotation")
adata <- addReproduciblePeakSet(
    ArchRProj = adata, cutOff=0.05,
    groupBy = "Annotation", #maxPeaks=250000,
    pathToMacs2 = pathToMacs2
)
adata <- addPeakMatrix(adata)


# Save Peak Matrix
Filter_by_proportion <- function(x, proportion){
    Binary_df = (x>0)
    peak_sum = rowSums(Binary_df)
    x[(peak_sum>proportion*dim(x)[2]),]
}
Peak_mtx = getMatrixFromProject(adata, useMatrix = 'PeakMatrix')
temp_df = Peak_mtx@assays@data$PeakMatrix
peak_meta = as.data.frame(Peak_mtx@rowRanges)
peaks = tidyr::unite(peak_meta, "peaks", seqnames, start, end)
rownames(temp_df) = peaks$peaks

xsp_Count = Filter_by_proportion(temp_df ,proportion = 0)
dim(xsp_Count)
Count_df = summary(xsp_Count)
colnames(Count_df) = c('Peaks', 'Cells', 'Counts')
Count_cells = colnames(xsp_Count)
Count_peaks = rownames(xsp_Count)

write.table(x=Count_df, file="Peak_df/ArchR/Count_df.txt", sep='\t', quote = FALSE, row.names = FALSE)
write.table(x=Count_cells, file="Peak_df/ArchR/Count_Cells.txt", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x=Count_peaks, file="Peak_df/ArchR/Count_Peaks.txt", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


# Save ArchR Project
saveArchRProject(adata)