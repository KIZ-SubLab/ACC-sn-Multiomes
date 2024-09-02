library(ArchR)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

# Load ArchR object
adata = loadArchRProject('HumanBrain')

# Add Union Peaks
peak_df = read.delim('../Union_Peak_Info.txt', sep='\t')
peak_df = peak_df[, c('hg38_chr', 'hg38_start', 'hg38_end')]
colnames(peak_df) = c('chr', 'start', 'end')

adata = addPeakSet(ArchRProj = adata, peakSet = Peak_GRanges, force=TRUE)
adata <- addPeakMatrix(adata, force=TRUE)



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

write.table(x=Count_df, file="Peak_df/Union/Count_df.txt", sep='\t', quote = FALSE, row.names = FALSE)
write.table(x=Count_cells, file="Peak_df/Union/Count_Cells.txt", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x=Count_peaks, file="Peak_df/Union/Count_Peaks.txt", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


# Save ArchR Project
saveArchRProject(adata)