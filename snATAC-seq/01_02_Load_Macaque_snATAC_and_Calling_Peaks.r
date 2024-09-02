# If you need the original fragment files for quality control analysis, please contact: wuhaixu@mail.kiz.ac.cn

library(ArchR)
addArchRThreads(threads = 8)
rheMac10_UCSC = createGenomeAnnotation(genome='BSgenome.Mmulatta.UCSC.rheMac10')


#Load rheMac10 Gene Annotation Information
mm10_annotation = read.delim(file = '../../Data/modify_mm10_used-UCSC.txt', sep='\t', header = TRUE)
gene_annotation = mm10_annotation[mm10_annotation['type']=='gene',]
exon_annotation = mm10_annotation[mm10_annotation['type']=='exon',]
mm10_annotation = mm10_annotation[mm10_annotation['type']!='gene',]

pos <- subset(mm10_annotation, strand == "+")
pos <- pos[order(pos$start),] 
pos <- pos[!duplicated(pos$transcript_id),] # remove all but the first exons per transcript
pos$end <- pos$start #+ 1 # make a 1 base pair marker of the TSS

neg <- subset(mm10_annotation, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
neg <- neg[!duplicated(neg$transcript_id),] # remove all but the first exons per transcript
neg$start <- neg$end #- 1

TSS_annotation <- rbind(pos, neg)
TSS_annotation = TSS_annotation[TSS_annotation$type=='transcript',]
TSS_annotation = GRanges(tx_id = TSS_annotation$transcript_id, tx_name=TSS_annotation$gene_name, TSS_annotation)
exon_annotation = GRanges(gene_id=exon_annotation$gene_id, symbol = exon_annotation$gene_name, exon_annotation)
gene_annotation = GRanges(gene_id=gene_annotation$gene_id, symbol=gene_annotation$gene_name, gene_annotation)
rheMac10_gene_annotation = createGeneAnnotation(TSS = TSS_annotation, exons = exon_annotation, genes = gene_annotation)


# Create ArchR Arrow Files
pathFragments = '../../Data/Macaque'
inputFiles =  c('../../Data/Macaque/Macaque_fragments.tsv.gz')
names(inputFiles) = "Macaque_ATAC"
ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    genomeAnnotation = rheMac10_UCSC,
    geneAnnotation=rheMac10_gene_annotation,
    sampleNames = names(inputFiles),
    minTSS = 3, #Dont set this too high because you can always increase later
    minFrags = 5000, maxFrags=1e7,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
)


# Preprocessing (Filtering Low-quality Cells)
ArrowFiles = 'Macaque_ATAC.arrow'
adata <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    genomeAnnotation = rheMac10_UCSC,
    geneAnnotation=rheMac10_gene_annotation,
    outputDirectory = "MacaqueBrain",
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
    ylabel = "TSS Enrichment",
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
    groupBy = "Annotation", 
    pathToMacs2 = pathToMacs2, genomeSize=2.2e9
)
adata <- addPeakMatrix(adata)


# Save Peak Matrix
Filter_by_proportion <- function(x, proportion){
    Binary_df = (x>0)
    peak_sum = rowSums(Binary_df)
    x[(peak_sum>proportion*dim(x)[2]),]
}
Peak_mtx = getMatrixFromProject(projHeme5, useMatrix = 'PeakMatrix')
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
