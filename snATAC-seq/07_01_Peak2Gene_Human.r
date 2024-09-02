library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRThreads(threads = 8) 
addArchRGenome("hg38")


# Load Data
adata = loadArchRProject('HumanBrain')


# Add scRNA Human
Human_scRNA = readRDS('../RDS_Seurat/scRNA_Human.rds')


adata <- addIterativeLSI(
    ArchRProj = adata,
    useMatrix = "PeakMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

adata <- addUMAP(
    ArchRProj = adata, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

#~5 minutes
adata <- addGeneIntegrationMatrix(
    ArchRProj = adata, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    groupATAC='Annotation',
    seRNA = Human_scRNA,
    addToArrow = TRUE,
    force= TRUE,
    #groupList = groupList,
    groupRNA = "BigCellType",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)


# Adjust corresponding cell information of 10x multi-ome
adata@cellColData[adata@cellColData$Technoloty=='MultiOme', 'predictedCell'] = adata@cellColData[adata@cellColData$Technoloty=='MultiOme', 'Barcode']


# Run Peak2Gene of ArchR package
adata <- addPeak2GeneLinks(
    ArchRProj = adata,
    reducedDims = "IterativeLSI"
)


# Save
saveArchRProject(adata)


# Save Peak2Gene result
p2g <- getPeak2GeneLinks(
    ArchRProj = adata,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)

write.table(p2g, 'Human_Peak2Gene.txt', sep='\t', quote=FALSE)
write.table(p2g@metadata$geneSet, 'Human_Peak2Gene_GeneSet.txt', sep='\t', quote=FALSE)
write.table(p2g@metadata$peakSet, 'Human_Peak2Gene_PeakSet.txt', sep='\t', quote=FALSE)