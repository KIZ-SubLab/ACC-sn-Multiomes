library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(JASPAR2020)
addArchRThreads(threads = 8) 
addArchRGenome("hg38")


# Load Data
adata = loadArchRProject('HumanBrain')


# Add JASPAR2020 motifs
PWM <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE, matrixtype="PWM")
)

adata <- addMotifAnnotations(ArchRProj = adata, motifPWMs = PWM,
                                 name = "Motif")


#######################################################
# Run Python Code
import pandas as pd
import numpy as np
import os
celltype_list = ['EX', 'IN', 'ASC', 'ODC', 'OPC', 'MG']
Union_df = pd.read_csv('../Union_Peak_Info.txt', sep='\t')

DiffPeak = pd.DataFrame(1, index=Union_df['hg38_Peak'], columns=celltype_list)
for temp_celltype in celltype_list:
    temp_diff = pd.read_csv('../DAcCRE/'+temp_celltype+'_Human.txt',
                            sep='\t', header=None)
    temp_diff.columns = ['chr', 'start', 'end']
    temp_diff['Peak'] = temp_diff.apply(lambda x: x['chr']+'_'+str(x['start'])+'_'+str(x['end']), axis=1)

    DiffPeak.loc[list(temp_diff['Peak']), temp_celltype] = 0
DiffPeak.index = [x.split('_')[0]+':'+x.split('_')[1]+'-'+x.split('_')[2] for x in DiffPeak.index]
DiffPeak.to_csv('DiffPeak_Human.txt', sep='\t')
#######################################################


# Add DAcCRE Information
DiffPeak_df = read.delim('DiffPeak_Human.txt', sep='\t')#, row.names=1)

Peak_Info = GRanges(DiffPeak_df$X)
Peak_Info = as.data.frame(Peak_Info)[, c('seqnames', 'start', 'end')]

Used_DE = DiffPeak_df[, c('EX', 'IN', 'ASC', 'ODC', 'OPC', 'MG')]
markersPeaks <- SummarizedExperiment(assays=list(FDR=Used_DE), rowData= Peak_Info)

markersPeaks@metadata = list()
markersPeaks@metadata$Params = list()
markersPeaks@metadata$Params$useMatrix = 'PeakMatrix'

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05")


enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = adata,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05"
  )

write.table(enrichMotifs@assays@data$mlog10Padj, 'DiffPeak_TF/mlog10Padj_BigCellType.txt', sep='\t', quote=FALSE)
write.table(enrichMotifs@assays@data$Enrichment, 'DiffPeak_TF/Enrichment_BigCellType.txt', sep='\t', quote=FALSE)


# Save ArchR Project
saveArchRProject(adata)