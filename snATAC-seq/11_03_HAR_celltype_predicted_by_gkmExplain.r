library(GenomicRanges)

# Get HAR SNC Overlap Information
SNC_df = read.delim('../../SNC_MAGGIE/SNCs_hg38_995-matched.txt', sep='\t')
SNC_df = SNC_df[, c('hg38_chr', 'hg38_start', 'hg38_end', 'hg38_name', 'hg38_Ref', 'Anc')]
colnames(SNC_df) = c('chr', 'start', 'end', 'Name', 'Ref', 'Anc')
SNCs = GRanges(SNC_df)
length(SNCs)

HAR_df = read.delim('../GSE180714_HARs.bed', sep='\t')
head(HAR_df)
HARs = GRanges(HAR_df)

res = findOverlaps(HARs, SNCs)

Overlap_HAR = HARs[res@from,]
Overlap_SNC = SNCs[res@to,]

write.table(Overlap_HAR, 'HAR_SNC_Overlap/HAR_Overlap.txt', sep='\t', row.names=FALSE, quote=FALSE)
write.table(Overlap_SNC, 'HAR_SNC_Overlap/SNC_Overlap.txt', sep='\t',  row.names=FALSE, quote=FALSE)