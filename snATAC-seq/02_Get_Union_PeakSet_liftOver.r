# We have published the information of the union peak set on https://ngdc.cncb.ac.cn/ (BioProject ID is PRJCA015229), and it can be directly accessed.

# Step1: liftOver Macaque Peaks (rheMac10) into hg38
system("liftOver ../Macaque_ATAC/Peak_df/ArchR/Count_Peaks.bed rheMac10ToHg38.over.chain.gz Macaque_Peaks_hg38.bed Macaque_Peaks_hg38_unmapped.bed")


# Step2: Merge Macaque Peaks (hg38) and Human Peaks (hg38)
library(GenomicRanges)
Macaque_Peaks_hg38 = read.delim('../Macaque_Peaks_hg38.bed', sep='\t', header=FALSE)
colnames(Macaque_Peaks_hg38) = c('chr', 'start', 'end')
Macaque_Peaks_hg38 = GRanges(Macaque_Peaks_hg38)

Human_Peaks_hg38 = read.delim('../Human_ATAC/Peak_df/ArchR/Count_Peaks.bed', sep='\t', header=FALSE)
colnames(Human_Peaks_hg38) = c('chr', 'start', 'end')
Human_Peaks_hg38 = GRanges(Human_Peaks_hg38)

Merged_Peaks_hg38 = educe(c(Human_peak, Monkey_peak))
Merged_Peaks_hg38 =  as.data.frame(resMerged_Peaks_hg38)[, c('seqnames', 'start', 'end')]

write.table(Merged_Peaks_hg38, 'Merged_Peaks_hg38.bed', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)


# Step3: liftOver Merged Peaks (hg38) into rheMac10
system("liftOver Merged_Peaks_hg38.bed hg38ToRheMac10.over.chain.gz Merged_Peaks_rheMac10.bed Merged_Peaks_rheMac10_unmapped.bed")


# Step4: Get liftOver relationship of Merged Peaks (Python Code)
import pandas as pd
import numpy as np
import os

peak_hg38 = pd.read_csv('Merged_Peaks_hg38.bed', sep='\t', header=None)
peak_hg38.columns = ['chr', 'start', 'end']
peak_hg38['Peak'] = peak_hg38.apply(lambda x: x['chr']+'_'+str(x['start'])+'_'+str(x['end']), axis=1)
peak_hg38.index = peak_hg38['Peak']

peak_rheMac10 = pd.read_csv('Merged_Peaks_rheMac10.bed', sep='\t', header=None)
peak_rheMac10.columns = ['chr', 'start', 'end']
peak_rheMac10['Peak'] = peak_rheMac10.apply(lambda x: x['chr']+'_'+str(x['start'])+'_'+str(x['end']), axis=1)
peak_rheMac10.index = peak_rheMac10['Peak']

unmap_df = pd.read_csv('Merged_Peaks_rheMac10_unmapped.bed',header=None)
unmap_df = unmap_df[1::2]
unmap_df[0] = unmap_df[0].map(lambda x: x.replace('\t', '_'))


peak_hg38.loc[list(unmap_df[0]), 'Used'] = 'F'
peak_hg38 = peak_hg38.loc[peak_hg38['Used']!='F',]
peak_hg38.loc[:, ['Peak', 'rheMac10_Peak']].to_csv('Merged_Peaks_liftOver.txt', sep='\t', index=None)


#Step5: Filtering Peaks Longer than 2000bp (Python Code)
import pandas as pd
import numpy as np
import os

Used_chr = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21',
            'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX']

temp_df = pd.read_csv('Merged_Peaks_liftOver.txt', sep='\t')
temp_df.columns = ['hg38_Peak', 'rheMac10_Peak']

temp_df['hg38_chr'] = temp_df['hg38_Peak'].map(lambda x: x.split('_')[0])
temp_df['hg38_start'] = temp_df['hg38_Peak'].map(lambda x: int(x.split('_')[1]))
temp_df['hg38_end'] = temp_df['hg38_Peak'].map(lambda x: int(x.split('_')[2]))

temp_df['rheMac10_chr'] = temp_df['rheMac10_Peak'].map(lambda x: x.split('_')[0])
temp_df['rheMac10_start'] = temp_df['rheMac10_Peak'].map(lambda x: int(x.split('_')[1]))
temp_df['rheMac10_end'] = temp_df['rheMac10_Peak'].map(lambda x: int(x.split('_')[2]))

temp_df['hg38_length'] = temp_df['hg38_end']-temp_df['hg38_start']
temp_df['rheMac_length'] = temp_df['rheMac10_end']-temp_df['rheMac10_start']

temp_df = temp_df.loc[temp_df['hg38_length']<2000,]
temp_df = temp_df.loc[temp_df['rheMac_length']<2000,]

Used_col = ['hg38_Peak',  'hg38_chr', 'hg38_start', 'hg38_end',
            'rheMac10_Peak', 'rheMac10_chr', 'rheMac10_start', 'rheMac10_end']
temp_df = temp_df.loc[:,Used_col]

temp_df.to_csv('../Union_Peak_Info.txt', sep='\t', index=None)
