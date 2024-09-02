import pandas as pd
import numpy as np
import os

import pysam
fa_file = os.path.join('../hg38.fa')
fa = pysam.FastaFile(fa_file)

# Get HAR sequence (Reference Genome and SNC Altered Genome)
HAR_df = pd.read_csv('../GSE180714_HARs.bed', sep='\t')
HAR_df.index = HAR_df['HAR_ID']
HAR_df = HAR_df.drop_duplicates()

Overlap_HAR = pd.read_csv('HAR_SNC_Overlap/HAR_Overlap.txt', sep='\t')
Overlap_SNC = pd.read_csv('HAR_SNC_Overlap/SNC_Overlap.txt', sep='\t')


HAR_SNC = pd.concat([Overlap_SNC.loc[:, ['Name', 'Ref', 'Anc']], Overlap_HAR['HAR_ID']], axis=1)
HAR_SNC.shape
HAR_SNC.columns = ['SNC_Name', 'Ref', 'Anc', 'HAR_ID']
HAR_SNC['SNC_chr'] = HAR_SNC['SNC_Name'].map(lambda x: x.split('-')[0])
HAR_SNC['SNC_pos'] = HAR_SNC['SNC_Name'].map(lambda x: int(x.split('-')[1]))

for it in HAR_SNC.index:
    temp_HAR = HAR_SNC.loc[it, 'HAR_ID']
    HAR_SNC.loc[it, 'HAR_chr'] = HAR_df.loc[temp_HAR, 'seqnames']
    HAR_SNC.loc[it, 'HAR_start'] = int(HAR_df.loc[temp_HAR, 'start_500bp'])
    HAR_SNC.loc[it, 'HAR_end'] = int(HAR_df.loc[temp_HAR, 'end_500bp'])
HAR_SNC['HAR_start'] = HAR_SNC['HAR_start'].map(int)
HAR_SNC['HAR_end'] = HAR_SNC['HAR_end'].map(int)

HAR_SNC = HAR_SNC.loc[HAR_SNC['Anc']!='-', ]
HAR_SNC['Anc'] = HAR_SNC['Anc'].map(str.upper)


fa_Ref = pd.Series()
fa_Anc = pd.Series()
for HAR_ID in HAR_df.index:
    fa_Ref.loc[fa_Ref.shape[0],] = '>'+HAR_ID
    fa_Anc.loc[fa_Anc.shape[0],] = '>'+HAR_ID
    
    #HAR_ID = HAR_df.loc[it, 'HAR_ID']
    HAR_chr = HAR_df.loc[HAR_ID, 'seqnames']
    HAR_start = HAR_df.loc[HAR_ID, 'start_500bp']
    HAR_end = HAR_df.loc[HAR_ID, 'end_500bp']
    
    chrs = fa.fetch(HAR_chr, HAR_start, HAR_end)
    chrs = str.upper(chrs)
    
    fa_Ref.loc[fa_Ref.shape[0],] = str.upper(chrs)
    
    temp_SNC_df = HAR_SNC.loc[HAR_SNC['HAR_ID']==HAR_ID,]
    for SNC_it in temp_SNC_df.index:
        SNC_chr = temp_SNC_df.loc[SNC_it, 'SNC_chr']
        assert(SNC_chr==HAR_chr)
        SNC_pos = temp_SNC_df.loc[SNC_it, 'SNC_pos']
        
        temp_Ref = temp_SNC_df.loc[SNC_it, 'Ref']
        temp_Anc = temp_SNC_df.loc[SNC_it, 'Anc']
        
        if HAR_start>SNC_pos:
            continue
        if SNC_pos-HAR_start>=500:
            continue
        assert(chrs[SNC_pos-HAR_start]==temp_Ref)
        
        chrs = chrs[:SNC_pos-HAR_start]+temp_Anc+chrs[SNC_pos-HAR_start+1:]
        
        assert(chrs[SNC_pos-HAR_start]==temp_Anc)
    fa_Anc.loc[fa_Anc.shape[0],] = str.upper(chrs)
    
fa_Ref.to_csv('HAR_fa_file/HAR_Ref.fa', sep='\t', header=None, index=None)
fa_Anc.to_csv('HAR_fa_file/HAR_Anc.fa', sep='\t', header=None, index=None)