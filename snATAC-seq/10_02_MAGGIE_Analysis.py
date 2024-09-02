import pandas as pd
import numpy as np
import os

import pysam

# Step2: Get hg38 Reference and Altered SNC fa files of DA-cCREs
fa = pysam.FastaFile("../hg38.fa")

SNCs_df = pd.read_csv('../SNCs_hg38_995-matched.txt', sep='\t')
SNCs_df['Ref'] = SNCs_df['Ref'].map(str.upper)
SNCs_df['Anc'] = SNCs_df['Anc'].map(str.upper)
SNCs_df = SNCs_df.loc[SNCs_df['Ref'].isin(['A', 'T', 'C', 'G']),]
SNCs_df = SNCs_df.loc[SNCs_df['Anc'].isin(['A', 'T', 'C', 'G']),]

SNCs_df.index = SNCs_df['hg38_name']

celltype_list = ['EX', 'IN', 'ASC', 'ODC', 'OPC', 'MG']
file_list = [x+'_Human.bed' for x in celltype_list]

for temp_file in file_list:
    Overlap_SNCs = pd.read_csv('../Overlap_SNCs/'+temp_file, sep='\t')
    Overlap_SNCs['name'] = Overlap_SNCs.apply(lambda x:x['seqnames']+'-'+str(x['start'])+'-'+str(x['end']), axis=1)
    Diff_SNCs_index = np.intersect1d(Overlap_SNCs['name'],SNCs_df.index)
    print(temp_file, Diff_SNCs_index.shape)
    print(Overlap_SNCs.shape)
    Diff_SNCs = SNCs_df.loc[Diff_SNCs_index, ['hg38_chr', 'hg38_start', 'Ref', 'Anc']]
    Diff_SNCs.columns = ['chr', 'start', 'Ref', 'Anc']
    non_SNCs_fa = pd.Series()
    SNCs_fa = pd.Series()
    for it in Diff_SNCs.index:
        temp_str = Diff_SNCs.loc[it, 'chr']
        temp_pos = Diff_SNCs.loc[it, 'start']+1
        temp_REF = str.upper(Diff_SNCs.loc[it, 'Ref'])
        temp_ALT = str.upper(Diff_SNCs.loc[it, 'Anc'])
        non_SNCs_fa.loc[non_SNCs_fa.shape[0],] = '>'+it
        SNCs_fa.loc[SNCs_fa.shape[0],] = '>'+it
        #non_SNCs_fa.loc[non_SNCs_fa.shape[0],]
        chrs = fa.fetch(temp_str, temp_pos-50, temp_pos+50)
        assert(str.upper(chrs[49])==str.upper(temp_REF))
        non_SNCs_fa.loc[non_SNCs_fa.shape[0],] = str.upper(chrs)
        chrs_SNC = chrs
        chrs_SNC = chrs_SNC[:49]+temp_ALT+chrs_SNC[50:]
        SNCs_fa.loc[SNCs_fa.shape[0],] = str.upper(chrs_SNC)
        #if str.upper(chrs)!=str.upper(temp_REF):
        #    zs+=1
    non_SNCs_fa.to_csv('fa_seq/'+temp_file.split('.')[0]+'_non_SNCs.fa', sep='\t', header=None, index=None)
    SNCs_fa.to_csv('fa_seq/'+temp_file.split('.')[0]+'_SNCs.fa', sep='\t', header=None, index=None)

