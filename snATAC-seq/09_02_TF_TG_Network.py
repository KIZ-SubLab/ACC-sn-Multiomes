import pandas as pd
import numpy as np
import os

celltype_list = ['ASC', 'EX', 'IN', 'ODC', 'OPC', 'MG']

Trans_df = pd.read_csv('../chromVAR/Motif_Trans.txt', sep='\t', header=None)
Trans_dict = dict(zip(Trans_df[0], Trans_df[1]))

# Peak Motif colocation Information (based on hg38)
Motif_Peak_df = pd.read_csv('Motif_Colocation/Peak_Motif_hg38.txt', sep='\t')
Motif_Peak_df['Peak'] = Motif_Peak_df['Peak'].map(lambda x: x.replace('-', '_'))
Motif_Peak_df['TF Name'] = Motif_Peak_df['TF'].map(Trans_dict)
Motif_Peak_df['TF Info'] = Motif_Peak_df['TF Name'].map(lambda x: x.split('(')[0])
Motif_Peak_df['TF Info'] = Motif_Peak_df['TF Info'].map(str.upper)
Motif_Peak_df.index = Motif_Peak_df['Peak']


# Peak2Gene Information (based on hg38)
Peak2Gene_df = pd.read_csv('../Human_ATAC/Human_Peak2Gene_Linkage.txt', sep='\t')
Peak2Gene_df.index = Peak2Gene_df['Peak']
Linked_Peak = np.unique(Peak2Gene_df['Peak'])


# Using Peak2Gene link Total PeakSet to Gene
input_dir = '../DAcCRE'
for temp_celltype in celltype_list:
    DiffPeak = pd.read_csv(os.path.join(input_dir, temp_celltype+'_Human.txt'), sep='\t', header=None)
    DiffPeak['Peak'] = DiffPeak.apply(lambda x: x[0]+'_'+str(x[1])+'_'+str(x[2]), axis=1)
    #DiffPeak['Peak'] = DiffPeak['Peak'].map(Peak_Trans_dict)
    
    Overlap_Peak = np.intersect1d(Linked_Peak, DiffPeak['Peak'])
    print(temp_celltype, DiffPeak.shape, Overlap_Peak.shape)
    temp_links = Peak2Gene_df.loc[Overlap_Peak,:]

    temp_Motif_Peak_df = Motif_Peak_df.loc[Overlap_Peak,]
    
    Link_dict = {}
    TF_dict = {}
    for it in Overlap_Peak:
        Link_dict[it] = np.unique(temp_links.loc[it,]['Gene'])
        TF_dict[it] = np.unique(temp_Motif_Peak_df.loc[it, 'TF Info'])

        
    Peak_Gene_TF = []
    for it in Overlap_Peak:
        for it1 in Link_dict[it]:
            for it2 in TF_dict[it]:
                Peak_Gene_TF.append(it+'**'+it1+'**'+it2)

                
    Peak_Gene_TF_df = pd.DataFrame(Peak_Gene_TF)
    Peak_Gene_TF_df.columns = ['Link']
    Peak_Gene_TF_df['Peak'] = Peak_Gene_TF_df['Link'].map(lambda x: x.split('**')[0])
    Peak_Gene_TF_df['Gene'] = Peak_Gene_TF_df['Link'].map(lambda x: x.split('**')[1])
    Peak_Gene_TF_df['TF'] = Peak_Gene_TF_df['Link'].map(lambda x: x.split('**')[2])
    
    output_file = os.path.join('Network/', temp_celltype+'.txt')
    Peak_Gene_TF_df.to_csv(output_file, sep='\t', index=None)