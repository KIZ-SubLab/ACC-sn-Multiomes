import pandas as pd
import numpy as np
import os

celltype_list = ['EX', 'IN', 'ASC', 'ODC', 'OPC', 'MG']

# Get DA-cCRE Set (Human)
for temp_celltype in celltype_list:
    Human_df = pd.read_csv('Res/'+temp_celltype+'_Human.txt', sep='\t')
    Human_df = Human_df.loc[Human_df['padj']<0.05,]
    Human_df = Human_df.loc[Human_df['logFC']<0,]

    Cross_df = pd.read_csv('Res/'+temp_celltype+'_CrossMvsH.txt', sep='\t')
    Cross_df = Cross_df.loc[Cross_df['padj']<0.05,]
    Cross_df = Cross_df.loc[Cross_df['logFC']<0,]

    print(temp_celltype, Human_df.shape, Cross_df.shape)
    Used_peak = np.intersect1d(Human_df.index, Cross_df.index)
    Used_df = Human_df.loc[Used_peak, ]
    Used_df = Used_df.sort_values(by='padj')

    print(Used_df.shape)
    Used_df['chr'] = [x.split('_')[0] for x in Used_df.index]
    Used_df['start'] = [int(x.split('_')[1]) for x in Used_df.index]
    Used_df['end'] = [int(x.split('_')[2]) for x in Used_df.index]

    Used_df.loc[:, ['chr','start','end']].to_csv('DAcCRE/'+temp_celltype+'_Human.txt', sep='\t', index=None, header=None)


# Get DA-cCRE Set (Macaque)
for temp_celltype in celltype_list:
    Human_df = pd.read_csv('Res/'+temp_celltype+'_Macaque.txt', sep='\t')
    Human_df = Human_df.loc[Human_df['padj']<0.05,]
    Human_df = Human_df.loc[Human_df['logFC']<0,]

    Cross_df = pd.read_csv('Res/'+temp_celltype+'_CrossMvsH.txt', sep='\t')
    Cross_df = Cross_df.loc[Cross_df['padj']<0.05,]
    Cross_df = Cross_df.loc[Cross_df['logFC']>0,]

    print(temp_celltype, Human_df.shape, Cross_df.shape)
    Used_peak = np.intersect1d(Human_df.index, Cross_df.index)
    Used_df = Human_df.loc[Used_peak, ]
    Used_df = Used_df.sort_values(by='padj')

    print(Used_df.shape)
    Used_df['chr'] = [x.split('_')[0] for x in Used_df.index]
    Used_df['start'] = [int(x.split('_')[1]) for x in Used_df.index]
    Used_df['end'] = [int(x.split('_')[2]) for x in Used_df.index]

    Used_df.loc[:, ['chr','start','end']].to_csv('DAcCRE/'+temp_celltype+'_Macaque.txt', sep='\t', index=None, header=None)