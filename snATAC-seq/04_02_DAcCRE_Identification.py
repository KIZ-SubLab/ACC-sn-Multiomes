import pandas as pd
import numpy as np
import os

celltype_list = ['EX', 'IN', 'ASC', 'ODC', 'OPC', 'MG']

Human_Sample = ['HM2013017', 'HM20200905', 'HM20200927', 'HM20201129',
                'HM20201213', 'HM20201222']
Monkey_Sample = ['RM06013', 'RM07071', 'RM07333', 'RM08403', 'RM1102009']

for temp_celltype in celltype_list:
    Used_df = pd.DataFrame()
    for temp_sample in Human_Sample:
        temp_df = pd.read_csv(temp_sample+'_'+temp_celltype+'.txt', sep='\t')
        Used_df[temp_sample+'_'+temp_celltype] = temp_df['x']

    Used_df.to_csv(temp_celltype+'_Human.txt', sep='\t')

for temp_celltype in celltype_list:
    Used_df = pd.DataFrame()
    for temp_sample in Monkey_Sample:
        temp_df = pd.read_csv(temp_sample+'_'+temp_celltype+'.txt', sep='\t')
        Used_df[temp_sample+'_'+temp_celltype] = temp_df['x']

    Used_df.to_csv(temp_celltype+'_Macaque.txt', sep='\t')


celltype_list = ['ODC', 'ASC', 'MG', 'EX', 'IN', 'OPC']
Human_df = pd.DataFrame()
for temp_celltype in celltype_list:
    temp_human = pd.read_csv(temp_celltype+'_Human.txt', sep='\t', index_col=0)
    temp_human = temp_human.loc[Used_peak, ]
    #temp_monkey = pd.read_csv('../Bulk_Celltype/'+temp_celltype+'_Monkey.txt', sep='\t', index_col=0)
    #temp_monkey = temp_monkey.loc[Used_peak, ]
    
    Human_df = pd.concat([Human_df, temp_human], axis=1)

Human_df.to_csv('Settings/Human_df.txt', sep='\t')

celltype_list = ['ODC', 'ASC', 'MG', 'EX', 'IN', 'OPC']
Human_df = pd.DataFrame()
for temp_celltype in celltype_list:
    temp_human = pd.read_csv(temp_celltype+'_Macaque.txt', sep='\t', index_col=0)
    temp_human = temp_human.loc[Used_peak, ]
    #temp_monkey = pd.read_csv('../Bulk_Celltype/'+temp_celltype+'_Monkey.txt', sep='\t', index_col=0)
    #temp_monkey = temp_monkey.loc[Used_peak, ]
    
    Human_df = pd.concat([Human_df, temp_human], axis=1)

Human_df.to_csv('Settings/Macaque_df.txt', sep='\t')


for temp_celltype in celltype_list:
    temp_human = pd.read_csv(temp_celltype+'_Human.txt', sep='\t', index_col=0)
    temp_human = temp_human.loc[Used_peak, ]
    temp_monkey = pd.read_csv(temp_celltype+'_Monkey.txt', sep='\t', index_col=0)
    temp_monkey = temp_monkey.loc[Used_peak, ]
    
    temp_df = pd.concat([temp_human, temp_monkey], axis=1)

    temp_df.to_csv('Settings/'+temp_celltype+'_Cross.txt', sep='\t')