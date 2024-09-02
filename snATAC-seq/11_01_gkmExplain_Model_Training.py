import pandas as pd
import numpy as np
import os

import pysam

chr_list = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
            'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21',
            'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX']

chrom_size = pd.read_csv('../hg38.chrom.sizes', sep='\t', header=None, index_col=0)

# Step1: Calculating Backgound GC rate (500bp window)
fa_file = os.path.join('../hg38.fa')
fa = pysam.FastaFile(fa_file)

for temp_chr in chr_list:
    max_length = chrom_size.loc[temp_chr, 1]
    
    Peak_list = []
    Num_list = [] 
    ATCG_dict = {}
    for temp_bp in ['A', 'T', 'C', 'G']:
        ATCG_dict[temp_bp] = []
    
    #chr_df = pd.DataFrame(index=range(int(max_length/200)))
    print(temp_chr)
    
    it = 501
    while it < max_length-1000:
        temp_peaks = temp_chr + '_' + str(it-500) + '_' + str(it+500)
        
        chrs = fa.fetch(temp_chr, it-500, it+500)
        chrs = chrs.upper()
        assert(len(chrs)==1000)
        
        it += 200
        
        temp_ATCG = pd.value_counts(list(chrs))
        if 'N' in temp_ATCG.index:
            continue
        if temp_ATCG.shape[0]<4:
            continue
        
        Peak_list.append(temp_peaks)
        for temp_bp in ['A', 'T', 'C', 'G']:
            ATCG_dict[temp_bp].append(int(temp_ATCG[temp_bp]))
        
        Num_list.append(int(temp_ATCG.shape[0]))
         
        
    chr_dict = {}
    chr_dict['Peak'] = Peak_list
    for temp_bp in ['A', 'T', 'C', 'G']:
        chr_dict[temp_bp] = ATCG_dict[temp_bp]
    chr_dict['Num'] = Num_list

    chr_df = pd.DataFrame(chr_dict)
    chr_df.to_csv('Back-GC/'+temp_chr+'.txt', sep='\t')


# Step2: Calculating GC rate of DA-cCREs
fa = pysam.FastaFile(fa_file)
celltype_list = ['ASC', 'EX', 'IN', 'ODC', 'OPC', 'MG']
file_list = [x+'_Human.txt' for x in celltype_list]

input_dir = '../DAcCRE'
output_dir = os.path.join('Diff_peaks-GC')


for temp_file in file_list:
    temp_df = pd.read_csv(os.path.join(input_dir, temp_file), sep='\t', header=None)
    temp_df.columns = ['chr', 'start', 'end']
    
    temp_df['mid'] = (temp_df['start']+temp_df['end'])/2
    temp_df['mid'] = temp_df['mid'].map(int)
    
    temp_df['start_1kb'] = temp_df['mid']-500
    temp_df['end_1kb'] = temp_df['mid']+500
    for it in temp_df.index:
        temp_peaks = temp_df.loc[it,]
        chrs = fa.fetch(temp_peaks['chr'], temp_peaks['start_1kb'], temp_peaks['end_1kb'])
        chrs = chrs.upper()
        assert(len(chrs)==1000)
                
        temp_ATCG = pd.value_counts(list(chrs))
        temp_df.loc[it, 'T'] = int(temp_ATCG['T'])
        temp_df.loc[it, 'A'] = int(temp_ATCG['A'])
        temp_df.loc[it, 'C'] = int(temp_ATCG['C'])
        temp_df.loc[it, 'G'] = int(temp_ATCG['G'])
        
        temp_df.loc[it, 'Num'] = int(temp_ATCG.shape[0])
    temp_df.to_csv(os.path.join(output_dir, temp_file), sep='\t', index=None)
        #assert(np.unique(list(chrs)).shape[0]==4)


for temp_file in file_list:
    temp_df = pd.read_csv(os.path.join(output_dir, temp_file), sep='\t')
    print(temp_file, '--', temp_df.shape)
    temp_df = temp_df.loc[temp_df['Num']<5,]
    print(temp_df.shape)
    temp_df['GC_rate'] = (temp_df['G']+temp_df['C'])/1000*100
    
    temp_df.to_csv(os.path.join(output_dir, temp_file), sep='\t', index=None)


# Step3: Get Backgound genomic regions with similar GC rate of DA-cCREs
Total_Back = pd.DataFrame()
for temp_chr in chr_list:
    temp_Back = pd.read_csv('Back-GC/'+temp_chr+'.txt', sep='\t', index_col=0)
    Total_Back = pd.concat([Total_Back, temp_Back])

Total_Back['GC_rate'] = (Total_Back['G']+Total_Back['C'])/10

celltype_list = ['EX', 'IN', 'ASC', 'ODC', 'OPC', 'MG']
file_list = [x+'_Human.txt' for x in celltype_list]

for temp_file in file_list:
    DiffPeak = pd.read_csv('Diff_peaks-GC/'+temp_file, sep='\t')
    DiffPeak['GC_group'] = 'R0'
    Total_Back['GC_group'] = 'R0'
    for it in range(1, 20):
        cutoff = np.percentile(DiffPeak['GC_rate'], it*5)
        DiffPeak.loc[DiffPeak['GC_rate']>cutoff, 'GC_group'] = 'R'+str(it)
        Total_Back.loc[Total_Back['GC_rate']>cutoff, 'GC_group'] = 'R'+str(it)
    for it in range(20):
        Sample_num = pd.value_counts(DiffPeak['GC_group'])['R'+str(it)]
        Sample_list = np.random.choice(list(Total_Back.loc[Total_Back['GC_group']=='R'+str(it), 'Peak']), Sample_num, replace=False)
        DiffPeak.loc[DiffPeak['GC_group']=='R'+str(it), 'Neg_Peak'] = Sample_list
    
    DiffPeak.to_csv('Diff_peaks_Neg/'+temp_file, sep='\t')



# Step4: Split Training set and test set for the 10-fold Cross-validation
chr_split = {}
chr_split['Fold_0'] = ['chr1']
chr_split['Fold_1'] = ['chr2', 'chr19']
chr_split['Fold_2'] = ['chr3', 'chr20']
chr_split['Fold_3'] = ['chr6', 'chr13', 'chr22']
chr_split['Fold_4'] = ['chr5', 'chr16']
chr_split['Fold_5'] = ['chr4', 'chr15', 'chr21']
chr_split['Fold_6'] = ['chr7', 'chr14', 'chr18']
chr_split['Fold_7'] = ['chr11', 'chr17', 'chrX']
chr_split['Fold_8'] = ['chr9', 'chr12']
chr_split['Fold_9'] = ['chr8', 'chr10']

celltype_list = ['EX', 'IN', 'ASC', 'ODC', 'OPC', 'MG']
file_list = [x+'_Human.txt' for x in celltype_list]

for temp_file in file_list:
    DiffPeak = pd.read_csv('Diff_peaks_Neg/'+temp_file, sep='\t', index_col=0)
    for it in range(10):
        Train_Set = DiffPeak.loc[(DiffPeak['chr'].isin(chr_split['Fold_'+str(it)]))==False,]
        Test_Set = DiffPeak.loc[DiffPeak['chr'].isin(chr_split['Fold_'+str(it)]),]
        
        # Train
        Train_Set_Pos = Train_Set.loc[:, ['chr', 'start_1kb', 'end_1kb']]
        Train_Set_Pos.to_csv('Final_Peaks/'+temp_file.split('.txt')[0]+'_Fold'+str(it)+'_Pos_Train.txt', sep='\t', header=None, index=None)
        Train_Set_Neg = Train_Set.loc[:, ['Neg_Peak']]
        Train_Set_Neg['chr'] = Train_Set_Neg['Neg_Peak'].map(lambda x: x.split('_')[0])
        Train_Set_Neg['start'] = Train_Set_Neg['Neg_Peak'].map(lambda x: int(x.split('_')[1]))
        Train_Set_Neg['end'] = Train_Set_Neg['Neg_Peak'].map(lambda x: int(x.split('_')[2]))
        Train_Set_Neg.loc[:, ['chr', 'start', 'end']].to_csv('Final_Peaks/'+temp_file.split('.txt')[0]+'_Fold'+str(it)+'_Neg_Train.txt', 
                                                             sep='\t', header=None, index=None)
        
        # Test
        Train_Set_Pos = Test_Set.loc[:, ['chr', 'start_1kb', 'end_1kb']]
        Train_Set_Pos.to_csv('Final_Peaks/'+temp_file.split('.txt')[0]+'_Fold'+str(it)+'_Pos_Test.txt', sep='\t', header=None, index=None)
        Train_Set_Neg = Test_Set.loc[:, ['Neg_Peak']]
        Train_Set_Neg['chr'] = Train_Set_Neg['Neg_Peak'].map(lambda x: x.split('_')[0])
        Train_Set_Neg['start'] = Train_Set_Neg['Neg_Peak'].map(lambda x: int(x.split('_')[1]))
        Train_Set_Neg['end'] = Train_Set_Neg['Neg_Peak'].map(lambda x: int(x.split('_')[2]))
        Train_Set_Neg.loc[:, ['chr', 'start', 'end']].to_csv('Final_Peaks/'+temp_file.split('.txt')[0]+'_Fold'+str(it)+'_Neg_Test.txt', 
                                                             sep='\t', header=None, index=None)


# Step5: Get .fa file
for temp_celltype in celltype_list:
    for temp_fold in Fold_list:
        input_file = os.path.join('Final_Peaks', temp_celltype+'_Human_'+temp_fold+'_Pos_Test.txt')
        output_file = os.path.join('FA_file', temp_celltype+'_Human_'+temp_fold+'_Pos_Test.fa')
        temp_df = pd.read_csv(input_file, sep='\t', header=None)
        temp_df.columns = ['chr', 'start', 'end']
        temp_df['Peak'] = temp_df.apply(lambda x: x['chr']+'_'+str(x['start'])+'_'+str(x['end']), axis=1)
        
        fa_peaks = pd.Series()
        for it in temp_df.index:
            temp_peaks = temp_df.loc[it, ['Peak']].values[0]
            temp_str = temp_df.loc[it, 'chr']
            temp_start = int(temp_df.loc[it, 'start'])
            temp_end = int(temp_df.loc[it, 'end'])
            
            fa_peaks.loc[fa_peaks.shape[0],] = '>'+temp_peaks
            
            chrs = fa.fetch(temp_str, temp_start, temp_end)
            
            fa_peaks.loc[fa_peaks.shape[0],] = str.upper(chrs)
        fa_peaks.to_csv(output_file, sep='\t', header=None, index=None)
for temp_celltype in celltype_list:
    for temp_fold in Fold_list:
        input_file = os.path.join('Final_Peaks', temp_celltype+'_Human_'+temp_fold+'_Neg_Test.txt')
        output_file = os.path.join('FA_file', temp_celltype+'_Human_'+temp_fold+'_Neg_Test.fa')
        temp_df = pd.read_csv(input_file, sep='\t', header=None)
        temp_df.columns = ['chr', 'start', 'end']
        temp_df['Peak'] = temp_df.apply(lambda x: x['chr']+'_'+str(x['start'])+'_'+str(x['end']), axis=1)
        
        fa_peaks = pd.Series()
        for it in temp_df.index:
            temp_peaks = temp_df.loc[it, ['Peak']].values[0]
            temp_str = temp_df.loc[it, 'chr']
            temp_start = int(temp_df.loc[it, 'start'])
            temp_end = int(temp_df.loc[it, 'end'])
            
            fa_peaks.loc[fa_peaks.shape[0],] = '>'+temp_peaks
            
            chrs = fa.fetch(temp_str, temp_start, temp_end)
            
            fa_peaks.loc[fa_peaks.shape[0],] = str.upper(chrs)
        fa_peaks.to_csv(output_file, sep='\t', header=None, index=None)
for temp_celltype in celltype_list:
    for temp_fold in Fold_list:
        input_file = os.path.join('Final_Peaks', temp_celltype+'_Human_'+temp_fold+'_Neg_Train.txt')
        output_file = os.path.join('FA_file', temp_celltype+'_Human_'+temp_fold+'_Neg_Train.fa')
        temp_df = pd.read_csv(input_file, sep='\t', header=None)
        temp_df.columns = ['chr', 'start', 'end']
        temp_df['Peak'] = temp_df.apply(lambda x: x['chr']+'_'+str(x['start'])+'_'+str(x['end']), axis=1)
        
        fa_peaks = pd.Series()
        for it in temp_df.index:
            temp_peaks = temp_df.loc[it, ['Peak']].values[0]
            temp_str = temp_df.loc[it, 'chr']
            temp_start = int(temp_df.loc[it, 'start'])
            temp_end = int(temp_df.loc[it, 'end'])
            
            fa_peaks.loc[fa_peaks.shape[0],] = '>'+temp_peaks
            
            chrs = fa.fetch(temp_str, temp_start, temp_end)
            
            fa_peaks.loc[fa_peaks.shape[0],] = str.upper(chrs)
        fa_peaks.to_csv(output_file, sep='\t', header=None, index=None)
for temp_celltype in celltype_list:
    for temp_fold in Fold_list:
        input_file = os.path.join('Final_Peaks', temp_celltype+'_Human_'+temp_fold+'_Pos_Train.txt')
        output_file = os.path.join('FA_file', temp_celltype+'_Human_'+temp_fold+'_Pos_Train.fa')
        temp_df = pd.read_csv(input_file, sep='\t', header=None)
        temp_df.columns = ['chr', 'start', 'end']
        temp_df['Peak'] = temp_df.apply(lambda x: x['chr']+'_'+str(x['start'])+'_'+str(x['end']), axis=1)
        
        fa_peaks = pd.Series()
        for it in temp_df.index:
            temp_peaks = temp_df.loc[it, ['Peak']].values[0]
            temp_str = temp_df.loc[it, 'chr']
            temp_start = int(temp_df.loc[it, 'start'])
            temp_end = int(temp_df.loc[it, 'end'])
            
            fa_peaks.loc[fa_peaks.shape[0],] = '>'+temp_peaks
            
            chrs = fa.fetch(temp_str, temp_start, temp_end)
            
            fa_peaks.loc[fa_peaks.shape[0],] = str.upper(chrs)
        fa_peaks.to_csv(output_file, sep='\t', header=None, index=None)