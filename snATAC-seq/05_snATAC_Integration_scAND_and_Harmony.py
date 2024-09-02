import pandas as pd
import numpy as np
import os
import sys
import time
import scanpy as sc
import matplotlib.pyplot as plt

# Step1: Concat Peak Matrix of Human and Macaque (Union Peakset)
import scipy.sparse

input_dir_Human = '../Human_ATAC/ArchR/Peak_df/Union'

Count_df = pd.read_csv(os.path.join(input_dir_Human, 'Count_df.txt'), sep='\t')
cells = pd.read_csv(os.path.join(input_dir_Human, 'Count_Cells.txt'), sep='\t', header=None)
cells = cells[0]
peaks = pd.read_csv(os.path.join(input_dir_Human, 'Count_Peaks.txt'), sep='\t', header=None)
peaks = peaks[0]
print('Number of peaks: %d' %peaks.shape[0])
print('Number of cells: %d' %cells.shape[0])

xsp_Human = scipy.sparse.coo_matrix((Count_df['Counts'], (Count_df['Peaks']-1, Count_df['Cells']-1)))

adata_Human = sc.AnnData(xsp_Human.T)
adata_Human.obs['Species'] = 'Human'
adata_Human.var_names = list(peaks)
adata_Human.obs_names = list(cells)

adata_Human.X = scipy.sparse.csr_matrix(adata_Human.X)


input_dir_Monkey = '../Macaque_ATAC/ArchR/Peak_df/Union'

Count_df = pd.read_csv(os.path.join(input_dir_Monkey, 'Count_df.txt'), sep='\t')
cells = pd.read_csv(os.path.join(input_dir_Monkey, 'Count_Cells.txt'), sep='\t', header=None)
cells = cells[0]
peaks = pd.read_csv(os.path.join(input_dir_Monkey, 'Count_Peaks_lifted.txt'), sep='\t', header=None)
peaks = peaks[0]
print('Number of peaks: %d' %peaks.shape[0])
print('Number of cells: %d' %cells.shape[0])

xsp_Monkey = scipy.sparse.coo_matrix((Count_df['Counts'], (Count_df['Peaks']-1, Count_df['Cells']-1)))
xsp_Monkey.shape

adata_Monkey = sc.AnnData(xsp_Monkey.T)
adata_Monkey.obs['Species'] = 'Macaque'
adata_Monkey.var_names = list(peaks)
adata_Monkey.obs_names = list(cells)

adata_Monkey.X = scipy.sparse.csr_matrix(adata_Monkey.X)

# Concat
adata = sc.concat([adata_Human, adata_Monkey])

Peak_Info = pd.read_csv('../Union_Peak_Info.txt', sep='\t')
Consensus_peaks = Peak_Info.loc[Peak_Info['Type']=='Consensus', 'hg38_Peak']

adata = adata[:,np.intersect1d(Consensus_peaks, adata.var_names)]

# Save
xsp = scipy.sparse.coo_matrix(adata.X)
Count_df = pd.DataFrame({'Peaks': xsp.col+1, 'Cells': xsp.row+1, 'Counts': xsp.data})
cells = pd.DataFrame(adata.obs_names)
peaks = pd.DataFrame(adata.var_names)

Count_df.to_csv('MergedPeakMatrix/Count_df.txt', sep='\t', index=None)
cells.to_csv('MergedPeakMatrix/Count_Cells.txt', sep='\t', index=None, header=None)
peaks.to_csv('MergedPeakMatrix/Count_Peaks.txt', sep='\t', index=None, header=None)




# Step2: Running scAND
sys.path.append('../scAND-code/')
import scAND
import scAND_utils

np.random.seed(2023)

input_dir = 'MergedPeakMatrix'
output_dir = 'scAND_Output'

beta_list = np.arange(20)/20
used_beta = 0.8

Count_df = pd.read_csv(os.path.join(input_dir, 'Count_df.txt'), sep='\t')
cells = pd.read_csv(os.path.join(input_dir, 'Count_Cells.txt'), sep='\t', header=None)
cells = cells[0]
peaks = pd.read_csv(os.path.join(input_dir, 'Count_Peaks.txt'), sep='\t', header=None)
peaks = peaks[0]
print('Number of peaks: %d' %peaks.shape[0])
print('Number of cells: %d' %cells.shape[0])

start_time = time.time()
Rep_cells = scAND.Run_scAND(Count_df, d=50, weights=beta_list, cells=cells, peaks=peaks, random_seed=2019)
end_time = time.time()
print('Time Costing: %f s' %(end_time-start_time))


print('The estimated beta is %.2f.' %used_beta)
print('The estimated dimension is %d.' %used_dim)
used_rep = scAND.Get_Result(Rep_cells, beta=used_beta, dim=50)

used_rep.to_csv(os.path.join(output_dir, 'scAND_50.txt'), sep='\t')



# Step3: Running Harmony based on scAND results
import harmonypy as hm
Used_Rep = pd.read_csv('scAND_Output/scAND_50.txt', sep='\t', index_col=0)
adata = sc.AnnData(Used_Rep)
data_mat = adata.X
meta_data = adata.obs

ho = hm.run_harmony(data_mat, meta_data, ['Sample'])#, sigma=0.5)

hm_df = pd.DataFrame(ho.Z_corr, columns=meta_data.index)
adata_hm = sc.AnnData(hm_df.T)

for it in adata.obs.columns:
    adata_hm.obs[it] = adata.obs.loc[adata_hm.obs_names, it]
sc.pp.neighbors(adata_hm, use_rep='scAND_hm')
sc.tl.umap(adata_hm)