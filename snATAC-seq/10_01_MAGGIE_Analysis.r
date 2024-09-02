library(Matrix)
library(GenomicRanges)
library(readr)


# Step1: Get Overlapped SNCs of DA-cCREs
celltype_list = c('ASC', 'EX', 'IN', 'MG', 'ODC', 'OPC')

input_dir = '../DAcCRE/'
output_dir = 'Overlap_SNCs/'

SNC_file = '../SNCs_hg38_995-matched.txt'
SNCs_df = read_delim(SNC_file, '\t')#, col_names=FALSE)

SNCs_df = SNCs_df[, c('hg38_chr', 'hg38_start', 'hg38_end')]
colnames(SNCs_df) = c('chrom', 'start', 'end')
SNCs = GRanges(SNCs_df)


for (temp_celltype in celltype_list){
    input_file = paste0(input_dir, temp_celltype, '_Human.txt')
    temp_df = read.delim(input_file, sep='\t', header=FALSE)
    colnames(temp_df) = c('chrom', 'start', 'end')
    
    diff_peaks = GRanges(temp_df)
    Overlap_id = findOverlaps(SNCs, diff_peaks)
    
    Overlap_SNCs = SNCs[Overlap_id@from,]
    write.table(Overlap_SNCs, paste0(output_dir, temp_celltype, '_Human.bed'), sep='\t',row.names = FALSE, quote=FALSE)
}


###########################################################################
#Run 10_02_MAGGIE_Analysis.py
###########################################################################

system("python ./maggie-master/bin/maggie_fasta_input.py ./fa_seq/EX_Human_non_SNCs.fa ./fa_seq/EX_Human_SNCs.fa -o ./maggie_Output/EX_Human -p 8")
system("python ./maggie-master/bin/maggie_fasta_input.py ./fa_seq/IN_Human_non_SNCs.fa ./fa_seq/IN_Human_SNCs.fa -o ./maggie_Output/IN_Human -p 8")
system("python ./maggie-master/bin/maggie_fasta_input.py ./fa_seq/ASC_Human_non_SNCs.fa ./fa_seq/ASC_Human_SNCs.fa -o ./maggie_Output/ASC_Human -p 8")
system("python ./maggie-master/bin/maggie_fasta_input.py ./fa_seq/ODC_Human_non_SNCs.fa ./fa_seq/ODC_Human_SNCs.fa -o ./maggie_Output/ODC_Human -p 8")
system("python ./maggie-master/bin/maggie_fasta_input.py ./fa_seq/OPC_Human_non_SNCs.fa ./fa_seq/OPC_Human_SNCs.fa -o ./maggie_Output/OPC_Human -p 8")
system("python ./maggie-master/bin/maggie_fasta_input.py ./fa_seq/MG_Human_non_SNCs.fa ./fa_seq/MG_Human_SNCs.fa -o ./maggie_Output/MG_Human -p 8")