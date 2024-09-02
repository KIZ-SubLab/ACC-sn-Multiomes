library(ArchR)
library(Signac)
library(Seurat)
library(Matrix)
library(readr)

celltype_list = c('EX', 'IN', 'ASC', 'ODC', 'OPC', 'MG')
Human_Sample = c('HM2013017', 'HM20200905', 'HM20200927', 'HM20201129',
                'HM20201213', 'HM20201222')
Monkey_Sample = c('RM06013', 'RM07071', 'RM07333', 'RM08403', 'RM1102009')


# Get pesudobulk accessibility (Human)
adata = loadArchRProject('../Human_ATAC/HumanBrain')
meta = adata@cellColData

temp_df <- read_delim("../Human_ATAC/Peak_df/Union/Count_df.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
peaknames <- read_delim("../Human_ATAC/Peak_df/Union/Count_Peaks.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
cellnames <- read_delim("../Human_ATAC/Peak_df/Union/Count_Cells.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)

datafr = sparseMatrix(i=temp_df$Peaks, j=temp_df$Cells, x=1)
colnames(datafr) = cellnames$X1
rownames(datafr) = peaknames$X1

for (temp_sample in Human_Sample){
    for (temp_celltype in celltype_list){
        output_file = paste0(temp_sample, '_', temp_celltype, '.txt')
        
        Used_meta = meta[meta$Annotaion==temp_celltype, ]
        Used_meta = Used_meta[Used_meta$Sample==temp_sample,]
        Used_cell = row.names(Used_meta)
        Bulk_df = rowSums(datafr[,Used_cell])
        write.table(Bulk_df, output_file, sep='\t', quote=FALSE) 
    }
}


# Get pesudobulk accessibility (Macaque)
adata = loadArchRProject('../Macaque_ATAC/MacaqueBrain')
meta = adata@cellColData

temp_df <- read_delim("../Monkey_ATAC/Peak_df/Union/Count_df.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
peaknames <- read_delim("../Monkey_ATAC/Peak_df/Union/Count_Peaks_lifted.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
cellnames <- read_delim("../Monkey_ATAC/Peak_df/Union/Count_Cells.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)

datafr = sparseMatrix(i=temp_df$Peaks, j=temp_df$Cells, x=1)
colnames(datafr) = cellnames$X1
rownames(datafr) = peaknames$X1

for (temp_sample in Monkey_Sample){
    for (temp_celltype in celltype_list){
        output_file = paste0(temp_sample, '_', temp_celltype, '.txt')
        
        Used_meta = meta[meta$Annotaion==temp_celltype, ]
        Used_meta = Used_meta[Used_meta$Sample==temp_sample,]
        Used_cell = row.names(Used_meta)
        Bulk_df = rowSums(datafr[,Used_cell])
        write.table(Bulk_df, output_file, sep='\t', quote=FALSE) 
    }
}



#####################################################################
# Run 04_02_DAcCRE_Identification.py
#####################################################################



# edgeR test for Human
library(edgeR)
library(stringr)

celltype_list = c('EX', 'IN', 'ODC', 'OPC', 'ASC', 'MG')
Human_Sample = c('HM2013017', 'HM20200905', 'HM20200927', 'HM20201129',
                'HM20201213', 'HM20201222')
Monkey_Sample = c('RM06013', 'RM07071', 'RM07333', 'RM08403', 'RM1102009')

Used_order = c()
for (temp_celltype in celltype_list){
    for (temp_sample in Human_Sample){
        Used_order = c(Used_order, paste0(temp_sample, '_', temp_celltype))
    }
}

bcv = 0.1
for (temp_celltype in celltype_list){
    Count_df = read.delim('Settings/Human_df.txt', sep='\t', row.names=1)
    Used_celltype = c()
    for (it in Used_order){
        if (strsplit(it, '_')[[1]][2]==temp_celltype){
            Used_celltype = c(Used_celltype, temp_celltype)
        }
        else{
            Used_celltype = c(Used_celltype, 'Other')
        }
    }
    Count_df = Count_df[, Used_order]
    meta = data.frame(Celltype=Used_celltype, row.names=Used_order)

    data.use = as.matrix(Count_df)
    group <- factor(meta$Celltype)
    design <- model.matrix(~group)
    y <- DGEList(counts=data.use, group=group)

    tb.pos <- exactTest(y, dispersion=bcv^2)$table;
    
    tb.pos$padj = p.adjust(tb.pos$PValue, method = 'hochberg')#, n = length(p))

    output_file = paste0('Res/', temp_celltype, '_Human.txt')
    write.table(tb.pos, output_file, sep='\t', quote=FALSE)
}


# edgeR test for Macaque
library(edgeR)
library(stringr)

celltype_list = c('EX', 'IN', 'ODC', 'OPC', 'ASC', 'MG')
Human_Sample = c('HM2013017', 'HM20200905', 'HM20200927', 'HM20201129',
                'HM20201213', 'HM20201222')
Monkey_Sample = c('RM06013', 'RM07071', 'RM07333', 'RM08403', 'RM1102009')

Used_order = c()
for (temp_celltype in celltype_list){
    for (temp_sample in Monkey_Sample){
        Used_order = c(Used_order, paste0(temp_sample, '_', temp_celltype))
    }
}

bcv = 0.1
for (temp_celltype in celltype_list){
    Count_df = read.delim('Settings/Macaque_df.txt', sep='\t', row.names=1)
    Used_celltype = c()
    for (it in Used_order){
        if (strsplit(it, '_')[[1]][2]==temp_celltype){
            Used_celltype = c(Used_celltype, temp_celltype)
        }
        else{
            Used_celltype = c(Used_celltype, 'Other')
        }
    }
    Count_df = Count_df[, Used_order]
    meta = data.frame(Celltype=Used_celltype, row.names=Used_order)

    data.use = as.matrix(Count_df)
    group <- factor(meta$Celltype)
    design <- model.matrix(~group)
    y <- DGEList(counts=data.use, group=group)

    tb.pos <- exactTest(y, dispersion=bcv^2)$table;
    
    tb.pos$padj = p.adjust(tb.pos$PValue, method = 'hochberg')#, n = length(p))

    output_file = paste0('Res/', temp_celltype, '_Macaque.txt')
    write.table(tb.pos, output_file, sep='\t', quote=FALSE)
}


# edgeR test for cross-species comparision
for (temp_celltype in celltype_list){
    Used_order = c()
    for (it in Human_Sample){
        Used_order = c(Used_order, paste0(it, '_', temp_celltype))
    }
    for (it in Monkey_Sample){
        Used_order = c(Used_order, paste0(it, '_', temp_celltype))
    }
    Used_celltype = c('Human', 'Human', 'Human', 'Human', 'Human', 'Human', 
                      'Monkey', 'Monkey', 'Monkey', 'Monkey', 'Monkey')
    meta = data.frame(Celltype=Used_celltype, row.names=Used_order)

    input_file = paste0('Settings/', temp_celltype, '_Cross.txt')
    output_file = paste0('Res/', temp_celltype, '_Cross_MvsH.txt')
    Count_df = read.delim(input_file, sep='\t', row.names=1)
    Count_df = Count_df[, Used_order]


    data.use = as.matrix(Count_df)
    group <- factor(meta$Celltype)
    design <- model.matrix(~group)
    y <- DGEList(counts=data.use, group=group)
    
    tb.pos <- exactTest(y, dispersion=bcv^2)$table;
    
    tb.pos$padj = p.adjust(tb.pos$PValue, method = 'hochberg')#, n = length(p))

    write.table(tb.pos, output_file, sep='\t', quote=FALSE)
}