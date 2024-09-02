# ACC-sn-Multiomes
Data Analysis Process of the Anterior Cingulate Cortex (ACC) in Humans and Macaque.

The project files include preprocessing of snRNA-seq and snATAC-seq data, differential accessibility peak analysis, peak functional analysis, and analysis of human genomic features (hsSNCs and HARs).

For citation or to learn more, please visit our article: xxxx

For data access, please visit: https://ngdc.cncb.ac.cn/, the BioProject ID is PRJCA015229. It includes filtered ExpressionMatrix files for snRNA-seq and PeakMatrix files for snATAC-seq. For other analysis files (such as fasta files or fragment files for snATAC-seq), please contact via email: sub@mail.kiz.ac.cn; wuhaixu@mail.kiz.ac.cn.

Below is the specific content introduction of the code file：

snRNA-seq:
  01: Preprocessing of snRNA-seq Data and Cell Type/Subtype/Cluster Annotation
  02: Data Analysis of snRNA-seq for VENs

snATAC-seq:
  01: Load Human and Macaque snATAC-seq data
  02: Get Human&Macaque Union Peak Set by LiftOver
  03: Add Union Peak information to Human and Macaque Peak Matrixs
  04: Idetification of DA-cCREs
  05: Intergration of snATAC-seq by scAND and Harmony
  06: TF enrichment analysis
  07: Peak2Gene analysis
  08-09: Construction of Gene Regulation Networks (TF to Genes)
  10: Predicting the Impact of hsSNCs on TF Binding Motifs Using MAGGIE
  11: Predicting the Role of HARs in Different Cell Types Using gkmExplain

snRNA_snATAC_Integration.R:  Integration of snRNA-seq and snATAC-seq data
