Single-cell TWAS of Lung Cancer

This repository contains summary statistics from a single-cell transcriptome-wide association study (scTWAS) of lung cancer, along with R scripts used to develop gene expression prediction models across 38 lung cell types in non-neoplastic lung tissues.

Below is a brief description of the files included in this repository.

Code
1. enet.r - R script for developing cell-type-specific gene expression prediction models using elastic net regression.

2. utmost.r - R script for developing cross-cell-type gene expression prediction models across 38 lung cell types using the UTMOST framework.

Input Data - The pseudobulk gene expression data required to run the scripts above are publicly available from the Gene Expression Omnibus (GEO: GSE227136). Genotype data used for model development were obtained from dbGaP under controlled access (dbGaP phs003521). Due to dbGaP data access policies, individual-level genotype data cannot be publicly redistributed in this repository. Both the pseudobulk gene expression and genotype datasets were originally reported in Natri HM et al., Nature Genetics, 2024 (PMID: 38548990). Researchers interested in reproducing the analyses should obtain the required datasets directly from GEO and dbGaP.

scTWAS Summary Statistics: The repository also contains summary statistics from scTWAS analyses for overall lung cancer and major histological subtypes. Each archive contains association statistics for all genes tested in each of the 38 cell types.

1. Overall_Lung_Cancer.zip - Summary statistics of scTWAS results for overall lung cancer risk.

2. LUAD.zip - Summary statistics of scTWAS results for lung adenocarcinoma (LUAD).

3. LSCC.zip - Summary statistics of scTWAS results for lung squamous cell carcinoma (LSCC).

4. SCLC.zip - Summary statistics of scTWAS results for small-cell lung cancer (SCLC).

Contact

For questions regarding the use of the code or replication of the analyses, please contact the corresponding author:

Dr. Yaohua Yang. Assistant Professor of Genome Sciences, University of Virginia School of Medicine. Email: vta8we@virginia.edu

We will do our best to respond within three working days.
