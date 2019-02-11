
# scRNAseqData
We have dowload 13 single cell RNA sequence (scRNA-seq) datasets from Gene Expression Omnibus (GEO), including blood, nerve, pancrease and so on. All the datasets are collected from the normal, rathre than patients. 

## data process
The code (data_process.R) is the procedure for us to process the datasets.

## data summary
| Data | Tissue | cell type |
| :------| :------ | :------ |
| GSE67835 | nerve | 9 |
| GSE70580 | tonsil | 4 |
| GSE73721 | nerve | 5 |
| GSE74310 | blood | 2 |
| GSE76381 | nerve | 27 |
| GSE81252 | liver | 15 |
| GSE81608 | pancreas | 8 |
| GSE83139 | pancreas | 8 |
| GSE84133 | pancreas | 14 |
| GSE89232 | blood | 3 |
| GSE94820 | blood | 28 |
| GSE102956 | nerve | 4 |
| GSE113197 | breast | 2 |

## differentially expressed gene from five methods
Identification of cell-subpopulation specific expression levels is the first step for integrative analysis. Intuitively, if a gene is expressed specifically in a cell-subpopulation, then the gene will likely have higher effects or causal probability for the GWAS trait. Here, we select five different methods, including zingeR, edgeR, MAST, t-statistics and high expression. Each method gives us P value for specific cell type. Then, we select top 10% genes, added a 100-kb window around their transcribed region. To show the performance of DE methods, we use LDSC to interactively analyze four neurological scRNA-seq datasets (GSE102956, GSE67835, GSE73721 and GSE76381) with summary data of 43 traits. The analysis includes 7,277 genes and 235,027 SNPs. The heatmaps and correlation plots of the five methods are shown in Figure 1 and Figure 2.
