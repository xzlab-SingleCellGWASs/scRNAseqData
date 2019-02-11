#!/usr/bin/env RScript
library(data.table)
library(optparse)
library(doParallel)

##Parameter setting
args_list = list(
  make_option("--data", type="numeric", default=NULL,
              help="INPUT: loop number", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

data_code <- read.table("/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/data_summ.txt")[opt$data, 1]

scSeq_path <- paste0("/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/", data_code, "/")
summ_path <- "/net/mulan/home/yasheng/summAnnot/analysis/gwas_data/range_100Kb/summary_data/"
de_path <- "/net/mulan/home/yasheng/summAnnot/analysis/result/single_cell_DE/"

###
setwd(scSeq_path)
cell_type <- read.table("cell_type.csv", header = T, stringsAsFactors = F)
cell_type_uni <- unique(cell_type[, 1])
scseq_mat <- data.frame(fread("inter_exp_mat.csv"))
summ_data <- read.table(paste0(summ_path, "pheno_1.summ.gz"), header = T, stringsAsFactors = F)
summ_gene <- unique(summ_data[, 2])
scseq_mat <- scseq_mat[scseq_mat[, 1] %in% summ_gene, ]
scseq_mat[, 1] <- factor(scseq_mat[, 1], levels = summ_gene)
scseq_mat <- scseq_mat[order(scseq_mat[, 1]), ]
gene_freq <- data.frame(table(scseq_mat[, 1]))
get_row <- vector()
for (i in 1: dim(gene_freq)[1]){
  if (gene_freq[i, 2] == 1){
    get_row <- c(get_row, which(scseq_mat[, 1] == gene_freq[i, 1]))
  }else{
    get_row <- c(get_row, which(scseq_mat[, 1] == gene_freq[i, 1])[1])
  }
}
scseq_mat <- scseq_mat[get_row, ]
row.names(scseq_mat) <- scseq_mat[, 1]
scseq_mat <- scseq_mat[, -1]


### code of scSeq data(GSEXXXXX)
scSeq_code <- substr(scSeq_path, 61, c(nchar(scSeq_path)-1))
setwd(de_path)
system(paste0("mkdir ", scSeq_code))
# rank
cutoff_rank <- floor(nrow(scseq_mat) * 0.1)

##ct = number of cell types
cl <- makeCluster(length(cell_type_uni))
registerDoParallel(cl)

ys <- foreach (cell=1: ct) %dopar% {
  
  source("/net/mulan/home/yasheng/summAnnot/analysis/code/DE_gene.R")
  
  setwd(paste0(de_path, scSeq_code))
  ### DE analysis (MAST)
  system(paste0("mkdir ", cell_type_uni[cell]))
  ind_info <- data.frame(cell_type = cell_type[, 1], total_counts = colSums(scseq_mat))
  MAST_mat <- MAST_ANAL(scseq_mat, ind_info, cell_type_uni[cell], ncol(scseq_mat)/10)
  
  ### DE analysis (zinger_DEseq)
  traits <- ifelse(grepl(cell_type_uni[cell], ind_info[, 1]), "cas", "con")
  zingeR_mat <- zingeR_DESeq2(scseq_mat, traits)
  
  ### DE analysis (edgeR)
  edgeR_mat <- edgeR_DE(scseq_mat, traits)
  
  ### DE analysis (t test)
  t_mat <- t_DE(scseq_mat, traits)
  
  ### result for LDSC
  MAST_DE_mat <- MAST_mat[order(MAST_mat$fdr)[c(1: cutoff_rank)], ]
  zingeR_DE_mat <- rownames(scseq_mat)[order(zingeR_mat@listData$pvalue)][c(1: cutoff_rank)]
  edgeR_DE_mat <- rownames(scseq_mat)[order(edgeR_mat$PValue)[c(1: cutoff_rank)]]
  exp_high_mat <- rownames(scseq_mat)[order(abs(edgeR_mat$logFC))[c(1: cutoff_rank)]]
  t_DE_mat <- names(t_mat)[order(t_mat)[c(1: cutoff_rank)]]
  
  ### result for TORUS
  MAST_anno <- cbind(summ_data[, 1],
                    ifelse(summ_data[, 2] %in% as.data.frame(MAST_DE_mat[, 1])[, 1], 1, 0))
  zingeR_anno <- cbind(summ_data[, 1],
                       ifelse(summ_data[, 2] %in% zingeR_DE_mat, 1, 0))
  edgeR_anno <- cbind(summ_data[, 1],
                     ifelse(summ_data[, 2] %in% edgeR_DE_mat, 1, 0))
  high_anno <- cbind(summ_data[, 1],
                     ifelse(summ_data[, 2] %in% exp_high_mat, 1, 0))
  t_anno <- cbind(summ_data[, 1],
                  ifelse(summ_data[, 2] %in% t_DE_mat, 1, 0))
  colnames(t_anno) <- colnames(high_anno) <- colnames(edgeR_anno) <- colnames(zingeR_anno) <- c("SNP", "feature_d")
  
  ### output for LDSC
  write.table(MAST_DE_mat, file = paste0(de_path, scSeq_code, "/", cell_type_uni[cell], "/MAST_sig_rank_10.txt"),
              quote = F, row.names = F, col.names = F)
  write.table(zingeR_DE_mat, file = paste0(de_path, scSeq_code, "/", cell_type_uni[cell], "/zingeR_sig_rank_10.txt"),
              quote = F, row.names = F, col.names = F)
  write.table(edgeR_DE_mat, file = paste0(de_path, scSeq_code, "/", cell_type_uni[cell], "/edgeR_sig_rank_10.txt"),
              quote = F, row.names = F, col.names = F)
  write.table(exp_high_mat, file = paste0(de_path, scSeq_code, "/", cell_type_uni[cell], "/high_sig_rank_10.txt"),
              quote = F, row.names = F, col.names = F)
  write.table(t_DE_mat, file = paste0(de_path, scSeq_code, "/", cell_type_uni[cell], "/t_sig_rank_10.txt"),
              quote = F, row.names = F, col.names = F)
  
  ### output for TORUS      
  MAST_res <- paste0(de_path, scSeq_code, "/", cell_type_uni[cell], "/MAST_anno_rank_10.txt")
  zingeR_res <- paste0(de_path, scSeq_code, "/", cell_type_uni[cell], "/zingeR_anno_rank_10.txt")
  edgeR_res <- paste0(de_path, scSeq_code, "/", cell_type_uni[cell], "/edgeR_anno_rank_10.txt")
  high_res <- paste0(de_path, scSeq_code, "/", cell_type_uni[cell], "/high_anno_rank_10.txt")
  t_res <- paste0(de_path, scSeq_code, "/", cell_type_uni[cell], "/t_anno_rank_10.txt")
  write.table(MAST_anno, file = MAST_res, quote = F, row.names = F)
  write.table(zingeR_anno, file = zingeR_res, quote = F, row.names = F)
  write.table(edgeR_anno, file = edgeR_res, quote = F, row.names = F)
  write.table(high_anno, file = high_res, quote = F, row.names = F)
  write.table(t_anno, file = t_res, quote = F, row.names = F)
  system(paste0("gzip ", MAST_res))
  system(paste0("gzip ", zingeR_res))
  system(paste0("gzip ", edgeR_res))
  system(paste0("gzip ", high_res))
  system(paste0("gzip ", t_res))
}  
