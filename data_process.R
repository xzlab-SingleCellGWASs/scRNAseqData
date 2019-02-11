library(data.table)
###############
###GSE102956
###############
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/GSE102956_neuron/"
setwd(path)
files <- list.files()
exp_mat_1 <- data.frame(fread(files[1]))
exp_mat_2 <- data.frame(fread(files[2]))
exp_mat_3 <- data.frame(fread(files[3]))
exp_mat_4 <- data.frame(fread(files[4]))
exp_mat <- cbind(exp_mat_1, exp_mat_2[, -1], exp_mat_3[, -1], exp_mat_4[, -1])
write.table(exp_mat, file = "expression_mat.csv", row.names = F, quote = F)
cell_type <- substr(colnames(exp_mat)[-1], 1, 3)
write.table(cell_type, file = "cell_type.csv", row.names = F, quote = F)
rm(list = ls())

###############
###GSE70580
###############
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/GSE70580_tonsil"
setwd(path)
files <- list.files()
gene_num <- nrow(read.table(gzfile(files[1]), header = T))
gene <- read.table(gzfile(files[1]), header = T, stringsAsFactors = F)[, 1]
exp_mat <- matrix(NA, ncol = length(files), nrow = gene_num)
for (i in 1: length(files)){
  exp_mat[, i] <- read.table(gzfile(files[i]), header = T)[, 4]
}
exp_mat <- cbind(gene, exp_mat)
colnames(exp_mat)[-1] <- paste0("S", c(1: length(files)))
cell_type <- laply(strsplit(files, "_"), function(a) a[5])
write.csv(cell_type, file = "cell_type.csv", row.names = F, quote = F)
write.csv(exp_mat, file = "expression_mat.csv", row.names = F, quote = F)
rm(list = ls())

###############
###GSE70580
###############
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/GSE73721_brain"
setwd(path)
exp_mat <- read.csv("GSE73721_Human_and_mouse_table.csv", header = T, stringsAsFactors = F)
exp_mat <- exp_mat[, -c(2:5, 39:46)]
cell_type <- laply(strsplit(colnames(exp_mat)[-1], "\\."), function(a) a[length(a)])
exp_mat[, 1] <- toupper(exp_mat[, 1])
write.csv(cell_type, file = "cell_type.csv", row.names = F, quote = F)
write.csv(exp_mat, file = "expression_mat.csv", row.names = F, quote = F)

##############
###
##############
rm(list = ls())
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/GSE75478_bone"
setwd(path)
files <- list.files()[-1]
gene_num <- nrow(read.csv(gzfile(files[1]), header = T))
gene <- data.frame(read.csv(gzfile(files[1]), header = T, stringsAsFactors = F)[, 1])
colnames(gene) <- "ensembl"
exp_mat <- matrix(NA, ncol = length(files), nrow = gene_num)
for (i in 1: length(files)){
  exp_mat[, i] <- read.csv(gzfile(files[i]), header = T)[, 2]
}
#GSM2289606_counts_I1_plate6_A_5.csv.gz is empty file, 2035 file
colnames(exp_mat)[-1] <- paste0("S", c(1: length(files)))
write.csv(exp_mat, file = "expression_mat.csv", row.names = F, quote = F)
cell_info <- read.table("GSE75478_series_matrix.txt", stringsAsFactors = F)
cell_type <- alply(strsplit(t(cell_info[1, -2035]), " "), function(a) a[1])
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id", "hgnc_symbol"),
                values=gene[, 1],mart= mart)
dim(G_list)
# ENSG00000187510 duplicate ensembl gene
G_list <- G_list[-15835, ]
gene_syb <- left_join(gene, G_list, by = c("ensembl" = "ensembl_gene_id"))
exp_mat <- cbind(gene_syb, exp_mat)
exp_mat <- exp_mat[-which(is.na(exp_mat[, 2])), ]
exp_mat <- exp_mat[, -1]
write.csv(exp_mat, file = "expression_mat.csv", row.names = F, quote = F)
write.csv(cell_type, file = "cell_type.csv", row.names = F, quote = F)

##############
###GSE81252
##############
rm(list = ls())
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/GSE81252_liver"
setwd(path)
cell_info <- t(read.table("GSE81252_series_matrix.txt", stringsAsFactors = F))
lineage <- data.frame(fread("GSE81252_data.cast.log2.lineage.csv", header = T, stringsAsFactors = F))
bud <- data.frame(fread("GSE81252_data.cast.log2.liverbud.csv", header = T, stringsAsFactors = F))
sample_info <- read.table("sample.csv", header = T, stringsAsFactors = F)
sample_info <- laply(strsplit(sample_info[, 1], "\\_"), function(a) paste0(a[2], "_", a[1]))

cell_inter <- intersect(lineage[, 1], bud[, 1])
bud_inter <- bud[-which(bud[, 1] %in% cell_inter), -3]
exp_mat <- rbind(lineage, bud_inter)

cell_inter2 <- intersect(exp_mat[, 1], cell_info[, 4])
cell_info2 <- cell_info[cell_info[, 4]%in%cell_inter2, ]
exp_mat2 <- exp_mat[exp_mat[, 1] %in% cell_inter2, ]
cell_info2 <- data.frame(cell_info2)
cell_info2[, 4] <- factor(cell_info2[, 4], levels = exp_mat2[, 1])
cell_info2 <- cell_info2[order(cell_info2[, 4]), ]

exp_mat2 <- t(exp_mat2)
cell_type <- gsub("cell type or stage: ", "", cell_info2[, 3])
cell_type <- gsub(" ", "_", cell_type)
write.csv(exp_mat2[-2, ], file = "expression_mat.csv", row.names = F, quote = F, col.names = F)
write.csv(cell_info2[, 4], file = "cell_type.csv", row.names = F, quote = F)

##############
###GSE81547
##############
rm(list = ls())
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/GSE81547_pancreas"
setwd(path)
files <- list.files()[-1]
gene_num <- nrow(read.csv(gzfile(files[1]), header = T))
gene <- data.frame(read.csv(gzfile(files[1]), header = T, stringsAsFactors = F)[, 1])
cell_info <- t(read.table("GSE81547_series_matrix.txt"))
exp_mat <- matrix(NA, nrow = gene_num, ncol = length(file))

##############
###GSE81608
##############
rm(list = ls())
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/GSE81608_pancreas"
setwd(path)
exp_mat <- data.frame(fread("GSE81608_log_counts.csv", header = T))
cell_type <- read.table("GSE81608_celltype.csv", header = T)


##############
###GSE85908
##############
rm(list = ls())
exp_mat <- read.csv("GSE85908_expression.csv", header = T)
exp_mat <- t(exp_mat)
ensembl <- data.frame(rownames(exp_mat))
colnames(ensembl) <- "ensembl"
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id", "hgnc_symbol"),
                values=ensembl[, 1], mart= mart)
#ENSG00000187510, ENSG00000230417 duplicate
G_list <- G_list[-c(21154,15961), ]
gene_syb <- left_join(ensembl, G_list, by = c("ensembl" = "ensembl_gene_id"))
exp_mat <- cbind(gene_syb, exp_mat)
exp_mat <- exp_mat[-which(is.na(exp_mat[, 2])), ]
exp_mat <- exp_mat[-which(exp_mat[, 2] == ""), ]
write.csv(exp_mat[, -1], file = "expression_mat.csv", row.names = F, quote = F)

##############
###GSE89232
##############
rm(list = ls())
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/GSE89232_blood"
setwd(path)
cell_info <- t(read.table("GSE89232_series_matrix.txt"))
exp_mat <- data.frame(fread("GSE89232_expMatrix.txt"))
write.csv(exp_mat, file = "expression_mat.csv", row.names = F, quote = F)

##############
###GSE94820
##############
rm(list = ls())
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/GSE94820_blood"
setwd(path)
mono <- data.frame(fread("GSE94820_raw.expMatrix_DCnMono.discovery.txt"))
deep <- data.frame(fread("GSE94820_raw.expMatrix_deeper.characterization.txt"))
exp_mat <- cbind(mono , deep[, -1])
cell_info1 <- t(read.table("GSE94820-GPL15520_series_matrix.txt", header = T))
cell_info2 <- t(read.table("GSE94820-GPL16791_series_matrix.txt", header = T))
cell_info1 <- cbind(rownames(cell_info1), cell_info1[, 1])
cell_info2 <- cbind(rownames(cell_info2), cell_info2[, 1])
cell_info <- rbind(cell_info1, cell_info2)
cell_info[, 2] <- laply(strsplit(cell_info[, 2], ": "), function(a) a[2])
cell_info <- data.frame(cell_info)
cell_info[, 1] <- factor(cell_info[, 1], levels = colnames(exp_mat)[-1])
cell_type <- cell_info[order(cell_info[, 1]), 2]
write.csv(exp_mat, file = "expression_mat.csv", row.names = F, quote = F)

##############
###GSE113197
##############
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/GSE113197_breast"
setwd(path)
cell_type1 <- t(read.table("GSE113099_series_matrix.txt"))
cell_type2 <- t(read.table("GSE113127_series_matrix.txt"))
cell_type3 <- t(read.table("GSE113198_series_matrix.txt"))
cell_type <- c(cell_type1[, 1], cell_type2[, 1], cell_type3[, 1])
files <- list.files()[grep("GSM", list.files())]
gene_num <- nrow(read.csv(gzfile(files[1]), header = T))
exp_mat1 <- matrix(NA, ncol = length(files), nrow = gene_num)
for (i in 1: length(files)){
  exp_mat1[, i] <- read.table(gzfile(files[i]), header = T)[, 2]
}
exp_mat2 <- matrix(NA, ncol = length(files), nrow = gene_num)
for (i in 1: length(files)){
  exp_mat2[, i] <- read.table(gzfile(files[i]), header = T)[, 2]
}
files <- list.files()[grep("GSM", list.files())]
exp_mat3 <- matrix(NA, ncol = length(files), nrow = gene_num)
for (i in 1: length(files)){
  exp_mat3[, i] <- read.table(gzfile(files[i]), header = T)[, 2]
}
gene <- read.table(gzfile(files[1]), header = T, stringsAsFactors = F)[, 1]
exp_mat <- cbind(gene, exp_mat)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id", "hgnc_symbol"),
                values=exp_mat[, 1], mart= mart)
#ENSG00000187510, ENSG00000230417, ENSG00000276085 duplicate
G_list <- G_list[-c(15919,20904, 32227), ]
ensembl <- data.frame(ensembl = exp_mat[, 1])
gene_syb <- left_join(ensembl, G_list, by = c("ensembl" = "ensembl_gene_id"))
exp_mat <- cbind(gene_syb, exp_mat)
exp_mat <- exp_mat[, -c(1, 3, 4)]
exp_mat <- exp_mat[-which(exp_mat[, 1] == ""), ]
write.csv(exp_mat, file = "expression_mat.csv", row.names = F, quote = F)
write.csv(cell_type, file = "cell_type.csv", row.names = F, quote = F)


##############
###GSE74310
##############
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/GSE74310_blood"
setwd(path)
exp_mat <- data.frame(fread("GSE74310_scATACseq_All_Counts.txt"))
exp_gene_mat <- exp_mat[, c(1:3)]
gene_loci <- data.frame(fread("/net/mulan/home/yasheng/summAnnot/real_data/summary_data/Gene_annotation/gencode_gene_V25.csv"))
gene_loci <- gene_loci[, c(1, 4, 5, 11)]
chr <- paste0("chr", c(1:22))
gene_label <- vector()
for (i in 2:22){
  print(i)
  gene_loci_temp <- gene_loci[gene_loci[, 1] == chr[i], ]
  exp_gene_mat_temp <- exp_gene_mat[exp_gene_mat[, 1] == chr[i], ]
  for (j in 1: nrow(exp_gene_mat_temp)){
    
    ro <- which(gene_loci_temp[, 3] > exp_gene_mat_temp[j, 2])
    if (length(ro) == 0){
      gene_label <- c(gene_label, NA)
    }else{
      temp <- gene_loci_temp[min(which(gene_loci_temp[, 3] > exp_gene_mat_temp[j, 2])), ]
      if(temp[, 3] > exp_gene_mat_temp[j, 3]){
        
        gene_label <- c(gene_label, as.character(temp[4]))
      }else{
        gene_label <- c(gene_label, NA)
      }
    }
  }
}
exp_mat <- exp_mat[which(exp_mat[, 1] == "chrX"), ]
gene_label_idx <- which(is.na(gene_label))
gene_label_na <- gene_label[-gene_label_idx]
exp_mat_na <- exp_mat[-gene_label_idx, -c(1:3)]

gene_label_uni <- unique(gene_label_na)
exp_mat_fin <- matrix(NA, nrow = length(gene_label_uni), ncol = ncol(exp_mat_na))
for (i in 1: length(gene_label_uni)){
  exp_mat_fin[i, ] <- colSums(exp_mat_na[which(gene_label_na == gene_label_uni[i]), ])
}
exp_mat1 <- cbind(gene_label_uni, exp_mat_fin)
colnames(exp_mat1) <- c("gene", colnames(exp_mat)[-c(1:3)])
cell_info <- t(read.table("GSE74310_series_matrix.txt"))
sample_info <- gsub("\\.", "-", colnames(exp_mat1)[-1])
cell_info <- data.frame(cell_info)
cell_info[, 1] <- factor(cell_info[, 1], levels = sample_info)
exp_mat2 <- exp_mat1[, c(1, grep("LMPP", colnames(exp_mat1)), grep("mono", colnames(exp_mat1)))]
write.csv(exp_mat2, file = "expression_mat.csv", row.names = F, quote = F)

###gene intersect
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/"
data_summ <- read.table(paste0(path, "data_summ.txt"), header = F, stringsAsFactors = F)[, 1]
gene <- vector()
for (i in 1: length(data_summ)){
  gene_temp <- data.frame(fread(paste0(path, data_summ[i], "/expression_mat.csv")))[, 1]
  gene_temp <- unique(gene_temp)
  gene <- c(gene, gene_temp)
}
gene_freq <- data.frame(table(gene))
gene_inter <- gene_freq[gene_freq[, 2] == 16, 1]
gene_inter <- gene_inter[-1]
gene_inter <- as.character(gene_inter)
write.table(gene_inter, file = paste0(path, "inter_gene.txt"), row.names = F, quote = F)

###expression intersect
path <- "/net/mulan/home/yasheng/summAnnot/analysis/single_cell_data/"
data_summ <- read.table(paste0(path, "data_summ.txt"), header = F, stringsAsFactors = F)[, 1]
for (i in 1: length(data_summ)){

  gene_temp <- data.frame(fread(paste0(path, data_summ[i], "/expression_mat.csv")))
  gene_temp <- gene_temp[gene_temp[, 1] %in% gene_inter, ]
  print(nrow(gene_temp))
  write.csv(gene_temp, file = paste0(path, data_summ[i], "/inter_exp_mat.csv"), 
            row.names = F, quote = F)
}
                        
