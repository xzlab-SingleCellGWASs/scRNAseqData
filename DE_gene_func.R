############
### MAST ###
############
library(MAST)
library(zingeR)
library(DESeq2)
library(edgeR)
library(plyr)

MAST_ANAL <- function(sc_exp, ind_info, cell_type, sc_cutoff){
  
  #cDat
  ind_label <- vector("character")
  ind_uni <- as.character(unique(ind_info[, 1]))
  for(i in 1: length(ind_uni)){
    
    temp <- as.character(ind_info[ind_info[, 1] == ind_uni[i], 1])
    temp <- paste0(temp, "_", c(1: length(temp)))
    ind_label <- c(ind_label, temp)
  }
  cDat <- data.frame(wellKey = ind_label,
                     condition = ifelse(ind_info[, 1] == cell_type, cell_type, paste0("Non_", cell_type)),
                     count = ind_info[, 2])
  rownames(cDat) <- ind_label
  
  #fDat
  fDat <- data.frame(primerid = rownames(sc_exp), entrez = rownames(sc_exp),
                     symbolid = rownames(sc_exp))
  rownames(fDat) <- rownames(sc_exp)

  #MAST
  res <- MAST_DE(fDat, cDat, sc_exp, cell_type, sc_cutoff)
  
  return(res)
}

MAST_DE <- function(fDat, cDat, sc_exp, cell_type, sc_cutoff){

  colnames(sc_exp) <- rownames(cDat)
  sc_mat <- FromMatrix(as.matrix(sc_exp), cDat, fDat, check_sanity = FALSE)
  filterCrit <- with(colData(sc_mat), count > sc_cutoff)
  sc_mat_sub <- subset(sc_mat, filterCrit)
  cond <- factor(colData(sc_mat_sub)$condition)
  cond <- relevel(cond, paste0("Non_", cell_type))
  colData(sc_mat_sub)$condition <- cond
  zlm_cond <- zlm( ~ condition, sc_mat_sub)
  cond_cell_type <- paste0("condition", cell_type)
  
  summary_cond <- summary(zlm_cond, doLRT = cond_cell_type) 
  summary_dat <- summary_cond$datatable
  fc_hurdle <- merge(summary_dat[contrast==cond_cell_type & component=='H',.(primerid, `Pr(>Chisq)`)], 
                     summary_dat[contrast==cond_cell_type & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], 
                     by='primerid') 
  fc_hurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  
  return(fc_hurdle)
}

#####################
### zingeR-DESeq2 ###
#####################
zingeR_DESeq2 <- function(counts, traits){
  
  iCounts <- apply(counts,  2,  function(x){
    storage.mode(x) <- 'integer'
    return(x) }
  )
  rm(counts)
  colData <- data.frame(traits = traits)

  design <- model.matrix(~traits)
  dse <- DESeqDataSetFromMatrix(countData = iCounts, colData = colData, design = ~traits)
  weights <- zeroWeightsLS(counts = iCounts, design = design, maxit = 500, normalization = "DESeq2_poscounts", colData = colData, designFormula = ~traits, verbose = TRUE)
  assays(dse)[["weights"]] <- weights
  dse <- estimateSizeFactors(dse, type="poscounts")
  dse <- estimateDispersions(dse)
  dse <- nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-2)
  res.ZingeR.DESeq2 = results(dse)
  return(res.ZingeR.DESeq2)
}

#############
### edgeR ###
#############
edgeR_DE <- function(counts, traits){

  cd <- DGEList(counts, group=factor(traits))
  y <- calcNormFactors(cd)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  et <- exactTest(y)$table
  return(et)
}
  
##############
### t test ###
##############
t_DE <- function(counts, traits){
  
  counts <- as.matrix(counts)
  t_stat <- aaply(counts, 1, function(a){
    x <- data.frame(exp = a, label = factor(traits))
    t <- try(abs(t.test(exp~traits, data = x)$statistic), silent = T)
    t <- ifelse(inherits(t, "try-error"), 0, t)
  })
  return(t_stat)
}