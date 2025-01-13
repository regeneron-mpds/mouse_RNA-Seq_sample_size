library(dplyr)
library(irr)
library(vcd)
# library(scales)

# INIT
NUM_TRIALS <- 25
N_RANGE <- 14:3
IGNORE_GENES <- c("Dchs1", "LacZ", "SV40")
PTYPE <- "Padj"
TISSUES_IN_ORDER <- c('heart', 'kidney', 'liver', 'lung')

# set up data frame for storing overlap-with-gold-standard stats
PTYPE <- "Padj"
pairwise.results <- data.frame(expand.grid(tissue = TISSUES_IN_ORDER, 
                                           tpm.cutoff = c(0, 1, 5),
                                           pval.cutoff = c(0.01, 0.05, 0.1, 1),
                                           trial = 1:NUM_TRIALS, N = N_RANGE,
                                           ovlap1 = NA, ovlap1.2=NA, ovlap1.5 = NA, ovlap2=NA, 
                                           uniqA1 = NA, uniqA1.2 = NA, uniqA1.5 = NA, uniqA2 = NA,
                                           uniqB1 = NA, uniqB1.2 = NA, uniqB1.5 = NA, uniqB2 = NA,
                                           ICCval1 = NA, ICCval1.2 = NA, ICCval1.5 = NA, ICCval2 = NA,
                                           cohenK1 = NA, cohenK1.2 = NA, cohenK1.5 = NA, cohenK2 = NA,
                                           cohenK_limTPM1 = NA, cohenK_limTPM1.2 = NA, cohenK_limTPM1.5 = NA, cohenK_limTPM2 = NA,
                                           all.cor = NA, stringsAsFactors = FALSE))
head(pairwise.results)


## obtain ovlap.results stats from DE_runs_subsample/Dchs1*.RDS files
prev.fm.fileA <- ""
for (i in 1:nrow(pairwise.results)) { 
  tissue <- pairwise.results[i, "tissue"]
  tpm.cutoff <- pairwise.results[i, 'tpm.cutoff']
  num.sample <- pairwise.results[i, "N"]
  pval.cutoff <- pairwise.results[i, 'pval.cutoff']
  
  if (i %% 1000 == 1) { message(paste(i, tissue, date())) }
  
  # read in formattedInfo (only once per tissue)
  fm.fileA <- paste0("DE_runs_pairwise/Dchs1_N", num.sample, "_trial", pairwise.results[i,"trial"], "_A.RDS")
  fm.fileB <- paste0("DE_runs_pairwise/Dchs1_N", num.sample, "_trial", pairwise.results[i,"trial"], "_B.RDS")
  
  if (fm.fileA != prev.fm.fileA) {
    fmA <- readRDS(fm.fileA)
    fmB <- readRDS(fm.fileB)
  }
  prev.fm.fileA <- fm.fileA
  
  #shorthand
  pvalsA <- fmA[[PTYPE]][,tissue]
  fcsA <- fmA$Fold.Change[,tissue]
  tpmA <- fmA$meanY.TPM[,tissue]
  downregA <- !is.na(fcsA) & fcsA < 1
  tpmA[downregA] <- fmA$meanX.TPM[downregA, tissue]
  #
  pvalsB <- fmB[[PTYPE]][,tissue]
  fcsB <- fmB$Fold.Change[,tissue]
  tpmB <- fmB$meanY.TPM[,tissue]
  downregB <- !is.na(fcsB) & fcsB < 1
  tpmB[downregB] <- fmB$meanX.TPM[downregB, tissue]
  
  # determine which genes are significant (prior to fold change filtering)
  passA <- !is.na(fcsA) & !is.na(pvalsA) 
  passA <- passA & pvalsA <= pval.cutoff   #pval or padj check
  passA <- passA & tpmA >= tpm.cutoff # tpm check
  #
  passB <- !is.na(fcsB) & !is.na(pvalsB) 
  passB <- passB & pvalsB <= pval.cutoff #pval or padj check
  passB <- passB & tpmB >= tpm.cutoff # tpm check
  
  # correlation between log2 fold changes
  pairwise.results[i, 'all.cor'] <- cor(log2(fcsA), log2(fcsB), use = 'complete.obs', method = 'spearman')

  # derive num significant, number overlapping, etc. Also include ICC & CohenK agreement measures
  for (this.FC in c(1, 1.2, 1.5, 2)) {
    tmpA <- passA & (fcsA >= this.FC | fcsA <= 1/this.FC) # fold change check
    tmpB <- passB & (fcsB >= this.FC | fcsB <= 1/this.FC) # fold change check
    this.signifA <- rownames(fmA[[1]])[tmpA]
    this.signifB <- rownames(fmB[[1]])[tmpB]
    # remove genes we don't care about e.g. LacZ (already done for gold.standards)
    this.signifA <- setdiff(this.signifA, IGNORE_GENES)
    this.signifB <- setdiff(this.signifB, IGNORE_GENES)
    
    # fill out pairwise.results
    n.ovlap <- length(intersect(this.signifA, this.signifB))
    n.uniqA <- length(setdiff(this.signifA, this.signifB))
    n.uniqB <- length(setdiff(this.signifB, this.signifA))
    pairwise.results[i, paste0("ovlap", this.FC)] <- n.ovlap
    pairwise.results[i, paste0("uniqA", this.FC)] <- n.uniqA
    pairwise.results[i, paste0("uniqB", this.FC)] <- n.uniqB

    # Cohen's kappa stats
    n.both_out <- length(fcsA) - (n.ovlap+n.uniqA+n.uniqB)
    n.both_out_tpmFilt <- max(0, sum(tpmA >= tpm.cutoff & tpmB >= tpm.cutoff) - (n.ovlap+n.uniqA+n.uniqB))            
    confusion_mat <- matrix(c(n.ovlap, n.uniqA, n.uniqB, n.both_out), 2, 2)
    confusion_mat_tpmFilt <- matrix(c(n.ovlap, n.uniqA, n.uniqB, n.both_out_tpmFilt), 2, 2)
    cohen_test <- Kappa(confusion_mat)
    pairwise.results[i, paste0("cohenK", this.FC)] <- cohen_test$Unweighted['value']
    pairwise.results[i, paste0("cohenK_limTPM", this.FC)] <- Kappa(confusion_mat_tpmFilt)$Unweighted['value']

	if (FALSE) { # quite slow, gives highly similar results to Cohen's K
		# ICC stats
		a <- rep(0, nrow(fmA[[1]]))
		b <- rep(0, nrow(fmA[[1]]))
		ovlap.max <- pairwise.results[i, paste0("ovlap", this.FC)]
		A.max <- ovlap.max + pairwise.results[i, paste0("uniqA", this.FC)]
		B.min <- A.max + 1 
		B.max <- A.max + pairwise.results[i, paste0("uniqB", this.FC)]
		a[1:A.max] <- 1
		b[1:ovlap.max] <- 1
		b[B.min:B.max] <- 1
		icc_out <- icc(data.frame(a,b), model = "oneway", type = "agreement", unit = "single")
		pairwise.results[i, paste0("ICC", "val", this.FC)] <- icc_out$value
	}
  } 
}
saveRDS(pairwise.results, file = 'result_objs/Dchs1_pairwise.results.RDS')
