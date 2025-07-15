# library(dplyr)

# INIT
NUM_TRIALS <- 25
N_RANGE <- 14:3
IGNORE_GENES <- c("Dchs1", "LacZ", "SV40")
PTYPE <- "Padj"
PADJ_CUTOFF <- 0.05 
TISSUES_IN_ORDER <- c('heart', 'kidney', 'liver', 'lung')
Dchs1.N30 <- readRDS("full_signature_refs/Dchs1_res.RDS")

# set up data frame for storing overlap-with-gold-standard stats
PTYPE <- "Padj"
WT.results <- data.frame(expand.grid(tissue = TISSUES_IN_ORDER, 
                                         tpm.cutoff = c(0, 1, 5),
                                         trial = 1:NUM_TRIALS, N = N_RANGE,
                                         FCtype = c("default", "ashr", "apeglm"),
                                         uniq1=NA, uniq1.2=NA, uniq1.5=NA, uniq2=NA,
										  wald.cor=NA, log2fc.cor=NA), stringsAsFactors = FALSE)
head(WT.results)


prev.DE.file <- ""
for (i in 1:nrow(WT.results)) {
  tissue <- WT.results[i, "tissue"]
  tpm.cutoff <- WT.results[i, 'tpm.cutoff']
  num.sample <- WT.results[i, "N"]
  FCtype <- WT.results[i, "FCtype"]

  if (i %% 500 == 1) { message(paste(i, tissue, date())) }

  # read in DE results (only once per tissue)
  DEres.file <- paste0("DE_runs_WTcompare/Dchs1_N", num.sample, "_trial", WT.results[i,"trial"], ".RDS")
  if (DEres.file != prev.DE.file) {
    DEres <- readRDS(DEres.file)
  }
  prev.DE.file <- DEres.file

  # decide which fold change to use
  if (FCtype == "default") {
      fc_str <- 'Fold.Change'            
  } else if (FCtype == "ashr") {
      fc_str <- 'FC.ashr'
  } else if (FCtype == 'apeglm') {
      fc_str <- 'FC.apeglm'
  }
  fcs <- DEres[[fc_str]][,tissue]
	
  #shorthand
  pvals <- DEres[[PTYPE]][,tissue]
  tpm <- DEres$meanY.TPM[,tissue]
  downreg <- !is.na(fcs) & fcs < 1
  tpm[downreg] <- DEres$meanX.TPM[downreg, tissue]
  
	
  # determine which genes are significant (prior to fold change filtering)
  tpm.pass <- tpm >= tpm.cutoff
  pass <- !is.na(fcs) & !is.na(pvals) 
  pass <- pass & pvals < PADJ_CUTOFF
  pass <- pass & tpm.pass 
  tmp2 <- pass & (fcs >= 2 | fcs <= 1/2) 
  tmp1.5 <- pass & (fcs >= 1.5 | fcs <= 1/1.5)
  tmp1.2 <- pass & (fcs >= 1.2 | fcs <= 1/1.2)
  tmp1 <- pass
  
  (this.signif.2 <- rownames(DEres[[1]])[tmp2])
  (this.signif.1.5 <- rownames(DEres[[1]])[tmp1.5])
  (this.signif.1.2 <- rownames(DEres[[1]])[tmp1.2])
  (this.signif.1 <- rownames(DEres[[1]])[tmp1])
  
  # remove genes we don't care about e.g. LacZ
  (this.signif.2 <- setdiff(this.signif.2, IGNORE_GENES))
  (this.signif.1.5 <- setdiff(this.signif.1.5, IGNORE_GENES))
  (this.signif.1.2 <- setdiff(this.signif.1.2, IGNORE_GENES))  
  (this.signif.1 <- setdiff(this.signif.1, IGNORE_GENES))

  # fill out WT.results
  WT.results[i,"uniq2"] <- length(this.signif.2)
  WT.results[i,"uniq1.5"] <- length(this.signif.1.5)
  WT.results[i,"uniq1.2"] <- length(this.signif.1.2)
  WT.results[i,"uniq1"] <- length(this.signif.1)

  # correlation between log2 fold changes
  WT.results[i, 'log2fc.cor'] <- cor(log2(fcs[tpm.pass]), log2(Dchs1.N30[[fc_str]][tpm.pass, tissue]), 
                                      use = 'complete.obs', method = 'spearman')

  WT.results[i, 'wald.cor'] <- cor(DEres$Wald[tpm.pass, tissue], Dchs1.N30$Wald[tpm.pass, tissue], 
                                      use = 'complete.obs', method = 'spearman')
}
saveRDS(WT.results, file = paste0('result_objs/Dchs1_WTcompare.results.RDS'))

