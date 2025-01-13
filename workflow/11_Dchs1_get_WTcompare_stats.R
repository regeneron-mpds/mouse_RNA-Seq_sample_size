val# library(dplyr)

# INIT
NUM_TRIALS <- 25
N_RANGE <- 14:3
IGNORE_GENES <- c("Dchs1", "LacZ", "SV40")
PTYPE <- "Padj"
PADJ_CUTOFF <- 0.05 
TISSUES_IN_ORDER <- c('heart', 'kidney', 'liver', 'lung')

# set up data frame for storing overlap-with-gold-standard stats
PTYPE <- "Padj"
pairwise.results <- data.frame(expand.grid(tissue = TISSUES_IN_ORDER, 
                                         tpm.cutoff = c(0, 1, 5),
                                         trial = 1:NUM_TRIALS, N = N_RANGE,
                                         uniq1=NA, uniq1.2=NA, uniq1.5=NA, uniq2=NA),
                             stringsAsFactors = FALSE)
head(pairwise.results)


prev.DE.file <- ""
for (i in 1:nrow(pairwise.results)) {
  tissue <- pairwise.results[i, "tissue"]
  tpm.cutoff <- pairwise.results[i, 'tpm.cutoff']
  num.sample <- pairwise.results[i, "N"]
  if (i %% 500 == 1) { message(paste(i, tissue, date())) }

  # read in DE results (only once per tissue)
  DEres.file <- paste0("DE_runs_WTcompare/Dchs1_N", num.sample, "_trial", pairwise.results[i,"trial"], ".RDS")
  if (DEres.file != prev.DE.file) {
    DEres <- readRDS(DEres.file)
  }
  prev.DE.file <- DEres.file
  
  #shorthand
  pvals <- DEres[[PTYPE]][,tissue]
  fcs <- DEres$Fold.Change[,tissue]
  tpm <- DEres$meanY.TPM[,tissue]
  downreg <- !is.na(fcs) & fcs < 1
  tpm[downreg] <- DEres$meanX.TPM[downreg, tissue]
  
  # determine which genes are significant (prior to fold change filtering)
  pass <- !is.na(fcs) & !is.na(pvals) 
  pass <- pass & pvals < PADJ_CUTOFF
  pass <- pass & tpm >= tpm.cutoff 
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

  # fill out pairwise.results
  pairwise.results[i,"uniq2"] <- length(this.signif.2)
  pairwise.results[i,"uniq1.5"] <- length(this.signif.1.5)
  pairwise.results[i,"uniq1.2"] <- length(this.signif.1.2)
  pairwise.results[i,"uniq1"] <- length(this.signif.1)
} 
saveRDS(pairwise.results, file = paste0('result_objs/Dchs1_WTcompare.results_', PTYPE, PADJ_CUTOFF, '.RDS'))