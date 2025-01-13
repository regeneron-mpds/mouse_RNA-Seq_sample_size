# library(dplyr)

# INIT
NUM_TRIALS <- 20
N_RANGE <- 3
IGNORE_GENES <- c("Dchs1", "LacZ", "SV40")
PTYPE <- "Padj"
TISSUES_IN_ORDER <- c('heart', 'kidney', 'liver', 'lung')


# read in gold standard (N30) signature info
dchs1.N30 <- readRDS('full_signature_refs/dchs1_res.RDS')
gold_standards <- readRDS('gold_standards.RDS')


# set up data frame for storing overlap-with-gold-standard stats
ovlap.results <- data.frame(expand.grid(tissue = TISSUES_IN_ORDER, 
                                        pval.cutoff = c(0.01, 0.05, 0.1, 1),
                                        tpm.cutoff = c(0, 1, 5),
                                        N = N_RANGE, trial = 1:NUM_TRIALS,
                                        signif1 = NA, signif1.2 = NA, signif1.5 = NA, signif2 = NA,
                                        ovlap1 = NA, ovlap1.2 = NA, ovlap1.5 = NA, ovlap2 = NA,
                                        detection1 = NA, detection1.2 = NA, detection1.5 = NA, detection2 = NA,
                                        all.cor = NA, stringsAsFactors = FALSE))
head(ovlap.results)


## obtain ovlap.results stats from DE_runs_subsample/Dchs1*.RDS files
prev.DE.file <- ""
for (i in 1:nrow(ovlap.results)) { 
  tissue <- ovlap.results[i, "tissue"]
  tpm.cutoff <- ovlap.results[i, 'tpm.cutoff']
  pval.cutoff <- ovlap.results[i, 'pval.cutoff']
  # if (i %% 500 == 1) { message(paste(i, tissue, date())) }
  
  # read in DE results (only once per tissue)
  this.DE.file <- paste0("DE_runs_subsample/Dchs1_N", ovlap.results[i,"N"],
                         "_trial", ovlap.results[i, "trial"], ".RDS")
  if (this.DE.file != prev.DE.file) {
    this.DE <- readRDS(this.DE.file)
  }
  prev.DE.file <- this.DE.file
  
  #shorthand
  this.pvals <- this.DE[[PTYPE]][,tissue]
  this.fcs <- this.DE$Fold.Change[,tissue]
  this.tpm <- this.DE$meanY.TPM[,tissue]
  downreg <- !is.na(this.fcs) & this.fcs < 1
  this.tpm[downreg] <- this.DE$meanX.TPM[downreg, tissue]
  
  # determine which genes are significant (prior to fold change filtering)
  tmp.pass <- !is.na(this.fcs) & !is.na(this.pvals) 
  tmp.pass <- tmp.pass & this.tpm >= tpm.cutoff 
  tmp.pass <- tmp.pass & this.pvals <= pval.cutoff
  
  # correlation between log2 fold changes
  ovlap.results[i, 'all.cor'] <- cor(log2(this.fcs), log2(dchs1.N30$Fold.Change[, tissue]), 
                                     use = 'complete.obs', method = 'spearman')

  # derive num significant, number overlapping, etc.
  for (this.FC in c(1, 1.2, 1.5, 2)) { 
    N30.genes <- gold_standards$dchs1[[paste(tissue, tpm.cutoff, PTYPE, pval.cutoff, this.FC, sep=',')]]

    tmp <- tmp.pass & (this.fcs >= this.FC | this.fcs <= 1/this.FC) # fold change check
    (this.signif <- rownames(this.DE[[1]])[tmp])
    # remove genes we don't care about e.g. LacZ (already done for gold_standards)
    (this.signif <- setdiff(this.signif, IGNORE_GENES))  
    
    # fill out ovlap.results
    ovlap.results[i, paste0("signif", this.FC)] <- length(this.signif)
    ovlap.results[i, paste0("ovlap", this.FC)] <- length(intersect(this.signif, N30.genes))
    ovlap.results[i, paste0("detection", this.FC)] <- round(100 * ovlap.results[i, paste0("ovlap", this.FC)] / length(N30.genes), 2)
  }
} 
saveRDS(ovlap.results, 'Dchs1_ovlap.results.RDS')
