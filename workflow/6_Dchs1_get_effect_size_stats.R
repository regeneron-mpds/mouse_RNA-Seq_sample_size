# INIT
NUM_TRIALS <- 40
N_RANGE <- 3:28
IGNORE_GENES <- c("Dchs1", "LacZ", "SV40")
PTYPE <- "Padj"
PADJ_CUTOFF <- 0.05 
TISSUES_IN_ORDER <- c('heart', 'kidney', 'liver', 'lung')
gold_standards <- readRDS('result_objs/gold_standards.RDS')
N30_DEres <- readRDS("full_signature_refs/Dchs1_res.RDS")

effect.size.results <- data.frame(expand.grid(tissue = TISSUES_IN_ORDER, 
                                              tpm.cutoff = c(0, 1, 5),
                                              N = N_RANGE, trial = 1:NUM_TRIALS,
                                              signif1 = NA, signif1.2 = NA, signif1.5 = NA, signif2 = NA,
                                              ovlap1 = NA, ovlap1.2 = NA, ovlap1.5 = NA, ovlap2 = NA,
                                              gt.orig1 = NA, gt.orig1.2 = NA, gt.orig1.5 = NA, gt.orig2 = NA,
                                              lt.orig1 = NA, lt.orig1.2 = NA, lt.orig1.5 = NA, lt.orig2 = NA,
                                              stringsAsFactors = FALSE))

# just doing this for P-value threshold of 0.05
prev.DEres.file <- ""
for (i in 1:nrow(effect.size.results)) { 
  tissue <- effect.size.results[i, "tissue"]
  tpm.cutoff <- effect.size.results[i, 'tpm.cutoff']
  if (i %% 500 == 1) { message(paste(i, tissue, date())) }
  
  # read in DE results (only once per tissue)
  this.DEres.file <- paste0("DE_runs_subsample/Dchs1_N", effect.size.results[i,"N"],
                         "_trial", effect.size.results[i, "trial"], ".RDS")
  if (this.DEres.file != prev.DEres.file) {
    this.DEres <- readRDS(this.DEres.file)
  }
  prev.DEres.file <- this.DEres.file
  
  #shorthand
  this.pvals <- this.DEres[[PTYPE]][,tissue]
  this.fcs <- this.DEres$Fold.Change[,tissue]
  this.tpm <- this.DEres$meanY.TPM[,tissue]
  downreg <- !is.na(this.fcs) & this.fcs < 1
  this.tpm[downreg] <- this.DEres$meanX.TPM[downreg, tissue]
  
  # determine which genes are significant (prior to fold change filtering)
  tmp.pass <- !is.na(this.fcs) & !is.na(this.pvals) 
  tmp.pass <- tmp.pass & this.pvals <=  PADJ_CUTOFF
  tmp.pass <- tmp.pass & this.tpm >= tpm.cutoff
  
  for (this.FC in c(1, 1.2, 1.5, 2)) { 
    N30.genes <- gold_standards$Dchs1[[paste(tissue, tpm.cutoff, PTYPE, PADJ_CUTOFF, this.FC, sep=',')]]

    tmp <- tmp.pass & (this.fcs >= this.FC | this.fcs <= 1/this.FC)
    this.signif <- rownames(this.DEres[[1]])[tmp]
    # remove genes we don't care about e.g. LacZ (already done for gold_standards)
    this.signif <- setdiff(this.signif, IGNORE_GENES)
    
    # fill out effect.size.results
    effect.size.results[i, paste0("signif", this.FC)] <- length(this.signif)    
    effect.size.results[i, paste0("ovlap", this.FC)] <- length(intersect(this.signif, N30.genes))                                        

    # how many of the ovlaps have greater/lesser fold changes than the N30 version?
    abs.log.N30.fcs <- abs(log2(N30_DEres$Fold.Change[this.signif, tissue]))
    abs.log.this.fcs <- abs(log2(this.fcs[this.signif]))
    
    effect.size.results[i, paste0("gt.orig", this.FC)] <- sum(abs.log.this.fcs > abs.log.N30.fcs)
    effect.size.results[i, paste0("lt.orig", this.FC)] <- sum(abs.log.this.fcs < abs.log.N30.fcs)
    stopifnot(effect.size.results[i, paste0("gt.orig", this.FC)] + 
              effect.size.results[i, paste0("lt.orig", this.FC)] == length(this.signif)) # sanity check: no genes have exactly equal FCs
  } 
}
saveRDS(effect.size.results, file='result_objs/Dchs1_effect.size.results.RDS')