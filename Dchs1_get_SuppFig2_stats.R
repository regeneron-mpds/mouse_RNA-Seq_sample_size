# INIT
# NUM_TRIALS <- 40
NUM_TRIALS <- 12
# N_RANGE <- 3:28
N_RANGE <- 3
IGNORE_GENES <- c("Dchs1", "LacZ", "SV40")
PTYPE <- "Padj"
PADJ_CUTOFF <- 0.05 
TPM_CUTOFF <- 0
TISSUES_IN_ORDER <- c('heart', 'kidney', 'liver', 'lung')
gold_standards <- readRDS('gold_standards.RDS')
goldstandard.FC <- 1.5
MIN_NUM_SIGNIF <- 10

prob_correct <- data.frame(expand.grid(tissue = TISSUES_IN_ORDER,
                                       N = N_RANGE, trial = 1:NUM_TRIALS,
                                       num.signif = NA,
                                       fc1.2 = NA, fc1.5 = NA, fc2 = NA, fc3 = NA,
                                       fc5 = NA, fc10 = NA, fcInf = NA, 
                                       stringsAsFactors = FALSE))
prob_detect = prob_detect_cumul = prob_correct

for (this.type in c("", "_uponly")) {
  prev.DEres.file <- ""
  for (i in 1:nrow(prob_detect)) { 
    tissue <- prob_detect[i, "tissue"]

    if (i %% 500 == 1) { message(paste(i, tissue, date())) }
    
    # read in formattedInfo (only once per tissue)
    this.DEres.file <- paste0("DE_runs_subsample/Dchs1_N", prob_detect[i,"N"],
                           "_trial", prob_detect[i, "trial"], ".RDS")
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
    # converting down-regulated fold changes to up-regulated, only care about magnitude
    if (this.type != "_uponly") {
      this.fcs[downreg] <- 1/this.fcs[downreg]
    }
    
    # determine which genes are significant 
    tmp.pass <- !is.na(this.fcs) & !is.na(this.pvals) 
    tmp.pass <- tmp.pass & this.tpm >= TPM_CUTOFF
    tmp.pass <- tmp.pass & this.pvals <= PADJ_CUTOFF
    
    N30.genes <- gold_standards$dchs1[[paste(tissue, TPM_CUTOFF, PTYPE, PADJ_CUTOFF, goldstandard.FC, sep=',')]]
    tmp <- tmp.pass & (this.fcs >= goldstandard.FC) 
    this.signif <- rownames(this.DEres[[1]])[tmp]
    # remove genes we don't care about e.g. LacZ (already done for gold_standards)
    this.signif <- setdiff(this.signif, IGNORE_GENES)

    # skip if we have too few genes to derive meaningful stats    
    if (length(this.signif) < MIN_NUM_SIGNIF) { next }

    # record num significant genes
    prob_correct[i, "num.signif"] <- length(this.signif)
    prob_detect[i, "num.signif"] <- length(this.signif)
    prob_detect_cumul[i, "num.signif"] <- length(this.signif)
    
    gaps <- as.numeric(c('0', "1.2", "1.5", "2", "3", "5", "10", "Inf"))
    for (k in 2:length(gaps)) {
      lower <- gaps[k-1]
      upper <- gaps[k]
      indchr <- paste0('fc', upper)

      # probability of a called gene within a given FC range belonging to GS
      x <- names(this.fcs)[tmp.pass & this.fcs > lower & this.fcs <= upper]
      prob_correct[i, indchr] <- sum(x %in% N30.genes) / length(x)
            
      # probability of detecting GS for a given FC range (non-cumulative)
      x <- names(this.fcs)[tmp.pass & this.fcs > lower & this.fcs <= upper]
      prob_detect[i, indchr] <- sum(N30.genes %in% x) / length(N30.genes)
      
      # probability of detecting GS FC for a given FC range (cumulative)
      x <- names(this.fcs)[tmp.pass & this.fcs > lower]
      prob_detect_cumul[i, indchr] <- sum(N30.genes %in% x) / length(N30.genes)
    }
    
    # x <- names(this.fcs)[tmp.pass & this.fcs <= 1.2]
    # FCranges[i, 'fc1.2'] <- sum(x %in% N30.genes) / length(x)
    # x <- names(this.fcs)[tmp.pass & this.fcs > 1.2 & this.fcs <= 1.5]
    # FCranges[i, 'fc1.5'] <- sum(x %in% N30.genes) / length(x)
    # x <- names(this.fcs)[tmp.pass & this.fcs > 1.5 & this.fcs <= 2]
    # FCranges[i, 'fc2'] <- sum(x %in% N30.genes) / length(x)
    # x <- names(this.fcs)[tmp.pass & this.fcs > 2 & this.fcs <= 3]
    # FCranges[i, 'fc3'] <- sum(x %in% N30.genes) / length(x)
    # x <- names(this.fcs)[tmp.pass & this.fcs > 3 & this.fcs <= 5]
    # FCranges[i, 'fc5'] <- sum(x %in% N30.genes) / length(x)            
    # x <- names(this.fcs)[tmp.pass & this.fcs > 5 & this.fcs <= 10]
    # FCranges[i, 'fc10'] <- sum(x %in% N30.genes) / length(x)
    # x <- names(this.fcs)[tmp.pass & this.fcs > 10]
    # FCranges[i, 'fcInf'] <- sum(x %in% N30.genes) / length(x)
  } 
  saveRDS(prob_correct, paste0('Dchs1_prob_correct', this.type, '_GSFC', goldstandard.FC, '_minSignif', MIN_NUM_SIGNIF, '.RDS'))
  saveRDS(prob_detect, paste0('Dchs1_prob_detect', this.type, '_GSFC', goldstandard.FC, '_minSignif', MIN_NUM_SIGNIF, '.RDS'))
  saveRDS(prob_detect_cumul, paste0('Dchs1_prob_detect_cumul', this.type, '_GSFC', goldstandard.FC, '_minSignif', MIN_NUM_SIGNIF, '.RDS'))
}
