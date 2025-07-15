library(dplyr)

# INIT
NUM_TRIALS <- 40
N_RANGE <- 3:28
PTYPE <- "Padj"
PADJ_CUTOFF <- 0.05
TPM_CUTOFF <- 0
TISSUES_IN_ORDER <- c('heart', 'kidney', 'liver', 'lung')
gold_standards <- readRDS("result_objs/gold_standards.RDS")$Dchs1

### SELECT SUBSET OF GENES
if (USE_ORIG_GENES_FOR_S4) {
	genes.for.effect.size.exploration <- readRDS('result_objs/Dchs1_genes.for.effect.size.exploration_ORIG_SUBMISSION.RDS')
	suffix <- "_ORIG_SUBMISSION"
} else {
	# Randomly select 20 genes per tissue, grouped into four fold change ranges
	genes.for.effect.size.exploration <- list()
	num.per.group <- 5  # 5 genes for each combo of tissue-by-fold.change.regime
	set.seed(17)
	for (this.tissue in TISSUES_IN_ORDER) {
	  str <- paste(this.tissue, TPM_CUTOFF, PTYPE, 'std.05', 'default', PADJ_CUTOFF, sep=',')
	  stopifnot(paste0(str, ',1.5') %in% names(gold_standards))

	  set2 <-   setdiff(unlist(gold_standards[[paste0(str, ",2")]]), "")  # genes with FC >= 2 or FC < 0.5
	  set1.5 <- setdiff(unlist(gold_standards[[paste0(str, ",1.5")]]), set2) # genes with (1.5 < FC < 2) OR (1/2 < FC < 2/3)
	  set1.2 <- setdiff(unlist(gold_standards[[paste0(str, ",1.2")]]), c(set2, set1.5)) # genes with (1.2 < FC < 1.5) OR (2/3 < FC < 5/6)
	  set1 <-   setdiff(unlist(gold_standards[[paste0(str, ",1")]]), c(set2, set1.5, set1.2)) # genes with (5/6 < FC < 1.2)
	  genes.for.effect.size.exploration[[this.tissue]] <- c(
		sample(set1, num.per.group, replace = FALSE),
		sample(set1.2, num.per.group, replace = FALSE),
		sample(set1.5, num.per.group, replace = FALSE),
		sample(set2, num.per.group, replace = FALSE))
	}
	saveRDS(genes.for.effect.size.exploration, 'result_objs/Dchs1_genes.for.effect.size.exploration.RDS')
	suffix <- ""
}

### GET EFFECT SIZE INFO FOR THESE GENES
uniq.explore.genes <- as.character(unlist(genes.for.effect.size.exploration))
genewise.effect.size.res <- data.frame(expand.grid(gene = uniq.explore.genes,
                                                   tissue = TISSUES_IN_ORDER,
                                                   tpm.cutoff = c(0, 1, 5),
                                                   N = N_RANGE, trial = 1:NUM_TRIALS,
                                                   FC = NA, FC.apeglm = NA, FC.ashr = NA,
												   is.signif = NA, meanX.TPM = NA,
                                                   stringsAsFactors = FALSE))
prev.DEres.file <- ""
block_size <- length(uniq.explore.genes) * length(TISSUES_IN_ORDER)

for (i in seq(1, nrow(genewise.effect.size.res), by = block_size)) { 
  tpm.cutoff <- genewise.effect.size.res[i, 'tpm.cutoff']
  this.N <- genewise.effect.size.res[i,"N"]
  this.trial <- genewise.effect.size.res[i,"trial"]
  this_range <- i:(i+block_size-1)
  stopifnot(all(genewise.effect.size.res[this_range, 'tpm.cutoff'] == tpm.cutoff))
  stopifnot(all(genewise.effect.size.res[this_range, 'trial'] == this.trial))
  stopifnot(all(genewise.effect.size.res[this_range, 'N'] == this.N))
  if (i %% 10000 < 50) { message(paste(i, date())) }
  
  # read in DE results
  this.DEres.file <- paste0("DE_runs_subsample//Dchs1_N", genewise.effect.size.res[i,"N"],
                            "_trial", genewise.effect.size.res[i, "trial"], ".RDS")
  if (this.DEres.file != prev.DEres.file) {
    this.DEres <- readRDS(this.DEres.file)
  }
  prev.DEres.file <- this.DEres.file
  
  for (j in seq(i, max(this_range), by = length(uniq.explore.genes))) {
    tissue <- genewise.effect.size.res[j, "tissue"]
    this_tissue_range <- j:(j+length(uniq.explore.genes)-1)
    stopifnot(genewise.effect.size.res[this_tissue_range, 'gene'] == uniq.explore.genes)
    stopifnot(all(genewise.effect.size.res[this_tissue_range, 'tissue'] == tissue))
    
    #shorthand
    this.pvals <- this.DEres[[PTYPE]][uniq.explore.genes, tissue]
	# this.tpms <- pmax(this.DEres$meanX.TPM[uniq.explore.genes, tissue], 
	# 				  this.DEres$meanY.TPM[uniq.explore.genes, tissue])
    this.tpms <- this.DEres$meanY.TPM[uniq.explore.genes, tissue]
    alt.tpm <- this.DEres$meanX.TPM[uniq.explore.genes, tissue]
    tmp.pass <- !is.na(this.pvals) & this.pvals < PADJ_CUTOFF

	# default
	this.fcs <- this.DEres$Fold.Change[uniq.explore.genes, tissue]
    downreg <- !is.na(this.fcs) & this.fcs < 1
    this.tpms[downreg] <- alt.tpm[downreg]
	this.pass <- tmp.pass & this.tpms >= tpm.cutoff
    this.pass <- this.pass & !is.na(this.fcs)
    genewise.effect.size.res[this_tissue_range, 'FC'] <- this.fcs
    genewise.effect.size.res[this_tissue_range, 'is.signif'] <- this.pass
    genewise.effect.size.res[this_tissue_range, 'meanX.TPM'] <- this.tpms
	  
	# apeglm
	this.fcs <- this.DEres$FC.apeglm[uniq.explore.genes, tissue]
    genewise.effect.size.res[this_tissue_range, 'FC.apeglm'] <- this.fcs

  	# ashr
	this.fcs <- this.DEres$FC.ashr[uniq.explore.genes, tissue]
    genewise.effect.size.res[this_tissue_range, 'FC.ashr'] <- this.fcs
  } # j
} # i
head(genewise.effect.size.res)
table(genewise.effect.size.res$is.signif, useNA = 'always')
saveRDS(genewise.effect.size.res, paste0('result_objs/Dchs1_genewise.effect.size.res', suffix, '.RDS'))
