library(dplyr)
library(parallel)

# INIT
N_RANGE <- 3:28


run_5_piece <- function(this.N) {
	NUM_TRIALS <- 40
	IGNORE_GENES <- c("Dchs1", "LacZ", "SV40")
	PTYPE <- "Padj"
	PADJ_CUTOFF <- 0.05 
	MIN_NUM_SIGNIF <- 10
	gold_standards <- readRDS('result_objs/gold_standards.RDS')


	prob_correct <- data.frame(expand.grid(tissue = c('heart', 'kidney', 'liver', 'lung'),
								   trial = 1:NUM_TRIALS, N = this.N, 
								   FCtype = c('default', 'apeglm', 'ashr'),
								   goldstandard.FC = c(1, 1.2, 1.5, 2),
									tpm.cutoff = c(0,1,5),
								   num.signif = NA,
								   fc1.2 = NA, fc1.5 = NA, fc2 = NA, fc3 = NA,
								   fc5 = NA, fc10 = NA, fcInf = NA, 
								   num.in.range1.2 = NA, num.in.range1.5 = NA, num.in.range2 = NA, num.in.range3 = NA,
								   num.in.range5 = NA, num.in.range10 = NA, num.in.rangeInf = NA, 
								   stringsAsFactors = FALSE))
	prob_detect = prob_detect_cumul = prob_correct %>% select(!contains("num.in.range"))

	for (this.type in c("")) { # , "_uponly")) {
	  prev.DEres.file <- ""
	  for (i in 1:nrow(prob_detect)) { 
		tissue <- prob_detect[i, "tissue"]
		this.FCtype <- prob_detect[i, "FCtype"]
		goldstandard.FC <- prob_detect[i, "goldstandard.FC"]
		tpm.cutoff <- prob_detect[i, "tpm.cutoff"]
		if (i %% 500 == 1) { message(paste(i, tissue, date())) }

		# read in formattedInfo (only once per tissue)
		this.DEres.file <- paste0("DE_runs_subsample/Dchs1_N", this.N,
							   "_trial", prob_detect[i, "trial"], ".RDS")
		if (this.DEres.file != prev.DEres.file) {
		  this.DEres <- readRDS(this.DEres.file)
		}
		prev.DEres.file <- this.DEres.file

	  # decide which fold change to use
	  if (this.FCtype == "default") { 
		  this.fcs <- this.DEres$Fold.Change[,tissue]
	  } else if (this.FCtype == "ashr") {
		   this.fcs <- this.DEres$FC.ashr[,tissue]
	  } else if (this.FCtype == 'apeglm') {
		   this.fcs <- this.DEres$FC.apeglm[,tissue]
	  }

		#shorthand
		this.pvals <- this.DEres[[PTYPE]][,tissue]
		this.tpm <- this.DEres$meanY.TPM[,tissue]
		downreg <- !is.na(this.fcs) & this.fcs < 1
		this.tpm[downreg] <- this.DEres$meanX.TPM[downreg, tissue]
		# converting down-regulated fold changes to up-regulated, only care about magnitude
		if (this.type != "_uponly") {
		  this.fcs[downreg] <- 1/this.fcs[downreg]
		}

		# determine which genes are significant 
		tmp.pass <- !is.na(this.fcs) & !is.na(this.pvals) 
		tmp.pass <- tmp.pass & this.tpm >= tpm.cutoff
		tmp.pass <- tmp.pass & this.pvals <= PADJ_CUTOFF

		str <- paste(tissue, tpm.cutoff, PTYPE, 'std.05', this.FCtype, PADJ_CUTOFF, goldstandard.FC, sep=',')
		N30.genes <- gold_standards$Dchs1[[str]]
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
		  indchr2 <- paste0('num.in.range', upper)

		  # probability of a called gene within a given FC range belonging to GS
		  x <- names(this.fcs)[tmp.pass & this.fcs > lower & this.fcs <= upper]
		  prob_correct[i, indchr] <- sum(x %in% N30.genes) / length(x)
  		  prob_correct[i, indchr2] <- length(x)

		  # probability of detecting GS for a given FC range (non-cumulative)
		  prob_detect[i, indchr] <- sum(N30.genes %in% x) / length(N30.genes)

		  # probability of detecting GS FC for a given FC range (cumulative)
		  x <- names(this.fcs)[tmp.pass & this.fcs > lower]
		  prob_detect_cumul[i, indchr] <- sum(N30.genes %in% x) / length(N30.genes)
		} # k
	  } # i
		suffix <- paste0(this.type, '_minSignif', MIN_NUM_SIGNIF, "_N", this.N, '.RDS')
		saveRDS(prob_correct, paste0('result_objs/Dchs1_prob_correct', suffix))
		saveRDS(prob_detect, paste0('result_objs/Dchs1_prob_detect', suffix))
		saveRDS(prob_detect_cumul, paste0('result_objs/Dchs1_prob_detect_cumul', suffix))
	} #this.type
}  # 

# Run the function in parallel using mclapply
mclapply(N_RANGE, run_5_piece, mc.cores = round(detectCores() / 3))


# concatenate results back together
MIN_NUM_SIGNIF <- 10
for (N in N_RANGE) {
	suffix <- paste0('_minSignif', MIN_NUM_SIGNIF, "_N", N, '.RDS')
	pc <- readRDS(paste0('result_objs/Dchs1_prob_correct', suffix))
	pd <- readRDS(paste0('result_objs/Dchs1_prob_detect', suffix))
	pdc <- readRDS(paste0('result_objs/Dchs1_prob_detect_cumul', suffix))

	if (N == min(as.numeric(N_RANGE))) {
		prob_correct <- pc
		prob_detect <- pd
		prob_detect_cumul <- pdc
	} else {
		prob_correct <- rbind(prob_correct, pc)
		prob_detect <- rbind(prob_detect, pd)
		prob_detect_cumul <- rbind(prob_detect_cumul, pdc)
	}
}

saveRDS(prob_correct, file = paste0('result_objs/Dchs1_prob_correct_minSignif', MIN_NUM_SIGNIF, '.RDS'))
saveRDS(prob_detect, file = paste0('result_objs/Dchs1_prob_detect_minSignif', MIN_NUM_SIGNIF, '.RDS'))
saveRDS(prob_detect_cumul, file = paste0('result_objs/Dchs1_prob_detect_cumul_minSignif', MIN_NUM_SIGNIF, '.RDS'))
