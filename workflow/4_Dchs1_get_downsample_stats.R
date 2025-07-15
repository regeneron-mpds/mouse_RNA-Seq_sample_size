library(parallel)


N_RANGE <- 3:28
# TPM.CUTOFFs <- c(0,1,5)

RUN_INFO <- expand.grid(N = N_RANGE) # , tpm.cutoff = TPM.CUTOFFs)
TOTAL_RUNS <- nrow(RUN_INFO)

run_4_piece <- function(this.N) { # tpm.cutoff, 

# INIT
NUM_TRIALS <- 40
TISSUES_IN_ORDER <- c('heart', 'kidney', 'liver', 'lung')
IGNORE_GENES <- c("Dchs1", "LacZ", "SV40")
PTYPE <- "Padj"

# read in gold standard (N30) signature info
Dchs1.N30 <- readRDS('full_signature_refs/Dchs1_res.RDS')
gold.standards <- readRDS('result_objs/gold_standards.RDS')
GS_passfail_results <- readRDS("result_objs/gold_standard_passfail_info.RDS")

# set up data frame for storing overlap-with-gold-standard stats
ovlap.results <- data.frame(expand.grid(tissue = TISSUES_IN_ORDER, 
                                        pval.cutoff = c(0.01, 0.05, 0.1, 1),
                                        tpm.cutoff = c(0, 1, 5),
                                        alphaType = c('std.05', 'match.padj.cutoff'),
                                        FCtype = c("default", "ashr", "apeglm"),
                                        N = this.N, trial = 1:NUM_TRIALS,
										# "row" vars correspond to rows of supp. table 1, commented out
										# row1_1 = NA, row1_1.2 = NA, row1_1.5 = NA, row1_2 = NA,
										# row2_1 = NA, row2_1.2 = NA, row2_1.5 = NA, row2_2 = NA,
										# row3_1 = NA, row3_1.2 = NA, row3_1.5 = NA, row3_2 = NA,
										# row4_1 = NA, row4_1.2 = NA, row4_1.5 = NA, row4_2 = NA,
										# row5_1 = NA, row5_1.2 = NA, row5_1.5 = NA, row5_2 = NA,
										# row6_1 = NA, row6_1.2 = NA, row6_1.5 = NA, row6_2 = NA,
										# row7_1 = NA, row7_1.2 = NA, row7_1.5 = NA, row7_2 = NA,
										# row8_1 = NA, row8_1.2 = NA, row8_1.5 = NA, row8_2 = NA,
										# row9_1 = NA, row9_1.2 = NA, row9_1.5 = NA, row9_2 = NA,
										# row10_1 = NA, row10_1.2 = NA, row10_1.5 = NA, row10_2 = NA,
										# row11_1 = NA, row11_1.2 = NA, row11_1.5 = NA, row11_2 = NA,
										# row12_1 = NA, row12_1.2 = NA, row12_1.5 = NA, row12_2 = NA,
										signif1 = NA, signif1.2 = NA, signif1.5 = NA, signif2 = NA,                                     
										ovlap1 = NA, ovlap1.2 = NA, ovlap1.5 = NA, ovlap2 = NA,
										detection1 = NA, detection1.2 = NA, detection1.5 = NA, detection2 = NA,
										signif_alt1 = NA, signif_alt1.2 = NA, signif_alt1.5 = NA, signif_alt2 = NA,
										ovlap_alt1 = NA, ovlap_alt1.2 = NA, ovlap_alt1.5 = NA, ovlap_alt2 = NA,
										detection_alt1 = NA, detection_alt1.2 = NA, detection_alt1.5 = NA, detection_alt2 = NA,
										detection_altb1 = NA, detection_altb1.2 = NA, detection_altb1.5 = NA, detection_altb2 = NA,
                                       wald.cor = NA, log2fc.cor = NA, stringsAsFactors = FALSE
									   ))
head(ovlap.results)

## obtain ovlap.results stats from DE_runs_subsample/Dchs1*.RDS files
prev.DE.file <- ""
for (i in 1:nrow(ovlap.results)) { 
    tissue <- ovlap.results[i, "tissue"]
    tpm.cutoff <- ovlap.results[i, 'tpm.cutoff']
    pval.cutoff <- ovlap.results[i, 'pval.cutoff']
    alphaType <- ovlap.results[i, 'alphaType']
    FCtype <- ovlap.results[i, 'FCtype']
    if (i %% 1000 == 1) { message(paste(i, tissue, date())) }

    # read in DE results (only once per tissue)
    this.DE.file <- paste0("DE_runs_subsample/Dchs1_N", ovlap.results[i,"N"],
                         "_trial", ovlap.results[i, "trial"], ".RDS")
    if (this.DE.file != prev.DE.file) {
    this.DE <- readRDS(this.DE.file)
    }
    prev.DE.file <- this.DE.file
    
    # decide which padj to use
    if (alphaType == 'std.05' || pval.cutoff == 0.05) {
        this.pvals <- this.DE[[PTYPE]][,tissue]
    } else if (alphaType == 'match.padj.cutoff' && pval.cutoff == .01) {
        this.pvals <- this.DE[['padj.01']][,tissue]
    } else if (alphaType == 'match.padj.cutoff' && pval.cutoff == .1) {
        this.pvals <- this.DE[['padj.1']][,tissue]
    } else if (alphaType == 'match.padj.cutoff' && pval.cutoff == 1) {
        this.pvals <- this.DE[['padj.999']][,tissue]
    }
    
    # decide which fold change to use
    if (FCtype == "default") { #  "ashr", "apeglm"),
        this.fcs <- this.DE$Fold.Change[,tissue]
    } else if (FCtype == "ashr") {
         this.fcs <- this.DE$FC.ashr[,tissue]
    } else if (FCtype == 'apeglm') {
         this.fcs <- this.DE$FC.apeglm[,tissue]
    }
        
    #shorthand
    this.tpm <- this.DE$meanY.TPM[,tissue]
    downreg <- !is.na(this.fcs) & this.fcs < 1
    this.tpm[downreg] <- this.DE$meanX.TPM[downreg, tissue]

	# determine which genes are significant (prior to fold change filtering)
	pval.pass <- !is.na(this.pvals) & this.pvals <= pval.cutoff
	tpm.pass <- this.tpm >= tpm.cutoff
    tpm.pass <- tpm.pass & !names(this.tpm) %in% IGNORE_GENES 

    # derive num significant, number overlapping, etc.
    for (this.FC in c(1, 1.2, 1.5, 2)) { 
		str <- paste(tissue, tpm.cutoff, PTYPE, alphaType, FCtype, pval.cutoff, this.FC, sep=',')

		N30_passfail_res <- GS_passfail_results$Dchs1[[str]]
		N30.genes <- sum(N30_passfail_res$grp1)
		stopifnot(N30.genes == length(gold.standards$Dchs1[[str]]))
        N30.genes_alt <- sum(N30_passfail_res$grp1 | N30_passfail_res$grp2)

        fc.pass <- this.fcs >= this.FC | this.fcs <= 1/this.FC # fold change check
        pos <- tpm.pass & pval.pass & fc.pass 
        low <- tpm.pass & pval.pass & !fc.pass
        neg <- tpm.pass & !pval.pass 
    
        # ovlap.results[i, paste0("row1_", this.FC)] <- sum(N30_passfail_res$grp1 & pos)
        # ovlap.results[i, paste0("row2_", this.FC)] <- sum(N30_passfail_res$grp1 & low)
        # ovlap.results[i, paste0("row3_", this.FC)] <- sum(N30_passfail_res$grp1 & neg)
        # ovlap.results[i, paste0("row4_", this.FC)] <- sum(N30_passfail_res$grp2 & pos)
        # ovlap.results[i, paste0("row5_", this.FC)] <- sum(N30_passfail_res$grp2 & low)
        # ovlap.results[i, paste0("row6_", this.FC)] <- sum(N30_passfail_res$grp2 & neg)
        # ovlap.results[i, paste0("row7_", this.FC)] <- sum(N30_passfail_res$grp3 & pos)
        # ovlap.results[i, paste0("row8_", this.FC)] <- sum(N30_passfail_res$grp3 & low)
        # ovlap.results[i, paste0("row9_", this.FC)] <- sum(N30_passfail_res$grp3 & neg)
        # ovlap.results[i, paste0("row10_", this.FC)] <- sum(N30_passfail_res$grp4 & pos)
        # ovlap.results[i, paste0("row11_", this.FC)] <- sum(N30_passfail_res$grp4 & low)
        # ovlap.results[i, paste0("row12_", this.FC)] <- sum(N30_passfail_res$grp4 & neg)
            
        ovlap.results[i, paste0("signif", this.FC)] <- sum(pos)
        ovlap.results[i, paste0("signif_alt", this.FC)] <- sum(pos | low)
        
        num.ovlap <- sum(N30_passfail_res$grp1 & pos)
        num.ovlap_alt <- sum(N30_passfail_res$grp1 & (pos | low)) +
                          sum(N30_passfail_res$grp2 & (pos | low))
        ovlap.results[i, paste0("ovlap", this.FC)] <- num.ovlap 
        ovlap.results[i, paste0("ovlap_alt", this.FC)] <- num.ovlap_alt
        
        ovlap.results[i, paste0("detection", this.FC)] <- round(100 * num.ovlap / N30.genes, 2)
		ovlap.results[i, paste0("detection_alt", this.FC)] <- round(100 * num.ovlap / N30.genes_alt, 2)
        ovlap.results[i, paste0("detection_altb", this.FC)] <- round(100 * num.ovlap_alt / N30.genes_alt, 2) # same for all FCs
    } # FCs

    
    # correlation between log2 fold changes
	ovlap.results[i, 'log2fc.cor'] <- cor(log2(this.fcs[tpm.pass]), log2(Dchs1.N30$Fold.Change[tpm.pass, tissue]), 
                                     use = 'complete.obs', method = 'spearman')

    ovlap.results[i, 'wald.cor'] <- cor(this.DE$Wald[tpm.pass,tissue], Dchs1.N30$Wald[tpm.pass, tissue], 
                                     use = 'complete.obs', method = 'spearman')
	
	} # i
	saveRDS(ovlap.results, file = paste0(paste('result_objs/Dchs1_ovlap.results_piece', this.N, sep="_"), ".RDS"))
}


# Run the function in parallel using mclapply
mclapply(1:TOTAL_RUNS, 
         function(ind) {
           run_4_piece(this.N = RUN_INFO[ind, 'N'])
         },
		 mc.cores = round(detectCores() / 2))

# concatenate everything (assuming this waits until everything is finished?)
for (i in 1:TOTAL_RUNS) {
	this_res <- readRDS(paste0(paste('result_objs/Dchs1_ovlap.results_piece', RUN_INFO[i, 'N'], sep="_"), ".RDS"))
	if (i==1) {
		ovlap.results <- this_res
	} else {
		ovlap.results <- rbind(ovlap.results, this_res)
	}
}
saveRDS(ovlap.results, 'result_objs/Dchs1_ovlap.results.RDS')
