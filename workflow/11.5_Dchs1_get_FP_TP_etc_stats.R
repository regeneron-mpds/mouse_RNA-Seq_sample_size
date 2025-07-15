library(dplyr)
library(tidyr)
library(parallel)

FCtypes <- c('default', 'apeglm', 'ashr')
TPM.CUTOFFs <- c(0,1,5)
FC.CUTOFFs <- c(1, 1.2, 1.5, 2, 4)

RUN_INFO <- expand.grid(FCtype = FCtypes, tpm.cutoff = TPM.CUTOFFs, fc.cutoff = FC.CUTOFFs)
TOTAL_RUNS <- nrow(RUN_INFO)

run_11.5_piece <- function(this.FCtype, tpm.cutoff, fc.cutoff) {
	Dchs1_res <- readRDS("full_signature_refs/Dchs1_res.RDS")
	NUM_TRIALS <- 40
	N_RANGE <- 3:28
	TISSUES_IN_ORDER <- colnames(Dchs1_res$Fold.Change)
	PADJ_CUTOFF <- 0.05 

	  # decide which fold change to use
	  if (this.FCtype == "default") {
		  fc_str <- 'Fold.Change'            
	  } else if (this.FCtype == "ashr") {
		  fc_str <- 'FC.ashr'
	  } else if (this.FCtype == 'apeglm') {
		  fc_str <- 'FC.apeglm'
	  }

    # order genes by FC in the N30 experiment (separately for each tissue)
    genes_in_N30.FC_order <- list()
    N30.fc.is_up <- list()
    N30.pval.pass <- list()
    for (tissue in TISSUES_IN_ORDER) {
		this_fcs <- abs(log2(Dchs1_res[[fc_str]][,tissue]))
        this_fcs[is.na(this_fcs)] <- 0
		this.tpms <- pmax(Dchs1_res$meanX.TPM, Dchs1_res$meanY.TPM)
		this_fcs <- this_fcs[this_fcs >= log2(fc.cutoff) & this.tpms >= tpm.cutoff]
        genes_in_order <- names(sort(this_fcs))
        genes_in_N30.FC_order[[tissue]] <- genes_in_order
        N30.fc.is_up[[tissue]] <- Dchs1_res[[fc_str]][genes_in_order, tissue] > 1

        N30.pvals <- Dchs1_res$Padj[genes_in_order, tissue]
        N30.pval.pass[[tissue]] <- !is.na(N30.pvals) & N30.pvals < PADJ_CUTOFF
    }

	# create N-by-tissue
    FPR = TNR = FNR = TPR = FPR_sd = TNR_sd = FNR_sd = TPR_sd = 
	FPR_med = TNR_med = FNR_med = TPR_med = FPR_N = TNR_N = FNR_N = TPR_N = 
	FDR = FDR_sd = FDR_med = FDR_N =
        matrix(NA, length(N_RANGE), length(TISSUES_IN_ORDER),
               dimnames=list(as.character(N_RANGE), TISSUES_IN_ORDER)) %>% as.data.frame()
	
    for (N in as.character(N_RANGE)) { 
        if (as.numeric(N) %% 4 == 1) { message(paste0("fc=", fc.cutoff, " N=", N)) }

        # create trial-by-tissue matrices
        TP = FP = TN = FN = 
            matrix(NA, NUM_TRIALS, length(TISSUES_IN_ORDER),
                   dimnames=list(as.character(1:NUM_TRIALS), TISSUES_IN_ORDER)) %>% as.data.frame()

        for (trial in 1:NUM_TRIALS) {        
            this_res <- readRDS(paste0("DE_runs_subsample/Dchs1_N", N, "_trial", trial, ".RDS"))
            for (tissue in TISSUES_IN_ORDER) {
                genes_in_order <- genes_in_N30.FC_order[[tissue]]

                this_fcs <- this_res[[fc_str]][genes_in_order, tissue]
                this_pval <- this_res$Padj[genes_in_order, tissue]
                signagree <- (this_fcs > 1) == N30.fc.is_up[[tissue]]
                this_pval_pass <- !is.na(this_pval) & this_pval < 0.05 & signagree                
                TP[trial, tissue] <- sum(this_pval_pass & N30.pval.pass[[tissue]])
                FP[trial, tissue] <- sum(this_pval_pass & ! N30.pval.pass[[tissue]])
                TN[trial, tissue] <- sum(! this_pval_pass & ! N30.pval.pass[[tissue]])
                FN[trial, tissue] <- sum(! this_pval_pass & N30.pval.pass[[tissue]])
            }
        }

        for (tissue in TISSUES_IN_ORDER) {
            tTPR <- TP[, tissue] / (TP[, tissue] + FN[, tissue])
            tFPR <- FP[, tissue] / (FP[, tissue] + TN[, tissue])            
            tTNR <- 1 - tFPR 
            tFNR <- 1 - tTPR 
            tTPR <- tTPR[!is.na(tTPR)]
            tFPR <- tFPR[!is.na(tFPR)]
            tTNR <- tTNR[!is.na(tTNR)]
            tFNR <- tFNR[!is.na(tFNR)]

            TPR[N, tissue] <- mean(tTPR)
            FPR[N, tissue] <- mean(tFPR)
            TNR[N, tissue] <- mean(tTNR)
            FNR[N, tissue] <- mean(tFNR)
            TPR_sd[N, tissue] <- sd(tTPR)
            FPR_sd[N, tissue] <- sd(tFPR)
            TNR_sd[N, tissue] <- sd(tTNR)
            FNR_sd[N, tissue] <- sd(tFNR)
			TPR_med[N, tissue] <- median(tTPR)
            FPR_med[N, tissue] <- median(tFPR)
            TNR_med[N, tissue] <- median(tTNR)
            FNR_med[N, tissue] <- median(tFNR)			
            TPR_N[N, tissue] <- sum(TP[,tissue])
            FPR_N[N, tissue] <- sum(FP[,tissue])
            TNR_N[N, tissue] <- sum(TN[,tissue])
            FNR_N[N, tissue] <- sum(FN[,tissue])
			
			tFDR <- FP[, tissue] / (FP[, tissue] + TP[, tissue])
            tFDR <- tFDR[!is.na(tFDR)]
            FDR[N, tissue] <- mean(tFDR)
            FDR_sd[N, tissue] <- sd(tFDR)
            FDR_med[N, tissue] <- median(tFDR)
        } # tissue
    } # N
    this_res <- rbind(
        TPR %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='mean_TPR'),
        TPR_med %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='median_TPR'),
        TPR_sd %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='sd_TPR'),
		TPR_N %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='TP'),
        FPR %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='mean_FPR'),
        FPR_med %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='median_FPR'),
        FPR_sd %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='sd_FPR'),
		FPR_N %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='FP'),
        TNR %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='mean_TNR'),
        TNR_med %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='median_TNR'),
        TNR_sd %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='sd_TNR'),
		TNR_N %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='TN'),
        FNR %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='mean_FNR'),
        FNR_med %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='median_FNR'),		
        FNR_sd %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='sd_FNR'),
		FNR_N %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='FN'),
        FDR %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='mean_FDR'),
        FDR_med %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='median_FDR'),		
        FDR_sd %>% tibble::rownames_to_column('N') %>% pivot_longer(-N, names_to='tissue') %>% mutate(type ='sd_FDR')) %>%
    mutate(fc = fc.cutoff,
		   FCtype = this.FCtype,
		  tpm = tpm.cutoff)
	
	saveRDS(this_res, file = paste0(paste('result_objs/Dchs1_FP_FN_etc_piece', fc.cutoff, this.FCtype, tpm.cutoff, sep="_"), ".RDS"))
}


# Run the function in parallel using mclapply
mclapply(1:TOTAL_RUNS, 
         function(ind) {
           run_11.5_piece(this.FCtype = RUN_INFO[ind, 'FCtype'],
						  tpm.cutoff = RUN_INFO[ind, 'tpm.cutoff'],
						  fc.cutoff = RUN_INFO[ind, 'fc.cutoff'])
         },
		 mc.cores = round(detectCores() / 2))

# concatenate everything (assuming this waits until everything is finished?)
for (i in 1:TOTAL_RUNS) {
	this_res <- readRDS(paste0(paste('result_objs/Dchs1_FP_FN_etc_piece',
								RUN_INFO[i, 'fc.cutoff'], 
								RUN_INFO[i, 'FCtype'],
								RUN_INFO[i, 'tpm.cutoff'], sep="_"), ".RDS"))
	if (i==1) {
		all_res <- this_res
	} else {
		all_res <- rbind(all_res, this_res)
	}
}
saveRDS(all_res, file = "result_objs/Dchs1_FP_FN_etc_stats_removeNAs.RDS")
