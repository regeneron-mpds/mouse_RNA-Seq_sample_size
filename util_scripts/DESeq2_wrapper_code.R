## DESEQ2 wrapper code
FC_CUTOFF = 1.5
PADJ_CUTOFF = 0.05
PVAL_CUTOFF = 0.01

DESeq2_wrapper <- function(counts, TPM, design, tissues_in_order) {
  # init
  design$genotype <- gsub("\\S+\\.", "", design$genotype) # Dchs1.WT -> WT
  design$genotype <- factor(design$genotype, levels = c("WT", "Het"))
  genes <- rownames(counts)
  DE_res <- list()
  DE_res[["Fold.Change"]] = matrix(1, length(genes), length(tissues_in_order), dimnames =list(genes, tissues_in_order))
  DE_res[["Pval"]] = DE_res[["Fold.Change"]]
  DE_res[["Padj"]] = DE_res[["Fold.Change"]]
  DE_res[["Selected"]] = DE_res[["Fold.Change"]] - 1
  DE_res[["Selected.pval"]] = DE_res[["Fold.Change"]] - 1
  DE_res[["Wald"]] = DE_res[["Fold.Change"]] * NA
  DE_res[["padj.01"]] = DE_res[["padj.1"]] = DE_res[["padj.999"]] = DE_res[["Fold.Change"]] * NA
  DE_res[["FC.ashr"]] = DE_res[["FC.apeglm"]] = DE_res[["Fold.Change"]] * NA

  if (! missing(TPM)) {
    DE_res[["meanX.TPM"]] <- matrix(NA, length(genes), length(tissues_in_order), dimnames = list(genes, tissues_in_order))
    DE_res[["meanY.TPM"]] <- matrix(NA, length(genes), length(tissues_in_order), dimnames = list(genes, tissues_in_order))
  }
  
  #perform DESeq2 to populate DE_res object
  for (this_tissue in tissues_in_order) { # pretty fast 1-2 minutes total for the 4 tissues
    this_design <- design %>% filter(tissue == this_tissue)
    this_counts <- counts[, rownames(this_design)]
    this_TPM <- TPM[, rownames(this_design)]

    print(paste("starting", this_tissue, ncol(this_counts)))
    this_dds <- DESeqDataSetFromMatrix(countData = round(this_counts),
                                       colData = this_design, 
                                       design =~ genotype)

    system.time(this_dds <- DESeq(this_dds))
    resultsNames(this_dds)
    
    this_res <- results(this_dds, name = "genotype_Het_vs_WT", alpha = 0.05)
    stopifnot(identical(rownames(this_res), genes))
    
    DE_res[["Fold.Change"]][, this_tissue] <- 2^this_res$log2FoldChange
    DE_res[["Pval"]][, this_tissue] <- this_res$pvalue
    DE_res[["Padj"]][, this_tissue] <- this_res$padj
    DE_res[["Wald"]][, this_tissue] <- this_res$stat
	DE_res[["Selected"]][, this_tissue] <- as.numeric(!is.na(this_res$padj) &
                                                          this_res$padj < PADJ_CUTOFF &
                                                          abs(this_res$log2FoldChange) >= log2(FC_CUTOFF))
    
    DE_res[["Selected.pval"]][, this_tissue] <- as.numeric(!is.na(this_res$pvalue) &
                                                               this_res$pvalue < PVAL_CUTOFF &
                                                               abs(this_res$log2FoldChange) >= log2(FC_CUTOFF))
    
    (WT.samples <- this_dds@colData@rownames[which(this_dds@colData@listData$genotype == "WT")])
    DE_res[["meanX.TPM"]][genes, this_tissue] <- apply(TPM[genes, WT.samples], 1, function(x) round(mean(x), 3))
    (Het.samples <- this_dds@colData@rownames[which(this_dds@colData@listData$genotype == "Het")])
    DE_res[["meanY.TPM"]][genes, this_tissue] <- apply(TPM[genes, Het.samples], 1, function(x) round(mean(x), 3))

	# reviewer1.alpha.plus.shrinkage
    DE_res[["padj.01"]][, this_tissue] <- results(this_dds, name = "genotype_Het_vs_WT", alpha = 0.01)$padj
    DE_res[["padj.1"]][, this_tissue] <- results(this_dds, name = "genotype_Het_vs_WT", alpha = 0.1)$padj
    DE_res[["padj.999"]][, this_tissue] <- results(this_dds, name = "genotype_Het_vs_WT", alpha = 0.999)$padj
    DE_res[["FC.ashr"]][, this_tissue] <- 2^(lfcShrink(this_dds, type = 'ashr', coef = 2)$log2FoldChange)
    DE_res[["FC.apeglm"]][, this_tissue] <- 2^(lfcShrink(this_dds, type = 'apeglm', coef = 2)$log2FoldChange)
  }
  
  return(DE_res)
}
