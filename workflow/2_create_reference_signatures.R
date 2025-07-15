set.seed(17)
library(DESeq2)
library(dplyr)
library(ashr)
library(apeglm)
packageVersion('DESeq2')
source("util_scripts/DESeq2_wrapper_code.R")

suppressWarnings(dir.create("full_signature_refs/"))
suppressWarnings(dir.create("result_objs/"))

# read in gene annotation information
gene_annots <- read.table("annotation/gene_annotation.tsv", row.names = 1, sep='\t', header=T, comment.char = "")
head(gene_annots)

#################
##### Dchs1 #####
#################
# read in design file
design <- read.table("combined_sample_meta.tsv", header=F, sep='\t', row.names = 1)
colnames(design) <-c('sample_id', "mouse_id", 'species', 'tissue', 'genotype', 'RNA_type', 
                     'read_type', 'sequencer_type', 'count_file', 'tpm_file', 'fastq_file')
design <- design %>% filter(grepl("Dchs1", genotype))
head(design)
tissues_in_order <- sort(unique(design$tissue))

# read in counts and TPM files
counts <- as.matrix(read.table("GSE272152_Dchs1_counts.txt", header=T, row.names = 1))
TPM <- as.matrix(read.table("GSE272152_Dchs1_TPM.txt", header=T, row.names = 1))
stopifnot(identical(rownames(counts), rownames(gene_annots)))
stopifnot(identical(dimnames(counts), dimnames(TPM)))
rownames(counts) = rownames(TPM) = gene_annots$GeneName

# remove mouse outlier
design <- design %>% filter(mouse_id != '5153265')
TPM <- TPM[, rownames(design)]
counts <- counts[, rownames(design)]

# run DESeq2
Dchs1_res <- DESeq2_wrapper(counts, TPM, design, tissues_in_order)
saveRDS(Dchs1_res, file = "full_signature_refs/Dchs1_res.RDS")


#################
##### Fat4 #####
#################
# read in design file
design <- read.table("combined_sample_meta.tsv", header=F, sep='\t', row.names = 1)
colnames(design) <-c('sample_id', "mouse_id", 'species', 'tissue', 'genotype', 'RNA_type', 
                     'read_type', 'sequencer_type', 'count_file', 'tpm_file', 'fastq_file')
design <- design %>% filter(grepl("Fat4", genotype))
head(design)
stopifnot(identical(tissues_in_order, sort(unique(design$tissue))))

# read in counts and TPM files
counts <- as.matrix(read.table("GSE272152_Fat4_counts.txt", header=T, row.names = 1))
TPM <- as.matrix(read.table("GSE272152_Fat4_TPM.txt", header=T, row.names = 1))
stopifnot(identical(rownames(counts), rownames(gene_annots)))
stopifnot(identical(dimnames(counts), dimnames(TPM)))
rownames(counts) = rownames(TPM) = gene_annots$GeneName

# run DESeq2
Fat4_res <- DESeq2_wrapper(counts, TPM, design, tissues_in_order)
saveRDS(Fat4_res, file = "full_signature_refs/Fat4_res.RDS")
##################################
### Create gold standard lists ###
##################################
gold_standards <- list()
GS_passfail_results <- list()

genes_in_order <- rownames(Dchs1_res[[1]])
stopifnot(identical(genes_in_order, rownames(Fat4_res[[1]])))

for (this_expt in c("Dchs1", "Fat4")) {
this.DE <- get(paste0(this_expt, '_res'))
gold_standards[[this_expt]] <- list() 
GS_passfail_results[[this_expt]] <- list()
for (tissue in tissues_in_order) {
	for (ptype in c('Pval', 'Padj')) {
		for (alphaType in c('std.05', 'match.padj.cutoff')) {
			if (ptype == 'Pval' && alphaType == 'match.padj.cutoff') { next }
			
			for (pval.cutoff in c(0.01, 0.05, 0.1, 1)) {
				# decide which padj to use
				if (alphaType == 'std.05' || pval.cutoff == 0.05 || ptype == 'Pval') {
					this.pvals <- this.DE[[ptype]][,tissue]
				} else if (alphaType == 'match.padj.cutoff' && pval.cutoff == .01) {
					this.pvals <- this.DE[['padj.01']][,tissue]
				} else if (alphaType == 'match.padj.cutoff' && pval.cutoff == .1) {
					this.pvals <- this.DE[['padj.1']][,tissue]
				} else if (alphaType == 'match.padj.cutoff' && pval.cutoff == 1) {
					this.pvals <- this.DE[['padj.999']][,tissue]
				}
				pval.pass <- !is.na(this.pvals) & this.pvals <= pval.cutoff
				
				for (FCtype in c("default", "ashr", "apeglm")) {
					# decide which fold change to use
					if (FCtype == "default") { 
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

					
					for (tpm.cutoff in c(0, 1, 5)) {
						tpm.pass <- this.tpm >= tpm.cutoff
						tpm.pass <- tpm.pass & !names(this.tpm) %in% c(this_expt, "LacZ", "SV40") 
						for (fc.cutoff in c(1, 1.2, 1.5, 2)) { 
							fc.pass <- !is.na(this.fcs) & (this.fcs >= fc.cutoff | this.fcs <= 1/fc.cutoff)

							str <- paste(tissue, tpm.cutoff, ptype, alphaType, FCtype, pval.cutoff, fc.cutoff, sep = ',')
							gold_standards[[this_expt]][[str]] <- genes_in_order[pval.pass & fc.pass & tpm.pass]

							pf <- matrix(NA, nrow(this.DE[[1]]), 5, 
										 dimnames = list(rownames(this.DE[[1]]), paste0('grp', 0:4)))
							pf[,'grp0'] <- !tpm.pass
							pf[,'grp1'] <- tpm.pass & pval.pass & fc.pass # same as gs
							pf[,'grp2'] <- tpm.pass & pval.pass & !fc.pass
							pf[,'grp3'] <- tpm.pass & !pval.pass & fc.pass
							pf[,'grp4'] <- tpm.pass & !pval.pass & !fc.pass
							GS_passfail_results[[this_expt]][[str]] <- as.data.frame(pf)
						} # tpm.cutoff
					} # fc.cutoff
				} # FCtype 
			 }#  pval.cutoff
			} # alphaType
		} # ptype
	} #tissue
} # transgenic
saveRDS(gold_standards, file = 'result_objs/gold_standards.RDS')
saveRDS(GS_passfail_results, file = "result_objs/gold_standard_passfail_info.RDS")
