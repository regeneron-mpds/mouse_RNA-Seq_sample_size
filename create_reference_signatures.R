set.seed(17)
library(DESeq2)
library(dplyr)
packageVersion('DESeq2')
source("extras/DESeq2_wrapper_code.R")

suppressWarnings(dir.create("full_signature_refs/"))

# read in gene annotation information
gene_annots <- read.table("gene_annotation.tsv", row.names = 1, sep='\t', header=T, comment.char = "")
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
dchs1_res <- DESeq2_wrapper(counts, TPM, design, tissues_in_order)
saveRDS(dchs1_res, file = "full_signature_refs/dchs1_res.RDS")


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
fat4_res <- DESeq2_wrapper(counts, TPM, design, tissues_in_order)
saveRDS(fat4_res, file = "full_signature_refs/fat4_res.RDS")


##################################
### Create gold standard lists ###
##################################
gold_standards <- list()
gold_standards$dchs1 <- list()  # genes significant in Dchs1
gold_standards$fat4 <- list()   # genes significant in Fat4

genes_in_order <- rownames(dchs1_res[[1]])
stopifnot(identical(genes_in_order, rownames(fat4_res[[1]])))

for (tissue in tissues_in_order) {
  for (tpm.cutoff in c(0, 1, 5)) {
    for (ptype in c('Pval', 'Padj')) {
      for (fc.cutoff in c(1, 1.2, 1.5, 2)) {                
        for (p.thresh in c(0.01, 0.01, 0.05, 0.1, 1)) {
          str <- paste(tissue, tpm.cutoff, ptype, p.thresh, fc.cutoff, sep = ',')
              
          gold_standards$dchs1[[str]] <- as.numeric(!is.na(dchs1_res[[ptype]][,tissue])) &
                                              as.numeric(!is.na(dchs1_res$Fold.Change[,tissue])) &
                                              dchs1_res[[ptype]][,tissue] <= p.thresh & 
                                              ( (dchs1_res$Fold.Change[,tissue] >= fc.cutoff & dchs1_res$meanY.TPM[,tissue] >= tpm.cutoff) |
                                                (dchs1_res$Fold.Change[,tissue] <=  1/fc.cutoff & dchs1_res$meanX.TPM[,tissue] >= tpm.cutoff))
          gold_standards$dchs1[[str]] <- genes_in_order[gold_standards$dchs1[[str]]]
          gold_standards$dchs1[[str]] <- setdiff(gold_standards$dchs1[[str]], c("Dchs1", "LacZ", "SV40"))
  
          
          gold_standards$fat4[[str]] <- as.numeric(!is.na(fat4_res[[ptype]][,tissue])) &
                                              as.numeric(!is.na(fat4_res$Fold.Change[,tissue])) &
                                              fat4_res[[ptype]][,tissue] <= p.thresh & 
                                              ( (fat4_res$Fold.Change[,tissue] >= fc.cutoff & fat4_res$meanY.TPM[,tissue] >= tpm.cutoff) |
                                                (fat4_res$Fold.Change[,tissue] <=  1/fc.cutoff & fat4_res$meanX.TPM[,tissue] >= tpm.cutoff))
          gold_standards$fat4[[str]] <- genes_in_order[gold_standards$fat4[[str]]]
          gold_standards$fat4[[str]] <- setdiff(gold_standards$fat4[[str]], c("Fat4", "LacZ", "SV40"))
        } # for pval.thresh
      } #fc.cutoff
    } # pvalue type
  } # tpm.cutoff
} #tissue
saveRDS(gold_standards, file = 'gold_standards.RDS')

