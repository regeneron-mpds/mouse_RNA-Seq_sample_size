library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
source("util_scripts/openPng.R")
source("util_scripts/embiggen.R")
source("util_scripts/plot_genes_by_genotype.R")
packageVersion('DESeq2')

#### Dchs1 ####
# read in design file
design <- read.table("combined_sample_meta.tsv", header=F, sep='\t', row.names = 1)
colnames(design) <-c('sample_id', "mouse_id", 'species', 'tissue', 'genotype', 'RNA_type', 
                     'read_type', 'sequencer_type', 'count_file', 'tpm_file', 'fastq_file')
design <- design %>% filter(grepl("Dchs1", genotype))
head(design)
tissues_in_order <- sort(unique(design$tissue))

# read in gene annotation information
gene_annots <- read.table("gene_annotation.tsv", row.names = 1, sep='\t', header=T, comment.char = "")
head(gene_annots)

# read in counts and TPM files
counts <- as.matrix(read.table("GSE272152_Dchs1_counts.txt", header=T, row.names = 1))
TPM <- as.matrix(read.table("GSE272152_Dchs1_TPM.txt", header=T, row.names = 1))
stopifnot(identical(rownames(counts), rownames(gene_annots)))
stopifnot(identical(dimnames(counts), dimnames(TPM)))
rownames(counts) = rownames(TPM) = gene_annots$GeneName

# barplots showing genotype-relevant (Dchs1 and LacZ)  expression
for (this_tissue in tissues_in_order) {
  p2 <- plot_genes_by_genotype(c("Dchs1", "LacZ"), 
                               this_tissue = this_tissue, 
                               this_design = design) + 
    geom_text(aes(x=-Inf, y=Inf, hjust=0, vjust=1, label=this_tissue), fontface = 'bold', size = 10, show.legend = FALSE)
  openPng(paste0("Dchs1_LacZ_barplots_", this_tissue), p = p2, out_dir = "figures_EDA/")
}
plot_genes_by_genotype(c("Dchs1", "LacZ"), this_tissue = 'kidney')
plot_genes_by_genotype(c("Dchs1", "LacZ", "SV40"), this_tissue = 'kidney')
plot_genes_by_genotype(c("Dchs1", "LacZ"), this_tissue = 'liver')
plot_genes_by_genotype(c("Dchs1", "LacZ", "SV40"), this_tissue = 'liver')
plot_genes_by_genotype(c("Dchs1", "LacZ"), this_tissue = 'lung')
plot_genes_by_genotype(c("Dchs1", "LacZ", "SV40"), this_tissue = 'lung')

# PCA using all (Dchs1 samples)
vst_cts <- vst(as.matrix(round(counts)))
vst_pca <- prcomp(t(vst_cts))
tmp_design <- design
tmp_design$PC1 <- vst_pca$x[,"PC1"]
tmp_design$PC2 <- vst_pca$x[,"PC2"]
(p <- tmp_design %>% ggplot(aes(x=PC1, y=PC2, color=tissue)) + geom_point(size=4) + theme(legend.position = 'bottom'))
openPng("PCA_Dchs1_allsamples", p = p + embiggen(), out_dir = "figures_EDA/")

# PCA per tissue
for (this_tissue in tissues_in_order) {
  tmp_design <- design %>% filter(tissue == this_tissue)
  vst_pca <- prcomp(t(vst_cts[,tmp_design$sample_id]))
  tmp_design$PC1 <- vst_pca$x[,"PC1"]
  tmp_design$PC2 <- vst_pca$x[,"PC2"]
  p1 <- tmp_design %>% 
    ggplot(aes(x=PC1, y=PC2, color=genotype, label = mouse_id)) + 
    geom_point(size=5) + 
    theme(legend.position = 'bottom') +
    embiggen() + 
    geom_text(aes(x=-Inf, y=Inf, hjust=0, vjust=1, label=this_tissue), 
              color = 'black', fontface = 'bold', size = 10, show.legend = FALSE)
  openPng(paste0("PCA_Dchs1_", this_tissue), p = p1, out_dir = "figures_EDA/")
  openPng(paste0("PCA_Dchs1_", this_tissue, "_labeled"), p = p1 + geom_text(vjust = 2, size = 5, show.legend = FALSE), out_dir = "figures_EDA/")
  
  # correlation with Dchs1 expression?
  stopifnot(identical(rownames(tmp_design), rownames(vst_pca$x)))
  pca_gene_cors <- apply(vst_pca$x, 2, function(pcaval) cor(pcaval, vst_cts["Dchs1",tmp_design$sample_id]))
  p3 <- ggplot(data.frame(PCA_Dchs1_crlxn = pca_gene_cors,
                          PC = 1:ncol(vst_pca$x),
                          correlation_sign = ifelse(pca_gene_cors > 0, "positive", "negative")),
               aes(x = PC, y = abs(PCA_Dchs1_crlxn), fill = correlation_sign)) +
    geom_bar(stat='identity') + 
    xlab("Principal Component") + ylab("Correlation with Dchs1 expression") + 
    embiggen() + 
    geom_text(aes(x=Inf, y=Inf, hjust=1, vjust=1, label=this_tissue), 
              color = 'black', fontface = 'bold', size = 10, show.legend = FALSE) + 
    theme(legend.position = c(0.85,0.85))
  openPng(paste0("PCA_Dchs1_crlxn_before_", this_tissue), p = p3, out_dir = "figures_EDA/")
}

# remove mouse 5153265, repeat PCA per tissue
design <- design %>% filter(mouse_id != '5153265')
TPM <- TPM[, rownames(design)]
counts <- counts[, rownames(design)]

for (this_tissue in tissues_in_order) {
  tmp_design <- design %>% filter(tissue == this_tissue)
  vst_pca <- prcomp(t(vst_cts[,tmp_design$sample_id]))
  tmp_design$PC1 <- vst_pca$x[,"PC1"]
  tmp_design$PC2 <- vst_pca$x[,"PC2"]
  p1 <- tmp_design %>% 
    ggplot(aes(x=PC1, y=PC2, color=genotype, label = mouse_id)) + 
    geom_point(size=5) + 
    theme(legend.position = 'bottom') +
    embiggen() + 
    geom_text(aes(x=-Inf, y=Inf, hjust=0, vjust=1, label=this_tissue), 
              color = 'black', fontface = 'bold', size = 10, show.legend = FALSE)
  openPng(paste0("PCA_Dchs1_OLremoved_", this_tissue), p = p1, out_dir = "figures_EDA/")
  openPng(paste0("PCA_Dchs1_OLremoved_", this_tissue, "_labeled"), p = p1 + geom_text(vjust = 2, size = 5, show.legend = FALSE), out_dir = "figures_EDA/")
  
  p3 <- ggplot(data.frame(PCA_Dchs1_crlxn = pca_gene_cors,
                          PC = 1:ncol(vst_pca$x),
                          correlation_sign = ifelse(pca_gene_cors > 0, "positive", "negative")),
               aes(x = PC, y = abs(PCA_Dchs1_crlxn), fill = correlation_sign)) +
    geom_bar(stat='identity') + 
    xlab("Principal Component") + ylab("Correlation with Dchs1 expression") + 
    embiggen() + 
    geom_text(aes(x=Inf, y=Inf, hjust=1, vjust=1, label=this_tissue), 
              color = 'black', fontface = 'bold', size = 10, show.legend = FALSE) + 
    theme(legend.position = c(0.85,0.85))
  openPng(paste0("PCA_Dchs1_crlxn_after_", this_tissue), p = p3, out_dir = "figures_EDA/")
}


#### SAME AS ABOVE, FOR Fat4 ####
# read in design file
design <- read.table("combined_sample_meta.tsv", header=F, sep='\t', row.names = 1)
colnames(design) <-c('sample_id', "mouse_id", 'species', 'tissue', 'genotype', 'RNA_type', 
                     'read_type', 'sequencer_type', 'count_file', 'tpm_file', 'fastq_file')
design <- design %>% filter(grepl("Fat4", genotype))
head(design)
tissues_in_order <- sort(unique(design$tissue))

# read in gene annotation information
gene_annots <- read.table("gene_annotation.tsv", row.names = 1, sep='\t', header=T, comment.char = "")
head(gene_annots)

# read in counts and TPM files
counts <- as.matrix(read.table("GSE272152_Fat4_counts.txt", header=T, row.names = 1))
TPM <- as.matrix(read.table("GSE272152_Fat4_TPM.txt", header=T, row.names = 1))
stopifnot(identical(rownames(counts), rownames(gene_annots)))
stopifnot(identical(dimnames(counts), dimnames(TPM)))
rownames(counts) = rownames(TPM) = gene_annots$GeneName

# barplots showing genotype-relevant (Fat4 and LacZ)  expression
for (this_tissue in tissues_in_order) {
  p2 <- plot_genes_by_genotype(c("Fat4", "LacZ"), 
                               this_tissue = this_tissue, 
                               this_design = design) + 
    geom_text(aes(x=-Inf, y=Inf, hjust=0, vjust=1, label=this_tissue), fontface = 'bold', size = 10, show.legend = FALSE)
  openPng(paste0("Fat4_LacZ_barplots_", this_tissue), p = p2, out_dir = "figures_EDA/")
}
plot_genes_by_genotype(c("Fat4", "LacZ"), this_tissue = 'kidney')
plot_genes_by_genotype(c("Fat4", "LacZ", "SV40"), this_tissue = 'kidney')
plot_genes_by_genotype(c("Fat4", "LacZ"), this_tissue = 'liver')
plot_genes_by_genotype(c("Fat4", "LacZ", "SV40"), this_tissue = 'liver')
plot_genes_by_genotype(c("Fat4", "LacZ"), this_tissue = 'lung')
plot_genes_by_genotype(c("Fat4", "LacZ", "SV40"), this_tissue = 'lung')

# PCA using all (Fat4 samples)
vst_cts <- vst(as.matrix(round(counts)))
vst_pca <- prcomp(t(vst_cts))
tmp_design <- design
tmp_design$PC1 <- vst_pca$x[,"PC1"]
tmp_design$PC2 <- vst_pca$x[,"PC2"]
(p <- tmp_design %>% ggplot(aes(x=PC1, y=PC2, color=tissue)) + geom_point(size=4) + theme(legend.position = 'bottom'))
openPng("PCA_Fat4_allsamples", p = p + embiggen(), out_dir = "figures_EDA/")

# PCA per tissue
for (this_tissue in tissues_in_order) {
  tmp_design <- design %>% filter(tissue == this_tissue)
  vst_pca <- prcomp(t(vst_cts[,tmp_design$sample_id]))
  tmp_design$PC1 <- vst_pca$x[,"PC1"]
  tmp_design$PC2 <- vst_pca$x[,"PC2"]
  p1 <- tmp_design %>% 
    ggplot(aes(x=PC1, y=PC2, color=genotype, label = mouse_id)) + 
    geom_point(size=5) + 
    theme(legend.position = 'bottom') +
    embiggen() + 
    geom_text(aes(x=-Inf, y=Inf, hjust=0, vjust=1, label=this_tissue), 
              color = 'black', fontface = 'bold', size = 10, show.legend = FALSE)
  openPng(paste0("PCA_Fat4_", this_tissue), p = p1, out_dir = "figures_EDA/")
  openPng(paste0("PCA_Fat4_", this_tissue, "_labeled"), p = p1 + geom_text(vjust = 2, size = 5, show.legend = FALSE), out_dir = "figures_EDA/")
  
  # correlation with Fat4 expression?
  stopifnot(identical(rownames(tmp_design), rownames(vst_pca$x)))
  pca_gene_cors <- apply(vst_pca$x, 2, function(pcaval) cor(pcaval, vst_cts["Fat4",tmp_design$sample_id]))
  p3 <- ggplot(data.frame(PCA_Fat4_crlxn = pca_gene_cors,
                    PC = 1:ncol(vst_pca$x),
                    #                  PC = factor(1:ncol(vst_pca$x), levels=1:ncol(vst_pca$x)), # factor(colnames(vst_pca$x), levels=colnames(vst_pca$x)),
                    correlation_sign = ifelse(pca_gene_cors > 0, "positive", "negative")),
         aes(x = PC, y = abs(PCA_Fat4_crlxn), fill = correlation_sign)) +
    geom_bar(stat='identity') + 
    xlab("Principal Component") + ylab("Correlation with Fat4 expression") + 
    embiggen() + 
    geom_text(aes(x=Inf, y=Inf, hjust=1, vjust=1, label=this_tissue), 
              color = 'black', fontface = 'bold', size = 10, show.legend = FALSE) + 
    theme(legend.position = c(0.85,0.85))
  openPng(paste0("PCA_Fat4_crlxn_", this_tissue), p = p3, out_dir = "figures_EDA/")
}



         
       

this_tissue <- 'liver'
tmp_design$PC1 <- vst_pca$x[,"PC2"]
tmp_design$PC2 <- vst_pca$x[,"PC3"]
p1 <- tmp_design %>% 
  ggplot(aes(x=PC1, y=PC2, color=genotype, label = mouse_id)) + 
  geom_point(size=5) + 
  theme(legend.position = 'bottom') +
  embiggen() + 
  geom_text(aes(x=-Inf, y=Inf, hjust=0, vjust=1, label=this_tissue), 
            color = 'black', fontface = 'bold', size = 10, show.legend = FALSE); p1
