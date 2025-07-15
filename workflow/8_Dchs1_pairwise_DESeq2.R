library(DESeq2)
packageVersion('DESeq2')
library(dplyr)
library(parallel)
source("util_scripts/DESeq2_wrapper_code.R")

# INIT
NUM_TRIALS <- 25
N_RANGE <- 14:3
RUN_INFO <- expand.grid(trial = 1:NUM_TRIALS, N = N_RANGE)[,c(2,1)]
TOTAL_RUNS <- nrow(RUN_INFO)

# For reproducibility: Generate unique RNG streams for each worker
RNGkind("L'Ecuyer-CMRG")
set.seed(17)
streams <- list() 
s <- .Random.seed
for (i in 1:TOTAL_RUNS) {
  streams[[i]] <- s
  s <- nextRNGStream(s)
}

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
genes <- rownames(counts)

## Perform pairwise comparisons
table(design$tissue, design$genotype)
unique.WT.animals <- unique(design$mouse_id[design$genotype == "Dchs1.WT"])
unique.Het.animals <- unique(design$mouse_id[design$genotype == "Dchs1.Het"])
suppressWarnings(dir.create("DE_runs_pairwise/"))

run_pairwise_DESeq2 <- function(seed, this_N, trial_num) {
  .Random.seed <<- seed
  OUTPUT_FILE1 <- paste0("DE_runs_pairwise/Dchs1_N", this_N, "_trial", trial_num, "_A.RDS")
  OUTPUT_FILE2 <- paste0("DE_runs_pairwise/Dchs1_N", this_N, "_trial", trial_num, "_B.RDS")
  # if (file.exists(OUTPUT_FILE1)) { print(paste(OUTPUT_FILE1, 'already exists, skipping...')); next }

  # randomly select animals
  this.WT.animals1 <- sample(unique.WT.animals, this_N, replace = FALSE)
  this.WT.animals2 <- sample(setdiff(unique.WT.animals, this.WT.animals1), this_N, replace = FALSE)
  this.Het.animals1 <- sample(unique.Het.animals, this_N, replace = FALSE)
  this.Het.animals2 <- sample(setdiff(unique.Het.animals, this.Het.animals1), this_N, replace = FALSE)

  for (grp in 1:2) {
    this.WT.animals <- get(paste0('this.WT.animals', grp))
    this.Het.animals <- get(paste0('this.Het.animals', grp))
    
    extra_liver_sample <- 'none'
    if ("5144552" %in% this.Het.animals) { # only mouse w/ less than 4 tissues; pick one more Het liver sample
      remaining_Het_animals <- setdiff(unique.Het.animals, c(this.Het.animals1, this.Het.animals2))
      extra_animal <- sample(remaining_Het_animals, 1)
      extra_liver_sample <- design %>% 
        filter(mouse_id == extra_animal, tissue == 'liver') %>%
        pull(sample_id)
    }
    this_animals <- c(this.WT.animals, this.Het.animals)
    this_design <- design %>% 
      filter(mouse_id %in% this_animals | sample_id == extra_liver_sample)
    
    # run DESeq2 and record results
    DE_res <- DESeq2_wrapper(counts = counts[, rownames(this_design)],
                             TPM = TPM[,rownames(this_design)],
                             design = this_design, 
                             tissues_in_order = tissues_in_order)
    DE_res$WT <- sort(this.WT.animals)
    DE_res$Het <- sort(this.Het.animals)
    DE_res$extra_liver_sample <- extra_liver_sample
    saveRDS(DE_res, file = get(paste0('OUTPUT_FILE', grp)))
  }
}

# Run the function in parallel using mclapply
mclapply(1:TOTAL_RUNS, 
         function(ind) {
           run_pairwise_DESeq2(streams[[ind]], 
                                this_N = RUN_INFO[ind, 'N'], 
                                trial_num = RUN_INFO[ind, 'trial'])
         },
         mc.cores = detectCores() / 2)
