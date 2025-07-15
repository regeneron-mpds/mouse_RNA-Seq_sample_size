library(ggplot2)
library(dplyr)
library(tidyr)
source("util_scripts/embiggen.R")
source("util_scripts/pub_figure_plotting_params.R")
source("util_scripts/openPng.R")


# INIT
PTYPE <- 'Padj'
TISSUES_IN_ORDER <- c('heart', 'kidney', 'liver', 'lung')
N30_DEres <- readRDS("full_signature_refs/Dchs1_res.RDS")
gold_standards <- readRDS("result_objs/gold_standards.RDS")
suppressWarnings(dir.create('figures/'))

################################
### FIGURE 2 pre-processing ####
################################
# read in statistics from down-sampling
ovlap.results <- readRDS('result_objs/Dchs1_ovlap.results.RDS')
ovlap.results$N <- factor(ovlap.results$N, 
                          levels = as.character(min(ovlap.results$N):max(ovlap.results$N)))

# rejigger columns 
rejigger_ovlap_res <- function(ovlap.results) {
  
ovlap.results.full <- ovlap.results %>%
  mutate(FDR1 = 100 * (signif1 - ovlap1) / signif1,
        FDR1.2 = 100 * (signif1.2 - ovlap1.2) / signif1.2,
        FDR1.5 = 100 * (signif1.5 - ovlap1.5) / signif1.5,
        FDR2 = 100 *  (signif2 - ovlap2) / signif2) %>%
        select(!matches("^row\\d+_")) %>% 
  replace(is.na(.), 0)
head(ovlap.results.full)

# re-rejigger columns
value_start_col <- which(colnames(ovlap.results.full) == 'signif1')
ovlap.results.long <- ovlap.results.full %>%
  tidyr::pivot_longer(cols = all_of(value_start_col:ncol(ovlap.results.full)), names_to = "type") %>%
  mutate(FC = gsub("^[^\\d\\.]+", "", type, perl=T),
         FC = ifelse(grepl("cor", FC), NA, FC),
         type = gsub("\\d(\\.\\d)?$", "", type))
return(ovlap.results.long)
}
ovlap.results.long <- rejigger_ovlap_res(ovlap.results)
table(ovlap.results.long$FC, ovlap.results.long$type, useNA = 'always')

###############################################
#### FIGURE 2A - FDR/sensitivity by tissue ####
###############################################
plot_2A <- function(ovlap.results.long, fig_str = 'Figure2A',
                    FCtypes = c('default'),  # 'apeglm', 'ashr'),
                    FCs = 1.5, # c(1, 1.2, 1.5, 2)) {
                    TPM.CUTOFFs = c(0,1)) {
for (this.FCtype in FCtypes) {
for (this.FC in FCs) {  
  for (this.TPM.CUTOFF in TPM.CUTOFFs) {  # c(0, 1, 5)) { 
    for (this.PVAL.CUTOFF in c(0.05)) {  # c(0.01, 0.05, 0.1, 1)
      p <- ovlap.results.long %>% 
        filter(pval.cutoff == this.PVAL.CUTOFF,
               type %in% c("detection", "FDR"),
               FC == this.FC,
               alphaType == 'std.05', # shouldn't matter for 0.05
               FCtype == this.FCtype,
               tpm.cutoff == this.TPM.CUTOFF) %>%
        mutate(type = gsub("detection", "Sensitivity", type),
               type = gsub("FDR", "False Discovery Rate", type)) %>%
        ggplot(aes(x = N, y = value)) + 
        facet_grid(tissue ~ type) + 
        ylab("") + 
        xlab(paste0("Sample Size")) +  
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"),
              legend.position = 'none') + 
        scale_x_discrete(labels=PRETTY_LABELS)
      
      p1 <- pub_boxplot_theme(p, 'png')
      p2 <- pub_boxplot_theme(p, 'tiff')
      
      outfile <- paste0(fig_str, '_sensitivity.FDR_FC', this.FC, '_TPM', this.TPM.CUTOFF, '_', PTYPE, this.PVAL.CUTOFF, '_', this.FCtype)
      #openPng(outfile, p = p1, height = 900*3, width = 1200*3, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    }
  }
} # FC
} # FCtype
}
plot_2A(ovlap.results.long, 'Figure2A')

###################################
### FIGURE 2B - vary FC cutoff ####
# & Supp Fig 1A ####
###################################
plot_2B <- function(ovlap.results.long, fig_str = 'Figure2B',
                    FCtypes = c('default'), #'apeglm', 'ashr'),
                    TPM.CUTOFFs = c(0,1)) {
  for (this.FCtype in FCtypes) {
    for (this.TISSUE in TISSUES_IN_ORDER) {
      for (this.TPM.CUTOFF in TPM.CUTOFFs) {  # c(0, 1, 5)) { 
        for (this.PVAL.CUTOFF in c(0.05)) { #c(0.01, 0.05, 0.1, 1)) {
      p <- ovlap.results.long %>% 
        filter(pval.cutoff == this.PVAL.CUTOFF,
               alphaType == 'std.05', # shouldn't matter
               FCtype == this.FCtype,
               type %in% c("detection", "FDR"),
               tissue == this.TISSUE,
               tpm.cutoff == this.TPM.CUTOFF) %>%
        mutate(type = gsub("detection", "Sensitivity", type),
               type = gsub("FDR", "False Discovery Rate", type)) %>%
        ggplot(aes(x = N, y = value)) + 
        facet_grid(FC ~ type, labeller = labeller(FC = FCfunc, tpm.cutoff = TPMfunc)) +  
        ylab("") + 
        xlab(paste0("Sample Size")) + 
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"),
              legend.position = 'none') + 
        scale_x_discrete(labels=PRETTY_LABELS)
      
      p1 <- pub_boxplot_theme(p, 'png')
      p2 <- pub_boxplot_theme(p, 'tiff')
      
      outfile <- paste0(fig_str, '_sensitivity.FDR_TPM', this.TPM.CUTOFF, '_', PTYPE, '_', this.PVAL.CUTOFF, '_', this.TISSUE, '_', this.FCtype)
      #openPng(outfile, p = p1, height = 900*3, width = 1200*3, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    } # P
  } # TPM
} # tissues
} # FCtype
}
plot_2B(ovlap.results.long, 'Figure2B')



#################################
### FIGURE 2C - vary P-value ####
# & Supp Fig 1B ####
#################################
plot_2C <- function(ovlap.results.long, fig_str = 'Figure2C',
                    FCtypes = c('default'), #'apeglm', 'ashr'),
                    FCs = 1.5, 
                    TPM.CUTOFFs = c(0,1)) {

  # arrange P-value from less to more specific
  ovlap.results.long$pval.cutoff <- factor(ovlap.results.long$pval.cutoff, 
                                           levels=c(1, 0.1, 0.05, 0.01))

  for (this.FCtype in FCtypes) {
    for (this.TISSUE in TISSUES_IN_ORDER) {
      for (this.FC.CUTOFF in FCs) { # c(1, 1.2, 1.5, 2)
        for (this.TPM.CUTOFF in TPM.CUTOFFs) {  # c(0, 1, 5)) { 
          p <- ovlap.results.long %>% 
            filter(type %in% c("detection", "FDR"),
                   tissue == this.TISSUE,
                   FC == this.FC.CUTOFF,
                   FCtype == this.FCtype,
                   alphaType == 'std.05', # affects things very slightly, other than for P=1
                   tpm.cutoff == this.TPM.CUTOFF) %>%
            mutate(type = gsub("detection", "Sensitivity", type)) %>%
            mutate(type = gsub("FDR", "False Discovery Rate", type)) %>%
            ggplot(aes(x = N, y = value)) + 
            facet_grid(pval.cutoff ~ type, labeller = labeller(pval.cutoff = PVALfunc, tpm.cutoff = TPMfunc)) +  
            ylab("") + 
            xlab(paste0("Sample Size")) + 
            theme(panel.grid.major = element_blank(), # leave blank
                  panel.grid.minor = element_blank(), # leave blank
                  axis.line = element_line(colour = "black"),
                  legend.position = 'none') + 
            scale_x_discrete(labels=PRETTY_LABELS)
          
          p1 <- pub_boxplot_theme(p, 'png')
          p2 <- pub_boxplot_theme(p, 'tiff')
          
          outfile <- paste0(fig_str, '_sensitivity.FDR_TPM', this.TPM.CUTOFF, '_FC', this.FC.CUTOFF, '_', PTYPE, '_', this.TISSUE, '_', this.FCtype)
          #openPng(outfile, p = p1, height = 900*3, width = 1200*3, res = thisres)
          tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
      }
    }
  } # tissues
  } # FCtype
}
plot_2C(ovlap.results.long, 'Figure2C')



#############################
### FIGURE 2D - vary TPM ####
# & Supp Fig 1C ####
#############################
plot_2D <- function(ovlap.results.long, fig_str = 'Figure2D',
                    FCtypes = c('default'), #'apeglm', 'ashr'),
                    FCs = 1.5, 
                    PVAL.CUTOFFs = c(0.05)) {  #c(0.01, 0.05, 0.1, 1)) {
  for (this.FCtype in FCtypes) {
    for (this.tissue in TISSUES_IN_ORDER) {
      for (this.FC.CUTOFF in FCs) { 
        for (this.PVAL.CUTOFF in PVAL.CUTOFFs) { 
          p <- ovlap.results.long %>% 
            filter(pval.cutoff == this.PVAL.CUTOFF,
                   alphaType == 'std.05', # shouldn't matter
                   FCtype == this.FCtype,
                   FC == this.FC.CUTOFF,
                   tissue == this.tissue,
                   type %in% c("detection", "FDR")) %>%
            mutate(type = gsub("detection", "Sensitivity", type),
                   type = gsub("FDR", "False Discovery Rate", type)) %>% 
            ggplot(aes(x = N, y = value)) + 
            facet_grid(tpm.cutoff ~ type, labeller = labeller(FC = FCfunc, tpm.cutoff = TPMfunc), scales = 'free') +  
            ylab("") + 
            xlab(paste0("Sample Size")) +
            theme(panel.grid.major = element_blank(), # leave blank
                  panel.grid.minor = element_blank(), # leave blank
                  axis.line = element_line(colour = "black"),
                  legend.position = 'none') + 
            scale_x_discrete(labels=PRETTY_LABELS)
            
          p1 <- pub_boxplot_theme(p, 'png')
          p2 <- pub_boxplot_theme(p, 'tiff')

          outfile <- paste0(fig_str, '_sensitivity.FDR_varyTPM_FC', this.FC.CUTOFF, '_', PTYPE, this.PVAL.CUTOFF, '_', this.tissue, '_', this.FCtype)
          #openPng(outfile, p = p1, height = 900*3, width = 1200*3, res = thisres)
          tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
        } # Pval
      } # FC
    } # tissue
  } # FCtype
}
plot_2D(ovlap.results.long, 'Figure2D')




################################
### FIGURE 3 pre-processing ####
################################
# read in statistics from pairwise comparisons
pairwise.results <- readRDS('result_objs/Dchs1_pairwise.results.RDS')


# re-format
rejigger_pwise_res <- function(pairwise.results) {
pw.res.final <- pairwise.results %>%
  mutate(percent.overlapA1 = 100 * ovlap1 / (ovlap1 + uniqA1),
         percent.overlapA1.2 = 100 * ovlap1.2 / (ovlap1.2 + uniqA1.2),
         percent.overlapA1.5 = 100 * ovlap1.5 / (ovlap1.5 + uniqA1.5),
         percent.overlapA2 = 100 * ovlap2 / (ovlap2 + uniqA2),
         percent.overlapB1 = 100 * ovlap1 / (ovlap1 + uniqB1),
         percent.overlapB1.2 = 100 * ovlap1.2 / (ovlap1.2 + uniqB1.2),
         percent.overlapB1.5 = 100 * ovlap1.5 / (ovlap1.5 + uniqB1.5),
         percent.overlapB2 = 100 * ovlap2 / (ovlap2 + uniqB2),
         percent.overlap.min1 = pmin(percent.overlapA1, percent.overlapB1),
         percent.overlap.min1.2 = pmin(percent.overlapA1.2, percent.overlapB1.2),
         percent.overlap.min1.5 = pmin(percent.overlapA1.5, percent.overlapB1.5),
         percent.overlap.min2 = pmin(percent.overlapA2, percent.overlapB2),
         percent.overlap.max1 = pmax(percent.overlapA1, percent.overlapB1),
         percent.overlap.max1.2 = pmax(percent.overlapA1.2, percent.overlapB1.2),
         percent.overlap.max1.5 = pmax(percent.overlapA1.5, percent.overlapB1.5),
         percent.overlap.max2 = pmax(percent.overlapA2, percent.overlapB2))
pw.res.final$N <- factor(pw.res.final$N,
                         levels = as.character(min(pw.res.final$N):max(pw.res.final$N)))

value_start_col <- which(colnames(pw.res.final) == 'ovlap1')
pw.res.final <- pw.res.final %>%
  tidyr::pivot_longer(cols = all_of(value_start_col:ncol(pw.res.final)), names_to = "type") %>%
  mutate(type2 = gsub("\\d+\\S*", "", type),
         type2 = gsub("^log$", "log2fc.cor", type2),
         FC = gsub("^\\.+", "", gsub('[^\\d\\.]', "", type, perl = TRUE)))
return(pw.res.final)
}
pw.res.final <- rejigger_pwise_res(pairwise.results)


######################################
### FIGURE 3 - pairwise/splitting ####
######################################
plot_3 <- function(pw.res.final, fig_str = 'Figure3') {
for (this.FCtype in c('default')) { # }, 'apeglm', 'ashr')) {
  for (this.pval in c(0.05)) { # 0.01, 0.05, 0.1, 1)) {
    for (this.TPM.CUTOFF in c(0,1)) { # c(0, 1, 5) {
      for (this.type2 in c('cohenK', 'cohenK_limTPM')) {
        p <- pw.res.final %>% 
          filter(tpm.cutoff == this.TPM.CUTOFF,
                 pval.cutoff == this.pval,
                 alphaType == 'std.05', # shouldn't matter
                 FCtype == this.FCtype,
                 type2 == this.type2) %>%
          ggplot(aes(x = N, y = value)) + 
          facet_grid(tissue ~ FC, labeller = labeller(FC = FCfunc)) +
          ylab(expression(bold(atop("Agreement between signatures" , "(Cohen's"~kappa*")")))) + 
          xlab(paste0("Sample Size")) +  # , " (TPM cutoff = ", this.TPM.CUTOFF, ")")) + 
          theme(panel.grid.major = element_blank(), # leave blank
                panel.grid.minor = element_blank(), # leave blank
                axis.line = element_line(colour = "black"), 
                legend.position = 'none') +
          scale_x_discrete(labels=PRETTY_LABELS)
        
        p1 <- pub_boxplot_theme(p, 'png', fixedcolor='blue')
        p2 <- pub_boxplot_theme(p, 'tiff', fixedcolor='blue')
        
        # outfile <- paste0(fig_str, '_rangeFCs_allTissues_TPM', this.TPM.CUTOFF, '_', PTYPE, this.pval, '_', this.type2)
        outfile <- paste0(fig_str, '_rangeFCs_allTissues_TPM', this.TPM.CUTOFF, '_', PTYPE, this.pval, '_', this.type2, '_', this.FCtype)
        #openPng(outfile, p = p1, height = 700*3, width = 1200*3, res = thisres)
        tiff(paste0("figures/", outfile, '.tiff'), height = 350*4, width = 500*4+200, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
      }
    } #TPM
  } #this.pval
} # FCtype
}
plot_3(pw.res.final, 'Figure3')

#############################
### FIGURE 4 pre-process ####
#############################
effect.size.results <- readRDS("result_objs/Dchs1_effect.size.results.RDS")
head(effect.size.results)

merged.ovlap.effect <- merge(ovlap.results %>% filter(pval.cutoff == 0.05) %>% select(!contains("ovlap")), 
                             effect.size.results) %>% 
  select(tissue, tpm.cutoff, N, trial, FCtype, matches("signif\\d"), contains("ovlap"), contains("gt.orig")) %>% 
  pivot_longer(cols = contains("gt.orig") | contains("signif") | contains("ovlap"), names_to = "type") %>% 
  mutate(FC = gsub("ovlap", "", gsub("signif", "", gsub("\\S*gt.orig", "", type))),
         type = gsub("\\d\\S*", "", type)) %>% 
  unique() %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(perc = 100 * gt.orig / signif,
         ov.perc = 100 * ov.gt.orig / ovlap,
         text = 'fraction\noverestimated')
nrow(merged.ovlap.effect)
head(merged.ovlap.effect)

################################
### FIGURE 4A - effect size ####
################################

for (gtype in c("perc")) { # "ov.perc"
  for (this.FCtype in c('default', 'apeglm', 'ashr')) {
    for (this.TPM.CUTOFF in c(0,1)) { # c(0, 1, 5)) {
        p <- merged.ovlap.effect %>%
        filter(tpm.cutoff == this.TPM.CUTOFF,
               FCtype == this.FCtype) %>%
        ggplot(aes(x = N, y = !!sym(gtype))) + 
        facet_grid(tissue ~ FC, labeller = labeller(FC = FCfunc)) + 
        xlab(paste0("Sample Size")) + 
        ylab("% of signature genes with fold change\ngreater than in N30 experiment") +
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"), 
              legend.position = 'none') +
        scale_x_discrete(labels=PRETTY_LABELS)
      
      p1 <- pub_boxplot_theme(p, 'png', fixedcolor='forestgreen') + 
        theme(plot.title = element_text(face = 'bold', colour = 'black', size = this4b.par$plot.title['png'])) + 
        geom_hline(yintercept = 50, col = "black", lwd = this4b.par$lwd['png'], linetype = 'dashed')
    
      p2 <- pub_boxplot_theme(p, 'tiff', fixedcolor='forestgreen') + 
        theme(plot.title = element_text(face = 'bold', colour = 'black', size = this4b.par$plot.title['tiff'])) + 
        geom_hline(yintercept = 50, col = "black", lwd = this4b.par$lwd['tiff'], linetype = 'dashed')
      
      outfile <- paste0("Figure4A_effect_size_proportion_alltissues_PVAL0.05_TPM", this.TPM.CUTOFF, '_', this.FCtype, '_', gtype)
      #openPng(outfile, p = p1, height = 900*3, width = 1200*3*1.6, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4*2, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()    
    } 
  } #FCtype
} # gtype


################################
### FIGURE 4Av2 - effect size ####
################################
effect.size.results <- readRDS("result_objs/Dchs1_effect.size.results_v2.RDS")
head(effect.size.results)
trial.FC <- 1    # NO THRESHOLD on mini-expt fold changes
stopifnot(identical(effect.size.results %>% mutate(x = gt.orig2+lt.orig2+ns.gt.orig2+ns.lt.orig2+num.NA) %>% pull(x), effect.size.results$N30.sig))
stopifnot(identical(effect.size.results %>% mutate(x = gt.orig1+lt.orig1+ns.gt.orig1+ns.lt.orig1+num.NA) %>% pull(x), effect.size.results$N30.sig))


merged.ovlap.effect <- effect.size.results %>% 
  select(tissue, goldstandard.FC, tpm.cutoff, N, trial, FCtype, contains("ovlap"), contains("gt.orig"), N30.sig, num.NA) %>%  # don't actually need: matches("signif\\d"), 
  pivot_longer(cols = contains("gt.orig") | contains("ovlap"), names_to = "type") %>%    # don't need: contains("signif"), 
  mutate(FC = gsub("ovlap", "", gsub("signif", "", gsub("\\S*gt.orig", "", type))),
         type = gsub("\\d\\S*", "", type)) %>% 
  unique() %>%
  filter(FC == trial.FC) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(perc.sig = 100 * gt.orig / ovlap,
         perc.ns = 100 * ns.gt.orig / (N30.sig - ovlap - num.NA),
         text = 'fraction\noverestimated')
nrow(merged.ovlap.effect)
head(merged.ovlap.effect)

for (gtype in c("perc.sig", "perc.ns", "both")) {
  for (this.FCtype in c('default')) { #}, 'apeglm', 'ashr')) {
    for (this.TPM.CUTOFF in c(0)) { # c(0, 1, 5)) {
      tmp <- merged.ovlap.effect %>%
        filter(tpm.cutoff == this.TPM.CUTOFF,
               FCtype == this.FCtype)
      
      if (gtype == 'perc.ns') {
        p <- tmp %>% ggplot(aes(x = N, y = perc.ns))
      } else {
        p <- tmp %>% ggplot(aes(x = N, y = perc.sig))
      }
      p <- p + 
        facet_grid(tissue ~ goldstandard.FC, labeller = labeller(goldstandard.FC = FCfunc)) + 
        xlab(paste0("Sample Size")) + 
        ylab("% of gold standard genes whose fold change\nis greater in sub-sampled experiment") +
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"), 
              legend.position = 'none') +
        scale_x_discrete(labels=PRETTY_LABELS)
      
      p1 <- pub_boxplot_theme(p, 'png', fixedcolor='forestgreen') + 
        theme(plot.title = element_text(face = 'bold', colour = 'black', size = this4b.par$plot.title['png'])) + 
        geom_hline(yintercept = 50, col = "black", lwd = this4b.par$lwd['png'], linetype = 'dashed')
      if (gtype == 'both') {
        p1 <- p1 + geom_boxplot(aes(y = perc.ns), lwd = usual.par$lwd['png'], outlier.size = usual.par$outlier.size['png'], color = 'honeydew3', alpha=0.1, fatten = 0.85)
      }
      
      p2 <- pub_boxplot_theme(p, 'tiff', fixedcolor='forestgreen') +
        theme(plot.title = element_text(face = 'bold', colour = 'black', size = this4b.par$plot.title['tiff'])) +
        geom_hline(yintercept = 50, col = "black", lwd = this4b.par$lwd['tiff'], linetype = 'dashed')
      if (gtype == 'both') {
        p2 <- p2 + geom_boxplot(aes(y = perc.ns), lwd = usual.par$lwd['tiff'], outlier.size = usual.par$outlier.size['tiff'], color = 'honeydew3', alpha=0.1, fatten = 0.85)
      }

      outfile <- paste0("Figure4Av2_effect_size_PVAL0.05_TPM", this.TPM.CUTOFF, '_', this.FCtype, "FC", trial.FC, '_', gtype)
      #openPng(outfile, p = p1, height = 900*3, width = 1200*3*1.6, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4*2, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()    
    }
  } #FCtype
} # gtype

###################


###################
### FIGURE 4B #######
### example gene ####
###################
for (USE_ORIG_GENES_FOR_S4 in c(TRUE)) {  
  if (USE_ORIG_GENES_FOR_S4) {
    genewise.effect.size.res <- readRDS("result_objs/Dchs1_genewise.effect.size.res_ORIG_SUBMISSION.RDS")
    this.gene <- 'Trex1'
    } else {
    genewise.effect.size.res <- readRDS("result_objs/Dchs1_genewise.effect.size.res.RDS") 
    this.gene <- genewise.effect.size.res %>% filter(tissue=='liver') %>% pull(gene) %>% head(1)
  }
  this.tissue = 'liver'
  this.TPM.CUTOFF = 0
  head(genewise.effect.size.res)
  
  for (this.FCtype in c("default", "apeglm", "ashr")) {
    gw_tmp <- genewise.effect.size.res
    if (this.FCtype == 'default') {
      gw_tmp$FC <- gw_tmp$FC
      N30_FCs <- N30_DEres$FC.apeglm
    } else if (this.FCtype == 'ashr') {
      gw_tmp$FC <- gw_tmp$FC.ashr
      N30_FCs <- N30_DEres$FC.ashr
    } else if (this.FCtype == 'ashr') {
      gw_tmp$FC <- gw_tmp[['FC.apeglm']]
      N30_FCs <- N30_DEres$FC.apeglm
    }
    
    str <- paste(this.tissue, this.TPM.CUTOFF, "Padj,std.05", this.FCtype, "0.05,1", sep=",")
    p <- gw_tmp %>% 
      filter(gene == this.gene, tissue == this.tissue, tpm.cutoff == this.TPM.CUTOFF) %>%
      mutate(N = factor(N, levels = min(N):max(N))) %>% 
      mutate(is.signif = ifelse(is.signif, "YES", "NO")) %>% 
      dplyr::rename(`significant trial` = is.signif) %>%
      ggplot(aes(x = N, y = log2(FC))) + 
      ylab(bquote(bold(.(this.gene)~~log[2]*FC))) + 
      xlab("Sample Size") +  
      theme(panel.grid.major = element_blank(), # leave blank
            panel.grid.minor = element_blank(), # leave blank
            axis.line = element_line(colour = "black"), 
            legend.position = c(0.8, 0.2)) + 
      geom_hline(yintercept = log2(N30_FCs[this.gene, this.tissue]),
                 linetype = ifelse(this.gene %in% gold_standards$Dchs1[[str]], "solid", "dashed")) +
      scale_x_discrete(labels=PRETTY_LABELS) 
    
    p1 <- pub_boxplot_theme_4b(p, 'png')
    p2 <- pub_boxplot_theme_4b(p, 'tiff')
    
    outfile <- paste0("Figure4B_", this.gene, "_", this.tissue, "_TPM", this.TPM.CUTOFF, "_", USE_ORIG_GENES_FOR_S4, "_", this.FCtype)
    #openPng(outfile, p = p1, height = 900, width = 1200*1.75, res = thisres)
    tiff(paste0("figures/", outfile, '.tiff'), height = 450*1.5, width = 500*2.5, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()    
  }
}
rm(gw_tmp, str, p1, p2, N30_FCs, USE_ORIG_GENES_FOR_S4, this.tissue, this.gene)

####################################
### SUPPLEMENTARY FIGURE 1A ########
### see Figure 2B section above ####
####################################

####################################
### SUPPLEMENTARY FIGURE 1B ########
### see Figure 2C section above ####
####################################

####################################
### SUPPLEMENTARY FIGURE 1C ########
### see Figure 2D section above ####
####################################

####################################
### SUPPLEMENTARY FIGURE 2A ########
### P(correct) given signature #####
####################################
MIN_NUM_SIGNIF <- 10
prob_correct <- readRDS('result_objs/Dchs1_prob_correct_minSignif10.RDS') 

for (this.tissue in TISSUES_IN_ORDER) {
  for (this.FCtype in c('default')) {    # , 'apeglm', 'ashr')) {
    for (this.goldstandard.FC in c(1, 1.5)) {  # , 1, 1.2, 1.5, 2)) {
      prob_correct %>%
        filter(tissue == this.tissue,
               FCtype == this.FCtype,
               goldstandard.FC == this.goldstandard.FC,
               num.signif > MIN_NUM_SIGNIF) %>%
        pivot_longer(cols = starts_with("fc", ignore.case = FALSE), names_to = 'fc_range', values_to = 'fraction_detected') -> tmp
      tmp$fc_range[tmp$fc_range == 'fc1.2'] <- '1-1.2'
      tmp$fc_range[tmp$fc_range == 'fc1.5'] <- '1.2-1.5'
      tmp$fc_range[tmp$fc_range == 'fc2'] <- '1.5-2'
      tmp$fc_range[tmp$fc_range == 'fc3'] <- '2-3'
      tmp$fc_range[tmp$fc_range == 'fc5'] <- '3-5'
      tmp$fc_range[tmp$fc_range == 'fc10'] <- '5-10'
      tmp$fc_range[tmp$fc_range == 'fcInf'] <- '>10'
      tmp$fc_range <- factor(tmp$fc_range, levels = c('1-1.2', '1.2-1.5', '1.5-2', '2-3', '3-5', '5-10', '>10'))

      # plot
      outfile <- paste0("SuppFig2A_", this.tissue, '_GSFC', this.goldstandard.FC, "_minSignif", MIN_NUM_SIGNIF, '_', this.FCtype)
      p <- tmp %>% 
        mutate(fraction_detected = 100 * fraction_detected,
               N = factor(N, levels = unique(sort(N))),
               trial = factor(trial, levels = unique(rev(sort(trial))))) %>%
        ggplot(aes(x = fc_range, y = fraction_detected)) + 
        xlab("Fold Change in down-sampled trial") + 
        ylab(paste("% of significant genes also detected in\ngold standard signature with fold change >", this.goldstandard.FC)) + 
        facet_wrap(~N) + 
        theme_bw() + embiggen() + 
        theme(legend.position = 'none')
      
      # PNG VERSION
      p1 <- p + 
        geom_boxplot() +
        theme(axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle=60, hjust = 1))
      #openPng(outfile, p = p1, width = 1400, height = 1200)
      
      # TIFF VERSION
      p2 <- p + 
        geom_boxplot(outlier.size = 0.7, lwd=0.35) +  
        theme(axis.text.x = element_text(colour = "black", size = 8, face = "bold", angle=60, hjust = 1),
              axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
              axis.title.x = element_text(face = "bold", size = 13.5), 
              axis.title.y = element_text(face = "bold", size = 14),
              strip.text = element_text(colour = "black", size = 10, face = "bold"))
       tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = 12, res= 300, compression = 'lzw'); plot(p2); dev.off()
    } # goldstandard.FC
  } # FCtype
} # tissue

#
if (FALSE) {
tpc <- prob_correct %>% 
  filter(goldstandard.FC==1.5, tpm.cutoff==0, FCtype=='default') %>% 
  select(-goldstandard.FC, -tpm.cutoff, -FCtype) %>%
  select(!contains('num.in.range')) %>%
  pivot_longer(cols = starts_with("fc", ignore.case = FALSE), names_to = "FC", values_to = "prob") %>%
  mutate(FC = as.numeric(gsub("fc", "", FC))) %>%
  # select(!starts_with('fc', ignore.case=FALSE)) %>% 
  # pivot_longer(cols = starts_with("num.in"), names_to = "num.in.range", values_to = "prob") %>% 
  # mutate(FC = as.numeric(gsub("num.in.range", "", num.in.range))) %>% 
  group_by(tissue, N, FC) %>% summarize(median = median(prob, na.rm=TRUE))
tpc %>% filter(FC == 1.2, N < 10) %>% pivot_wider(names_from='tissue', values_from = 'median')
tpc %>% filter(FC == 1.5, N < 10) %>% pivot_wider(names_from='tissue', values_from = 'median')

tpc %>% filter(FC == 2, N < 20) %>% pivot_wider(names_from='tissue', values_from = 'median')
rm(tpc)
}

###################################
### SUPPLEMENTARY FIGURE 2B #######
### Prob. Detection given sig. ####
################################
MIN_NUM_SIGNIF <- 10
prob_detect_cumul <- readRDS('result_objs/Dchs1_prob_detect_cumul_minSignif10.RDS')

for (this.tissue in TISSUES_IN_ORDER) {
  for (this.FCtype in c('default')) { # 'apeglm', 'ashr')) {
    for (this.goldstandard.FC in c(1, 1.5)) { # 1, 1.2, 1.5, 2)) {
      prob_detect_cumul %>%
        filter(tissue == this.tissue,
               FCtype == this.FCtype,
               goldstandard.FC == this.goldstandard.FC,
               num.signif > MIN_NUM_SIGNIF) %>%
        pivot_longer(cols = starts_with("fc", ignore.case = FALSE), names_to = 'fc_range', values_to = 'fraction_detected') -> tmp
      tmp$fc_range[tmp$fc_range == 'fc1.2'] <- '1-1.2'
      tmp$fc_range[tmp$fc_range == 'fc1.5'] <- '1.2-1.5'
      tmp$fc_range[tmp$fc_range == 'fc2'] <- '1.5-2'
      tmp$fc_range[tmp$fc_range == 'fc3'] <- '2-3'
      tmp$fc_range[tmp$fc_range == 'fc5'] <- '3-5'
      tmp$fc_range[tmp$fc_range == 'fc10'] <- '5-10'
      tmp$fc_range[tmp$fc_range == 'fcInf'] <- '>10'
      tmp$fc_range <- factor(tmp$fc_range, levels = c('1-1.2', '1.2-1.5', '1.5-2', '2-3', '3-5', '5-10', '>10'))

      # plot
      outfile <- paste0("SuppFig2B_", this.tissue, '_GSFC', this.goldstandard.FC, "_minSignif", MIN_NUM_SIGNIF, '_', this.FCtype)
      p <- tmp %>%
        mutate(fraction_detected = 100 * fraction_detected,
               N = factor(N, levels = unique(sort(N))),
               trial = factor(trial, levels = unique(rev(sort(trial))))) %>%
        ggplot(aes(x = fc_range, y = fraction_detected)) +
        xlab("Fold Change in down-sampled trial") + 
        ylab(paste("% of gold standard genes with fold change >", this.goldstandard.FC, "\nalso detected in sub-sampled signature")) + 
        facet_wrap(~N) + 
        theme_bw() + embiggen() + 
        theme(legend.position = 'none')
      
      # PNG VERSION
      p1 <- p + 
        geom_boxplot() +
        theme(axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle=60, hjust = 1))
      #openPng(outfile, p = p1, width = 1400, height = 1200)    
      
      # TIFF VERSION
      p2 <- p + 
        geom_boxplot(outlier.size = 0.7, lwd=0.35) +  #geom_boxplot(stat='identity') + 
        theme(axis.text.x = element_text(colour = "black", size = 8, face = "bold", angle=60, hjust = 1),
              axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
              axis.title.x = element_text(face = "bold", size = 13.5), 
              axis.title.y = element_text(face = "bold", size = 14),
              strip.text = element_text(colour = "black", size = 10, face = "bold"))
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = 12, res= 300, compression = 'lzw'); plot(p2); dev.off()
    }
  }
}


################################
### SUPPLEMENTARY FIGURE 3 ####
### WT-vs-WT comparisons #######
################################
pairwise_WT_comp <- readRDS('result_objs/Dchs1_WTcompare.results.RDS') # padj cutoff = 0.05

start_col <- which(colnames(pairwise_WT_comp) == 'uniq1')
pairwise_WT_comp_long <- pairwise_WT_comp %>%
  tidyr::pivot_longer(cols = all_of(start_col:ncol(pairwise_WT_comp)), names_to = "FC_cutoff") %>%
  mutate(FC_cutoff = gsub("uniq", "", FC_cutoff),
         N = factor(N))
head(pairwise_WT_comp_long)

for (this.FCtype in c('default')) {    # }, 'apeglm', 'ashr')) {
  for (this.FC in c(1, 1.2, 1.5, 2)) {
    p <- pairwise_WT_comp_long %>% 
      filter(FC_cutoff == this.FC,
             FCtype == this.FCtype) %>%
      mutate(tpm.cutoff = paste0('TPM > ', tpm.cutoff)) %>% 
      ggplot(aes(x = N, y = value)) + 
      geom_boxplot(lwd = 1.6, outlier.size = 2) + 
      facet_wrap(tpm.cutoff ~ tissue, scales = "free", ncol = 4, labeller = label_wrap_gen(multi_line=FALSE, width=30)) +
      ylab("Number of genes in WT-vs-WT signature") + 
      xlab("Sample Size") + 
      scale_x_discrete(labels=PRETTY_LABELS) + 
      theme(axis.text = element_text(size = 20, face = "bold"),
            axis.title = element_text(size = 25, face = 'bold'),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            strip.text = element_text(colour="black", size = 22, face = "bold"),
            panel.background = element_rect(fill = "white", colour = NA)); p
    
    outfile <- paste0('SuppFig3_FC', this.FC, "_", PTYPE, '_', this.FCtype)
    #openPng(outfile, p = p, height = 900*3, width = 1200*4, res = thisres)
    tiff(paste0("figures/", outfile, '.tiff'), height = 450*6, width = 500*8, pointsize = ps, res= thisres, compression = 'lzw'); plot(p); dev.off()
    
    p <-p + ylim(c(0,30))
    outfile <- paste0(outfile, "_ylim30")
    #openPng(outfile, p = p, height = 900*3, width = 1200*4, res = thisres)
    tiff(paste0("figures/", outfile, '.tiff'), height = 450*6, width = 500*8, pointsize = ps, res= thisres, compression = 'lzw'); plot(p); dev.off()
  } # FC
} # FCtype




##########################################
### SUPPLEMENTARY FIGURE 4 ###############
### pairwise agreement, vary TPM & FC ####
#######################################
for (this.FCtype in c('default')) { # , 'apeglm', 'ashr')) {
  for (this.tissue in TISSUES_IN_ORDER) {
    for (this.pval in c(0.05)) { # c(0.01, 0.05, 0.1, 1)) {
      p <- pw.res.final %>% 
        filter(tissue == this.tissue, 
               FCtype == this.FCtype,
               pval.cutoff == this.pval,
               alphaType == 'std.05', # shouldn't matter
               type2 == 'cohenK') %>%
        ggplot(aes(x = N, y = value)) + 
        facet_grid(tpm.cutoff ~ FC, labeller = labeller(FC = FCfunc, tpm.cutoff = TPMfunc)) + 
        ylab(expression(bold(atop("Agreement between signatures" , "(Cohen's"~kappa*")")))) + 
        xlab("Sample Size") +
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"), 
              legend.position = 'none') +
        scale_x_discrete(labels=PRETTY_LABELS)
      
      p1 <- pub_boxplot_theme(p, 'png', fixedcolor='blue')
      p2 <- pub_boxplot_theme(p, 'tiff', fixedcolor='blue')
      outfile <- paste0("SuppFig4_rangeFC.TPM_", this.tissue, "_PVAL", this.pval, "_cohenK_", this.FCtype)
      
      #openPng(outfile, p = p1, height = 700*3*6/7, width = 1200*3, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 350*4*6/7, width = 500*4+200, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    }  # pval
  } #tissue
} # FCtype


######################################
### SUPPLEMENTARY FIGURE 5 ###########
### multiple effect size examples ####
####################################
for (USE_ORIG_GENES_FOR_S4 in c(TRUE)) { 
  if (USE_ORIG_GENES_FOR_S4) {
    genes.for.effect.size.exploration <- readRDS("result_objs/Dchs1_genes.for.effect.size.exploration_ORIG_SUBMISSION.RDS")
    genewise.effect.size.res <- readRDS("result_objs/Dchs1_genewise.effect.size.res_ORIG_SUBMISSION.RDS")
  } else {
    genes.for.effect.size.exploration <- readRDS("result_objs/Dchs1_genes.for.effect.size.exploration.RDS")
    genewise.effect.size.res <- readRDS("result_objs/Dchs1_genewise.effect.size.res.RDS")
  }
  this.TPM.CUTOFF = 0
  
  for (this.FCtype in c("default", "apeglm", "ashr")) {
    gw_tmp <- genewise.effect.size.res
    if (this.FCtype == 'default') {
      gw_tmp$FC <- gw_tmp$FC
      N30_FCs <- N30_DEres$FC.apeglm
    } else if (this.FCtype == 'ashr') {
      gw_tmp$FC <- gw_tmp$FC.ashr
      N30_FCs <- N30_DEres$FC.ashr
    } else if (this.FCtype == 'ashr') {
      gw_tmp$FC <- gw_tmp[['FC.apeglm']]
      N30_FCs <- N30_DEres$FC.apeglm
    }

    for (this.tissue in TISSUES_IN_ORDER) {
      # get fold changes from N30 experiment
      this_genes <- genes.for.effect.size.exploration[[this.tissue]]
      tmp.N30.fcs <- data.frame(
        gene = this_genes,
        fc.N30 = log2(N30_FCs[this_genes, this.tissue])) %>%
        mutate(abs.fc.N30 = abs(fc.N30),
               tissue = this.tissue) %>% 
        arrange(desc(abs.fc.N30))
      genes_in_order <- tmp.N30.fcs$gene
      
      # prepare effect size table (including above FCs), filter to tissue, etc.
      tmp0 <- left_join(tmp.N30.fcs, gw_tmp) %>%
        filter(tpm.cutoff == this.TPM.CUTOFF)
      
      tmp1 <- tmp0 %>% group_by(gene) %>% summarize(low = min(meanX.TPM, na.rm=T), 
                                                    median = median(meanX.TPM, na.rm=T), 
                                                    high = max(meanX.TPM, na.rm=T))
      
      tmp2 <- left_join(tmp0, tmp1) %>% 
        mutate(gene = factor(gene, levels = genes_in_order)) %>%
        mutate(titlelabel = paste0(gene, " (", round(median, 1), ")")) %>%  # CHANGE TO only show genes!!
        mutate(titlelabel = factor(titlelabel, levels = unique(titlelabel))) %>%
        mutate(`sample size` = factor(N))
      tmp2$is.signif[is.na(tmp2$is.signif)] <- FALSE # NOTE: treating NA values as NOT significant
      
      p <- tmp2 %>% 
        filter(!is.na(is.signif)) %>%        
        mutate(is.signif = gsub("TRUE", "P < 0.05", gsub("FALSE", "P \u2265 0.05", is.signif))) %>%
        mutate(is.signif = factor(is.signif, levels = c("P \u2265 0.05", "P < 0.05"))) %>%
        dplyr::rename(`P < 0.05` = is.signif) %>%
        ggplot(aes(x = `sample size`, y = log2(FC))) + 
        xlab('Sample Size') + 
        ylab(bquote(bold(log[2]*FC))) + 
        theme(axis.text.x = element_blank(),
              panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              legend.title = element_blank()) + 
        facet_wrap(~titlelabel, scales = "free_y", ncol = 5) + 
        scale_x_discrete(labels=c("4"="", "5"="", "7"="", "8"="", "9"="", "12"="", "14"="", "11" = "", "13" = "", "16" = "", "17" = "", "18"="", "19" = "", "21" = "", "22" = "", "23" = "", "24"="", "26"="",  "27" = "", "28" = ""))
      
      
      p1 <- p + geom_boxplot(lwd = S4.par$lwd['png'], outlier.size = S4.par$outlier.size['png'], aes(color = `P < 0.05`)) +
        geom_hline(data = tmp.N30.fcs %>% left_join(tmp2), aes(yintercept = fc.N30), colour = "black", size = S4.par$lwd['png'], linetype = 'dashed') + 
        theme(axis.text.x = element_text(colour = "black", size = S4.par$axis.text.x['png'], face = "bold"), 
              axis.text.y = element_text(colour = "black", size = S4.par$axis.text.y['png'], face = "bold"), 
              axis.title.x = element_text(face = "bold", size = S4.par$axis.title.x['png']), 
              axis.title.y = element_text(face = "bold", size = S4.par$axis.title.y['png']), 
              strip.text = element_text(colour = "black", size = S4.par$strip.text['png'], face = "bold"))
      
      p2 <- p + geom_boxplot(lwd = S4.par$lwd['tiff'], outlier.size = S4.par$outlier.size['tiff'], aes(color = `P < 0.05`)) +
        geom_hline(data = tmp.N30.fcs %>% left_join(tmp2), aes(yintercept = fc.N30), colour = "black", size = S4.par$lwd['tiff'], linetype = 'dashed') + 
        theme(axis.text.x = element_text(colour = "black", size = S4.par$axis.text.x['tiff'], face = "bold"), 
              axis.text.y = element_text(colour = "black", size = S4.par$axis.text.y['tiff'], face = "bold"), 
              axis.title.x = element_text(face = "bold", size = S4.par$axis.title.x['tiff']), 
              axis.title.y = element_text(face = "bold", size = S4.par$axis.title.y['tiff']), 
              strip.text = element_text(colour = "black", size = S4.par$strip.text['tiff'], face = "bold"))
      
      outfile <- paste0("SuppFig5_manygenes_", this.tissue, "_tpm", this.TPM.CUTOFF, "_", USE_ORIG_GENES_FOR_S4, "_", this.FCtype)
      #openPng(outfile, p = p1, height = 900*3, width = 1200*3*1.6, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4*2, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    } # tissue
  } # FCtype
} # USE_ORIG

###########################################
### SUPPLEMENTARY FIGURE 6 ################
# See Figure 4 - Effect size  #############
# (adding TPM filter (A) or shrinkage (B) #
# from pairwise experiments ###############
###########################################
################################
#### ALTERNATIVE SUPP FIG 7 ####
################################
# log2FC
for (this.FCtype in c('default')) { # , 'apeglm', 'ashr')) {
  for (this.TPM.CUTOFF in c(0)) { # c(0,1,5)) { 
    for (this.PVAL.CUTOFF in c(0.05)) {  #  c(0.01, 0.05, 0.1, 1)) {
      tmp <- ovlap.results.long %>% 
        filter(pval.cutoff == this.PVAL.CUTOFF,
               tpm.cutoff == this.TPM.CUTOFF,
               FCtype == this.FCtype,
               type == 'log2fc.cor')
      tmp2 <- pairwise_WT_comp_long %>% 
        mutate(type = FC_cutoff) %>%
        filter(tpm.cutoff == this.TPM.CUTOFF,
               FCtype == this.FCtype,
               type == 'log2fc.cor') %>%
        select(-FC_cutoff)
      
      p <- tmp[,colnames(tmp2)] %>% 
        ggplot(aes(x = N, y = value)) + 
        facet_grid(tissue ~ type) + 
        ylab(expression(bold("Correlation ("~log[2]~"FC) with gold standard"))) + 
        xlab(paste0("Sample Size")) +
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"),
              legend.position = 'none') + 
        scale_x_discrete(labels=PRETTY_LABELS) + 
        geom_hline(yintercept=0, linetype='dashed')
      
      outfile <- paste0("SuppFig7_FC_crlxns_TPM", this.TPM.CUTOFF, "_", PTYPE, this.PVAL.CUTOFF, '_', this.FCtype)
      
      p1 <- pub_boxplot_theme(p, 'png') + theme(strip.text.x = element_blank()) + geom_boxplot(data = tmp2, color='#00BFC4')
      p2 <- pub_boxplot_theme(p, 'tiff') + theme(strip.text.x = element_blank()) + geom_boxplot(data = tmp2, color='#00BFC4')
      #openPng(outfile, p = p1, height = 900*3, width = 1200*1.5, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*2.4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    }
  }
} # FCtype

# Wald
for (this.FCtype in c('default')) { # , 'apeglm', 'ashr')) {
  for (this.TPM.CUTOFF in c(0)) { # c(0,1,5)) { 
    for (this.PVAL.CUTOFF in c(0.05)) {  #  c(0.01, 0.05, 0.1, 1)) {
      tmp <- ovlap.results.long %>% 
        filter(pval.cutoff == this.PVAL.CUTOFF,
               tpm.cutoff == this.TPM.CUTOFF,
               FCtype == this.FCtype,
               type == 'wald.cor')
      tmp2 <- pairwise_WT_comp_long %>% 
        mutate(type = FC_cutoff) %>%
        filter(tpm.cutoff == this.TPM.CUTOFF,
               FCtype == this.FCtype,
               type == 'wald.cor') %>%
        select(-FC_cutoff)
      
      p <- tmp[,colnames(tmp2)] %>% 
        ggplot(aes(x = N, y = value)) + 
        facet_grid(tissue ~ type) + 
        ylab(expression(bold("Correlation (Wald Stat) with gold standard"))) + 
        xlab(paste0("Sample Size")) +
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"),
              legend.position = 'none') + 
        scale_x_discrete(labels=PRETTY_LABELS) + 
        geom_hline(yintercept=0, linetype='dashed')

      outfile <- paste0("SuppFig7_Wald_crlxns_TPM", this.TPM.CUTOFF, "_", PTYPE, this.PVAL.CUTOFF, '_', this.FCtype)
      
      p1 <- pub_boxplot_theme(p, 'png') + theme(strip.text.x = element_blank()) + geom_boxplot(data = tmp2, color='#00BFC4')
      p2 <- pub_boxplot_theme(p, 'tiff') + theme(strip.text.x = element_blank()) + geom_boxplot(data = tmp2, color='#00BFC4')
      #openPng(outfile, p = p1, height = 900*3, width = 1200*1.5, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*2.4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    }
  }
} # FCtype


# LEGEND for SUPP FIG. 7
if (FALSE) {
  this.FCtype <- 'default'
  this.TPM.CUTOFF <- 0
  this.PVAL.CUTOFF <- 0.05
  tmp <- ovlap.results.long %>% 
    filter(pval.cutoff == this.PVAL.CUTOFF,
           tpm.cutoff == this.TPM.CUTOFF,
           FCtype == this.FCtype,
           type == 'log2fc.cor') %>%
    mutate(type = 'sub-sampled')
  tmp2 <- pairwise_WT_comp_long %>% 
    filter(tpm.cutoff == this.TPM.CUTOFF,
           FCtype == this.FCtype) %>%
    mutate(type = 'WT-vs-WT') %>%
    select(-FC_cutoff)
  p <- rbind(tmp[,colnames(tmp2)], tmp2) %>% ggplot(aes(x = N, y = value)) + 
    facet_grid(tissue ~ type) + 
    ylab(expression(bold(atop("Correlation ("~log[2]~"fold changes) with gold standard")))) + 
    xlab(paste0("Sample Size")) +
    theme(panel.grid.major = element_blank(), # leave blank
          panel.grid.minor = element_blank(), # leave blank
          axis.line = element_line(colour = "black")) + 
    scale_x_discrete(labels=PRETTY_LABELS) + 
    geom_boxplot(data = tmp2, color='#00BFC4')
  p1 <- pub_boxplot_theme(p, 'png') + theme(strip.text.x = element_blank()) + 
    theme(legend.title = element_blank(), legend.text = element_text(face = 'bold', size=20), legend.key.size = unit(1.2, units='cm')); plot(p1)
}

##################################
### SUPPLEMENTARY FIGURE 7(A) ####
# Crlxn between log2fcs ##########
# from down-sampled vs N30 #######
##################################
for (this.FCtype in c('default')) { # , 'apeglm', 'ashr')) {
  for (this.TPM.CUTOFF in c(0,1)) { # c(0,1,5)) { 
    for (this.PVAL.CUTOFF in c(0.05)) {  #  c(0.01, 0.05, 0.1, 1)) {
      p <- ovlap.results.long %>% 
        filter(pval.cutoff == this.PVAL.CUTOFF,
               tpm.cutoff == this.TPM.CUTOFF,
               FCtype == this.FCtype,
               type == 'log2fc.cor') %>% 
        ggplot(aes(x = N, y = value)) + 
        facet_grid(tissue ~ type) + 
        ylab(expression(bold(atop("Correlation between"~log[2]~"fold changes from",
                                  "gold standard vs sub-sampled experiment")))) + 
        xlab(paste0("Sample Size")) +
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"),
              legend.position = 'none') + 
        scale_x_discrete(labels=PRETTY_LABELS)
      
      outfile <- paste0("SuppFig7A_FC_crlxns_TPM", this.TPM.CUTOFF, "_", PTYPE, this.PVAL.CUTOFF, '_', this.FCtype)
  
      p1 <- pub_boxplot_theme(p, 'png') + theme(strip.text.x = element_blank())
      p2 <- pub_boxplot_theme(p, 'tiff') + theme(strip.text.x = element_blank())
      
      #openPng(outfile, p = p1, height = 900*3, width = 1200*1.5, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*2.4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    }
  }
} # FCtype


################################
### SUPPLEMENTARY FIGURE 7B ####
# Crlxn between Wald Stats #####
# from down-sampled vs N30 #####
################################
for (this.FCtype in c('default')) { # }, 'apeglm', 'ashr')) {
  for (this.TPM.CUTOFF in c(0,1)) { # c(0,1,5)) { 
    for (this.PVAL.CUTOFF in c(0.05)) {  #  c(0.01, 0.05, 0.1, 1)) {
      p <- ovlap.results.long %>% 
        filter(pval.cutoff == this.PVAL.CUTOFF,
               tpm.cutoff == this.TPM.CUTOFF,
               FCtype == this.FCtype,
               type == 'wald.cor') %>% 
        ggplot(aes(x = N, y = value)) + 
        facet_grid(tissue ~ type) + 
        ylab(expression(bold(atop("Correlation between Wald statistics from",
                                  "gold standard vs sub-sampled experiment")))) + 
        xlab(paste0("Sample Size")) +
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"),
              legend.position = 'none') + 
        scale_x_discrete(labels=PRETTY_LABELS)
      
      outfile <- paste0("SuppFig7B_Wald_crlxns_TPM", this.TPM.CUTOFF, "_", PTYPE, this.PVAL.CUTOFF, '_', this.FCtype)
      
      p1 <- pub_boxplot_theme(p, 'png') + theme(strip.text.x = element_blank())
      p2 <- pub_boxplot_theme(p, 'tiff') + theme(strip.text.x = element_blank())
      
      #openPng(outfile, p = p1, height = 900*3, width = 1200*1.5, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*2.4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    }
  }
} # FCtype
#####################################
### SUPPLEMENTARY FIGURE 7C/D  #####
### Crlxn between WT-vs-WT expts ####
### log2fc only
#####################################
for (this.FCtype in c('default')) { # }, 'apeglm', 'ashr')) {
  for (this.TPM.CUTOFF in c(0,1)) { # c(0,1,5)) { 
    p <- pairwise_WT_comp_long %>% 
      mutate(type = FC_cutoff) %>%
      filter(tpm.cutoff == this.TPM.CUTOFF,
             FCtype == this.FCtype,
             type == 'wald.cor') %>% 
      ggplot(aes(x = N, y = value)) + 
      facet_grid(tissue ~ type) + 
      ylab(expression(bold(atop("Correlation between Wald statistics from",
                                "gold standard vs WT-vs-WT experiment")))) + 
      xlab(paste0("Sample Size")) +
      theme(panel.grid.major = element_blank(), # leave blank
            panel.grid.minor = element_blank(), # leave blank
            axis.line = element_line(colour = "black"),
            legend.position = 'none') + 
      scale_x_discrete(labels=PRETTY_LABELS)
    
    outfile <- paste0("SuppFig7D_Wald_crlxns_TPM", this.TPM.CUTOFF, "_Padj0.05_", this.FCtype)
    p1 <- pub_boxplot_theme(p, 'png') + theme(strip.text.x = element_blank())
    p2 <- pub_boxplot_theme(p, 'tiff') + theme(strip.text.x = element_blank())
    
    #openPng(outfile, p = p1, height = 900*3, width = 1200*1.5, res = thisres)
    tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*2.4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
  }
} # FCtype
#
for (this.FCtype in c('default')) {  # }, 'apeglm', 'ashr')) {
  for (this.TPM.CUTOFF in c(0,1)) {  # c(0,1,5)) { 
    p <- pairwise_WT_comp_long %>% 
      mutate(type = FC_cutoff) %>%
      filter(tpm.cutoff == this.TPM.CUTOFF,
             FCtype == this.FCtype,
             type == 'log2fc.cor') %>% 
      ggplot(aes(x = N, y = value)) + 
      facet_grid(tissue ~ type) + 
      ylab(expression(bold(atop("Correlation between"~log[2]~"fold changes from",
                                "gold standard vs WT-vs-WT experiment")))) + 
      xlab(paste0("Sample Size")) +
      theme(panel.grid.major = element_blank(), # leave blank
            panel.grid.minor = element_blank(), # leave blank
            axis.line = element_line(colour = "black"),
            legend.position = 'none') + 
      scale_x_discrete(labels=PRETTY_LABELS)
    
    outfile <- paste0("SuppFig7C_log2FC_crlxns_TPM", this.TPM.CUTOFF, "_Padj0.05_", this.FCtype)
    p1 <- pub_boxplot_theme(p, 'png') + theme(strip.text.x = element_blank())
    p2 <- pub_boxplot_theme(p, 'tiff') + theme(strip.text.x = element_blank())
    
    #openPng(outfile, p = p1, height = 900*3, width = 1200*1.5, res = thisres)
    tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*2.4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
  }
} # FCtype


###################################
### SUPPLEMENTARY FIGURE 8(A) ####
# Crlxn between log2fcs     #######
# from pairwise experiments #######
###################################
for (this.FCtype in c('default')) { # , 'apeglm', 'ashr')) {
  for (this.TPM.CUTOFF in c(0,1)) { # c(0, 1, 5)) {
    for (this.PVAL.CUTOFF in c(0.05)) { # c(0.01, 0.05, 0.1, 1)) {
      p <- pw.res.final %>% 
        filter(tpm.cutoff == this.TPM.CUTOFF,
               pval.cutoff == this.PVAL.CUTOFF,
               FCtype == this.FCtype,
               type2 == 'log2fc.cor') %>%
        ggplot(aes(x = N, y = value)) + 
        facet_grid(tissue ~ type) + 
        ylab(expression(bold(log[2]~"fold change correlations"))) + 
        xlab("Sample Size") +
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"), 
              legend.position = 'none') +
        scale_x_discrete(labels=PRETTY_LABELS) + 
        geom_hline(yintercept=0, linetype='dashed')
      
      outfile <- paste0("SuppFig8A_pairwise_crlxns_TPM", this.TPM.CUTOFF, "_", PTYPE, this.PVAL.CUTOFF, '_', this.FCtype)
      
      p1 <- pub_boxplot_theme(p, 'png', fixedcolor = 'blue') + 
        theme(strip.text.x = element_blank(),
              plot.title = element_text(face = 'bold', colour = 'black', size = usual.par$plot.title['png']))
      
      p2 <- pub_boxplot_theme(p, 'tiff', fixedcolor = 'blue') + 
        theme(strip.text.x = element_blank(),
              plot.title = element_text(face = 'bold', colour = 'black', size = usual.par$plot.title['tiff']))
      
      #openPng(outfile, p = p1, height = 900*3, width = 1200*1.5, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*2.4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    } #pval
  } # TPM
} # FCtype


###################################
### SUPPLEMENTARY FIGURE 8(B) ####
# Crlxn between Wald stats  #######
# from pairwise experiments #######
###################################
for (this.FCtype in c('default')) { # }, 'apeglm', 'ashr')) {
  for (this.TPM.CUTOFF in c(0,1)) { # c(0, 1, 5)) {
    for (this.PVAL.CUTOFF in c(0.05)) { # c(0.01, 0.05, 0.1, 1)) {
      p <- pw.res.final %>% 
        filter(tpm.cutoff == this.TPM.CUTOFF,
               pval.cutoff == this.PVAL.CUTOFF,
               FCtype == this.FCtype,
               type2 == 'wald.cor') %>%
        ggplot(aes(x = N, y = value)) + 
        facet_grid(tissue ~ type) + 
        ylab(expression(bold("Wald statistic correlations"))) + 
        xlab("Sample Size") +
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"), 
              legend.position = 'none') +
        scale_x_discrete(labels=PRETTY_LABELS) + 
        geom_hline(yintercept=0, linetype='dashed')
      
      outfile <- paste0("SuppFig8B_wald_pairwise_crlxns_TPM", this.TPM.CUTOFF, "_", PTYPE, this.PVAL.CUTOFF, "_", this.FCtype)
      
      p1 <- pub_boxplot_theme(p, 'png', fixedcolor = 'blue') + 
        theme(strip.text.x = element_blank(),
              plot.title = element_text(face = 'bold', colour = 'black', size = usual.par$plot.title['png']))
      
      p2 <- pub_boxplot_theme(p, 'tiff', fixedcolor = 'blue') + 
        theme(strip.text.x = element_blank(),
              plot.title = element_text(face = 'bold', colour = 'black', size = usual.par$plot.title['tiff']))
      
      #openPng(outfile, p = p1, height = 900*3, width = 1200*1.5, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*2.4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    } # pval
  } # TPM
} # FCtype


##################################
### SUPPLEMENTARY FIGURE 9  #####
# analogous to Figure 2, but #####
# using Pval-based thresholds ####
# (not FC)                    #### 
###############################
# create alternative ovlap.results, which uses Schurch criteria
# ovlap.results_ALT <- ovlap.results %>% 
#   select(tissue, pval.cutoff, tpm.cutoff, alphaType, FCtype, N, trial, contains("_alt")) %>%
#   select(!contains("altb"))
# # re-format
# colnames(ovlap.results_ALT) <- gsub("_alt", "", colnames(ovlap.results_ALT))
# ovlap.results_ALT.long <- rejigger_ovlap_res(ovlap.results_ALT)
# table(ovlap.results_ALT.long$FC, ovlap.results_ALT.long$type, useNA = 'always')
# # plot 2A/B/C
# plot_2A(ovlap.results_ALT.long, 'SuppFig9A')
# plot_2B(ovlap.results_ALT.long, 'SuppFig9B')
# plot_2C(ovlap.results_ALT.long, 'SuppFig9C')

ovlap.results_ALT <- ovlap.results %>% 
  select(tissue, pval.cutoff, tpm.cutoff, alphaType, FCtype, N, trial, contains("_alt")) %>%
  select(!matches("detection_alt\\d"))
# re-format
colnames(ovlap.results_ALT) <- gsub("_altb?", "", colnames(ovlap.results_ALT))
ovlap.results_ALT.long <- rejigger_ovlap_res(ovlap.results_ALT)
table(ovlap.results_ALT.long$FC, ovlap.results_ALT.long$type, useNA = 'always')

# plot 2A/B/C
plot_2A(ovlap.results_ALT.long, 'SuppFig9A', FCs = 1.5, TPM.CUTOFFs = c(0))
# plot_2B(ovlap.results_ALT.long, 'SuppFig9B', TPM.CUTOFFs = 1) # no reason to vary fold changes, as they are ignored
# plot_2C(ovlap.results_ALT.long, 'SuppFig9B', FCs = 1.5, TPM.CUTOFFs = c(0)) # minor diffs
# plot_2D(ovlap.results_ALT.long, 'SuppFig9C', FCs = 1.5) # minor diffs

###################################
### SUPPLEMENTARY FIGURE 10A #######
### (Schurch et al. Fig 1B) #######
###################################
# replace_res <- readRDS(file = "result_objs/Dchs1_FP_FN_etc_stats_replaceNAs.RDS")
remove_res <- readRDS(file = "result_objs/Dchs1_FP_FN_etc_stats_removeNAs.RDS")

for (this.FCtype in 'default') { #  unique(remove_res$FCtype)) {
  for (tpm.cutoff in c(0,1)) { # unique(remove_res$tpm)) {
    tmp <- remove_res %>% 
      filter(FCtype == this.FCtype, tpm == tpm.cutoff) %>%
      mutate(N = factor(N, levels = sort(unique(as.numeric(N))))) %>% 
      filter(grepl("_TPR", type) | grepl("_FPR", type)) %>%
      pivot_wider(names_from = 'type', values_from = 'value')
    
    tmp %>%
      ggplot(aes(x = N)) + 
      geom_line(data = tmp %>% filter(fc==1), aes(y = mean_TPR, group = 1), color = "slategrey", size = 1) + # grey
      geom_line(data = tmp %>% filter(fc==1.2), aes(y = mean_TPR, group = 1), color = "darkgreen", size = 1) + # green4
      geom_line(data = tmp %>% filter(fc==1.5), aes(y = mean_TPR, group = 1), color = "blue", size = 1) +  # lightblue
      geom_line(data = tmp %>% filter(fc==2), aes(y = mean_TPR, group = 1), color = "red", size = 1) +  # tomato
      geom_line(data = tmp %>% filter(fc==4), aes(y = mean_TPR, group = 1), color = "cyan", size = 1) +  # lightcyan3
      geom_ribbon(data = tmp %>% filter(fc==1), aes(ymin = mean_FPR - sd_FPR, ymax = mean_FPR + sd_FPR, group=1), fill = "grey", alpha = 0.4) +  # Shading for ±1 SD
      facet_wrap(~tissue) +
      theme_bw() + 
      embiggen() + 
      scale_x_discrete(labels=PRETTY_LABELS) +
    #  ylim(0,1.1) + 
      xlab(paste0("Sample Size")) +  
      ylab("TPR (FPR)") -> p; p
    
    outfile <- paste0("SuppFig10A_Schurch1B_FCtype_", this.FCtype, "_TPM", tpm.cutoff)
    img_type <- 'png'
    p1 <- p + theme(axis.text.x = element_text(colour = "black", size = usual.par$axis.text.x[img_type], face = "bold"), 
                    axis.text.y = element_text(colour = "black", size = usual.par$axis.text.y[img_type], face = "bold"), 
                    axis.title.x = element_text(face = "bold", size = usual.par$axis.title.x[img_type]), 
                    axis.title.y = element_text(face = "bold", size = usual.par$axis.title.y[img_type]), 
                    strip.text = element_text(colour = "black", size = usual.par$strip.text[img_type], face = "bold"))
    img_type <- 'tiff'
    p2 <- p + theme(axis.text.x = element_text(colour = "black", size = usual.par$axis.text.x[img_type], face = "bold"), 
                    axis.text.y = element_text(colour = "black", size = usual.par$axis.text.y[img_type], face = "bold"), 
                    axis.title.x = element_text(face = "bold", size = usual.par$axis.title.x[img_type]), 
                    axis.title.y = element_text(face = "bold", size = usual.par$axis.title.y[img_type]), 
                    strip.text = element_text(colour = "black", size = usual.par$strip.text[img_type], face = "bold"))
    
   #openPng(outfile, p = p1, res = thisres, height = 700*2, width = 1200*2 + 200)
   tiff(paste0("figures/", outfile, '.tiff'), height = 350*4*3/3, width = 1200*1.8, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
  }
}

###################################
### SUPPLEMENTARY FIGURE 10B ######
### (Schurch et al. but for FDR) ##
###################################
for (this.FCtype in 'default') { #  unique(remove_res$FCtype)) {
  for (tpm.cutoff in c(0,1)) { # unique(remove_res$tpm)) {
    tmp <- remove_res %>% 
      filter(FCtype == this.FCtype, tpm == tpm.cutoff) %>%
      mutate(N = factor(N, levels = sort(unique(as.numeric(N))))) %>% 
      filter(grepl("_TPR", type) | grepl("_FDR", type)) %>% 
      pivot_wider(names_from = 'type', values_from = 'value')
    
    tmp %>%
      ggplot(aes(x = N)) + 
      geom_line(data = tmp %>% filter(fc==1), aes(y = mean_FDR, group = 1), color = "slategrey", size = 1) +
      geom_line(data = tmp %>% filter(fc==1.2), aes(y = mean_FDR, group = 1), color = "darkgreen", size = 1) +
      geom_line(data = tmp %>% filter(fc==1.5), aes(y = mean_FDR, group = 1), color = "blue", size = 1) +
      geom_line(data = tmp %>% filter(fc==2), aes(y = mean_FDR, group = 1), color = "red", size = 1) +
      geom_line(data = tmp %>% filter(fc==4), aes(y = mean_FDR, group = 1), color = "cyan", size = 1) +
      facet_wrap(~tissue) +
      theme_bw() + 
      embiggen() + 
      scale_x_discrete(labels=PRETTY_LABELS) +
      #  ylim(0,1.1) + 
      xlab(paste0("Sample Size")) +  
      ylab("FDR") -> p; p
    
    outfile <- paste0("SuppFig10B_Schurch1B_FCtype_", this.FCtype, "_TPM", tpm.cutoff)
    img_type <- 'png'
    p1 <- p + theme(axis.text.x = element_text(colour = "black", size = usual.par$axis.text.x[img_type], face = "bold"), 
                    axis.text.y = element_text(colour = "black", size = usual.par$axis.text.y[img_type], face = "bold"), 
                    axis.title.x = element_text(face = "bold", size = usual.par$axis.title.x[img_type]), 
                    axis.title.y = element_text(face = "bold", size = usual.par$axis.title.y[img_type]), 
                    strip.text = element_text(colour = "black", size = usual.par$strip.text[img_type], face = "bold"))
    img_type <- 'tiff'
    p2 <- p + theme(axis.text.x = element_text(colour = "black", size = usual.par$axis.text.x[img_type], face = "bold"), 
                    axis.text.y = element_text(colour = "black", size = usual.par$axis.text.y[img_type], face = "bold"), 
                    axis.title.x = element_text(face = "bold", size = usual.par$axis.title.x[img_type]), 
                    axis.title.y = element_text(face = "bold", size = usual.par$axis.title.y[img_type]), 
                    strip.text = element_text(colour = "black", size = usual.par$strip.text[img_type], face = "bold"))
    
    #openPng(outfile, p = p1, res = thisres, height = 700*2, width = 1200*2 + 200)
    tiff(paste0("figures/", outfile, '.tiff'), height = 350*4*3/3, width = 1200*1.8, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
  }
}



###################################
### SUPPLEMENTARY FIGURE 10C  #####
### (Schurch et al. Fig 1A, not included) #######
###################################
# this.PVAL.CUTOFF <- 0.05
# this.TPM.CUTOFF <- 0
# Dchs1_res <- readRDS('full_signature_refs/Dchs1_res.RDS')
# 
# for (this.FCtype in c('default')) { # }, 'apeglm', 'ashr')) { no fc cutoff so doesn't matter
#   p <- ovlap.results %>% 
#     filter(pval.cutoff == this.PVAL.CUTOFF, 
#            tpm.cutoff == this.TPM.CUTOFF,
#            FCtype == this.FCtype) %>%
#     ggplot(aes(x = N, y = signif1 / nrow(Dchs1_res[[1]]) * 100)) + 
#     # geom_boxplot(outlier.size = 0.7, lwd=0.35) +  
#     facet_wrap(~tissue) + ylab("% transcriptome DE") + 
#     embiggen() + 
#     xlab(paste0("Sample Size")) +
#     theme(panel.grid.major = element_blank(), # leave blank
#           panel.grid.minor = element_blank(), # leave blank
#           axis.line = element_line(colour = "black"),
#           legend.position = 'none') + 
#     scale_x_discrete(labels=PRETTY_LABELS); p
#   
#   outfile <- paste0("SuppFig10B_unused_Schurch.Fig1A_TPM", this.TPM.CUTOFF, "_Padj", this.PVAL.CUTOFF, "_", this.FCtype)
#   p1 <- pub_boxplot_theme(p, 'png', fixedcolor = 'slategray3') + theme(strip.text.x = element_blank()) # grey80
#   p2 <- pub_boxplot_theme(p, 'tiff', fixedcolor = 'slategray3') + theme(strip.text.x = element_blank())
#   
#   #openPng(outfile, p = p1, height = 700*2, width = 1200*2 + 200, res = thisres)
# tiff(paste0("figures/", outfile, '.tiff'), height = 350*4*2/3, width = 1200*1.8, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
# } # FCtype
#######################################
### SUPPLEMENTARY FIGURES 11-14 #######
### (Fat4 version of main figures) ####
#######################################



#### Code supporting various statements in the manuscript ####
#### compare_runs_lineplot function ####
if (FALSE) {
  compare_runs_lineplot <- function(str1, str2, smry_stat_type, stats = fig2A_stats) {
    FIELDS <- c('FCtype', 'out_type', 'FC', 'TPM', 'Padj')
    
    one <- stats[[str1]] %>% select(N, contains(smry_stat_type))
    two <- stats[[str2]] %>% select(N, contains(smry_stat_type))
    fields1 <- unlist(strsplit(str1, ","))
    fields2 <- unlist(strsplit(str2, ","))
    stopifnot(length(fields1) == length(fields2))
    discrepant <- fields1 != fields2
    if (sum(discrepant) != 1) {
      print("MULTIPLE (or no) DISCREPANT FIELDS, ARE YOU SURE YOU INTENDED THIS?")
    }
    out_str1 <- gsub("out_type", "", gsub("FCtype", "", paste(paste(FIELDS[discrepant], fields1[discrepant], sep=""), collapse="_")))
    out_str2 <- gsub("out_type", "", gsub("FCtype", "", paste(paste(FIELDS[discrepant], fields2[discrepant], sep=""), collapse="_")))
    shared <- !discrepant & FIELDS != 'out_type'
    shared_str <- paste0("  (", gsub("FCtype=", "", paste(paste(FIELDS[shared], fields2[shared], sep="="), collapse=", ")), ")")
    
    rbind(one %>% pivot_longer(cols= -N, values_to = 'value', names_to='tissue') %>% mutate(tissue=gsub("\\S+_", "", tissue)) %>% mutate(type=out_str1),
          two %>% pivot_longer(cols= -N, values_to = 'value', names_to='tissue') %>% mutate(tissue=gsub("\\S+_", "", tissue)) %>% mutate(type=out_str2)) %>% as.data.frame() %>%
      ggplot(aes(x = N, y = value, color = type)) + geom_line(aes(group=type)) + geom_point() + facet_wrap(~tissue) + 
      ylab(paste(smry_stat_type, fields1[2], shared_str)) + embiggen() + 
      theme(legend.title = element_blank())
  }
}

#### reading exact numbers from Figure 2A and figure 3 ####
if (TRUE) {
  unique(ovlap.results.long$type)
  fig2A_stats <- list()
  for (this.FCtype in c('default', 'ashr', 'apeglm')) {
    for (this.metric in c('FDR', 'detection', "cohenK", "cohenK_limTPM")) {
      for (this.FC in c(1, 1.2, 1.5, 2)) {
        for (this.TPM.CUTOFF in c(0,1)) {
          for (this.PVAL.CUTOFF in c(0.01, 0.05)) {
            if (grepl("cohen", this.metric)) {
              tmp <- pw.res.final %>% filter(type2 == this.metric)
            } else {
              tmp <- ovlap.results.long %>% filter(type == this.metric)
            }
            str <- paste(this.FCtype, this.metric, this.FC, this.TPM.CUTOFF, this.PVAL.CUTOFF, sep=",")
            
            fig2A_stats[[str]] <- tmp %>%
              filter(pval.cutoff==this.PVAL.CUTOFF, tpm.cutoff==this.TPM.CUTOFF, alphaType=='std.05', FCtype==this.FCtype, FC==this.FC) %>%
              select(tissue, N, trial, value) %>%
              group_by(N, tissue) %>% summarize(lower=quantile(value, 0.25), med = median(value), mean=mean(value), upper=quantile(value, 0.75), sd=sd(value)) %>% 
              pivot_wider(names_from='tissue', values_from = c('lower', 'med', 'mean', 'upper', 'sd')) %>% as.data.frame()
          }
        }
      }
    }
  }
  # head(fig2A_stats[['default,cohenK_limTPM,1.5,0,0.05']])
  
  # For a sample size of 3, over a third (38%) of  genes found to be perturbed in the heart represent false discoveries...
  fig2A_stats$`default,FDR,1.5,0,0.05` %>% select(N, contains('med'))
  
  # For kidney, and liver, a median sensitivity of 50% is attained by N=8...
  fig2A_stats$`default,detection,1.5,0,0.05` %>% select(N, contains('med'))
  fig2A_stats$`default,detection,1.5,0,0.05` %>% select(N, contains('mean'))
  
  # In all tissues, this variability drops markedly by N=6...    [sd (variability) of FDR in Fig2A]
  fig2A_stats$`default,FDR,1.5,0,0.05` %>% select(N, contains('sd')) %>% pivot_longer(cols= -N, values_to = 'value', names_to='tissue') %>% mutate(tissue=gsub("\\S+_", "", tissue)) %>%
    ggplot(aes(x = N, y = value, colour = tissue)) + geom_line(aes(group=tissue)) + geom_point()
  
  # The variation in sensitivity across trials also inversely falls more gradually...  [sd (variability) of detection in Fig2A]
  fig2A_stats$`default,detection,1.5,0,0.05` %>% select(N, contains('sd')) %>% pivot_longer(cols= -N, values_to = 'value', names_to='tissue') %>% mutate(tissue=gsub("\\S+_", "", tissue)) %>%
    ggplot(aes(x = N, y = value, colour = tissue)) + geom_line(aes(group=tissue)) + geom_point()
  
  # the differences between using an alpha of 0.05 vs 0.01 are minor...
  compare_runs_lineplot("default,FDR,1.5,0,0.05",
                        "default,FDR,1.5,0,0.01", "med")
  compare_runs_lineplot("default,detection,1.5,0,0.05",
                        "default,detection,1.5,0,0.01", "med")
  
  # By contrast, applying a minimum abundance threshold leaves the FDR increase largely unaffected...
  compare_runs_lineplot("default,FDR,1.5,0,0.05",
                        "default,FDR,1.5,1,0.05", "med")
  
  # ... while greatly increasing sensitivity to detect gold standard genes also passing this expression filter 
  compare_runs_lineplot("default,detection,1.5,0,0.05",
                        "default,detection,1.5,1,0.05", "med")
  
  # ashr
  compare_runs_lineplot("default,FDR,1.5,0,0.05",
                        "ashr,FDR,1.5,0,0.05", "med")
  compare_runs_lineplot("default,detection,1.5,0,0.05",
                        "ashr,detection,1.5,0,0.05", "med")
  # apeglm
  compare_runs_lineplot("default,FDR,1.5,0,0.05",
                        "apeglm,FDR,1.5,0,0.05", "med")
  compare_runs_lineplot("default,detection,1.5,0,0.05",
                        "apeglm,detection,1.5,0,0.05", "med")
  
  # ashr, at TPM=1 (asks: does ashr reduce FDR because of eliminating low-expressed genes?)
  compare_runs_lineplot("default,FDR,1.5,1,0.05",
                        "ashr,FDR,1.5,1,0.05", "med")
  compare_runs_lineplot("default,detection,1.5,1,0.05",
                        "ashr,detection,1.5,1,0.05", "med")
  
  compare_runs_lineplot("default,FDR,1.5,1,0.05",
                        "ashr,FDR,1.5,0,0.05", "med")
  compare_runs_lineplot("default,detection,1.5,1,0.05",
                        "ashr,detection,1.5,0,0.05", "med")
  
  ### Compare Figure 3 stuff ##
  compare_runs_lineplot("default,cohenK,1.5,0,0.05",
                        "default,cohenK_limTPM,1.5,0,0.05", "med")
  compare_runs_lineplot("default,cohenK,1.5,1,0.05",
                        "default,cohenK_limTPM,1.5,1,0.05", "med")
  compare_runs_lineplot("default,cohenK,1.5,1,0.01",
                        "default,cohenK,1.5,1,0.05", "med")
}


#### how many genes do we lose going from TPM0->TPM1 (FC=1)? ####
if (TRUE) {
  for (tissue in TISSUES_IN_ORDER) {
    (tot <- length(gold_standards$Dchs1[[paste(tissue, 0, "Padj", "std.05", "default", "0.05", "1", sep=",")]]))
    (lim <- length(gold_standards$Dchs1[[paste(tissue, 1, "Padj", "std.05", "default", "0.05", "1", sep=",")]]))
    print(paste(tissue, tot, lim, lim/tot))
  }
  # how many genes do we lose going from TPM0->TPM1 (FC=1.5)?
  for (tissue in TISSUES_IN_ORDER) {
    (tot <- length(gold_standards$Dchs1[[paste(tissue, 0, "Padj", "std.05", "default", "0.05", "1.5", sep=",")]]))
    (lim <- length(gold_standards$Dchs1[[paste(tissue, 1, "Padj", "std.05", "default", "0.05", "1.5", sep=",")]]))
    print(paste(tissue, tot, lim, lim/tot))
  }
}
# 50% detection (mean): H=9,K=8/9,Li=8,Lu=11/12
# 50% detection (median): H=8,K=9,Li=7,Lu=11


#################################
### FIGURE 2Cmatch - vary match.padj.cutoff P-value ####
#################################
plot_2Cmatch <- function(ovlap.results.long, fig_str = 'Figure2Cmatch',
                    FCtypes = c('default'), #'apeglm', 'ashr'),
                    FCs = 1.5, 
                    TPM.CUTOFFs = c(0,1)) {
  
  # arrange P-value from less to more specific
  ovlap.results.long$pval.cutoff <- factor(ovlap.results.long$pval.cutoff, 
                                           levels=c(1, 0.1, 0.05, 0.01))
  
  for (this.FCtype in FCtypes) {
    for (this.TISSUE in TISSUES_IN_ORDER) {
      for (this.FC.CUTOFF in FCs) { # c(1, 1.2, 1.5, 2)
        for (this.TPM.CUTOFF in TPM.CUTOFFs) {  # c(0, 1, 5)) { 
          p <- ovlap.results.long %>% 
            filter(type %in% c("detection", "FDR"),
                   tissue == this.TISSUE,
                   FC == this.FC.CUTOFF,
                   FCtype == this.FCtype,
                   alphaType == 'match.padj.cutoff', # affects things very slightly, other than for P=1
                   #alphaType == 'std.05', # affects things very slightly, other than for P=1
                   tpm.cutoff == this.TPM.CUTOFF) %>%
            mutate(type = gsub("detection", "Sensitivity", type)) %>%
            mutate(type = gsub("FDR", "False Discovery Rate", type)) %>%
            ggplot(aes(x = N, y = value)) + 
            facet_grid(pval.cutoff ~ type, labeller = labeller(pval.cutoff = PVALfunc, tpm.cutoff = TPMfunc)) +  
            ylab("") + 
            xlab(paste0("Sample Size")) + 
            theme(panel.grid.major = element_blank(), # leave blank
                  panel.grid.minor = element_blank(), # leave blank
                  axis.line = element_line(colour = "black"),
                  legend.position = 'none') + 
            scale_x_discrete(labels=PRETTY_LABELS)
          
          p1 <- pub_boxplot_theme(p, 'png')
          p2 <- pub_boxplot_theme(p, 'tiff')
          
          outfile <- paste0(fig_str, '_sensitivity.FDR_TPM', this.TPM.CUTOFF, '_FC', this.FC.CUTOFF, '_', PTYPE, '_', this.TISSUE, '_', this.FCtype)
          #openPng(outfile, p = p1, height = 900*3, width = 1200*3, res = thisres)
          tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
        }
      }
    } # tissues
  } # FCtype
}
plot_2Cmatch(ovlap.results.long, FCtypes = c('default'), FCs = 1.5, TPM.CUTOFFs = c(0) )



#### std.05 - match.padj.cutoff comparison scatterplots ####
# ovlap.results %>% head
# ovlap.results.long %>% head
# # ovlap.results %>% filter(pval.cutoff == 0.1, tpm.cutoff == 0, FCtype == 'default', tissue=='heart') %>% 
# for (this.pval in c(0.01, 0.05, 0.1, 1)) {
#   for (this.type in c("FDR", "detection")) {
#     for (this.tissue in TISSUES_IN_ORDER) {
#       for (this.FC in c(1, 1.2, 1.5, 2)) {
#         for (this.TPM.CUTOFF in c(0,1)) {
#         ovlap.results.long %>% 
#           filter(pval.cutoff == this.pval, tpm.cutoff == this.TPM.CUTOFF, FCtype == 'default', tissue==this.tissue, FC==this.FC, type==this.type) %>%
#           #select(alphaType, N, trial, FDR2) %>%
#           select(alphaType, N, trial, type, value) %>% 
#           pivot_wider(names_from = 'alphaType', values_from = 'value') %>% 
#           #filter(std.05 != match.padj.cutoff) %>% nrow()
#           ggplot(aes(x = std.05, y = match.padj.cutoff)) +
#           geom_point()  + geom_abline(slope = 1, intercept = 0) + #  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
#           facet_wrap(~N) + 
#           xlab("Using alpha=0.05") + 
#           ylab(paste0("Using alpha=", this.pval)) +
#             ggtitle(paste(this.tissue, this.type, paste0("FC=", this.FC), paste0("TPM=", this.TPM.CUTOFF))) + 
#           embiggen() -> p
#           #openPng(paste(this.tissue, this.type, this.pval, paste0("FC", this.FC), paste0("TPM", this.TPM.CUTOFF), sep="_"),
#                   p = p)
#         }
#       }
#     }
#   }
# }
