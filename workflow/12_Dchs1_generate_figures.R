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


# read in statistics from down-sampling
ovlap.results <- readRDS('result_objs/Dchs1_ovlap.results.RDS')
ovlap.results$N <- factor(ovlap.results$N, 
                          levels = as.character(min(ovlap.results$N):max(ovlap.results$N)))

# rejigger columns 
ovlap.results.full <- ovlap.results %>%
  mutate(false.pos1 = signif1 - ovlap1,
        false.pos1.2 = signif1.2 - ovlap1.2,
        false.pos1.5 = signif1.5 - ovlap1.5,
        false.pos2 = signif2 - ovlap2,
        precision_FC1 = 100 * ovlap1 / signif1,
        precision_FC1.2 = 100 * ovlap1.2 / signif1.2,
        precision_FC1.5 = 100 * ovlap1.5 / signif1.5,
        precision_FC2 = 100 * ovlap2 / signif2,
        FDR_FC1 = 100 * false.pos1 / signif1,
        FDR_FC1.2 = 100 * false.pos1.2 / signif1.2,
        FDR_FC1.5 = 100 * false.pos1.5 / signif1.5,
        FDR_FC2 = 100 * false.pos2 / signif2) %>%
  replace(is.na(.), 0)
head(ovlap.results.full)

# re-rejigger columns 
ovlap.results.long <- ovlap.results.full %>%
  tidyr::pivot_longer(cols = 6:ncol(ovlap.results.full), names_to = "type") %>%
  mutate(FC = gsub("detection", "", type),
        FC = gsub("false.pos", "", FC),
        FC = gsub("FDR_FC", "", FC),
        FC = gsub("ovlap", "", FC),
        FC = gsub("precision_FC", "", FC),
        FC = gsub("signif", "", FC),
        FC = gsub("\\S*cor\\S*", NA, FC),
        type = gsub("1.5", "", type),
        type = gsub("1.2", "", type),
        type = gsub("1", "", type),
        type = gsub("2", "", type),
        type = gsub("_FC", "", type))
table(ovlap.results.long$FC, ovlap.results.long$type, useNA = 'always')


#################
### FIGURE 2A ###
#################
for (this.FC in 1.5) {  # c(1, 1.2, 1.5, 2)) {
  for (this.TPM.CUTOFF in c(0)) {  # c(0, 1, 5)) { 
    for (this.PVAL.CUTOFF in c(0.05)) {  # c(0.01, 0.05, 0.1, 1)) {
      p <- ovlap.results.long %>% 
        filter(pval.cutoff == this.PVAL.CUTOFF,
               type %in% c("detection", "FDR"),
               FC == this.FC,
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
      
      outfile <- paste0('Figure2A_sensitivity.FDR_FC', this.FC, '_TPM', this.TPM.CUTOFF, '_', PTYPE, this.PVAL.CUTOFF)
      openPng(outfile, p = p1, height = 900*3, width = 1200*3, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    }
  }
}

#################
### FIGURE 2B ###
# & Supp Fig 1A #
#################
for (this.TISSUE in TISSUES_IN_ORDER) {
  for (this.TPM.CUTOFF in c(0)) { # c(0, 1, 5)) {
    for (this.PVAL.CUTOFF in c(0.05)) { #c(0.01, 0.05, 0.1, 1)) {
      p <- ovlap.results.long %>% 
        filter(pval.cutoff == this.PVAL.CUTOFF,
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
      
      outfile <- paste0('Figure2B_sensitivity.FDR_TPM', this.TPM.CUTOFF, '_', PTYPE, '_', this.PVAL.CUTOFF, '_', this.TISSUE)
      openPng(outfile, p = p1, height = 900*3, width = 1200*3, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    }
  }
}


#################
### FIGURE 2C ###
# & Supp Fig 1B #
#################
for (this.TISSUE in TISSUES_IN_ORDER) {
  for (this.FC.CUTOFF in 1.5) { # c(1, 1.2, 1.5, 2)
    for (this.TPM.CUTOFF in c(0)) { # c(0, 1, 5)) { 
      p <- ovlap.results.long %>% 
        filter(FC == this.FC.CUTOFF,
               type %in% c("detection", "FDR"),
               tissue == this.TISSUE,
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
      
      outfile <- paste0('Figure2C_sensitivity.FDR_TPM', this.TPM.CUTOFF, '_FC', this.FC.CUTOFF, '_', PTYPE, '_', this.TISSUE)
      openPng(outfile, p = p1, height = 900*3, width = 1200*3, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    }
  }
}


###############################
### FIGURE 3 pre-processing ###
###############################
# read in statistics from pairwise comparisons
pairwise.results <- readRDS('result_objs/Dchs1_pairwise.results.RDS')

# re-format
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
#
pw.res.final <- pw.res.final %>%
  tidyr::pivot_longer(cols = 6:ncol(pw.res.final), names_to = "type") %>%
  mutate(type2 = gsub("\\d+\\S*", "", type),
         FC = gsub("^\\.+", "", gsub('[^\\d\\.]', "", type, perl = TRUE)))


#################
### FIGURE 3 ####
#################
for (this.pval in c(0.05)) { # 0.01, 0.05, 0.1, 1)) {
  for (this.TPM.CUTOFF in c(0)) { # c(0, 1, 5) {
    for (this.type2 in c('cohenK', 'cohenK_limTPM')) {
      p <- pw.res.final %>% 
        filter(tpm.cutoff == this.TPM.CUTOFF,
               pval.cutoff == this.pval,
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
      
      outfile <- paste0('Figure3_rangeFCs_allTissues_TPM', this.TPM.CUTOFF, '_', PTYPE, this.pval, '_', this.type2)
      openPng(outfile, p = p1, height = 700*3, width = 1200*3, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 350*4, width = 500*4+200, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    }
  }
}


############################
### FIGURE 4 pre-process ###
############################
effect.size.results <- readRDS("result_objs/Dchs1_effect.size.results.RDS")

merged.ovlap.effect <- merge(ovlap.results %>% filter(pval.cutoff == 0.05),
                             effect.size.results) %>%
  select(tissue, tpm.cutoff, N, trial, contains("signif"), contains("gt.orig")) %>%
  pivot_longer(cols = contains("gt.orig") | contains("signif"), names_to = "type") %>% 
  mutate(FC = gsub("signif", "", gsub("gt.orig", "", type)),
         type = gsub("\\d\\S*", "", type)) %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  mutate(perc = 100 * gt.orig / signif,
         text = 'fraction\noverestimated')
nrow(merged.ovlap.effect)
head(merged.ovlap.effect)


#################
### FIGURE 4A ###
#################
for (this.TPM.CUTOFF in c(0)) { # c(0, 1, 5)) {
  p <- merged.ovlap.effect %>%
    filter(tpm.cutoff == this.TPM.CUTOFF) %>%
    ggplot(aes(x = N, y = perc)) + 
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
  
  outfile <- paste0("Figure4A_effect_size_proportion_alltissues_PVAL0.05_TPM", this.TPM.CUTOFF)
  openPng(outfile, p = p1, height = 900*3, width = 1200*3*1.6, res = thisres)
  tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4*2, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()    
}

#################
### FIGURE 4B ###
#################
genewise.effect.size.res <- readRDS("result_objs/Dchs1_genewise.effect.size.res.RDS")
this.gene <- 'Trex1'
  
this.tissue = 'liver'
this.TPM.CUTOFF = 0

p <- genewise.effect.size.res %>% 
  filter(gene == this.gene, tissue == this.tissue, tpm.cutoff == this.TPM.CUTOFF) %>%
  mutate(N = factor(N)) %>% 
  mutate(is.signif = ifelse(is.signif, "YES", "NO")) %>% 
  rename(`significant trial` = is.signif) %>%
  ggplot(aes(x = N, y = log2(FC))) + 
  ylab(paste(this.gene, "log2(FC)")) + 
  xlab("Sample Size") +  
  theme(panel.grid.major = element_blank(), # leave blank
        panel.grid.minor = element_blank(), # leave blank
        axis.line = element_line(colour = "black"), 
        legend.position = c(0.8, 0.2)) + 
  geom_hline(yintercept = log2(N30_DEres[['Fold.Change']][this.gene, this.tissue]),
             linetype = ifelse(this.gene %in% gold_standards$Dchs1[[paste(this.tissue, this.TPM.CUTOFF, 'Padj,1', sep=",")]], "solid", "dashed")) + 
  scale_x_discrete(labels=PRETTY_LABELS)

p1 <- pub_boxplot_theme_4b(p, 'png')
p2 <- pub_boxplot_theme_4b(p, 'tiff')

outfile <- paste0("Figure4B_", this.gene, "_", this.tissue, "_TPM", this.TPM.CUTOFF)
openPng(outfile, p = p1, height = 900, width = 1200*1.75, res = thisres)
tiff(paste0("figures/", outfile, '.tiff'), height = 450*1.5, width = 500*2.5, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()    


####################################
### SUPPLEMENTARY FIGURE 1A ########
### see Figure 2B section above ####
####################################

####################################
### SUPPLEMENTARY FIGURE 1B ########
### see Figure 2C section above ####
####################################

###############################
### SUPPLEMENTARY FIGURE 1C ##
###############################
for (this_TISSUE in TISSUES_IN_ORDER) {
  for (this.FC.CUTOFF in c(1.5)) { # c(1, 1.2, 1.5, 2)) {
    for (this.PVAL.CUTOFF in c(0.05)) { #c(0.01, 0.05, 0.1, 1)) {
      p <- ovlap.results.long %>% 
        filter(pval.cutoff == this.PVAL.CUTOFF,
               FC == this.FC.CUTOFF,
               tissue == this_TISSUE,
               type %in% c("detection", "FDR")) %>%
        mutate(type = gsub("detection", "Sensitivity", type),
               type = gsub("FDR", "False Discovery Rate", type)) %>%
        ggplot(aes(x = N, y = value)) + 
        facet_grid(tpm.cutoff ~ type, labeller = labeller(FC = FCfunc, tpm.cutoff = TPMfunc)) +  
        ylab("") + 
        xlab(paste0("Sample Size")) +
        theme(panel.grid.major = element_blank(), # leave blank
              panel.grid.minor = element_blank(), # leave blank
              axis.line = element_line(colour = "black"),
              legend.position = 'none') + 
        scale_x_discrete(labels=PRETTY_LABELS)
      
      p1 <- pub_boxplot_theme(p, 'png')
      p2 <- pub_boxplot_theme(p, 'tiff')
      
      outfile <- paste0('SuppFig1C_sensitivity.FDR_varyTPM_FC', this.FC.CUTOFF, '_', PTYPE, this.PVAL.CUTOFF, '_', this_TISSUE)
      openPng(outfile, p = p1, height = 900*3, width = 1200*3, res = thisres)
      tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
    }
  }
}



##############################
### SUPPLEMENTARY FIGURE 2A ##
##############################
MIN_NUM_SIGNIF <- 10
goldstandard.FC <- 1.5
prob_correct <- readRDS(file = paste0('result_objs/Dchs1_prob_correct_GSFC', goldstandard.FC, '_minSignif', MIN_NUM_SIGNIF, '.RDS'))

for (THIS_TISSUE in TISSUES_IN_ORDER) {
  prob_correct %>%
    filter(tissue == THIS_TISSUE, num.signif > MIN_NUM_SIGNIF) %>%
    pivot_longer(cols = starts_with("fc"), names_to = 'fc_range', values_to = 'fraction_detected') -> tmp
  tmp$fc_range[tmp$fc_range == 'fc1.2'] <- '1-1.2'
  tmp$fc_range[tmp$fc_range == 'fc1.5'] <- '1.2-1.5'
  tmp$fc_range[tmp$fc_range == 'fc2'] <- '1.5-2'
  tmp$fc_range[tmp$fc_range == 'fc3'] <- '2-3'
  tmp$fc_range[tmp$fc_range == 'fc5'] <- '3-5'
  tmp$fc_range[tmp$fc_range == 'fc10'] <- '5-10'
  tmp$fc_range[tmp$fc_range == 'fcInf'] <- '>10'
  tmp$fc_range <- factor(tmp$fc_range, levels = c('1-1.2', '1.2-1.5', '1.5-2', '2-3', '3-5', '5-10', '>10'))
  
  # plot
  outfile <- paste0("SuppFig2A_", THIS_TISSUE, '_GSFC', goldstandard.FC, "_minSignif", MIN_NUM_SIGNIF)
  p <- tmp %>% 
    mutate(fraction_detected = 100 * fraction_detected,
           N = factor(N, levels = unique(sort(N))),
           trial = factor(trial, levels = unique(rev(sort(trial))))) %>%
    ggplot(aes(x = fc_range, y = fraction_detected)) + 
    xlab("Fold Change in down-sampled trial") + 
    ylab(paste("% of significant genes also detected in\ngold standard signature with fold change >", goldstandard.FC)) + 
    facet_wrap(~N) + 
    theme_bw() + embiggen() + 
    theme(legend.position = 'none')
  
  # PNG VERSION
  p1 <- p + 
    geom_boxplot() +
    theme(axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle=60, hjust = 1))
  openPng(outfile, p = p1, width = 1400, height = 1200)
  
  # TIFF VERSION
  p2 <- p + 
    geom_boxplot(outlier.size = 0.7, lwd=0.35) +  
    theme(axis.text.x = element_text(colour = "black", size = 8, face = "bold", angle=60, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
          axis.title.x = element_text(face = "bold", size = 13.5), 
          axis.title.y = element_text(face = "bold", size = 14),
          strip.text = element_text(colour = "black", size = 10, face = "bold"))
  tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = 12, res= 300, compression = 'lzw'); plot(p2); dev.off()
}


##############################
### SUPPLEMENTARY FIGURE 2B ##
##############################
goldstandard.FC <- 1.5
MIN_NUM_SIGNIF <- 10
prob_detect_cumul <- readRDS(paste0('result_objs/Dchs1_prob_detect_cumul_GSFC', goldstandard.FC, '_minSignif', MIN_NUM_SIGNIF, '.RDS'))

for (THIS_TISSUE in TISSUES_IN_ORDER) {
  prob_detect_cumul %>%
    filter(tissue == THIS_TISSUE, num.signif > MIN_NUM_SIGNIF) %>%
    pivot_longer(cols = starts_with("fc"), names_to = 'fc_range', values_to = 'fraction_detected') -> tmp
  tmp$fc_range[tmp$fc_range == 'fc1.2'] <- '>1'
  tmp$fc_range[tmp$fc_range == 'fc1.5'] <- '>1.2'
  tmp$fc_range[tmp$fc_range == 'fc2'] <- '>1.5'
  tmp$fc_range[tmp$fc_range == 'fc3'] <- '>2'
  tmp$fc_range[tmp$fc_range == 'fc5'] <- '>3'
  tmp$fc_range[tmp$fc_range == 'fc10'] <- '>5'
  tmp$fc_range[tmp$fc_range == 'fcInf'] <- '>10'
  tmp$fc_range <- factor(tmp$fc_range, levels = c('>1', '>1.2', '>1.5', '>2', '>3', '>5', '>10')) # '>5'))
  
  # plot
  outfile <- paste0("SuppFig2B_", THIS_TISSUE, '_GSFC', goldstandard.FC, "_minSignif", MIN_NUM_SIGNIF)
  p <- tmp %>%
    mutate(fraction_detected = 100 * fraction_detected,
           N = factor(N, levels = unique(sort(N))),
           trial = factor(trial, levels = unique(rev(sort(trial))))) %>%
    ggplot(aes(x = fc_range, y = fraction_detected)) +
    xlab("Fold Change in down-sampled trial") + 
    ylab(paste("% of gold standard genes with fold change >", goldstandard.FC, "\nalso detected in sub-sampled signature")) + 
    facet_wrap(~N) + 
    theme_bw() + embiggen() + 
    theme(legend.position = 'none')
  
  # PNG VERSION
  p1 <- p + 
    geom_boxplot() +
    theme(axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle=60, hjust = 1))
  openPng(outfile, p = p1, width = 1400, height = 1200)    
  
  # TIFF VERSION
  p2 <- p + 
    geom_boxplot(outlier.size = 0.7, lwd=0.35) +  #geom_boxplot(stat='identity') + 
    theme(axis.text.x = element_text(colour = "black", size = 8, face = "bold", angle=60, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
          axis.title.x = element_text(face = "bold", size = 13.5), 
          axis.title.y = element_text(face = "bold", size = 14),
          strip.text = element_text(colour = "black", size = 10, face = "bold")) +
  tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4, pointsize = 12, res= 300, compression = 'lzw'); plot(p2); dev.off()
}


##############################
### SUPPLEMENTARY FIGURE 3 ###
##############################
for (this.tissue in TISSUES_IN_ORDER) {
  for (this.pval in c(0.05)) { # c(0.01, 0.05, 0.1, 1)) {
    p <- pw.res.final %>% 
      filter(tissue == this.tissue, 
             pval.cutoff == this.pval,
             type2 == 'cohenK') %>%
      ggplot(aes(x = N, y = value)) + 
      facet_grid(tpm.cutoff ~ FC, labeller = labeller(FC = FCfunc, tpm.cutoff = TPMfunc)) + 
      ylab(paste0("Agreement between signatures\n(", gsub("\\.", " ", this.type2), ")")) +  
      xlab("Sample Size") +
      theme(panel.grid.major = element_blank(), # leave blank
            panel.grid.minor = element_blank(), # leave blank
            axis.line = element_line(colour = "black"), 
            legend.position = 'none') +
      scale_x_discrete(labels=PRETTY_LABELS)
    
    p1 <- pub_boxplot_theme(p, 'png', fixedcolor='blue')
    p2 <- pub_boxplot_theme(p, 'tiff', fixedcolor='blue')
    outfile <- paste0("SuppFig3_rangeFC.TPM_", this.tissue, "_PVAL", this.pval, "_cohenK")
    
    openPng(outfile, p = p1, height = 700*3*6/7, width = 1200*3, res = thisres)
    tiff(paste0("figures/", outfile, '.tiff'), height = 350*4*6/7, width = 500*4+200, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
  }
}



##############################
### SUPPLEMENTARY FIGURE 4 ###
##############################
if (USE_ORIG_GENES_FOR_S4) {
  genes.for.effect.size.exploration <- readRDS("result_objs/Dchs1_genes.for.effect.size.explorationORIGgenes.RDS")  
} else {
  genes.for.effect.size.exploration <- readRDS("result_objs/Dchs1_genes.for.effect.size.exploration.RDS")
}
this.TPM.CUTOFF = 0

for (this.tissue in TISSUES_IN_ORDER) {
  # get fold changes from N30 experiment
  this_genes <- genes.for.effect.size.exploration[[this.tissue]]
  tmp.N30.fcs <- data.frame(
    gene = this_genes,
    fc = log2(N30_DEres$Fold.Change[this_genes, this.tissue])) %>%
    mutate(abs.fc = abs(fc),
           tissue = this.tissue) %>% 
    arrange(desc(abs.fc))
  genes_in_order <- tmp.N30.fcs$gene
  
  # prepare effect size table (including above FCs), filter to tissue, etc.
  tmp0 <- left_join(tmp.N30.fcs, genewise.effect.size.res) %>%
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
    rename(`P < 0.05` = is.signif) %>%
    ggplot(aes(x = `sample size`, y = log2(FC))) + 
    xlab('Sample Size') + 
    ylab(bquote(bold(log[2] ~FC))) + 
    theme(axis.text.x = element_blank(),
          panel.grid.major = element_blank(), # leave blank
          panel.grid.minor = element_blank(), # leave blank
          legend.title = element_blank()) + 
    facet_wrap(~titlelabel, scales = "free_y", ncol = 5) + 
    scale_x_discrete(labels=c("4"="", "5"="", "7"="", "8"="", "9"="", "12"="", "14"="", "11" = "", "13" = "", "16" = "", "17" = "", "18"="", "19" = "", "21" = "", "22" = "", "23" = "", "24"="", "26"="",  "27" = "", "28" = ""))
  
  
  p1 <- p + geom_boxplot(lwd = S4.par$lwd['png'], outlier.size = S4.par$outlier.size['png'], aes(color = `P < 0.05`)) +
    geom_hline(data = tmp.N30.fcs %>% left_join(tmp2), aes(yintercept = fc), colour = "black", size = S4.par$lwd['png'], linetype = 'dashed') + 
    theme(axis.text.x = element_text(colour = "black", size = S4.par$axis.text.x['png'], face = "bold"), 
          axis.text.y = element_text(colour = "black", size = S4.par$axis.text.y['png'], face = "bold"), 
          axis.title.x = element_text(face = "bold", size = S4.par$axis.title.x['png']), 
          axis.title.y = element_text(face = "bold", size = S4.par$axis.title.y['png']), 
          strip.text = element_text(colour = "black", size = S4.par$strip.text['png'], face = "bold"))
  
  p2 <- p + geom_boxplot(lwd = S4.par$lwd['tiff'], outlier.size = S4.par$outlier.size['tiff'], aes(color = `P < 0.05`)) +
    geom_hline(data = tmp.N30.fcs %>% left_join(tmp2), aes(yintercept = fc), colour = "black", size = S4.par$lwd['tiff'], linetype = 'dashed') + 
    theme(axis.text.x = element_text(colour = "black", size = S4.par$axis.text.x['tiff'], face = "bold"), 
          axis.text.y = element_text(colour = "black", size = S4.par$axis.text.y['tiff'], face = "bold"), 
          axis.title.x = element_text(face = "bold", size = S4.par$axis.title.x['tiff']), 
          axis.title.y = element_text(face = "bold", size = S4.par$axis.title.y['tiff']), 
          strip.text = element_text(colour = "black", size = S4.par$strip.text['tiff'], face = "bold"))
  
  outfile <- ifelse(USE_ORIG_GENES_FOR_S4, 
                    paste0("SuppFig4_manyORIGgenes_", this.tissue, "_tpm", this.TPM.CUTOFF),
                    paste0("SuppFig4_manygenes_", this.tissue, "_tpm", this.TPM.CUTOFF))
  message(outfile)
  
  openPng(outfile, p = p1, height = 900*3, width = 1200*3*1.6, res = thisres)
  tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*4*2, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()    
}


####################################
### SUPPLEMENTARY FIGURES 5-8 ######
### (Fat4 version of main figures) #
####################################


##############################
### SUPPLEMENTARY FIGURE 9 ###
##############################
for (this.TPM.CUTOFF in c(0)) { # c(0,1,5)) { 
  for (this.PVAL.CUTOFF in c(0.05)) {  #  c(0.01, 0.05, 0.1, 1)) {
    p <- ovlap.results.long %>% 
      filter(pval.cutoff == this.PVAL.CUTOFF,
             tpm.cutoff == this.TPM.CUTOFF,
             type == 'all.cor') %>% 
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
    
    outfile <- paste0("SuppFig9_crlxns_TPM", this.TPM.CUTOFF, "_", PTYPE, this.PVAL.CUTOFF)

    p1 <- pub_boxplot_theme(p, 'png') + theme(strip.text.x = element_blank())
    p2 <- pub_boxplot_theme(p, 'tiff') + theme(strip.text.x = element_blank())
    
    openPng(outfile, p = p1, height = 900*3, width = 1200*1.5, res = thisres)
    tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*2.4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
  }
}


##############################
### SUPPLEMENTARY FIGURE 10 ##
##############################
for (this.TPM.CUTOFF in c(0)) { # c(0, 1, 5)) {
  for (this.PVAL.CUTOFF in c(0.05)) { # c(0.01, 0.05, 0.1, 1)) {
    p <- pw.res.final %>% 
      filter(tpm.cutoff == this.TPM.CUTOFF,
             pval.cutoff == this.PVAL.CUTOFF,
             type2 == 'all.cor') %>%
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
    
    outfile <- paste0("SuppFig10_pairwise_crlxns_TPM", this.TPM.CUTOFF, "_", PTYPE, this.PVAL.CUTOFF)
    
    p1 <- pub_boxplot_theme(p, 'png', fixedcolor = 'blue') + 
      theme(strip.text.x = element_blank(),
            plot.title = element_text(face = 'bold', colour = 'black', size = usual.par$plot.title['png']))
    
    p2 <- pub_boxplot_theme(p, 'tiff', fixedcolor = 'blue') + 
      theme(strip.text.x = element_blank(),
            plot.title = element_text(face = 'bold', colour = 'black', size = usual.par$plot.title['tiff']))
    
    openPng(outfile, p = p1, height = 900*3, width = 1200*1.5, res = thisres)
    tiff(paste0("figures/", outfile, '.tiff'), height = 450*4, width = 500*2.4, pointsize = ps, res= thisres, compression = 'lzw'); plot(p2); dev.off()
  }
}



##############################
### SUPPLEMENTARY FIGURE 11 ##
##############################
this.PVAL.CUTOFF <- 0.05
pairwise_WT_comp <- readRDS(paste0('result_objs/Dchs1_WTcompare.results_', PTYPE, this.PVAL.CUTOFF, '.RDS'))

pairwise_WT_comp_long <- pairwise_WT_comp %>%
  tidyr::pivot_longer(cols = 5:ncol(pairwise_WT_comp), names_to = "FC_cutoff") %>%
  mutate(FC_cutoff = gsub("uniq", "", FC_cutoff),
         N = factor(N))
head(pairwise_WT_comp_long)

for (this.FC in c(1, 1.2, 1.5, 2)) {
  p <- pairwise_WT_comp_long %>% 
    filter(FC_cutoff == this.FC) %>%
    mutate(tpm.cutoff = paste0('TPM=', tpm.cutoff)) %>%
    ggplot(aes(x = N, y = value)) + 
    geom_boxplot(lwd = 1.6, outlier.size = 2) + 
    facet_wrap(tpm.cutoff ~ tissue, scales = "free", ncol = 4, labeller = label_wrap_gen(multi_line=FALSE, width=30)) +
    ylab("Number of genes in WT-vs-WT signature") + 
    xlab("Sample Size") + 
    theme(axis.text = element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 25, face = 'bold'),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.text = element_text(colour="black", size = 22, face = "bold"),
          panel.background = element_rect(fill = "white", colour = NA))
  openPng(paste0("wrap_12panel_numFPs_FC", this.FC, "_", PTYPE), p = p, height = 900, width = 1850) #, "_TPM", this.tpm, "_", this.ptype
  openPng(paste0("wrap_12panel_numFPs_FC", this.FC, "_", PTYPE, "_limy30"), p = p + ylim(c(0,30)), height = 900, width = 1850) #, "_TPM", this.tpm, "_", this.ptype    
}
p
