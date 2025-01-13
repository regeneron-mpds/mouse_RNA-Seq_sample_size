# Create a barplot of expression across samples, grouped by genotype
#
# @param  genes
# @param  this_tissue
# @param  this_expr
# @param  this_design
#
# @return a ggplot object 
plot_genes_by_genotype <- function(genes, 
                                   this_tissue, 
                                   this_expr = TPM,
                                   this_design = design) {
  exprType = deparse(substitute(this_expr))
  stopifnot(identical(colnames(this_expr), rownames(this_design)))
  
  tmp <- as.data.frame(t(this_expr[genes, , drop = FALSE])) %>% 
    tibble::rownames_to_column("sample_id") %>% 
    pivot_longer(cols = -1, names_to = "gene", values_to = 'expr') %>%
    inner_join(this_design) %>% 
    filter(tissue == this_tissue)
  stopifnot(all(table(tmp %>% filter(gene == genes[1]) %>% pull(mouse_id)) == 1))
  tmp$mouse_id <- factor(tmp$mouse_id, 
                         levels=tmp %>% filter(gene == genes[1]) %>% arrange(genotype, expr) %>% pull(mouse_id))

  tmp %>% 
    ggplot(aes(x = mouse_id, y = expr, fill = genotype)) + 
    geom_bar(aes(fill = genotype), stat = 'identity', position = 'dodge') + 
    scale_fill_discrete(name = 'genotype') +
    scale_color_discrete(name = 'genotype') +
    ylab(exprType) + 
    embiggen() + 
    theme(axis.text.x = element_text(colour = "black", size = 15, face = "bold", angle=45, hjust=1),
          axis.title.x = element_blank(),
          legend.position='bottom') + 
    facet_wrap(~gene, scales = "free_y", ncol = 1)
}