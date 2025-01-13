# PLOTTING DEFAULTS
ps <- 12
thisres <- 300

usual.par <- list(lwd = setNames(c(1, 0.7), nm = c("png", "tiff")),
                 outlier.size = setNames(c(2, 1), nm = c("png", "tiff")),
                 axis.text.x = setNames(c(14, 8), nm = c("png", "tiff")),
                 axis.text.y = setNames(c(18, 10), nm = c("png", "tiff")),
                 axis.title.x = setNames(c(25, 13.5), nm = c("png", "tiff")),
                 axis.title.y = setNames(c(25, 14), nm = c("png", "tiff")),
                 strip.text = setNames(c(22, 15), nm = c("png", "tiff")))
PRETTY_LABELS <- c("9"="", "11" = "", "13" = "", "15" = "", "17" = "", "19" = "", "21" = "", "23" = "", "25" = "", "27" = "")
# PRETTY_LABELS <- c("8"="", "10" = "", "12" = "", "14" = "", "16" = "", "18" = "", "20" = "", "22" = "", "24" = "", "26" = "", "28" = "")
FCfunc <- function(x) { paste("FC >", x) }
TPMfunc <- function(x) { paste("TPM >", x) }
PVALfunc <- function(x) { paste("P \u2264", x) }   # greater than equal to: \u2265

pub_boxplot_theme <- function(ggplot_obj, img_type, fixedcolor, this.par = usual.par) {
  if (missing(fixedcolor)) {
    ggplot_obj <- ggplot_obj + 
      geom_boxplot(lwd = this.par$lwd[img_type], outlier.size = this.par$outlier.size[img_type], aes(color = type), fatten = 0.85)
  } else {
    ggplot_obj <- ggplot_obj + 
      geom_boxplot(lwd = this.par$lwd[img_type], outlier.size = this.par$outlier.size[img_type], color = fixedcolor, fatten = 0.85)
  }
  
  ggplot_obj + 
    theme(axis.text.x = element_text(colour = "black", size = this.par$axis.text.x[img_type], face = "bold"), 
          axis.text.y = element_text(colour = "black", size = this.par$axis.text.y[img_type], face = "bold"), 
          axis.title.x = element_text(face = "bold", size = this.par$axis.title.x[img_type]), 
          axis.title.y = element_text(face = "bold", size = this.par$axis.title.y[img_type]), 
          strip.text = element_text(colour = "black", size = this.par$strip.text[img_type], face = "bold"))
}



this4b.par <- list(lwd = setNames(c(1, 0.7), nm = c("png", "tiff")),
                   outlier.size = setNames(c(2, 1), nm = c("png", "tiff")),
                   axis.text.x = setNames(c(14, 8), nm = c("png", "tiff")),
                   axis.text.y = setNames(c(14, 10), nm = c("png", "tiff")),
                   axis.title.x = setNames(c(19, 13.5), nm = c("png", "tiff")),
                   axis.title.y = setNames(c(19, 14), nm = c("png", "tiff")),
                   legend.title = setNames(c(13, 6.5), nm = c("png", "tiff")),
                   legend.text = setNames(c(13, 6.5), nm = c("png", "tiff")),                 
                   legend.key.height = setNames(c(0.5, 0.25), nm = c("png", "tiff")),
                   legend.key.width = setNames(c(1.5, 0.75), nm = c("png", "tiff")),
                   strip.text = setNames(c(18, 15), nm = c("png", "tiff")))

pub_boxplot_theme_4b <- function(ggplot_obj, img_type, this.par = this4b.par) {
  ggplot_obj + geom_boxplot(lwd = this.par$lwd[img_type], outlier.size = this.par$outlier.size[img_type], aes(color = `significant trial`), fatten = 0.8) + 
  theme(axis.text.x = element_text(colour = "black", size = this.par$axis.text.x[img_type], face = "bold"), 
        axis.text.y = element_text(colour = "black", size = this.par$axis.text.y[img_type], face = "bold"), 
        axis.title.x = element_text(face = "bold", size = this.par$axis.title.x[img_type]), 
        axis.title.y = element_text(face = "bold", size = this.par$axis.title.y[img_type]), 
        legend.title = element_text(colour = "black", size = this.par$legend.title[img_type], face = "bold"),
        legend.text = element_text(colour = "black", size = this.par$legend.text[img_type], face = "bold"),
        legend.key.height= unit(this.par$legend.key.height[img_type], 'cm'),
        legend.key.width= unit(this.par$legend.key.width[img_type], 'cm'))  
}

S4.par <- list(lwd = setNames(c(1, 0.7), nm = c("png", "tiff")),
               outlier.size = setNames(c(2, 1), nm = c("png", "tiff")),
               axis.text.x = setNames(c(14, 8), nm = c("png", "tiff")),
               axis.text.y = setNames(c(14, 10), nm = c("png", "tiff")),
               axis.title.x = setNames(c(19, 13.5), nm = c("png", "tiff")),
               axis.title.y = setNames(c(19, 14), nm = c("png", "tiff")),
               strip.text = setNames(c(18, 15), nm = c("png", "tiff")))

