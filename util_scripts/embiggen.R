# convenience function: better default sizing of ggplots
embiggen <- function() {
  theme(legend.background = element_rect(linetype = "solid", color = "black"),
        legend.title = element_text(size = 25, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20, face = "bold"),
        axis.text.x=element_text(colour="black", size = 20, face = "bold"),
        axis.text.y=element_text(colour="black", size = 20, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 25),
        axis.title.y = element_text(face = "bold", size = 25),
        strip.text = element_text(colour="black", size = 22, face = "bold"))
}
