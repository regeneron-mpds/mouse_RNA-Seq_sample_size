# Wrapper around png with sensible defaults
#
# @param  uniqstr basename for output png file
# @param  out_dir
# @param  p a (gg)plot object. If not NULL, it is displayed with plot(p), and connection is closed
# @param  height    height of png in pixels
# @param  width     width of png in pixels
# @param  pointsize default pointsize for plotted text
#
# @return opens a png graphics device if p is NULL, or creates a png file otherwise
openPng <- function(uniqstr, out_dir = "figures/", p = NULL,
                    height = 960, width = 1200, pointsize = 36, ...) {
  if(!dir.exists(out_dir)) {
    dir.create(out_dir, showWarnings=F)
  }

  filename <- paste0(out_dir, uniqstr, ".png")
  png(filename, height=height, width=width, pointsize=pointsize, bg='transparent', ...)

  if (!is.null(p)) {
    print(p)
    dev.off()
  }
}
