# function to extract the legend (borrowed from: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs )
extract_g_legend_1L_chr <- function(a.gplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
  }
