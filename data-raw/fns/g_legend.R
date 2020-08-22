# function to extract the legend (borrowed from: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs )
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
