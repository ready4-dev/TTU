#' Extract g legend 1L
#' @description extract_g_legend_1L_chr() is an Extract function that extracts data from an object. Specifically, this function implements an algorithm to extract g legend 1L. The function is called for its side effects and does not return a value.
#' @param a.gplot PARAM_DESCRIPTION
#' @return NULL
#' @rdname extract_g_legend_1L_chr
#' @export 
#' @importFrom ggplot2 ggplot_gtable ggplot_build
#' @keywords internal
extract_g_legend_1L_chr <- function (a.gplot) 
{
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}
