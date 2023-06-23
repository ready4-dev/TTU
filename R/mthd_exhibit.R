#' 
#' Exhibit features of a dataset by printing them to the R console
#' @name exhibit-TTUProject
#' @description exhibit method applied to TTUProject
#' @param x An object of class TTUProject
#' @param what_1L_chr What (a character vector of length one), Default: 'predictors'
#' @param ... Additional arguments
#' @return NULL
#' @rdname exhibit-methods
#' @aliases exhibit,TTUProject-method
#' @export 
#' @importFrom ready4 exhibit
methods::setMethod("exhibit", "TTUProject", function (x, what_1L_chr = "predictors", ...) 
{
    if (what_1L_chr == "predictors") {
        exhibitSlot(x, "b_SpecificParameters@predictors_lup", 
            ... = ...)
    }
})
