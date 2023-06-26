#' 
#' Manufacture a new object
#' @name manufacture-TTUProject
#' @description manufacture method applied to TTUProject
#' @param x An object of class TTUProject
#' @param type_1L_chr Type (a character vector of length one), Default: 'dummys'
#' @param what_1L_chr What (a character vector of length one), Default: 'factors'
#' @param ... Additional arguments
#' @return Object (an output object of multiple potential types)
#' @rdname manufacture-methods
#' @aliases manufacture,TTUProject-method
#' @export 
#' @importFrom ready4 manufacture
methods::setMethod("manufacture", "TTUProject", function (x, type_1L_chr = "dummys", what_1L_chr = "factors", 
    ...) 
{
    object_xx <- manufacture(x@c_SpecificProject@a_YouthvarsProfile@a_Ready4useDyad, 
        type_1L_chr = type_1L_chr, what_1L_chr = what_1L_chr, 
        restrict_to_chr = x@c_SpecificProject@b_SpecificParameters@candidate_covars_chr)
    return(object_xx)
})
