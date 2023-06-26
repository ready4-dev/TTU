#' 
#' Procure items from a dataset
#' @name procure-TTUProject
#' @description procure method applied to TTUProject
#' @param x An object of class TTUProject
#' @param type_1L_chr Type (a character vector of length one), Default: 'default'
#' @param variables_chr Variables (a character vector), Default: character(0)
#' @param what_1L_chr What (a character vector of length one), Default: 'records'
#' @param ... Additional arguments
#' @return Object (an output object of multiple potential types)
#' @rdname procure-methods
#' @aliases procure,TTUProject-method
#' @export 
#' @importFrom dplyr select
#' @importFrom ready4 procure
methods::setMethod("procure", "TTUProject", function (x, type_1L_chr = "default", variables_chr = character(0), 
    what_1L_chr = "records", ...) 
{
    if (what_1L_chr == "parameters") {
        if (type_1L_chr == "models_lup") {
            object_xx <- x@b_SpecificParameters@candidate_mdls_lup
        }
    }
    if (what_1L_chr %in% c("project")) {
        if (type_1L_chr == "models") {
            object_xx <- procureSlot(x, "c_SpecificProject", 
                use_procure_mthd_1L_lgl = T, what_1L_chr = "prefd_mdls")
        }
    }
    if (what_1L_chr %in% c("profile", "records")) {
        object_xx <- x@a_ScorzProfile@a_YouthvarsProfile
        if (!identical(variables_chr, character(0))) {
            object_xx@a_Ready4useDyad@ds_tb <- object_xx@a_Ready4useDyad@ds_tb %>% 
                dplyr::select(variables_chr)
        }
        if (what_1L_chr == "records") {
            if (type_1L_chr %in% c("default")) {
                object_xx <- object_xx@a_Ready4useDyad@ds_tb
            }
        }
    }
    return(object_xx)
})
