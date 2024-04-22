#' 
#' Exhibit features of model module data by printing them to the R console
#' @name exhibit-TTUProject
#' @description exhibit method applied to TTUProject
#' @param x An object of class TTUProject
#' @param captions_chr Captions (a character vector), Default: NULL
#' @param display_1L_chr Display (a character vector of length one), Default: 'all'
#' @param header_1L_chr Header (a character vector of length one), Default: ''
#' @param header_col_nms_chr Header column names (a character vector), Default: ' '
#' @param mkdn_tbl_refs_chr Markdown table references (a character vector), Default: NULL
#' @param profile_idx_int Profile index (an integer vector), Default: NA
#' @param output_type_1L_chr Output type (a character vector of length one), Default: 'HTML'
#' @param type_1L_chr Type (a character vector of length one), Default: 'default'
#' @param use_lbls_as_col_nms_1L_lgl Use labels as column names (a logical vector of length one), Default: T
#' @param use_rdocx_1L_lgl Use rdocx (a logical vector of length one), Default: F
#' @param variables_chr Variables (a character vector), Default: character(0)
#' @param what_1L_chr What (a character vector of length one), Default: 'predictors'
#' @param ... Additional arguments
#' @return No return value, called for side effects.
#' @rdname exhibit-methods
#' @aliases exhibit,TTUProject-method
#' @export 
#' @importFrom ready4 exhibit
methods::setMethod("exhibit", "TTUProject", function (x, captions_chr = NULL, display_1L_chr = "all", header_1L_chr = "", 
    header_col_nms_chr = " ", mkdn_tbl_refs_chr = NULL, profile_idx_int = NA_integer_, 
    output_type_1L_chr = "HTML", type_1L_chr = "default", use_lbls_as_col_nms_1L_lgl = T, 
    use_rdocx_1L_lgl = F, variables_chr = character(0), what_1L_chr = "predictors", 
    ...) 
{
    if (what_1L_chr == "predictors") {
        exhibitSlot(x, "b_SpecificParameters@predictors_lup", 
            ... = ...)
    }
    if (what_1L_chr == "profile") {
        exhibit(procure(x, variables_chr = variables_chr, type_1L_chr == 
            "default", what_1L_chr = what_1L_chr, ... = ...), 
            captions_chr = captions_chr, header_1L_chr = header_1L_chr, 
            header_col_nms_chr = header_col_nms_chr, mkdn_tbl_refs_chr = mkdn_tbl_refs_chr, 
            profile_idx_int = profile_idx_int, output_type_1L_chr = output_type_1L_chr, 
            what_1L_chr = type_1L_chr)
    }
    if (what_1L_chr == "records") {
        if (type_1L_chr %in% c("ds", "dict")) {
            exhibit(x@a_ScorzProfile@a_YouthvarsProfile@a_Ready4useDyad, 
                caption_1L_chr = {
                  if (is.null(captions_chr)) {
                    NA_character_
                  }
                  else {
                    captions_chr[1]
                  }
                }, display_1L_chr = display_1L_chr, mkdn_tbl_ref_1L_chr = {
                  if (is.null(mkdn_tbl_refs_chr)) {
                    ""
                  }
                  else {
                    mkdn_tbl_refs_chr[1]
                  }
                }, output_type_1L_chr = "HTML", type_1L_chr = type_1L_chr, 
                use_lbls_as_col_nms_1L_lgl = T, use_rdocx_1L_lgl = F, 
                ... = ...)
        }
    }
})
