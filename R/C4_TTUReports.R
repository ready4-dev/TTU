#' TTUReports
#' 
#' Metadata to produce utility mapping study reports.
#' 
#' @include 
#' @slot a_SpecificSynopsis  (an instance of the SpecificSynopsis class)
#' @slot catalogue_tmpl_chr Catalogue template (a character vector)
#' @slot catalogue_fl_nms_ls Catalogue file names (a list)
#' @slot manuscript_tmpl_chr Manuscript template (a character vector)
#' @slot manuscript_fl_nms_ls Manuscript file names (a list)
#' @slot dissemination_1L_chr Dissemination (a character vector of length one)
#' @import ready4
#' @name TTUReports-class
#' @rdname TTUReports-class
#' @export TTUReports
#' @exportClass TTUReports
TTUReports <- methods::setClass("TTUReports",
contains = "Ready4Module",
slots = c(a_SpecificSynopsis = "SpecificSynopsis",catalogue_tmpl_chr = "character",catalogue_fl_nms_ls = "list",manuscript_tmpl_chr = "character",manuscript_fl_nms_ls = "list",dissemination_1L_chr = "character"),
prototype =  list(a_SpecificSynopsis = specific::SpecificSynopsis(),catalogue_tmpl_chr = c("https://github.com/ready4-dev/ttu_mdl_ctlg","0.0.9.3"),catalogue_fl_nms_ls = ready4show::make_rmd_fl_nms_ls("Lngl_Mdls_HTML",
                                                                                                                           pdf_fl_nm_1L_chr = "Lngl_Mdls_PDF",
                                                                                                                           word_fl_nm_1L_chr = "Lngl_Mdls_Word"),manuscript_tmpl_chr = c("https://github.com/ready4-dev/ttu_lng_ss","0.6"),manuscript_fl_nms_ls = ready4show::make_rmd_fl_nms_ls(pdf_fl_nm_1L_chr = "Main_PDF",
                                                                                                                            word_fl_nm_1L_chr = "Main_Word")))


methods::setValidity(methods::className("TTUReports"),
function(object){
msg <- NULL
if (is.null(msg)) TRUE else msg
})
