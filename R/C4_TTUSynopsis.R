#' TTUSynopsis
#' 
#' Input, Output and Authorship Data For Generating Utility Mapping Study Reports.
#' 
#' @slot a_Ready4showPaths  (an instance of the Ready4showPaths class)
#' @slot b_SpecificResults  (an instance of the SpecificResults class)
#' @slot c_SpecificParameters  (an instance of the SpecificParameters class)
#' @slot d_YouthvarsProfile  (an instance of the YouthvarsProfile class)
#' @slot e_Ready4useRepos  (an instance of the Ready4useRepos class)
#' @slot abstract_args_ls Abstract arguments (a list)
#' @slot authors_r3 Authors (a ready4 submodule)
#' @slot background_1L_chr Background (a character vector of length one)
#' @slot coi_1L_chr Conflict of interest (a character vector of length one)
#' @slot conclusion_1L_chr Conclusion (a character vector of length one)
#' @slot correspondences_r3 Correspondences (a ready4 submodule)
#' @slot digits_int Digits (an integer vector)
#' @slot ethics_1L_chr Ethics (a character vector of length one)
#' @slot fl_nm_1L_chr File name (a character vector of length one)
#' @slot figures_in_body_lgl Figures in body (a logical vector)
#' @slot funding_1L_chr Funding (a character vector of length one)
#' @slot institutes_r3 Institutes (a ready4 submodule)
#' @slot interval_chr Interval (a character vector)
#' @slot keywords_chr Keywords (a character vector)
#' @slot outp_formats_chr Output formats (a character vector)
#' @slot rmd_fl_nms_ls R Markdown file names (a list)
#' @slot sample_desc_1L_chr Sample description (a character vector of length one)
#' @slot tables_in_body_lgl Tables in body (a logical vector)
#' @slot title_1L_chr Title (a character vector of length one)
#' @slot dissemination_1L_chr Dissemination (a character vector of length one)
#' @import specific
#' @name TTUSynopsis-class
#' @rdname TTUSynopsis-class
#' @export TTUSynopsis
#' @exportClass TTUSynopsis
TTUSynopsis <- methods::setClass("TTUSynopsis",
contains = "SpecificSynopsis",
slots = c(a_Ready4showPaths = "Ready4showPaths",b_SpecificResults = "SpecificResults",c_SpecificParameters = "SpecificParameters",d_YouthvarsProfile = "YouthvarsProfile",e_Ready4useRepos = "Ready4useRepos",abstract_args_ls = "list",authors_r3 = "ready4show_authors",background_1L_chr = "character",coi_1L_chr = "character",conclusion_1L_chr = "character",correspondences_r3 = "ready4show_correspondences",digits_int = "integer",ethics_1L_chr = "character",fl_nm_1L_chr = "character",figures_in_body_lgl = "logical",funding_1L_chr = "character",institutes_r3 = "ready4show_institutes",interval_chr = "character",keywords_chr = "character",outp_formats_chr = "character",rmd_fl_nms_ls = "list",sample_desc_1L_chr = "character",tables_in_body_lgl = "logical",title_1L_chr = "character",dissemination_1L_chr = "character"),
prototype =  list(a_Ready4showPaths = ready4show::Ready4showPaths(),b_SpecificResults = specific::SpecificResults(),c_SpecificParameters = specific::SpecificParameters(),d_YouthvarsProfile = youthvars::YouthvarsProfile(),e_Ready4useRepos = ready4use::Ready4useRepos()))


methods::setValidity(methods::className("TTUSynopsis"),
function(object){
msg <- NULL
if (is.null(msg)) TRUE else msg
})
