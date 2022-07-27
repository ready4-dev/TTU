#' TTUParameters
#' 
#' Input Parameters For Generating Utility Mapping Models.
#' 
#' @slot a_ScorzProfile  (an instance of the ScorzProfile class)
#' @slot b_SpecificParameters  (an instance of the SpecificParameters class)
#' @slot paths_chr Paths (a character vector)
#' @slot prefd_mdls_chr Preferred models (a character vector)
#' @slot prefd_covars_chr Preferred covariates (a character vector)
#' @slot scndry_anlys_params_ls Secondary analysis parameters (a list)
#' @slot dissemination_1L_chr Dissemination (a character vector of length one)
#' @import ready4
#' @name TTUParameters-class
#' @rdname TTUParameters-class
#' @export TTUParameters
#' @exportClass TTUParameters
TTUParameters <- methods::setClass("TTUParameters",
contains = "Ready4Module",
slots = c(a_ScorzProfile = "ScorzProfile",b_SpecificParameters = "SpecificParameters",paths_chr = "character",prefd_mdls_chr = "character",prefd_covars_chr = "character",scndry_anlys_params_ls = "list",dissemination_1L_chr = "character"),
prototype =  list(a_ScorzProfile = scorz::ScorzProfile(),b_SpecificParameters = specific::SpecificParameters(),paths_chr = NA_character_,prefd_mdls_chr = NA_character_,prefd_covars_chr = NA_character_,scndry_anlys_params_ls = list(list())))


methods::setValidity(methods::className("TTUParameters"),
function(object){
msg <- NULL
if (is.null(msg)) TRUE else msg
})
