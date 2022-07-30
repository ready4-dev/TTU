#' TTUProject
#' 
#' Input And Output Data For Undertaking and Reporting Utility Mapping Studies.
#' 
#' @include C4_TTUReports.R
#' @slot a_ScorzProfile  (an instance of the ScorzProfile class)
#' @slot b_SpecificParameters  (an instance of the SpecificParameters class)
#' @slot c_SpecificProject  (an instance of the SpecificProject class)
#' @slot d_TTUReports  (an instance of the TTUReports class)
#' @slot dissemination_1L_chr Dissemination (a character vector of length one)
#' @import ready4
#' @name TTUProject-class
#' @rdname TTUProject-class
#' @export TTUProject
#' @exportClass TTUProject
TTUProject <- methods::setClass("TTUProject",
contains = "Ready4Module",
slots = c(a_ScorzProfile = "ScorzProfile",b_SpecificParameters = "SpecificParameters",c_SpecificProject = "SpecificProject",d_TTUReports = "TTUReports",dissemination_1L_chr = "character"),
prototype =  list(a_ScorzProfile = scorz::ScorzProfile(),b_SpecificParameters = specific::SpecificParameters(),c_SpecificProject = specific::SpecificProject(),d_TTUReports = TTUReports()))


methods::setValidity(methods::className("TTUProject"),
function(object){
msg <- NULL
if (is.null(msg)) TRUE else msg
})
