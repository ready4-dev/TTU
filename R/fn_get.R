#' Get link from transformation
#' @description get_link_from_tfmn() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get link from transformation. Function argument tfmn_1L_chr specifies the where to look for the required object. The function returns Link (a character vector of length one).
#' @param tfmn_1L_chr Transformation (a character vector of length one)
#' @param is_OLS_1L_lgl Is OLS (a logical vector of length one), Default: F
#' @return Link (a character vector of length one)
#' @rdname get_link_from_tfmn
#' @export 

#' @keywords internal
get_link_from_tfmn <- function (tfmn_1L_chr, is_OLS_1L_lgl = F) 
{
    link_1L_chr <- ifelse(is_OLS_1L_lgl, "identity", ifelse(tfmn_1L_chr == 
        "LOG", "log", ifelse(tfmn_1L_chr == "LGT", "logit", ifelse(tfmn_1L_chr == 
        "CLL", "cloglog", ifelse(tfmn_1L_chr == "LOGLOG", "loglog", 
        ifelse(tfmn_1L_chr == "NTF", "identity", "ERROR"))))))
    if (link_1L_chr == "ERROR") 
        stop("Link cannot be identified - incorrect transformation argument tfmn_1L_chr")
    return(link_1L_chr)
}
#' Get significant covariates
#' @description get_signft_covars() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get significant covariates. Function argument mdls_with_covars_smry_tb specifies the where to look for the required object. The function returns Signt covariates (a character vector).
#' @param mdls_with_covars_smry_tb Models with covariates summary (a tibble)
#' @param covar_var_nms_chr Covariate variable names (a character vector)
#' @return Signt covariates (a character vector)
#' @rdname get_signft_covars
#' @export 
#' @importFrom purrr map flatten flatten_chr
get_signft_covars <- function (mdls_with_covars_smry_tb, covar_var_nms_chr) 
{
    signif_vars_chr <- mdls_with_covars_smry_tb$Significant %>% 
        purrr::map(~strsplit(.x, " ")) %>% purrr::flatten() %>% 
        purrr::flatten_chr() %>% unique()
    signt_covars_chr <- covar_var_nms_chr[covar_var_nms_chr %in% 
        signif_vars_chr]
    if (identical(signt_covars_chr, character(0))) 
        signt_covars_chr <- NA_character_
    return(signt_covars_chr)
}
