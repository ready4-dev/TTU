#' Get candidates for mxd models
#' @description get_cndts_for_mxd_mdls() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get candidates for mxd models. Function argument mdl_types_lup specifies the where to look for the required object. The function returns Candidates for mxd models (a lookup table).
#' @param mdl_types_lup Model types (a lookup table), Default: NULL
#' @return Candidates for mxd models (a lookup table)
#' @rdname get_cndts_for_mxd_mdls
#' @export 
#' @importFrom utils data
#' @importFrom dplyr filter
get_cndts_for_mxd_mdls <- function (mdl_types_lup = NULL) 
{
    if (is.null(mdl_types_lup)) 
        utils::data("mdl_types_lup", package = "TTU", envir = environment())
    cndts_for_mxd_mdls_lup <- mdl_types_lup %>% dplyr::filter(!tfmn_for_bnml_lgl, 
        short_name_chr != "BET_LOG")
    return(cndts_for_mxd_mdls_lup)
}
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
#' Get random intercept
#' @description get_random_intercept() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get random intercept. Function argument mdls_smry_tb specifies the where to look for the required object. The function returns Standard deviation (a double vector).
#' @param mdls_smry_tb Models summary (a tibble)
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param deterministic_1L_lgl Deterministic (a logical vector of length one), Default: T
#' @return Standard deviation (a double vector)
#' @rdname get_random_intercept
#' @export 
#' @importFrom dplyr filter
#' @importFrom ready4fun get_from_lup_obj
#' @keywords internal
get_random_intercept <- function (mdls_smry_tb, mdl_nm_1L_chr, deterministic_1L_lgl = T) 
{
    mdl_smry_tb <- mdls_smry_tb %>% dplyr::filter(Model == mdl_nm_1L_chr)
    sd_dbl <- c(mdl_smry_tb %>% ready4fun::get_from_lup_obj(match_value_xx = "SD (Intercept)", 
        match_var_nm_1L_chr = "Parameter", target_var_nm_1L_chr = "Estimate", 
        evaluate_lgl = F), ifelse(deterministic_1L_lgl, 0, mdl_smry_tb %>% 
        ready4fun::get_from_lup_obj(match_value_xx = "SD (Intercept)", 
            match_var_nm_1L_chr = "Parameter", target_var_nm_1L_chr = "SE", 
            evaluate_lgl = F)))
    return(sd_dbl)
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
