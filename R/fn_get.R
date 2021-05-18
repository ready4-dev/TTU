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
