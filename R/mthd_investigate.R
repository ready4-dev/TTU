#' 
#' Investigate solutions to an inverse problem
#' @name investigate-TTUProject
#' @description investigate method applied to TTUProject
#' @param x An object of class TTUProject
#' @param backend_1L_chr Backend (a character vector of length one), Default: 'cmdstanr'
#' @param combinations_1L_lgl Combinations (a logical vector of length one), Default: F
#' @param consent_1L_chr Consent (a character vector of length one), Default: ''
#' @param cores_1L_int Cores (an integer vector of length one), Default: 1
#' @param depnt_var_max_val_1L_dbl Dependent variable maximum value (a double vector of length one), Default: numeric(0)
#' @param depnt_var_min_val_1L_dbl Dependent variable minimum value (a double vector of length one), Default: numeric(0)
#' @param existing_predrs_ls Existing predictors (a list), Default: NULL
#' @param max_nbr_of_covars_1L_int Maximum number of covariates (an integer vector of length one), Default: integer(0)
#' @param new_dir_nm_1L_chr New directory name (a character vector of length one), Default: 'F_TS_Mdls'
#' @param scndry_anlys_params_ls Secondary analysis parameters (a list), Default: NULL
#' @param session_ls Session (a list), Default: NULL
#' @param signft_covars_cdn_1L_chr Significant covariates condition (a character vector of length one), Default: 'any'
#' @param ... Additional arguments
#' @return x (An object of class TTUProject)
#' @rdname investigate-methods
#' @aliases investigate,TTUProject-method
#' @export 
#' @importFrom rlang exec
#' @importFrom ready4 investigate
methods::setMethod("investigate", "TTUProject", function (x, backend_1L_chr = "cmdstanr", combinations_1L_lgl = F, 
    consent_1L_chr = "", cores_1L_int = 1L, depnt_var_max_val_1L_dbl = numeric(0), 
    depnt_var_min_val_1L_dbl = numeric(0), existing_predrs_ls = NULL, 
    max_nbr_of_covars_1L_int = integer(0), new_dir_nm_1L_chr = "F_TS_Mdls", 
    scndry_anlys_params_ls = NULL, session_ls = NULL, signft_covars_cdn_1L_chr = "any", 
    ...) 
{
    args_ls <- list(...)
    args_ls <- append(list(slot_nm_1L_chr = "c_SpecificProject", 
        consent_1L_chr = consent_1L_chr), args_ls)
    if (inherits(x@c_SpecificProject, what = "SpecificModels") & 
        !(inherits(x@c_SpecificProject, what = "SpecificFixed") | 
            inherits(x@c_SpecificProject, what = "SpecificMixed") | 
            inherits(x@c_SpecificProject, what = "SpecificPredictors"))) {
        args_ls <- append(list(depnt_var_max_val_1L_dbl = ifelse(identical(depnt_var_max_val_1L_dbl, 
            numeric(0)), 1, depnt_var_max_val_1L_dbl), depnt_var_min_val_1L_dbl = ifelse(identical(depnt_var_max_val_1L_dbl, 
            numeric(0)), -1, depnt_var_max_val_1L_dbl), session_ls = session_ls), 
            args_ls)
    }
    if (inherits(x@c_SpecificProject, what = "SpecificPredictors") & 
        !(inherits(x@c_SpecificProject, what = "SpecificFixed") | 
            inherits(x@c_SpecificProject, what = "SpecificMixed"))) {
        args_ls <- append(list(depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
            signft_covars_cdn_1L_chr = signft_covars_cdn_1L_chr), 
            args_ls)
    }
    if (inherits(x@c_SpecificProject, what = "SpecificFixed") & 
        !(inherits(x@c_SpecificProject, what = "SpecificMixed"))) {
        args_ls <- append(list(combinations_1L_lgl = combinations_1L_lgl, 
            depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
            existing_predrs_ls = existing_predrs_ls, max_nbr_of_covars_1L_int = max_nbr_of_covars_1L_int), 
            args_ls)
    }
    if (inherits(x@c_SpecificProject, what = "SpecificMixed")) {
        args_ls <- append(list(backend_1L_chr = backend_1L_chr, 
            combinations_1L_lgl = combinations_1L_lgl, cores_1L_int = cores_1L_int, 
            depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
            existing_predrs_ls = existing_predrs_ls, max_nbr_of_covars_1L_int = max_nbr_of_covars_1L_int, 
            new_dir_nm_1L_chr = new_dir_nm_1L_chr, scndry_anlys_params_ls = scndry_anlys_params_ls), 
            args_ls)
    }
    new_val_xx <- rlang::exec(investigateSlot, x, !!!args_ls)
    x <- renewSlot(x, "c_SpecificProject", new_val_xx)
    return(x)
})
