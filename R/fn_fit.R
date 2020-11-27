#' Fit clg log tfmn
#' @description fit_clg_log_tfmn() is a Fit function that fits a model of a specified type to a dataset Specifically, this function implements an algorithm to fit clg log tfmn. The function returns Mdl (a list).
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w_cloglog'
#' @param predr_vars_nms_chr Predr vars names (a character vector)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param iters_1L_int Iters (an integer vector of length one), Default: 4000
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @return Mdl (a list)
#' @rdname fit_clg_log_tfmn
#' @export 

#' @keywords internal
fit_clg_log_tfmn <- function (data_tb, dep_var_nm_1L_chr = "aqol6d_total_w_cloglog", 
    predr_vars_nms_chr, id_var_nm_1L_chr = "fkClientID", backend_1L_chr = getOption("brms.backend", 
        "rstan"), iters_1L_int = 4000L, seed_1L_int = 1000L) 
{
    mdl_ls <- fit_ts_model_with_brm(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, link_1L_chr = "identity", 
        id_var_nm_1L_chr = id_var_nm_1L_chr, backend_1L_chr = backend_1L_chr, 
        iters_1L_int = iters_1L_int, seed_1L_int = seed_1L_int)
    return(mdl_ls)
}
#' Fit gsn log lnk
#' @description fit_gsn_log_lnk() is a Fit function that fits a model of a specified type to a dataset Specifically, this function implements an algorithm to fit gsn log lnk. The function returns Mdl (a list).
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param predr_vars_nms_chr Predr vars names (a character vector)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param iters_1L_int Iters (an integer vector of length one), Default: 4000
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @return Mdl (a list)
#' @rdname fit_gsn_log_lnk
#' @export 

#' @keywords internal
fit_gsn_log_lnk <- function (data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", predr_vars_nms_chr, 
    id_var_nm_1L_chr = "fkClientID", backend_1L_chr = getOption("brms.backend", 
        "rstan"), iters_1L_int = 4000L, seed_1L_int = 1000L) 
{
    mdl_ls <- fit_ts_model_with_brm(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, link_1L_chr = "log", 
        id_var_nm_1L_chr = id_var_nm_1L_chr, backend_1L_chr = backend_1L_chr, 
        iters_1L_int = iters_1L_int, seed_1L_int = seed_1L_int)
    return(mdl_ls)
}
