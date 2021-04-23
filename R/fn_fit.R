#' Fit complementary log log transformation
#' @description fit_clg_log_tfmn() is a Fit function that fits a model of a specified type to a dataset Specifically, this function implements an algorithm to fit complementary log log transformation. The function returns Model list (a list of models).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w_cloglog'
#' @param predr_vars_nms_chr Predictor variables names (a character vector)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param iters_1L_int Iterations (an integer vector of length one), Default: 4000
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @return Model list (a list of models)
#' @rdname fit_clg_log_tfmn
#' @export 

#' @keywords internal
fit_clg_log_tfmn <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w_cloglog", 
    predr_vars_nms_chr, id_var_nm_1L_chr = "fkClientID", backend_1L_chr = getOption("brms.backend", 
        "rstan"), iters_1L_int = 4000L, seed_1L_int = 1000L) 
{
    mdl_ls <- fit_ts_model_with_brm(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, link_1L_chr = "identity", 
        id_var_nm_1L_chr = id_var_nm_1L_chr, backend_1L_chr = backend_1L_chr, 
        iters_1L_int = iters_1L_int, seed_1L_int = seed_1L_int)
    return(mdl_ls)
}
#' Fit gaussian log lnk
#' @description fit_gsn_log_lnk() is a Fit function that fits a model of a specified type to a dataset Specifically, this function implements an algorithm to fit gaussian log lnk. The function returns Model list (a list of models).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one), Default: 'utl_total_w'
#' @param predr_vars_nms_chr Predictor variables names (a character vector)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param iters_1L_int Iterations (an integer vector of length one), Default: 4000
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @return Model list (a list of models)
#' @rdname fit_gsn_log_lnk
#' @export 

#' @keywords internal
fit_gsn_log_lnk <- function (data_tb, depnt_var_nm_1L_chr = "utl_total_w", predr_vars_nms_chr, 
    id_var_nm_1L_chr = "fkClientID", backend_1L_chr = getOption("brms.backend", 
        "rstan"), iters_1L_int = 4000L, seed_1L_int = 1000L) 
{
    mdl_ls <- fit_ts_model_with_brm(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, link_1L_chr = "log", 
        id_var_nm_1L_chr = id_var_nm_1L_chr, backend_1L_chr = backend_1L_chr, 
        iters_1L_int = iters_1L_int, seed_1L_int = seed_1L_int)
    return(mdl_ls)
}
#' Fit time series model with bayesian regression model
#' @description fit_ts_model_with_brm() is a Fit function that fits a model of a specified type to a dataset Specifically, this function implements an algorithm to fit time series model with bayesian regression model. The function returns Model list (a list of models).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one)
#' @param predr_vars_nms_chr Predictor variables names (a character vector)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one)
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param link_1L_chr Link (a character vector of length one), Default: 'identity'
#' @param iters_1L_int Iterations (an integer vector of length one), Default: 4000
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @return Model list (a list of models)
#' @rdname fit_ts_model_with_brm
#' @export 
#' @importFrom brms brm
#' @importFrom stats as.formula gaussian
#' @importFrom purrr map_chr
#' @keywords internal
fit_ts_model_with_brm <- function (data_tb, depnt_var_nm_1L_chr, predr_vars_nms_chr, id_var_nm_1L_chr, 
    backend_1L_chr = getOption("brms.backend", "rstan"), link_1L_chr = "identity", 
    iters_1L_int = 4000L, seed_1L_int = 1000L) 
{
    mdl_ls <- brms::brm(formula = stats::as.formula(paste0(depnt_var_nm_1L_chr, 
        " ~ ", purrr::map_chr(predr_vars_nms_chr, ~paste0(.x, 
            "_baseline + ", .x, "_change + ")) %>% paste0(collapse = ""), 
        "(1|", id_var_nm_1L_chr, ")")), backend = backend_1L_chr, 
        data = data_tb, family = stats::gaussian(link = link_1L_chr), 
        iter = iters_1L_int, seed = seed_1L_int)
    return(mdl_ls)
}
