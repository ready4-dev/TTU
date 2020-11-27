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
