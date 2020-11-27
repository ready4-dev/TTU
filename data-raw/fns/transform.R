transform_data_tb_for_cmprsn <- function (data_tb, mdl, dep_var_nm_1L_chr = "aqol6d_total_w", 
    source_data_nm_1L_chr = "Original", tf_type_1L_chr = "Predicted", 
    pred_type_1L_chr = NULL, tfmn_for_bnml_1L_lgl = F, family_1L_chr = NA_character_) 
{
    if (tf_type_1L_chr == "Predicted") 
        new_data_dbl <- predict(mdl, type = pred_type_1L_chr)
    if (tf_type_1L_chr == "Simulated" & !tfmn_for_bnml_1L_lgl) 
        new_data_dbl <- simulate(mdl)$sim_1
    if (tf_type_1L_chr == "Simulated" & tfmn_for_bnml_1L_lgl) 
        new_data_dbl <- (predict(mdl) + rnorm(nrow(data_tb), 
            0, sigma(mdl))) %>% calculate_dep_var_tfmn(tfmn_1L_chr = ifelse(family_1L_chr == 
            "quasibinomial(log)", "LOG", ifelse(family_1L_chr == 
            "quasibinomial(logit)", "LOGIT", ifelse(family_1L_chr == 
            "quasibinomial(cloglog)", "CLL", "NTF"))), tfmn_is_outp_1L_lgl = T)
    tfd_data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(tf_type_1L_chr), 
        new_data_dbl), `:=`(!!rlang::sym(source_data_nm_1L_chr), 
        !!rlang::sym(dep_var_nm_1L_chr)))
    return(tfd_data_tb)
}
transform_dep_var_nm <- function (dep_var_nm_1L_chr, tfmn_1L_chr = "NTF") 
{
    tfd_dep_var_nm_1L_chr <- paste0(dep_var_nm_1L_chr, ifelse(tfmn_1L_chr == 
        "NTF", "", paste0("_", tfmn_1L_chr)))
    return(tfd_dep_var_nm_1L_chr)
}
transform_dep_var_nm_for_cll <- function (dep_var_nm_1L_chr) 
{
    tfd_dep_var_nm_1L_chr <- paste0(dep_var_nm_1L_chr, "_cloglog")
    return(tfd_dep_var_nm_1L_chr)
}
transform_ds_for_mdlng <- function (data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", predr_var_nm_1L_chr, 
    covar_var_nms_chr = NA_character_) 
{
    mdl_vars_chr <- c(names(data_tb)[names(data_tb) %>% startsWith(dep_var_nm_1L_chr)], 
        predr_var_nm_1L_chr, covar_var_nms_chr) %>% purrr::discard(is.na)
    tfd_data_tb <- data_tb %>% tidyr::drop_na(!!!rlang::syms(mdl_vars_chr)) %>% 
        dplyr::select(!!!rlang::syms(mdl_vars_chr))
    return(tfd_data_tb)
}
transform_tb_to_mdl_inp <- function (data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", predr_vars_nms_chr, 
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
    round_bl_val_1L_chr = "Baseline") 
{
    data_tb <- data.frame(data_tb)
    tfd_for_gsn_log_mdl_tb <- data_tb %>% dplyr::select(dplyr::all_of(dep_var_nm_1L_chr), 
        dplyr::all_of(id_var_nm_1L_chr), dplyr::all_of(round_var_nm_1L_chr), 
        dplyr::all_of(predr_vars_nms_chr)) %>% dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>% 
        dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr), !!rlang::sym(round_var_nm_1L_chr))
    tfd_for_gsn_log_mdl_tb <- tfd_for_gsn_log_mdl_tb %>% dplyr::mutate(dplyr::across(dplyr::all_of(predr_vars_nms_chr), 
        .fns = list(baseline = ~dplyr::first(.)/100, change = ~ifelse(round == 
            "Baseline", 0, (. - dplyr::lag(.))/100)))) %>% dplyr::mutate(`:=`(!!rlang::sym(transform_dep_var_nm_for_cll(dep_var_nm_1L_chr)), 
        log(-log(1 - ifelse(!!rlang::sym(dep_var_nm_1L_chr) == 
            1, 0.999, !!rlang::sym(dep_var_nm_1L_chr)))))) %>% 
        na.omit()
    return(tfd_for_gsn_log_mdl_tb)
}
