#' Transform data tibble for cmprsn
#' @description transform_data_tb_for_cmprsn() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform data tibble for cmprsn. Function argument data_tb specifies the object to be updated. Argument model_mdl provides the object to be updated. The function returns Transformed data (a tibble).
#' @param data_tb Data (a tibble)
#' @param model_mdl PARAM_DESCRIPTION
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'utl_total_w'
#' @param source_data_nm_1L_chr Source data name (a character vector of length one), Default: 'Original'
#' @param tf_type_1L_chr Transform type (a character vector of length one), Default: 'Predicted'
#' @param pred_type_1L_chr Pred type (a character vector of length one), Default: NULL
#' @param tfmn_for_bnml_1L_lgl Transformation for bnml (a logical vector of length one), Default: F
#' @param family_1L_chr Family (a character vector of length one), Default: 'NA'
#' @return Transformed data (a tibble)
#' @rdname transform_data_tb_for_cmprsn
#' @export 
#' @importFrom stats predict simulate rnorm sigma
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @keywords internal
transform_data_tb_for_cmprsn <- function (data_tb, model_mdl, dep_var_nm_1L_chr = "utl_total_w", 
    source_data_nm_1L_chr = "Original", tf_type_1L_chr = "Predicted", 
    pred_type_1L_chr = NULL, tfmn_for_bnml_1L_lgl = F, family_1L_chr = NA_character_) 
{
    if (tf_type_1L_chr == "Predicted") 
        new_data_dbl <- stats::predict(model_mdl, type = pred_type_1L_chr)
    if (tf_type_1L_chr == "Simulated" & !tfmn_for_bnml_1L_lgl) 
        new_data_dbl <- stats::simulate(model_mdl)$sim_1
    if (tf_type_1L_chr == "Simulated" & tfmn_for_bnml_1L_lgl) 
        new_data_dbl <- (stats::predict(model_mdl) + stats::rnorm(nrow(data_tb), 
            0, stats::sigma(model_mdl))) %>% calculate_dep_var_tfmn(tfmn_1L_chr = ifelse(family_1L_chr == 
            "quasibinomial(log)", "LOG", ifelse(family_1L_chr == 
            "quasibinomial(logit)", "LOGIT", ifelse(family_1L_chr == 
            "quasibinomial(cloglog)", "CLL", "NTF"))), tfmn_is_outp_1L_lgl = T)
    tfd_data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(tf_type_1L_chr), 
        new_data_dbl), `:=`(!!rlang::sym(source_data_nm_1L_chr), 
        !!rlang::sym(dep_var_nm_1L_chr)))
    return(tfd_data_tb)
}
#' Transform dep var name
#' @description transform_dep_var_nm() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dep var name. Function argument dep_var_nm_1L_chr specifies the object to be updated. Argument tfmn_1L_chr provides the object to be updated. The function returns Transformed dep var name (a character vector of length one).
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one)
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @return Transformed dep var name (a character vector of length one)
#' @rdname transform_dep_var_nm
#' @export 

#' @keywords internal
transform_dep_var_nm <- function (dep_var_nm_1L_chr, tfmn_1L_chr = "NTF") 
{
    tfd_dep_var_nm_1L_chr <- paste0(dep_var_nm_1L_chr, ifelse(tfmn_1L_chr == 
        "NTF", "", paste0("_", tfmn_1L_chr)))
    return(tfd_dep_var_nm_1L_chr)
}
#' Transform dep var name for cll
#' @description transform_dep_var_nm_for_cll() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dep var name for cll. Function argument dep_var_nm_1L_chr specifies the object to be updated. The function returns Transformed dep var name (a character vector of length one).
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one)
#' @return Transformed dep var name (a character vector of length one)
#' @rdname transform_dep_var_nm_for_cll
#' @export 

#' @keywords internal
transform_dep_var_nm_for_cll <- function (dep_var_nm_1L_chr) 
{
    tfd_dep_var_nm_1L_chr <- paste0(dep_var_nm_1L_chr, "_cloglog")
    return(tfd_dep_var_nm_1L_chr)
}
#' Transform dataset for mdlng
#' @description transform_ds_for_mdlng() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dataset for mdlng. Function argument data_tb specifies the object to be updated. Argument dep_var_nm_1L_chr provides the object to be updated. The function returns Transformed data (a tibble).
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'utl_total_w'
#' @param predr_var_nm_1L_chr Predr var name (a character vector of length one)
#' @param covar_var_nms_chr Covar var names (a character vector), Default: 'NA'
#' @return Transformed data (a tibble)
#' @rdname transform_ds_for_mdlng
#' @export 
#' @importFrom purrr discard
#' @importFrom tidyr drop_na
#' @importFrom rlang syms
#' @importFrom dplyr select
#' @keywords internal
transform_ds_for_mdlng <- function (data_tb, dep_var_nm_1L_chr = "utl_total_w", predr_var_nm_1L_chr, 
    covar_var_nms_chr = NA_character_) 
{
    mdl_vars_chr <- c(names(data_tb)[names(data_tb) %>% startsWith(dep_var_nm_1L_chr)], 
        predr_var_nm_1L_chr, covar_var_nms_chr) %>% purrr::discard(is.na)
    tfd_data_tb <- data_tb %>% tidyr::drop_na(!!!rlang::syms(mdl_vars_chr)) %>% 
        dplyr::select(!!!rlang::syms(mdl_vars_chr))
    return(tfd_data_tb)
}
#' Transform tibble to model input
#' @description transform_tb_to_mdl_inp() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform tibble to model input. Function argument data_tb specifies the object to be updated. Argument dep_var_nm_1L_chr provides the object to be updated. The function returns Transformed for gsn log model (a tibble).
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'utl_total_w'
#' @param predr_vars_nms_chr Predr vars names (a character vector)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round bl value (a character vector of length one), Default: 'Baseline'
#' @return Transformed for gsn log model (a tibble)
#' @rdname transform_tb_to_mdl_inp
#' @export 
#' @importFrom dplyr select all_of group_by arrange mutate across first lag
#' @importFrom rlang sym
#' @importFrom stats na.omit
#' @keywords internal
transform_tb_to_mdl_inp <- function (data_tb, dep_var_nm_1L_chr = "utl_total_w", predr_vars_nms_chr, 
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
        stats::na.omit()
    return(tfd_for_gsn_log_mdl_tb)
}
