#' Add preferred predictor variable to model summary
#' @description add_prefd_predr_var_to_mdl_smry_ls() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add preferred predictor variable to model summary list. Function argument mdl_smry_ls specifies the object to be updated. The function returns Model summary (a list).
#' @param mdl_smry_ls Model summary (a list)
#' @param ds_smry_ls Dataset summary (a list)
#' @return Model summary (a list)
#' @rdname add_prefd_predr_var_to_mdl_smry_ls
#' @export 
#' @importFrom ready4fun get_from_lup_obj
#' @keywords internal
add_prefd_predr_var_to_mdl_smry_ls <- function (mdl_smry_ls, ds_smry_ls) 
{
    mdl_smry_ls$predr_var_nm_1L_chr <- ds_smry_ls$candidate_predrs_chr[1]
    mdl_smry_ls$predr_var_desc_1L_chr <- ds_smry_ls$predictors_lup %>% 
        ready4fun::get_from_lup_obj(match_value_xx = mdl_smry_ls$predr_var_nm_1L_chr, 
            match_var_nm_1L_chr = "short_name_chr", target_var_nm_1L_chr = "long_name_chr", 
            evaluate_lgl = F)
    mdl_smry_ls$predr_vals_dbl <- make_predr_vals(mdl_smry_ls$predr_var_nm_1L_chr, 
        candidate_predrs_lup = ds_smry_ls$predictors_lup)
    return(mdl_smry_ls)
}
#' Add tfmd variable to dataset
#' @description add_tfmd_var_to_ds() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add tfmd variable to dataset. Function argument data_tb specifies the object to be updated. The function returns Data (a tibble).
#' @param data_tb Data (a tibble)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one)
#' @param tfmn_1L_chr Transformation (a character vector of length one)
#' @return Data (a tibble)
#' @rdname add_tfmd_var_to_ds
#' @export 
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @keywords internal
add_tfmd_var_to_ds <- function (data_tb, depnt_var_nm_1L_chr, tfmn_1L_chr) 
{
    data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(transform_depnt_var_nm(depnt_var_nm_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr)), !!rlang::sym(depnt_var_nm_1L_chr) %>% 
        calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr)))
    return(data_tb)
}
#' Add unique identifiers to tibbles
#' @description add_uids_to_tbs_ls() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add unique identifiers to tibbles list. Function argument tbs_ls specifies the object to be updated. The function returns Tibbles (a list).
#' @param tbs_ls Tibbles (a list)
#' @param prefix_1L_chr Prefix (a character vector of length one)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @return Tibbles (a list)
#' @rdname add_uids_to_tbs_ls
#' @export 
#' @importFrom purrr map map_chr
#' @importFrom dplyr mutate arrange
#' @importFrom rlang sym
#' @importFrom tidyselect all_of
#' @importFrom stringr str_replace
#' @importFrom stats setNames
#' @keywords internal
add_uids_to_tbs_ls <- function (tbs_ls, prefix_1L_chr, id_var_nm_1L_chr = "fkClientID") 
{
    participant_ids <- paste0(prefix_1L_chr, 1:nrow(tbs_ls[[1]])) %>% 
        sample(nrow(tbs_ls[[1]]))
    tbs_ls <- purrr::map(tbs_ls, ~{
        .x %>% dplyr::mutate(`:=`(!!rlang::sym(id_var_nm_1L_chr), 
            tidyselect::all_of(participant_ids[1:nrow(.x)]))) %>% 
            dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr) %>% 
                purrr::map_chr(~stringr::str_replace(.x, prefix_1L_chr, 
                  "")) %>% as.numeric())
    }) %>% stats::setNames(names(tbs_ls))
    return(tbs_ls)
}
#' Add utility prediction to dataset
#' @description add_utility_predn_to_ds() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add utility prediction to dataset. Function argument data_tb specifies the object to be updated. The function returns Data (a tibble).
#' @param data_tb Data (a tibble)
#' @param model_mdl Model (a model)
#' @param tfmn_1L_chr Transformation (a character vector of length one)
#' @param depnt_var_nm_1L_chr Dependent variable name (a character vector of length one)
#' @param predr_vars_nms_chr Predictor variables names (a character vector), Default: NULL
#' @param force_min_max_1L_lgl Force minimum maximum (a logical vector of length one), Default: T
#' @param impute_1L_lgl Impute (a logical vector of length one), Default: T
#' @param new_data_is_1L_chr New data is (a character vector of length one), Default: 'Predicted'
#' @param rmv_tfd_depnt_var_1L_lgl Remove transformed dependent variable (a logical vector of length one), Default: F
#' @param utl_cls_fn Utility class (a function), Default: NULL
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: 0.03
#' @param force_new_data_1L_lgl Force new data (a logical vector of length one), Default: F
#' @param is_brms_mdl_1L_lgl Is bayesian regression models model (a logical vector of length one), Default: T
#' @param sd_dbl Standard deviation (a double vector), Default: NA
#' @return Data (a tibble)
#' @rdname add_utility_predn_to_ds
#' @export 
#' @importFrom purrr reduce map flatten_chr
#' @importFrom dplyr mutate select
#' @importFrom rlang sym
#' @importFrom tidyselect all_of
add_utility_predn_to_ds <- function (data_tb, model_mdl, tfmn_1L_chr, depnt_var_nm_1L_chr, 
    predr_vars_nms_chr = NULL, force_min_max_1L_lgl = T, impute_1L_lgl = T, 
    new_data_is_1L_chr = "Predicted", rmv_tfd_depnt_var_1L_lgl = F, 
    utl_cls_fn = NULL, utl_min_val_1L_dbl = 0.03, force_new_data_1L_lgl = F, 
    is_brms_mdl_1L_lgl = T, sd_dbl = NA_real_) 
{
    dep_vars_chr <- c(depnt_var_nm_1L_chr, transform_depnt_var_nm(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr)) %>% unique()
    data_tb <- purrr::reduce(dep_vars_chr, .init = data_tb, ~dplyr::mutate(.x, 
        `:=`(!!rlang::sym(.y), NA_real_)))
    predictions_dbl <- predict_utility(data_tb = data_tb, tfmn_1L_chr = tfmn_1L_chr, 
        model_mdl = model_mdl, force_min_max_1L_lgl = force_min_max_1L_lgl, 
        utl_min_val_1L_dbl = utl_min_val_1L_dbl, impute_1L_lgl = impute_1L_lgl, 
        new_data_is_1L_chr = new_data_is_1L_chr, utl_cls_fn = utl_cls_fn, 
        is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl, force_new_data_1L_lgl = force_new_data_1L_lgl, 
        sd_dbl = NA_real_)
    data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(depnt_var_nm_1L_chr), 
        predictions_dbl))
    if (!is.null(predr_vars_nms_chr)) {
        data_tb <- data_tb %>% dplyr::select(-tidyselect::all_of(purrr::map(predr_vars_nms_chr, 
            ~paste0(.x, c("_baseline", "_change"))) %>% purrr::flatten_chr()))
    }
    if (rmv_tfd_depnt_var_1L_lgl) {
        data_tb <- data_tb %>% dplyr::select(-tidyselect::all_of(dep_vars_chr[dep_vars_chr != 
            depnt_var_nm_1L_chr]))
    }
    return(data_tb)
}
