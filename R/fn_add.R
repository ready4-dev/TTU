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
#' @param dep_var_max_val_1L_dbl Dep variable maximum value (a double vector of length one), Default: NULL
#' @return Data (a tibble)
#' @rdname add_tfmd_var_to_ds
#' @export 
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @keywords internal
add_tfmd_var_to_ds <- function (data_tb, depnt_var_nm_1L_chr, tfmn_1L_chr, dep_var_max_val_1L_dbl = NULL) 
{
    data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(transform_depnt_var_nm(depnt_var_nm_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr)), !!rlang::sym(depnt_var_nm_1L_chr) %>% 
        calculate_dpnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, tfmn_is_outp_1L_lgl = F, 
            dep_var_max_val_1L_dbl = dep_var_max_val_1L_dbl)))
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
#' @param force_min_max_1L_lgl Force minimum maximum (a logical vector of length one), Default: T
#' @param force_new_data_1L_lgl Force new data (a logical vector of length one), Default: F
#' @param impute_1L_lgl Impute (a logical vector of length one), Default: T
#' @param is_brms_mdl_1L_lgl Is bayesian regression models model (a logical vector of length one), Default: T
#' @param new_data_is_1L_chr New data is (a character vector of length one), Default: 'Predicted'
#' @param predn_type_1L_chr Prediction type (a character vector of length one), Default: NULL
#' @param predr_vars_nms_chr Predictor variables names (a character vector), Default: NULL
#' @param rmv_tfd_depnt_var_1L_lgl Remove transformed dependent variable (a logical vector of length one), Default: F
#' @param sd_dbl Standard deviation (a double vector), Default: NA
#' @param utl_cls_fn Utility class (a function), Default: NULL
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: -1
#' @return Data (a tibble)
#' @rdname add_utility_predn_to_ds
#' @export 
#' @importFrom purrr reduce map flatten_chr
#' @importFrom dplyr mutate select
#' @importFrom rlang sym
#' @importFrom tidyselect all_of
add_utility_predn_to_ds <- function (data_tb, model_mdl, tfmn_1L_chr, depnt_var_nm_1L_chr, 
    force_min_max_1L_lgl = T, force_new_data_1L_lgl = F, impute_1L_lgl = T, 
    is_brms_mdl_1L_lgl = T, new_data_is_1L_chr = "Predicted", 
    predn_type_1L_chr = NULL, predr_vars_nms_chr = NULL, rmv_tfd_depnt_var_1L_lgl = F, 
    sd_dbl = NA_real_, utl_cls_fn = NULL, utl_min_val_1L_dbl = -1) 
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
        predn_type_1L_chr = predn_type_1L_chr, sd_dbl = sd_dbl)
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
#' Add utility prediction to new dataset
#' @description add_utl_predn_to_new_ds() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add utility prediction to new dataset. Function argument data_tb specifies the object to be updated. The function returns Updated (a tibble).
#' @param data_tb Data (a tibble)
#' @param ingredients_ls Ingredients (a list)
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param analysis_1L_chr Analysis (a character vector of length one), Default: NULL
#' @param deterministic_1L_lgl Deterministic (a logical vector of length one), Default: T
#' @param force_min_max_1L_lgl Force minimum maximum (a logical vector of length one), Default: T
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: NULL
#' @param model_mdl Model (a model), Default: NULL
#' @param new_data_is_1L_chr New data is (a character vector of length one), Default: 'Simulated'
#' @param predr_vars_nms_chr Predictor variables names (a character vector), Default: NULL
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'Timepoint'
#' @param round_bl_val_1L_chr Round baseline value (a character vector of length one), Default: 'BL'
#' @param utl_cls_fn Utility class (a function), Default: NULL
#' @param utl_var_nm_1L_chr Utility variable name (a character vector of length one), Default: NULL
#' @return Updated (a tibble)
#' @rdname add_utl_predn_to_new_ds
#' @export 
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom dplyr filter pull rename left_join select
#' @importFrom purrr flatten_chr
#' @importFrom rlang sym
#' @importFrom tidyselect all_of
#' @keywords internal
add_utl_predn_to_new_ds <- function (data_tb, ingredients_ls, mdl_nm_1L_chr, analysis_1L_chr = NULL, 
    deterministic_1L_lgl = T, force_min_max_1L_lgl = T, id_var_nm_1L_chr = NULL, 
    model_mdl = NULL, new_data_is_1L_chr = "Simulated", predr_vars_nms_chr = NULL, 
    round_var_nm_1L_chr = "Timepoint", round_bl_val_1L_chr = "BL", 
    utl_cls_fn = NULL, utl_var_nm_1L_chr = NULL) 
{
    if (is.null(model_mdl)) 
        model_mdl <- get_table_predn_mdl(mdl_nm_1L_chr, ingredients_ls = ingredients_ls, 
            analysis_1L_chr = analysis_1L_chr)
    mdl_type_1L_chr <- get_mdl_type_from_nm(mdl_nm_1L_chr)
    tfmn_1L_chr <- ready4fun::get_from_lup_obj(ingredients_ls$mdl_types_lup, 
        match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
        target_var_nm_1L_chr = "tfmn_chr", evaluate_lgl = F)
    predn_type_1L_chr <- ready4fun::get_from_lup_obj(ingredients_ls$mdl_types_lup, 
        match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
        target_var_nm_1L_chr = "predn_type_chr", evaluate_lgl = F)
    if (is.na(predn_type_1L_chr)) 
        predn_type_1L_chr <- NULL
    id_var_nm_1L_chr <- ifelse(is.null(id_var_nm_1L_chr), ingredients_ls$id_var_nm_1L_chr, 
        id_var_nm_1L_chr)
    if (!is.null(predr_vars_nms_chr)) {
        data_tb <- rename_from_nmd_vec(data_tb, nmd_vec_chr = predr_vars_nms_chr, 
            vec_nms_as_new_1L_lgl = T)
    }
    mdl_predr_terms_chr <- ingredients_ls$mdls_lup %>% dplyr::filter(mdl_nms_chr == 
        mdl_nm_1L_chr) %>% dplyr::pull(predrs_ls) %>% purrr::flatten_chr()
    original_ds_vars_chr <- names(data_tb)[!names(data_tb) %in% 
        c(mdl_predr_terms_chr, ifelse(!is.null(utl_var_nm_1L_chr), 
            utl_var_nm_1L_chr, ingredients_ls$depnt_var_nm_1L_chr))]
    updated_tb <- data_tb %>% transform_ds_to_predn_ds(predr_vars_nms_chr = mdl_predr_terms_chr, 
        tfmn_1L_chr = tfmn_1L_chr, depnt_var_nm_1L_chr = ingredients_ls$depnt_var_nm_1L_chr, 
        id_var_nm_1L_chr = id_var_nm_1L_chr, round_var_nm_1L_chr = round_var_nm_1L_chr, 
        round_bl_val_1L_chr = round_bl_val_1L_chr, predictors_lup = ingredients_ls$predictors_lup) %>% 
        add_utility_predn_to_ds(model_mdl = model_mdl, tfmn_1L_chr = tfmn_1L_chr, 
            depnt_var_nm_1L_chr = ingredients_ls$depnt_var_nm_1L_chr, 
            predr_vars_nms_chr = mdl_predr_terms_chr, force_min_max_1L_lgl = force_min_max_1L_lgl, 
            force_new_data_1L_lgl = T, impute_1L_lgl = T, is_brms_mdl_1L_lgl = inherits(model_mdl, 
                "brmsfit"), new_data_is_1L_chr = new_data_is_1L_chr, 
            predn_type_1L_chr = NULL, rmv_tfd_depnt_var_1L_lgl = T, 
            utl_cls_fn = utl_cls_fn, utl_min_val_1L_dbl = ingredients_ls$utl_min_val_1L_dbl, 
            sd_dbl = get_random_intercept(ingredients_ls$mdls_smry_tb, 
                mdl_nm_1L_chr = mdl_nm_1L_chr, deterministic_1L_lgl = deterministic_1L_lgl))
    if (!is.null(utl_var_nm_1L_chr)) {
        updated_tb <- updated_tb %>% dplyr::rename(`:=`(!!rlang::sym(utl_var_nm_1L_chr), 
            tidyselect::all_of(ingredients_ls$depnt_var_nm_1L_chr)))
    }
    if (!is.null(names(predr_vars_nms_chr))) {
        updated_tb <- rename_from_nmd_vec(updated_tb, nmd_vec_chr = predr_vars_nms_chr, 
            vec_nms_as_new_1L_lgl = F)
    }
    names_to_incl_chr <- c(names(updated_tb), setdiff(names(data_tb), 
        names(updated_tb)))
    rename_tb <- make_uid_rename_lup(data_tb, id_var_nm_1L_chr = id_var_nm_1L_chr)
    updated_tb <- dplyr::left_join(data_tb %>% dplyr::select(tidyselect::all_of(original_ds_vars_chr)) %>% 
        transform_uid_var(id_var_nm_1L_chr = id_var_nm_1L_chr, 
            rename_tb = rename_tb), updated_tb) %>% transform_uid_var(id_var_nm_1L_chr = id_var_nm_1L_chr, 
        rename_tb = rename_tb, old_new_chr = c("new_id_int", 
            "old_id_xx"))
    return(updated_tb)
}
