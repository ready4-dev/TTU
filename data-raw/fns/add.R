add_utility_predn_to_ds <- function (data_tb, model_mdl, tfmn_1L_chr, dep_var_nm_1L_chr, predr_vars_nms_chr = NULL,
                                     force_min_max_1L_lgl = T, utl_min_val_1L_dbl = 0.03, impute_1L_lgl = T, utl_cls_fn = NULL,
                                     rmv_tfmd_dep_var_1L_lgl = F)
{
    dep_vars_chr <- c(dep_var_nm_1L_chr, transform_dep_var_nm(dep_var_nm_1L_chr = dep_var_nm_1L_chr,
        tfmn_1L_chr = tfmn_1L_chr)) %>% unique()
    data_tb <- purrr::reduce(dep_vars_chr, .init = data_tb, ~dplyr::mutate(.x,
                                                                           !!rlang::sym(.y):= NA_real_))
    predictions_dbl <- predict_utility(data_tb = data_tb, tfmn_1L_chr = tfmn_1L_chr,
                                       model_mdl = model_mdl,
                                       force_min_max_1L_lgl = force_min_max_1L_lgl,
                                       utl_min_val_1L_dbl = utl_min_val_1L_dbl,
                                       impute_1L_lgl = impute_1L_lgl,
                                       utl_cls_fn = utl_cls_fn)
    data_tb <- data_tb %>% dplyr::mutate(!!rlang::sym(dep_var_nm_1L_chr):=predictions_dbl)
    if(!is.null(predr_vars_nms_chr)){
      data_tb <- data_tb %>% dplyr::select(-tidyselect::all_of(purrr::map(predr_vars_nms_chr, ~ paste0(.x,c("_baseline","_change"))) %>%
                                             purrr::flatten_chr()))
    }
    if(rmv_tfmd_dep_var_1L_lgl){
      data_tb <- data_tb %>% dplyr::select(-tidyselect::all_of(dep_vars_chr[dep_vars_chr!=dep_var_nm_1L_chr]))
    }
    return(data_tb)
}
